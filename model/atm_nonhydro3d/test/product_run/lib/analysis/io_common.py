"""Reading SCALE-DG output (per-rank netCDF files) into xarray.

SCALE-DG writes one file per MPI rank, ``<basename>.pe%06d.nc``, with
variables on ``(time, z, y, x)`` (history) or ``(z, y, x)`` (restart).
The coordinates are the physical node coordinates; LGL nodes include the
element end points, so the coordinate values on element boundaries are
duplicated.  Ranks are ordered with x fastest: ``rank = ix + NprcX * iy``.

Main entry points
-----------------
* :class:`MeshInfo`            mesh/decomposition parameters (from conf or files)
* :func:`open_dg_dataset`      read all or part of the ranks of one output
* :func:`open_dg_mfdataset`    concatenate several outputs along time
* :func:`to_element_form` / :func:`to_node_form`  conversion between the two forms
"""
from __future__ import annotations

import glob
import os
import re
import warnings

# Allow reading history files of runs that are still being written (HDF5 file locking
# would make netCDF4 fail with "HDF error").  HDF5 reads this variable when it is
# initialised, i.e. at the first file open in the process, so import this module
# before opening any netCDF file.  The last output step of a running job may be
# incomplete.  Set HDF5_USE_FILE_LOCKING explicitly beforehand to override.
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import xarray as xr

from . import DG_common as dgc

__all__ = [
    "MeshInfo", "read_namelist_group", "list_rank_files", "rank_extents",
    "open_dg_dataset", "open_dg_mfdataset", "find_runs", "open_dg_runs",
    "to_element_form", "to_node_form",
]

_DIR_INFO = (("x", "kx", "px"), ("y", "ky", "py"), ("z", "kz", "pz"))


# ---------------------------------------------------------------------------
# Namelist parsing
# ---------------------------------------------------------------------------

def _strip_comment(line: str) -> str:
    s = line.strip()
    if s.startswith("#"):
        return ""
    out, q = [], None
    for ch in line:
        if q is None and ch == "!":
            break
        if ch in ("'", '"'):
            q = None if q == ch else (ch if q is None else q)
        out.append(ch)
    return "".join(out)


def _split_values(s: str) -> List[str]:
    toks, cur, q = [], [], None
    for ch in s:
        if ch in ("'", '"'):
            q = None if q == ch else (ch if q is None else q)
        if ch == "," and q is None:
            toks.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    toks.append("".join(cur))
    return [t.strip() for t in toks if t.strip()]


def _parse_value(tok: str):
    t = tok.strip()
    if len(t) >= 2 and t[0] == t[-1] and t[0] in ("'", '"'):
        return t[1:-1]
    tl = t.lower()
    if tl in (".true.", "t", ".t."):
        return True
    if tl in (".false.", "f", ".f."):
        return False
    if re.fullmatch(r"[+-]?\d+", t):
        return int(t)
    try:
        return float(re.sub(r"[dD]", "e", t))
    except ValueError:
        return t


def read_namelist_group(path: str, group: str) -> Dict[str, object]:
    """Read one namelist group (``&GROUP ... /``) of a SCALE conf file into a dict.

    Keys are returned as written in the file.  Lists of values become Python lists.
    Lines starting with ``#`` and ``!`` comments are ignored.
    """
    with open(path) as f:
        text = "\n".join(_strip_comment(l) for l in f)
    m = re.search(r"&%s\b(.*?)^\s*/" % re.escape(group), text, re.S | re.M | re.I)
    if m is None:
        raise KeyError(f"namelist group &{group} not found in {path}")
    body = m.group(1)
    keys = list(re.finditer(r"([A-Za-z_]\w*)\s*=", body))
    out: Dict[str, object] = {}
    for i, km in enumerate(keys):
        end = keys[i + 1].start() if i + 1 < len(keys) else len(body)
        vals = [_parse_value(t) for t in _split_values(body[km.end():end])]
        out[km.group(1)] = vals[0] if len(vals) == 1 else vals
    return out


# ---------------------------------------------------------------------------
# Mesh information
# ---------------------------------------------------------------------------

@dataclass
class MeshInfo:
    """Global mesh and decomposition of a SCALE-DG regional (cube domain) run.

    ``NeX, NeY, NeZ`` are the *global* numbers of elements.
    """
    NeX: int
    NeY: int
    NeZ: int
    PolyOrder_h: int
    PolyOrder_v: int
    NprcX: int = 1
    NprcY: int = 1
    node_type: str = "LGL"
    isPeriodicX: bool = True
    isPeriodicY: bool = True
    dom_x: Optional[Tuple[float, float]] = None
    dom_y: Optional[Tuple[float, float]] = None
    dom_z: Optional[Tuple[float, float]] = None
    Fz: Optional[List[float]] = field(default=None)

    # --- derived quantities
    @property
    def Np_h(self) -> int:
        return self.PolyOrder_h + 1

    @property
    def Np_v(self) -> int:
        return self.PolyOrder_v + 1

    @property
    def NeX_lc(self) -> int:
        return self.NeX // self.NprcX

    @property
    def NeY_lc(self) -> int:
        return self.NeY // self.NprcY

    @property
    def nrank(self) -> int:
        return self.NprcX * self.NprcY

    def rank_index(self, rank: int) -> Tuple[int, int]:
        """(ix, iy) of a rank (x fastest)."""
        return rank % self.NprcX, rank // self.NprcX

    def rank_of(self, ix: int, iy: int) -> int:
        return ix + self.NprcX * iy

    def rank_layout(self) -> np.ndarray:
        """Array of rank numbers with shape (NprcY, NprcX)."""
        return np.arange(self.nrank).reshape(self.NprcY, self.NprcX)

    # --- constructors
    @classmethod
    def from_conf(cls, path: str, group: str = "PARAM_ATMOS_MESH", node_type: str = "LGL") -> "MeshInfo":
        """Build from the mesh namelist of a run/init conf.

        In the conf, ``NeX``/``NeY`` are numbers of elements *per process*; they are
        multiplied by ``NprcX``/``NprcY`` here.
        """
        nml = {k.lower(): v for k, v in read_namelist_group(path, group).items()}
        NprcX = int(nml.get("nprcx", 1))
        NprcY = int(nml.get("nprcy", 1))
        p_h = int(nml.get("polyorder_h"))
        p_v = int(nml.get("polyorder_v", p_h))

        def _rng(a, b):
            return (float(nml[a]), float(nml[b])) if a in nml and b in nml else None
        Fz = nml.get("fz")
        return cls(NeX=int(nml["nex"]) * NprcX, NeY=int(nml["ney"]) * NprcY, NeZ=int(nml["nez"]),
                   PolyOrder_h=p_h, PolyOrder_v=p_v, NprcX=NprcX, NprcY=NprcY, node_type=node_type,
                   isPeriodicX=bool(nml.get("isperiodicx", False)),
                   isPeriodicY=bool(nml.get("isperiodicy", False)),
                   dom_x=_rng("dom_xmin", "dom_xmax"), dom_y=_rng("dom_ymin", "dom_ymax"),
                   dom_z=_rng("dom_zmin", "dom_zmax"),
                   Fz=[float(v) for v in Fz] if isinstance(Fz, list) else None)

    @classmethod
    def from_files(cls, basename: str, isPeriodicX: bool = True, isPeriodicY: bool = True) -> "MeshInfo":
        """Infer the mesh from the coordinates of the rank files (LGL output only).

        The number of nodes per element is detected from the duplicated element-boundary
        coordinates.  Periodicity cannot be inferred and must be given.
        """
        files = list_rank_files(basename)
        ext = rank_extents(files)
        x0s = np.unique(np.round([e["x"][0] for e in ext.values()], 6))
        y0s = np.unique(np.round([e["y"][0] for e in ext.values()], 6))
        with xr.open_dataset(files[0], decode_times=False) as ds:
            cx, cy, cz = ds["x"].values, ds["y"].values, ds["z"].values
        Np_h = _detect_np(cx)
        Np_v = _detect_np(cz)
        return cls(NeX=len(cx) // Np_h * len(x0s), NeY=len(cy) // Np_h * len(y0s),
                   NeZ=len(cz) // Np_v, PolyOrder_h=Np_h - 1, PolyOrder_v=Np_v - 1,
                   NprcX=len(x0s), NprcY=len(y0s), node_type="LGL",
                   isPeriodicX=isPeriodicX, isPeriodicY=isPeriodicY,
                   dom_x=(min(e["x"][0] for e in ext.values()), max(e["x"][1] for e in ext.values())),
                   dom_y=(min(e["y"][0] for e in ext.values()), max(e["y"][1] for e in ext.values())),
                   dom_z=(float(cz[0]), float(cz[-1])))


def _detect_np(c: np.ndarray) -> int:
    """Nodes per element from duplicated coordinates (LGL); whole length if none."""
    dup = np.nonzero(np.isclose(np.diff(c), 0.0, atol=1e-9 * max(1.0, np.abs(c).max())))[0]
    return int(dup[0] + 1) if dup.size else len(c)


# ---------------------------------------------------------------------------
# Files
# ---------------------------------------------------------------------------

def list_rank_files(basename: str) -> List[str]:
    """Rank files ``<basename>.pe%06d.nc`` sorted by rank."""
    files = glob.glob(basename + ".pe*.nc")
    pat = re.compile(re.escape(basename) + r"\.pe(\d+)\.nc$")
    ranked = sorted((int(m.group(1)), f) for f in files for m in [pat.search(f)] if m)
    if not ranked:
        raise FileNotFoundError(f"no files match {basename}.pe*.nc")
    ranks = [r for r, _ in ranked]
    if ranks != list(range(len(ranks))):
        warnings.warn(f"rank numbers are not contiguous: {ranks[:5]}...")
    return [f for _, f in ranked]


def rank_extents(files: Sequence[str]) -> Dict[int, Dict[str, Tuple[float, float]]]:
    """{rank: {"x": (xmin, xmax), "y": (ymin, ymax)}} read from the coordinates only."""
    out = {}
    for r, f in enumerate(files):
        with xr.open_dataset(f, decode_times=False) as ds:
            out[r] = {d: (float(ds[d].values[0]), float(ds[d].values[-1])) for d in ("x", "y")}
    return out


# ---------------------------------------------------------------------------
# Element form <-> node form
# ---------------------------------------------------------------------------

def _split_variable(var: xr.Variable, dim: str, ne: int, np_: int, kname: str, pname: str) -> xr.Variable:
    if dim not in var.dims:
        return var
    ax = var.get_axis_num(dim)
    if var.shape[ax] != ne * np_:
        raise ValueError(f"size of '{dim}' ({var.shape[ax]}) != Ne*Np ({ne}*{np_})")
    shape = var.shape[:ax] + (ne, np_) + var.shape[ax + 1:]
    dims = var.dims[:ax] + (kname, pname) + var.dims[ax + 1:]
    return xr.Variable(dims, var.data.reshape(shape), attrs=var.attrs, encoding=var.encoding)


def _merge_variable(var: xr.Variable, kname: str, pname: str, dim: str) -> xr.Variable:
    if kname not in var.dims or pname not in var.dims:
        return var
    order = []
    for d in var.dims:
        if d == pname:
            continue
        order.append(d)
        if d == kname:
            order.append(pname)
    var = var.transpose(*order)
    ax = var.get_axis_num(kname)
    shape = var.shape[:ax] + (var.shape[ax] * var.shape[ax + 1],) + var.shape[ax + 2:]
    dims = var.dims[:ax] + (dim,) + var.dims[ax + 2:]
    return xr.Variable(dims, var.data.reshape(shape), attrs=var.attrs)


def to_element_form(ds: xr.Dataset, mesh: MeshInfo, k_offset: Optional[Dict[str, int]] = None,
                    check: bool = True) -> xr.Dataset:
    """Convert a node-form Dataset (one rank, or ranks already glued in node form) to element form.

    Parameters
    ----------
    ds : Dataset with dims ``x, y, z`` whose sizes are multiples of the nodes per element.
    mesh : mesh information (Np, node type).
    k_offset : global element index of the first element, e.g. ``{"x": 2, "y": 0}``.
    check : verify that element end points coincide (LGL only).
    """
    k_offset = k_offset or {}
    nps = {"x": mesh.Np_h, "y": mesh.Np_h, "z": mesh.Np_v}
    ref = {d: dgc.ref_nodes(nps[d], mesh.node_type) for d in nps}

    ne = {}
    for d, _, _ in _DIR_INFO:
        if d in ds.dims:
            if ds.sizes[d] % nps[d]:
                raise ValueError(f"size of '{d}' ({ds.sizes[d]}) is not a multiple of Np={nps[d]}")
            ne[d] = ds.sizes[d] // nps[d]

    new_vars = {}
    for name, v in ds.data_vars.items():
        var = v.variable
        for d, kn, pn in _DIR_INFO:
            if d in ne:
                var = _split_variable(var, d, ne[d], nps[d], kn, pn)
        new_vars[name] = var

    coords = {}
    for cname, c in ds.coords.items():
        if not (set(c.dims) & set(ne)):
            coords[cname] = c.variable
    for d, kn, pn in _DIR_INFO:
        if d not in ne:
            continue
        c = np.asarray(ds[d].values, dtype=float).reshape(ne[d], nps[d])
        xi = ref[d]
        width = (c[:, -1] - c[:, 0]) * 2.0 / (xi[-1] - xi[0])
        if check and mesh.node_type.upper() == "LGL" and ne[d] > 1:
            if not np.allclose(c[1:, 0], c[:-1, -1], rtol=1e-6, atol=1e-6 * np.abs(c).max()):
                raise ValueError(f"element end points along '{d}' do not coincide; check Np/node_type")
        if check:
            xi_file = 2.0 * (c - c[:, :1]) / (c[:, -1:] - c[:, :1]) - 1.0 if mesh.node_type.upper() == "LGL" else None
            if xi_file is not None and not np.allclose(xi_file, xi[None, :], atol=1e-5):
                warnings.warn(f"node positions along '{d}' differ from {mesh.node_type} nodes")
        k0 = int(k_offset.get(d, 0))
        coords[kn] = (kn, np.arange(k0, k0 + ne[d]))
        coords[pn] = (pn, np.arange(nps[d]))
        coords[dgc.REF_COORD[d]] = (pn, np.asarray(xi))
        coords[d] = ((kn, pn), c, dict(ds[d].attrs))
        coords["d" + d] = (kn, width, {"long_name": f"element width ({d})", "units": ds[d].attrs.get("units", "")})

    out = xr.Dataset(new_vars, coords=coords, attrs=dict(ds.attrs))
    out.attrs.update({"dg_form": "element", "dg_node_type": mesh.node_type,
                      "dg_PolyOrder_h": mesh.PolyOrder_h, "dg_PolyOrder_v": mesh.PolyOrder_v})
    return out


def to_node_form(obj: Union[xr.Dataset, xr.DataArray]) -> Union[xr.Dataset, xr.DataArray]:
    """Convert element form back to node form (dims ``x, y, z``; duplicated boundary nodes kept)."""
    if isinstance(obj, xr.DataArray):
        name = obj.name if obj.name is not None else "__da__"
        return to_node_form(obj.to_dataset(name=name))[name].rename(obj.name)
    ds = obj
    new_vars = {}
    for name, v in ds.data_vars.items():
        var = v.variable
        for d, kn, pn in _DIR_INFO:
            var = _merge_variable(var, kn, pn, d)
        new_vars[name] = var
    # only the directions that have both k<d> and p<d> are merged; coordinates of other
    # directions (e.g. xi(px) of a composite) are kept
    merged = {kn for _, kn, pn in _DIR_INFO if kn in ds.dims and pn in ds.dims} | \
             {pn for _, kn, pn in _DIR_INFO if kn in ds.dims and pn in ds.dims}
    coords = {c: ds[c].variable for c in ds.coords if not (set(ds[c].dims) & merged)}
    for d, kn, pn in _DIR_INFO:
        if d in ds.coords and set(ds[d].dims) == {kn, pn}:
            coords[d] = (d, ds[d].transpose(kn, pn).values.ravel(), dict(ds[d].attrs))
    out = xr.Dataset(new_vars, coords=coords, attrs=dict(ds.attrs))
    out.attrs["dg_form"] = "node"
    return out


# ---------------------------------------------------------------------------
# Opening outputs
# ---------------------------------------------------------------------------

def _select_ranks(mesh: MeshInfo, files: Sequence[str], ranks, region) -> List[int]:
    if ranks is not None:
        sel = sorted(int(r) for r in ranks)
    else:
        sel = list(range(len(files)))
    if region:
        ext = rank_extents(files)
        def _ok(r):
            for d, (lo, hi) in region.items():
                if d not in ("x", "y"):
                    continue
                a, b = ext[r][d]
                if not (a < hi and b > lo):
                    return False
            return True
        sel = [r for r in sel if _ok(r)]
    if not sel:
        raise ValueError("no rank selected")
    return sel


def _rank_grid(mesh: MeshInfo, sel: Sequence[int]) -> List[List[int]]:
    ixs = sorted({mesh.rank_index(r)[0] for r in sel})
    iys = sorted({mesh.rank_index(r)[1] for r in sel})
    grid = [[mesh.rank_of(ix, iy) for ix in ixs] for iy in iys]
    if sorted(r for row in grid for r in row) != sorted(sel):
        raise ValueError(f"selected ranks {list(sel)} do not form a rectangle in the (ix, iy) layout")
    return grid


def open_dg_dataset(basename: str, mesh: Optional[MeshInfo] = None,
                    ranks: Optional[Sequence[int]] = None,
                    region: Optional[Dict[str, Tuple[float, float]]] = None,
                    variables: Optional[Sequence[str]] = None,
                    chunks: Optional[dict] = {},
                    form: str = "element",
                    as_float64: bool = True,
                    decode_times: bool = True,
                    check: bool = True) -> xr.Dataset:
    """Open all or part of the rank files of one SCALE-DG output.

    Parameters
    ----------
    basename : path without ``.peXXXXXX.nc``, e.g. ``"run_spinup1/history"``.
    mesh : mesh information; inferred from the files (LGL) if None.
    ranks : ranks to read (must form a rectangle in the rank layout).  Default: all.
    region : physical range, e.g. ``{"x": (0, 2e6), "y": (1e6, 3e6)}``.  Only the ranks
        overlapping it are opened, and the result is cropped to the overlapping elements.
    variables : data variables to keep (``time_bnds`` is kept when present).
    chunks : passed to ``xr.open_dataset``.  The default ``{}`` keeps data lazy (dask);
        use ``None`` to load into memory.
    form : ``"element"`` (default) or ``"node"``.
    as_float64 : cast float32 variables to float64.
    decode_times : decode ``time`` into cftime objects (SCALE-DG runs start in year 0001,
        which numpy datetime64 cannot represent).

    Notes
    -----
    Surface fluxes such as ``RAIN`` are written with the downward direction negative
    (RAIN <= 0); flip the sign for precipitation rates.
    check : verify element boundaries when converting to element form.

    Returns
    -------
    Dataset whose element indices ``kx, ky, kz`` are global.
    """
    files = list_rank_files(basename)
    if mesh is None:
        mesh = MeshInfo.from_files(basename)
    if len(files) != mesh.nrank:
        raise ValueError(f"{len(files)} files found but the mesh has {mesh.nrank} ranks")
    sel = _select_ranks(mesh, files, ranks, region)
    grid = _rank_grid(mesh, sel)

    def _open(r):
        with warnings.catch_warnings():
            # year-0001 dates cannot be numpy datetime64: use cftime silently
            warnings.filterwarnings("ignore", message=".*year zero convention.*")
            warnings.filterwarnings("ignore", message=".*Unable to decode time axis.*")
            try:
                ds = xr.open_dataset(files[r], chunks=chunks, decode_times=decode_times,
                                     use_cftime=True if decode_times else None)
            except OSError as e:
                raise OSError(f"cannot open {files[r]} ({e}). If the run is still being written, "
                              "set HDF5_USE_FILE_LOCKING=FALSE before the first netCDF file is opened "
                              "(importing analysis.io_common first does this).") from e
        if variables is not None:
            keep = [v for v in variables if v in ds.data_vars]
            missing = set(variables) - set(keep)
            if missing:
                raise KeyError(f"variables not found: {sorted(missing)}")
            if "time_bnds" in ds.data_vars and "time_bnds" not in keep:
                keep.append("time_bnds")
            ds = ds[keep]
        if as_float64:
            for v in ds.data_vars:
                if ds[v].dtype == np.float32:
                    ds[v] = ds[v].astype(np.float64)
        ix, iy = mesh.rank_index(r)
        return to_element_form(ds, mesh, k_offset={"x": ix * mesh.NeX_lc, "y": iy * mesh.NeY_lc},
                               check=check)

    dsgrid = [[_open(r) for r in row] for row in grid]
    # ranks of a running job may have written different numbers of steps: keep the common ones
    if "time" in dsgrid[0][0].dims:
        nts = {d.sizes["time"] for row in dsgrid for d in row}
        if len(nts) > 1:
            nt = min(nts)
            warnings.warn(f"{basename}: ranks have different numbers of time steps {sorted(nts)}; "
                          f"using the first {nt} (run still in progress?)")
            dsgrid = [[d.isel(time=slice(0, nt)) for d in row] for row in dsgrid]
    if len(sel) == 1:
        ds = dsgrid[0][0]
    else:
        ds = xr.combine_nested(dsgrid, concat_dim=["ky", "kx"], data_vars="minimal",
                               coords="minimal", compat="override", combine_attrs="override")
    ds.attrs.pop("rankid", None)
    ds.attrs.update({"dg_ranks": ",".join(str(r) for r in sel), "dg_basename": os.path.abspath(basename),
                     "dg_NeX": mesh.NeX, "dg_NeY": mesh.NeY, "dg_NeZ": mesh.NeZ,
                     "dg_isPeriodicX": int(mesh.isPeriodicX), "dg_isPeriodicY": int(mesh.isPeriodicY)})

    if region:
        isel = {}
        for d, (lo, hi) in region.items():
            if d not in ("x", "y", "z"):
                raise ValueError(f"unknown region key '{d}'")
            kn, pn = dgc.ELEM_DIMS[d]
            c = ds[d].transpose(kn, pn).values
            isel[kn] = np.nonzero((c[:, 0] < hi) & (c[:, -1] > lo))[0]
        ds = ds.isel(isel)

    if form == "node":
        return to_node_form(ds)
    if form != "element":
        raise ValueError("form must be 'element' or 'node'")
    return ds


def open_dg_mfdataset(basenames: Sequence[str], mesh: Optional[MeshInfo] = None,
                      overlap: str = "later", labels: Optional[Sequence[str]] = None,
                      **kwargs) -> xr.Dataset:
    """Open several outputs (e.g. chained runs) and concatenate them along ``time``.

    Parameters
    ----------
    basenames : outputs in chronological order of their start times.
    overlap : how to treat overlapping periods, e.g. a run restarted from day 25 of a
        previous run that went on until day 30.

        * ``"later"`` (default): the later run wins; each run is cut just before the
          first time of the next run.
        * ``"earlier"``: the earlier run wins; times of a later run that are not after
          the end of the previous runs are dropped.
    labels : optional name of each output; stored as the coordinate ``run(time)``.
    Other keyword arguments go to :func:`open_dg_dataset`.
    """
    if overlap not in ("later", "earlier"):
        raise ValueError("overlap must be 'later' or 'earlier'")
    if mesh is None:
        mesh = MeshInfo.from_files(basenames[0])
    dss = [open_dg_dataset(b, mesh=mesh, **kwargs) for b in basenames]
    if labels is not None:
        if len(labels) != len(dss):
            raise ValueError("len(labels) != len(basenames)")
        dss = [d.assign_coords(run=("time", [str(l)] * d.sizes["time"])) for d, l in zip(dss, labels)]
    if overlap == "later":
        for i in range(len(dss) - 1):
            t_next = dss[i + 1]["time"].values[0]
            dss[i] = dss[i].isel(time=np.nonzero(dss[i]["time"].values < t_next)[0])
    else:
        t_end = dss[0]["time"].values[-1]
        for i in range(1, len(dss)):
            dss[i] = dss[i].isel(time=np.nonzero(dss[i]["time"].values > t_end)[0])
            if dss[i].sizes["time"]:
                t_end = dss[i]["time"].values[-1]
    dss = [d for d in dss if d.sizes["time"]]
    ds = xr.concat(dss, dim="time", data_vars="minimal", coords="minimal",
                   compat="override", combine_attrs="override")
    if not ds.get_index("time").is_monotonic_increasing:
        raise ValueError("concatenated time axis is not increasing; check the order of basenames")
    ds.attrs["dg_basename"] = ",".join(os.path.abspath(b) for b in basenames)
    return ds


def _run_sort_key(name: str):
    m = re.search(r"(\d+)$", name)
    return (name[:m.start()] if m else name, int(m.group(1)) if m else -1)


def find_runs(case_dir: str, runs=None, history: str = "history", pattern: str = "run*",
              skip_missing: bool = True) -> List[str]:
    """Run directories of a case that contain history output, in run order.

    Parameters
    ----------
    case_dir : case directory containing ``run1, run2, ...`` (e.g. ``"../asis"``).
    runs : None (all directories matching ``pattern``), or a list of run numbers
        (``[1, 2]``) / names (``["run3", "run4"]``) / paths.
    history : basename of the history files in each run directory.
    skip_missing : skip runs without history files (not started yet) with a warning;
        if False, raise FileNotFoundError.

    Returns
    -------
    list of run directories (paths).  Runs are sorted by the trailing number of their
    name (run2 < run10) when ``runs`` is None; otherwise the given order is kept.
    """
    if runs is None:
        dirs = [d for d in glob.glob(os.path.join(case_dir, pattern)) if os.path.isdir(d)]
        dirs.sort(key=lambda d: _run_sort_key(os.path.basename(d)))
    else:
        dirs = []
        for r in runs:
            if isinstance(r, (int, np.integer)):
                dirs.append(os.path.join(case_dir, f"run{int(r)}"))
            elif os.path.isdir(str(r)):
                dirs.append(str(r))
            else:
                dirs.append(os.path.join(case_dir, str(r)))
    out, missing = [], []
    for d in dirs:
        (out if glob.glob(os.path.join(d, history + ".pe*.nc")) else missing).append(d)
    if missing:
        msg = f"no {history} output in: {', '.join(os.path.basename(m) for m in missing)}"
        if not skip_missing:
            raise FileNotFoundError(msg)
        warnings.warn(msg + " (skipped)")
    if not out:
        raise FileNotFoundError(f"no run with {history} output under {case_dir}")
    return out


def open_dg_runs(case_dir: Optional[str] = None, runs=None, mesh: Optional[MeshInfo] = None,
                 history: str = "history", run_dirs: Optional[Sequence[str]] = None,
                 overlap: str = "later", skip_missing: bool = True, **kwargs) -> xr.Dataset:
    """Open the history of several chained runs (e.g. ``asis/run1 .. run4``) as one Dataset.

    Parameters
    ----------
    case_dir, runs : see :func:`find_runs`.
    mesh : mesh information.  Default: from ``run.conf`` of the first run if present,
        otherwise inferred from the files.
    run_dirs : explicit list of run directories in chronological order (instead of
        ``case_dir``/``runs``), e.g. ``["../run_spinup1", "../run_spinup2"]``.
    overlap : see :func:`open_dg_mfdataset`.
    skip_missing : skip runs that have no output yet.
    Other keyword arguments (``variables``, ``region``, ``chunks``, ...) go to
    :func:`open_dg_dataset`.

    Returns
    -------
    Dataset concatenated along ``time`` with the coordinate ``run(time)`` (run directory
    name of each time) and the attribute ``dg_runs``.
    """
    if run_dirs is None:
        if case_dir is None:
            raise ValueError("give case_dir or run_dirs")
        run_dirs = find_runs(case_dir, runs, history=history, skip_missing=skip_missing)
    else:
        run_dirs = find_runs("", list(run_dirs), history=history, skip_missing=skip_missing)
    if mesh is None:
        conf = os.path.join(run_dirs[0], "run.conf")
        mesh = MeshInfo.from_conf(conf) if os.path.exists(conf) else \
            MeshInfo.from_files(os.path.join(run_dirs[0], history))
    labels = [os.path.basename(os.path.normpath(d)) for d in run_dirs]
    ds = open_dg_mfdataset([os.path.join(d, history) for d in run_dirs], mesh,
                           overlap=overlap, labels=labels, **kwargs)
    ds.attrs["dg_runs"] = ",".join(labels)
    if case_dir is not None:
        ds.attrs["dg_case"] = os.path.basename(os.path.normpath(os.path.abspath(case_dir)))
    return ds


# ---------------------------------------------------------------------------
# Self test (uses moist_conv/Eh8Ez4P7 output if available)
# ---------------------------------------------------------------------------

def _selftest(case_dir: Optional[str] = None, verbose: bool = True) -> None:
    here = os.path.dirname(os.path.abspath(__file__))
    case_dir = case_dir or os.path.join(here, "..", "..", "moist_conv", "Eh8Ez4P7")
    conf = os.path.join(case_dir, "init.conf")
    hist = os.path.join(case_dir, "run_spinup1", "history")
    if not (os.path.exists(conf) and glob.glob(hist + ".pe*.nc")):
        print("io_common: test data not found, skipped")
        return
    mesh = MeshInfo.from_conf(conf)
    mesh_f = MeshInfo.from_files(hist)
    for a in ("NeX", "NeY", "NeZ", "PolyOrder_h", "PolyOrder_v", "NprcX", "NprcY"):
        assert getattr(mesh, a) == getattr(mesh_f, a), a
    assert mesh.isPeriodicX and mesh.Fz is not None and len(mesh.Fz) == mesh.NeZ + 1

    ds = open_dg_dataset(hist, mesh, variables=["RAIN", "T"])
    assert ds.sizes["kx"] == mesh.NeX and ds.sizes["ky"] == mesh.NeY and ds.sizes["kz"] == mesh.NeZ
    assert ds["T"].dims == ("time", "kz", "pz", "ky", "py", "kx", "px")
    assert ds["RAIN"].dims == ("time", "ky", "py", "kx", "px")
    assert np.allclose(ds["dx"].values, 1.0e6) and np.allclose(ds["x"].values[1:, 0], ds["x"].values[:-1, -1])
    # node form round trip, compared with a direct read of one rank
    nd = to_node_form(ds.isel(time=-1))
    r5 = xr.open_dataset(hist + ".pe000005.nc", decode_times=False).isel(time=-1)
    ix, iy = mesh.rank_index(5)
    sub = nd["T"].isel(x=slice(ix * 16, (ix + 1) * 16), y=slice(iy * 32, (iy + 1) * 32))
    assert np.array_equal(sub.values, r5["T"].values.astype(np.float64))
    # partial read by region and by ranks
    dr = open_dg_dataset(hist, mesh, region={"x": (2.5e6, 3.5e6), "y": (0.0, 1.0e6)}, variables=["RAIN"])
    assert list(dr["kx"].values) == [2, 3] and list(dr["ky"].values) == [0]
    dk = open_dg_dataset(hist, mesh, ranks=[1, 2, 5, 6], variables=["RAIN"])
    assert list(dk["kx"].values) == [2, 3, 4, 5] and dk.sizes["ky"] == 8
    assert np.array_equal(dk["RAIN"].isel(time=-1).values,
                          ds["RAIN"].isel(time=-1, kx=slice(2, 6)).values)
    # chained runs: spinup2 restarts at day 25 of spinup1 (which ran to day 30)
    hist2 = os.path.join(case_dir, "run_spinup2", "history")
    if glob.glob(hist2 + ".pe*.nc"):
        mf = open_dg_mfdataset([hist, hist2], mesh, variables=["RAIN"])
        t = mf.get_index("time")
        assert t.is_monotonic_increasing and not t.duplicated().any()
        mfe = open_dg_mfdataset([hist, hist2], mesh, variables=["RAIN"], overlap="earlier")
        assert mfe.get_index("time").is_monotonic_increasing
    # several runs through open_dg_runs (explicit run directories)
    if glob.glob(hist2 + ".pe*.nc"):
        rr = open_dg_runs(run_dirs=[os.path.dirname(hist), os.path.dirname(hist2)], variables=["RAIN"])
        assert rr.get_index("time").is_monotonic_increasing
        assert list(np.unique(rr["run"].values)) == ["run_spinup1", "run_spinup2"]
        assert rr.sizes["time"] == mf.sizes["time"]
    # domain mean through DG_common
    m = dgc.mean(ds["RAIN"].isel(time=-1), "xy")
    assert np.isfinite(float(m))
    if verbose:
        print(f"io_common: selftest passed (domain-mean RAIN at last step = {float(m):.3e} kg/m2/s)")


if __name__ == "__main__":
    _selftest()
