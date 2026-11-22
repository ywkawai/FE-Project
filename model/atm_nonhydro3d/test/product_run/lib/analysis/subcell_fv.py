"""Conversion of DG nodal data (LGL) to sub-cell finite-volume data.

Definition (toy-model PreFV8; moist_conv/NOTES.md Sec. 7.1):
each element is divided into Nc sub-cells and the sub-cell mean of the element
polynomial is computed with exact quadrature.  The inverse operation is the L2
projection of the piecewise-constant function back onto the DG space.

* uniform=True  : equal-width sub-cells (``SubCellFV_uniform_flag = .true.``, inferred)
* uniform=False : sub-cell widths proportional to the Nc-point LGL weights
* :func:`subcell_to_dg` corresponds to ``PostSubCellFV_L2Proj_flag`` (inferred)

Element form of the sub-cell data uses the dimension ``c<d>`` (``cx``, ``cy``,
``cz``) in place of ``p<d>``; the coordinates are the cell centres ``<d>(k<d>, c<d>)``,
the reference cell centres ``<ref>_c(c<d>)`` and the cell bounds ``<d>_lbnd(k<d>, c<d>)``,
``<d>_ubnd(k<d>, c<d>)``.  :func:`to_fv_grid` flattens ``(k<d>, c<d>)`` into a plain FV
grid ``<d>_fv`` (coordinates ``<d>_fv``, ``<d>_fv_lbnd``, ``<d>_fv_ubnd``).
"""
from __future__ import annotations

from functools import lru_cache
from math import ceil
from typing import Callable, Iterable, Tuple, Union

import numpy as np
import xarray as xr

from . import DG_common as dgc

__all__ = ["subcell_edges", "subcell_average_matrix", "subcell_l2proj_matrix",
           "to_subcell", "subcell_to_dg", "to_fv_grid", "apply_pointwise"]

CELL_DIMS = {"x": "cx", "y": "cy", "z": "cz"}


@lru_cache(maxsize=None)
def subcell_edges(Nc: int, uniform: bool = True) -> np.ndarray:
    """Sub-cell edges (Nc+1,) in reference coordinates [-1, 1]."""
    if Nc < 1:
        raise ValueError("Nc >= 1")
    if uniform or Nc == 1:
        e = np.linspace(-1.0, 1.0, Nc + 1)
    else:
        w = dgc.lgl_nodes(Nc)[1]
        e = np.concatenate(([-1.0], -1.0 + np.cumsum(w)))
        e[-1] = 1.0
    e.setflags(write=False)
    return e


def _cell_quadrature(edges: np.ndarray, nq: int) -> Tuple[np.ndarray, np.ndarray]:
    """GL points (Nc, nq) and weights (Nc, nq) of each cell, weights normalised to the cell mean."""
    r, w = dgc.gl_nodes(nq)
    a, b = edges[:-1, None], edges[1:, None]
    x = a + 0.5 * (b - a) * (1.0 + r[None, :])
    return x, np.broadcast_to(0.5 * w[None, :], x.shape)


@lru_cache(maxsize=None)
def subcell_average_matrix(Np: int, Nc: int, uniform: bool = True, node_type: str = "LGL") -> np.ndarray:
    """Matrix S (Nc x Np): sub-cell means of the degree-(Np-1) element polynomial (exact)."""
    xi = dgc.ref_nodes(Np, node_type)
    edges = subcell_edges(Nc, uniform)
    nq = int(ceil(Np / 2.0)) + 1
    xq, wq = _cell_quadrature(edges, nq)
    S = np.einsum("cq,cqj->cj", wq, dgc.lagrange_matrix(xi, xq.ravel()).reshape(Nc, nq, Np))
    S.setflags(write=False)
    return S


@lru_cache(maxsize=None)
def subcell_l2proj_matrix(Np: int, Nc: int, uniform: bool = True, node_type: str = "LGL") -> np.ndarray:
    """Matrix P (Np x Nc): L2 projection of the piecewise-constant sub-cell data onto the
    degree-(Np-1) polynomials, returned as nodal values.  Conserves the element integral."""
    xi = dgc.ref_nodes(Np, node_type)
    edges = subcell_edges(Nc, uniform)
    nq = int(ceil(Np / 2.0)) + 1
    xq, wq = _cell_quadrature(edges, nq)
    h = np.diff(edges)
    # B[m, c] = integral over cell c of the orthonormal P_m
    Vq = dgc.legendre_vandermonde(xq.ravel(), Np).reshape(Nc, nq, Np)
    B = np.einsum("cq,cqm->mc", wq * h[:, None], Vq)
    P = dgc.legendre_vandermonde(xi, Np) @ B
    P.setflags(write=False)
    return P


def to_subcell(obj, Nc: int, dims: Union[str, Iterable[str]] = "xy", uniform: bool = True):
    """DG nodal data (element form) -> sub-cell means, dims ``p<d>`` -> ``c<d>``."""
    dirs = dgc._parse_dirs(dims)

    def _f(da):
        for d in dirs:
            kdim, pdim = dgc.ELEM_DIMS[d]
            Np = da.sizes[pdim]
            S = subcell_average_matrix(Np, Nc, uniform)
            edges = subcell_edges(Nc, uniform)
            has_phys = d in da.coords and set(da[d].dims) == {kdim, pdim}
            if has_phys and "d" + d not in da.coords:
                da = da.assign_coords({"d" + d: dgc.element_width(da, d)})
            if has_phys:
                xi = dgc.ref_nodes_of(da, d)
                c = da[d].transpose(kdim, pdim).values
                left = c[:, 0] - (xi[0] + 1.0) * 0.5 * da["d" + d].values
            da = dgc.apply_along(da, S, d, out_pdim=CELL_DIMS[d],
                                 out_ref=(dgc.REF_COORD[d] + "_c", 0.5 * (edges[:-1] + edges[1:])))
            if has_phys:
                w = da["d" + d].values
                lb = left[:, None] + 0.5 * (edges[None, :-1] + 1.0) * w[:, None]
                ub = left[:, None] + 0.5 * (edges[None, 1:] + 1.0) * w[:, None]
                da = da.assign_coords({d + "_lbnd": ((kdim, CELL_DIMS[d]), lb),
                                       d + "_ubnd": ((kdim, CELL_DIMS[d]), ub)})
            da = da.assign_coords({CELL_DIMS[d]: np.arange(Nc)})
        da.attrs["dg_op"] = (da.attrs.get("dg_op", "") +
                             f" subcellFV_{Nc}{'' if uniform else 'lgl'}_{''.join(dirs)}").strip()
        da.attrs["subcell_uniform"] = int(uniform)
        return da
    return dgc.map_dataset(_f, obj, dirs)


def subcell_to_dg(obj, Np: int, dims: Union[str, Iterable[str]] = "xy", uniform: bool = True,
                  node_type: str = "LGL"):
    """Sub-cell data (dims ``c<d>``) -> DG nodal data (dims ``p<d>``) by L2 projection."""
    dirs = dgc._parse_dirs(dims)
    xi = dgc.ref_nodes(Np, node_type)

    def _f(da):
        for d in dirs:
            kdim, pdim = dgc.ELEM_DIMS[d]
            cdim = CELL_DIMS[d]
            Nc = da.sizes[cdim]
            P = subcell_l2proj_matrix(Np, Nc, uniform, node_type)
            drop = [c for c in da.coords if cdim in da[c].dims]
            phys_left = None
            if d + "_lbnd" in da.coords and "d" + d in da.coords:
                phys_left = da[d + "_lbnd"].transpose(kdim, cdim).values[:, 0]
            base = da.drop_vars(drop)
            # rename the cell dim to the node dim name expected by apply_along
            base = base.rename({cdim: pdim})
            out = dgc.apply_along(base, P, d, transform_coord=False, out_ref=(dgc.REF_COORD[d], xi))
            if phys_left is not None:
                w = da["d" + d].values
                out = out.assign_coords({d: ((kdim, pdim), phys_left[:, None] + 0.5 * (xi[None, :] + 1.0) * w[:, None])})
            out = out.assign_coords({pdim: np.arange(Np)})
            da = out
        return da
    need = {CELL_DIMS[d] for d in dirs}

    if isinstance(obj, xr.DataArray):
        return _f(obj)
    out = {}
    for name, var in obj.data_vars.items():
        if need <= set(var.dims):
            out[name] = _f(var)
        elif not (set(var.dims) & (need | {dgc.ELEM_DIMS[d][0] for d in dirs})):
            out[name] = var
    return xr.Dataset(out, attrs=dict(obj.attrs))


def to_fv_grid(obj, dims: Union[str, Iterable[str]] = "xy"):
    """Flatten sub-cell element form ``(k<d>, c<d>)`` into an FV grid dimension ``<d>_fv``
    with coordinates ``<d>_fv`` (cell centres), ``<d>_fv_lbnd`` and ``<d>_fv_ubnd``.
    No duplicated points."""
    dirs = dgc._parse_dirs(dims)

    def _f(da):
        for d in dirs:
            kdim, cdim = dgc.ELEM_DIMS[d][0], CELL_DIMS[d]
            fv = d + "_fv"
            centres = da[d].transpose(kdim, cdim).values.ravel() if d in da.coords else None
            bnds = ({b: da[d + "_" + b].transpose(kdim, cdim).values.ravel() for b in ("lbnd", "ubnd")}
                    if d + "_lbnd" in da.coords else None)
            drop = [c for c in da.coords if (kdim in da[c].dims or cdim in da[c].dims)]
            base = da.drop_vars(drop)
            order = [dd for dd in base.dims if dd != cdim]
            ax = order.index(kdim)
            order.insert(ax + 1, cdim)
            var = base.transpose(*order).variable
            data = var.data.reshape(var.shape[:ax] + (var.shape[ax] * var.shape[ax + 1],) + var.shape[ax + 2:])
            st = xr.DataArray(data, dims=tuple(order[:ax]) + (fv,) + tuple(order[ax + 2:]),
                              coords={c: base.coords[c] for c in base.coords},
                              attrs=dict(da.attrs), name=da.name)
            if centres is not None:
                st = st.assign_coords({fv: centres})
            if bnds is not None:
                st = st.assign_coords({fv + "_" + b: (fv, v) for b, v in bnds.items()})
            da = st
        return da

    if isinstance(obj, xr.DataArray):
        return _f(obj)
    need = {CELL_DIMS[d] for d in dirs}
    out = {n: _f(v) for n, v in obj.data_vars.items() if need <= set(v.dims)}
    return xr.Dataset(out, attrs=dict(obj.attrs))


def apply_pointwise(func: Callable[..., xr.DataArray], *das: xr.DataArray, Nc: int,
                    dims: Union[str, Iterable[str]] = "xy", uniform: bool = True,
                    back_to_dg: bool = True) -> xr.DataArray:
    """Evaluate a (nonlinear) pointwise function on sub-cell means.

    ``func(*subcell_means)`` is evaluated on the sub-cell means of ``das``; the result
    is projected back to DG nodes (L2) when ``back_to_dg`` is True, as in PreFV8.
    """
    dirs = dgc._parse_dirs(dims)
    Np = {d: das[0].sizes[dgc.ELEM_DIMS[d][1]] for d in dirs}
    qbar = [to_subcell(da, Nc, dirs, uniform) for da in das]
    res = func(*qbar)
    if not back_to_dg:
        return res
    out = res
    for d in dirs:
        out = subcell_to_dg(out, Np[d], d, uniform)
    return out


# ---------------------------------------------------------------------------
# Self test
# ---------------------------------------------------------------------------

def _selftest(verbose: bool = True) -> None:
    Np = 8
    xi = dgc.lgl_nodes(Np)[0]
    for uniform in (True, False):
        for Nc in (1, 4, 8, 12):
            S = subcell_average_matrix(Np, Nc, uniform)
            P = subcell_l2proj_matrix(Np, Nc, uniform)
            e = subcell_edges(Nc, uniform)
            h = np.diff(e)
            # exact cell means of x^7
            exact = (e[1:] ** 8 - e[:-1] ** 8) / 8.0 / h
            assert np.allclose(S @ xi ** 7, exact, atol=1e-12)
            # conservation: element integral of the projection = sum of h * qbar
            qbar = np.random.default_rng(0).normal(size=Nc)
            w = dgc.quad_weights(xi)
            assert abs(w @ (P @ qbar) - h @ qbar) < 1e-12
            # constants are reproduced by the round trip (P S is not the identity in general)
            assert np.allclose(P @ (S @ np.ones(Np)), 1.0, atol=1e-12)
            if Nc == 1:
                assert np.allclose(S @ P, 1.0)

    # xarray path
    f = lambda X, Y: X ** 3 + Y
    da = dgc._make_test_array(NeX=4, NeY=3, Np=8, Lx=4.0, Ly=3.0, func=f)
    sc = to_subcell(da, 8, "xy")
    assert sc.dims == ("ky", "cy", "kx", "cx")
    assert np.allclose(sc["x"].values, np.arange(4)[:, None] + (np.arange(8)[None, :] + 0.5) / 8.0)
    lb, ub = sc["x_lbnd"].values, sc["x_ubnd"].values
    exact = (((ub ** 4 - lb ** 4) / 4.0 / (ub - lb))[None, None] + sc["y"].values[:, :, None, None])
    assert np.allclose(sc.values, exact, atol=1e-12)
    # element means / domain integral conserved
    assert np.allclose(sc.mean(["cx", "cy"]).values, dgc.element_mean(da, "xy").values, atol=1e-12)
    back = subcell_to_dg(sc, 8, "xy")
    assert back.dims == da.dims
    assert np.allclose(dgc.element_mean(back, "xy").values, dgc.element_mean(da, "xy").values, atol=1e-12)
    # the round trip is a smoothing (not the identity) but converges with Nc
    err8 = float(np.abs(back - da).max())
    err64 = float(np.abs(subcell_to_dg(to_subcell(da, 64, "xy"), 8, "xy") - da).max())
    assert err64 < 0.05 * err8, (err8, err64)
    assert np.allclose(back["x"].values, da["x"].values)
    fv = to_fv_grid(sc)
    assert fv.dims == ("y_fv", "x_fv") and fv.sizes["x_fv"] == 32
    assert np.all(np.diff(fv["x_fv"].values) > 0)
    # nonlinear pointwise evaluation conserves the integral of func(qbar)
    r = apply_pointwise(lambda q: np.maximum(q - 1.0, 0.0), da, Nc=8)
    ref = np.maximum(sc - 1.0, 0.0).mean(["cx", "cy"])
    assert np.allclose(dgc.element_mean(r, "xy").values, ref.values, atol=1e-12)
    if verbose:
        print("subcell_fv: selftest passed")


if __name__ == "__main__":
    _selftest()
