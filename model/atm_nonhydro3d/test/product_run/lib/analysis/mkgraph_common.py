"""Generic plotting helpers for element-form DG data (matplotlib).

The functions do not compute any physics-specific quantity: reduce the data
first with :mod:`DG_common` (``mean``, ``integrate``, ``composite``, ...) so that
only the dimensions to be drawn remain, then call

* :func:`plot_tone`          2-D tone (colour) plot over two directions, e.g. x-z or x-y
* :func:`plot_line`          1-D line plot along one direction (vertical profile for z)
* :func:`plot_composite2d`   element composite on the reference coordinates (xi, eta, ...)
* :func:`plot_composite1d`   1-D element composite (other node dims averaged in the element)
* :func:`plot_composite_summary`  2-D composite + the two 1-D profiles in one figure

Common conventions
------------------
* Directions are ``"x"``, ``"y"``, ``"z"``; element-form dims ``(k<d>, p<d>)``.
* ``scale`` converts physical coordinates for the axes (default 1e-3: m -> km);
  a float or a dict per direction.
* Labels are taken from ``long_name`` / ``units`` attributes unless given.
* Every function takes ``ax`` (created if None) and returns the Axes (and the
  mappable / line objects where relevant), so figures can be composed freely.
"""
from __future__ import annotations

from typing import Dict, Iterable, Optional, Sequence, Tuple, Union

import numpy as np
import xarray as xr

from . import DG_common as dgc
from . import io_common as io

__all__ = ["element_edges", "add_element_edges", "auto_label",
           "plot_tone", "plot_line", "plot_composite2d", "plot_composite1d",
           "plot_composite_summary"]

Scale = Union[float, Dict[str, float]]
_UNIT_PREFIX = {1.0: "", 1e-3: "k", 1e-6: "M", 1e3: "m"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _plt():
    import matplotlib.pyplot as plt
    return plt


def _get_ax(ax):
    if ax is None:
        _, ax = _plt().subplots()
    return ax


def _scale_of(scale: Scale, d: str) -> float:
    if isinstance(scale, dict):
        return float(scale.get(d, 1.0))
    return float(scale)


def _coord_label(obj, d: str, scale: Scale) -> str:
    s = _scale_of(scale, d)
    units = obj[d].attrs.get("units", "") if d in obj.coords else ""
    if units == "m" and s in _UNIT_PREFIX:
        return f"{d} [{_UNIT_PREFIX[s]}m]"
    if s == 1.0:
        return f"{d} [{units}]" if units else d
    return f"{d} [{units} × {s:g}]" if units else f"{d} (× {s:g})"


def auto_label(da: xr.DataArray, label: Optional[str] = None) -> str:
    """Label from ``long_name`` (or name) and ``units`` attributes."""
    if label is not None:
        return label
    name = da.attrs.get("long_name", da.name if da.name is not None else "")
    units = da.attrs.get("units", "")
    return f"{name} [{units}]" if units else str(name)


def _check_dims(da: xr.DataArray, keep: Sequence[str], hue: Optional[str] = None) -> xr.DataArray:
    allowed = set(keep) | ({hue} if hue else set())
    extra = [d for d in da.dims if d not in allowed]
    squeeze = [d for d in extra if da.sizes[d] == 1]
    da = da.squeeze(squeeze, drop=False) if squeeze else da
    extra = [d for d in da.dims if d not in allowed]
    if extra:
        raise ValueError(f"dimensions {extra} remain; reduce them first (e.g. .mean / .isel or "
                         f"DG_common.mean / integrate). Expected only {sorted(allowed)}")
    missing = [d for d in keep if d not in da.dims]
    if missing:
        raise ValueError(f"dimensions {missing} are missing in {da.dims}")
    return da


def _color_range(values: np.ndarray, vmin, vmax, symmetric, cmap):
    finite = values[np.isfinite(values)]
    lo = float(finite.min()) if finite.size else 0.0
    hi = float(finite.max()) if finite.size else 1.0
    if symmetric is None:
        symmetric = (vmin is None and vmax is None and lo < 0.0 < hi
                     and min(-lo, hi) > 0.05 * max(-lo, hi))
    if symmetric:
        m = max(abs(lo), abs(hi)) if (vmin is None and vmax is None) else max(abs(vmin or 0), abs(vmax or 0))
        vmin, vmax = -m, m
        cmap = cmap or "RdBu_r"
    else:
        vmin = lo if vmin is None else vmin
        vmax = hi if vmax is None else vmax
        cmap = cmap or "viridis"
    return vmin, vmax, cmap


def _unique_nodes(da: xr.DataArray, d: str) -> xr.DataArray:
    """Node-form array with duplicated element-boundary nodes averaged (monotonic coordinate)."""
    c = np.round(da[d].values, 9)
    cu, inv = np.unique(c, return_inverse=True)
    cnt = np.bincount(inv, minlength=cu.size)
    ax = da.get_axis_num(d)
    vals = np.moveaxis(np.asarray(da.values, dtype=float), ax, -1)
    acc = np.zeros(vals.shape[:-1] + (cu.size,))
    for j in range(vals.shape[-1]):
        acc[..., inv[j]] += vals[..., j]
    acc = np.moveaxis(acc / cnt, -1, ax)
    coords = {k: v for k, v in da.coords.items() if d not in v.dims}
    coords[d] = (d, cu, dict(da[d].attrs))
    return xr.DataArray(acc, dims=da.dims, coords=coords, attrs=da.attrs, name=da.name)


def element_edges(obj, d: str) -> np.ndarray:
    """Element boundary coordinates (Ne+1,) along ``d`` from ``<d>(k<d>, p<d>)`` and ``d<d>``."""
    kdim, pdim = dgc.ELEM_DIMS[d]
    if d in obj.coords and set(obj[d].dims) == {kdim, pdim}:
        c = obj[d].transpose(kdim, pdim).values
        xi = dgc.ref_nodes_of(obj, d)
        w = dgc.element_width(obj, d).values
        left = c[:, 0] - (xi[0] + 1.0) * 0.5 * w
        return np.r_[left, left[-1] + w[-1]]
    raise KeyError(f"coordinate '{d}({kdim},{pdim})' not found")


def add_element_edges(ax, obj, d: str, axis: str = "x", scale: Scale = 1e-3, **kw):
    """Draw element boundaries along direction ``d`` as lines on the given plot axis."""
    style = dict(color="k", lw=0.5, ls="--", alpha=0.6)
    style.update(kw)
    s = _scale_of(scale, d)
    for e in element_edges(obj, d) * s:
        (ax.axvline if axis == "x" else ax.axhline)(e, **style)
    return ax


# ---------------------------------------------------------------------------
# Tone plot
# ---------------------------------------------------------------------------

def plot_tone(da: xr.DataArray, dims: Tuple[str, str] = ("x", "z"), ax=None,
              method: str = "pcolormesh", cmap=None, vmin=None, vmax=None,
              symmetric: Optional[bool] = None, levels=None, scale: Scale = 1e-3,
              edges: bool = True, edge_kw: Optional[dict] = None,
              colorbar: bool = True, cbar_label: Optional[str] = None, cbar_kw: Optional[dict] = None,
              title: Optional[str] = None, xlim=None, ylim=None, **kw):
    """2-D tone plot of an element-form field over two directions.

    Parameters
    ----------
    da : DataArray whose only dims are the element dims of the two directions
        (size-1 dims are squeezed).  For each direction either ``(k<d>, p<d>)`` (drawn on
        the physical coordinate) or ``p<d>`` only (an element composite along that
        direction, drawn on the reference coordinate xi/eta/zeta).  Example: the
        x-composite of an x-z section, dims ``(kz, pz, px)``.
    dims : (horizontal-axis direction, vertical-axis direction), e.g. ("x", "z"), ("x", "y").
    method : ``"pcolormesh"`` (gouraud shading of the nodal values; element boundaries
        keep both one-sided values) or ``"contourf"`` (duplicated boundary nodes averaged).
    symmetric : colour range symmetric about 0 (default: automatic when the data change sign).
    levels : contour levels for ``contourf`` (int or sequence).
    edges : draw element boundaries (physical directions only).
    Other keyword arguments go to ``pcolormesh`` / ``contourf``.

    Returns
    -------
    ax, mappable
    """
    d1, d2 = dims
    phys = {d: dgc.ELEM_DIMS[d][0] in da.dims for d in dims}
    keep = []
    for d in dims:
        keep += list(dgc.ELEM_DIMS[d]) if phys[d] else [dgc.ELEM_DIMS[d][1]]
    da = _check_dims(da, keep).compute()
    ax = _get_ax(ax)
    nd = io.to_node_form(da)
    axdim = {d: d if phys[d] else dgc.ELEM_DIMS[d][1] for d in dims}
    if method == "contourf":
        for d in dims:
            if phys[d]:
                nd = _unique_nodes(nd, d)
    elif method != "pcolormesh":
        raise ValueError("method must be 'pcolormesh' or 'contourf'")
    Z = nd.transpose(axdim[d2], axdim[d1]).values

    def _axis_values(d):
        if phys[d]:
            return nd[d].values * _scale_of(scale, d)
        return dgc.ref_nodes_of(da, d)
    X, Y = _axis_values(d1), _axis_values(d2)
    vmin, vmax, cmap = _color_range(Z, vmin, vmax, symmetric, cmap)
    if method == "pcolormesh":
        m = ax.pcolormesh(X, Y, Z, shading="gouraud", cmap=cmap, vmin=vmin, vmax=vmax, **kw)
    else:
        if levels is None or np.isscalar(levels):
            levels = np.linspace(vmin, vmax, (levels or 20) + 1)
        m = ax.contourf(X, Y, Z, levels=levels, cmap=cmap, extend="both", **kw)
    if edges:
        for d, a in ((d1, "x"), (d2, "y")):
            if phys[d]:
                add_element_edges(ax, da, d, a, scale, **(edge_kw or {}))
    lab = {d: _coord_label(da, d, scale) if phys[d] else dgc.REF_COORD[d] for d in dims}
    ax.set_xlabel(lab[d1])
    ax.set_ylabel(lab[d2])
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if phys[d1] and phys[d2] and d1 in ("x", "y") and d2 in ("x", "y"):
        ax.set_aspect("equal")
    if colorbar:
        ax.figure.colorbar(m, ax=ax, label=auto_label(da, cbar_label), **(cbar_kw or {}))
    if title is not None:
        ax.set_title(title)
    return ax, m


# ---------------------------------------------------------------------------
# Line plot
# ---------------------------------------------------------------------------

def plot_line(da: xr.DataArray, dim: str = "x", ax=None, hue: Optional[str] = None,
              refine: Optional[int] = None, nodes: bool = True, vertical: Optional[bool] = None,
              scale: Scale = 1e-3, edges: bool = True, edge_kw: Optional[dict] = None,
              label: Optional[str] = None, value_label: Optional[str] = None,
              title: Optional[str] = None, node_kw: Optional[dict] = None, **plot_kw):
    """1-D line plot of an element-form field along one direction.

    Parameters
    ----------
    da : DataArray with dims ``(k<dim>, p<dim>)`` (+ ``hue``).
    dim : direction.
    hue : optional extra dimension; one line per value.
    refine : if given, draw each element polynomial on ``refine`` equispaced points per
        element (and the nodal values as markers when ``nodes``).  If None, the nodal
        values are joined (element-boundary jumps appear as vertical segments).
    vertical : put the coordinate on the vertical axis (default: True for ``"z"``).
    label : legend label (for ``hue``, the hue value is appended).

    Returns
    -------
    ax, list of Line2D
    """
    keep = list(dgc.ELEM_DIMS[dim])
    da = _check_dims(da, keep, hue).compute()
    ax = _get_ax(ax)
    vertical = (dim == "z") if vertical is None else vertical
    s = _scale_of(scale, dim)
    members = [(None, da)] if hue is None else [(v, da.sel({hue: v})) for v in da[hue].values]
    lines = []
    for i, (hv, d) in enumerate(members):
        lab = label if hv is None else (f"{label} {hv}" if label else f"{hue}={hv}")
        kw = dict(plot_kw)
        if hue is not None:
            kw.setdefault("color", f"C{i % 10}")
        if refine:
            df = dgc.interp_ref(d, np.linspace(-1.0, 1.0, int(refine)), dim)
            nd_f = io.to_node_form(df)
            c, v = nd_f[dim].values * s, nd_f.values
            ln, = (ax.plot(v, c, label=lab, **kw) if vertical else ax.plot(c, v, label=lab, **kw))
            lines.append(ln)
            if nodes:
                nd = io.to_node_form(d)
                mk = dict(ls="none", marker="o", ms=3, color=ln.get_color())
                mk.update(node_kw or {})
                c0, v0 = nd[dim].values * s, nd.values
                lines.append((ax.plot(v0, c0, **mk) if vertical else ax.plot(c0, v0, **mk))[0])
        else:
            nd = io.to_node_form(d)
            kw.setdefault("marker", "." if nodes else None)
            c, v = nd[dim].values * s, nd.values
            ln, = (ax.plot(v, c, label=lab, **kw) if vertical else ax.plot(c, v, label=lab, **kw))
            lines.append(ln)
    if edges:
        add_element_edges(ax, da, dim, "y" if vertical else "x", scale, **(edge_kw or {}))
    clab, vlab = _coord_label(da, dim, scale), auto_label(da, value_label)
    ax.set_xlabel(vlab if vertical else clab)
    ax.set_ylabel(clab if vertical else vlab)
    if label is not None or hue is not None:
        ax.legend()
    if title is not None:
        ax.set_title(title)
    return ax, lines


# ---------------------------------------------------------------------------
# Element composites
# ---------------------------------------------------------------------------

def plot_composite2d(comp: xr.DataArray, dims: Tuple[str, str] = ("x", "y"), ax=None,
                     refine: Optional[int] = None, nodes: bool = True,
                     cmap=None, vmin=None, vmax=None, symmetric: Optional[bool] = None,
                     levels=None, colorbar: bool = True, cbar_label: Optional[str] = None,
                     title: Optional[str] = None, xlim=(-1, 1), ylim=(-1, 1),
                     node_kw: Optional[dict] = None, **kw):
    """Element composite on the reference coordinates of two directions.

    Parameters
    ----------
    comp : DataArray with dims ``(p<d1>, p<d2>)``, e.g. from ``DG_common.composite_xy``.
    refine : None -> gouraud shading of the nodal values; int -> the element polynomial is
        interpolated to ``refine`` x ``refine`` points and drawn with ``contourf``.
    nodes : mark the node positions.
    vmin, vmax : colour range (default: data range, symmetric about 0 if the sign changes).
    xlim, ylim : axis ranges in reference coordinates (default (-1, 1)).

    Returns
    -------
    ax, mappable
    """
    d1, d2 = dims
    p1, p2 = dgc.ELEM_DIMS[d1][1], dgc.ELEM_DIMS[d2][1]
    comp = _check_dims(comp, [p1, p2]).compute()
    ax = _get_ax(ax)
    r1, r2 = dgc.ref_nodes_of(comp, d1), dgc.ref_nodes_of(comp, d2)
    if refine:
        xf = np.linspace(-1.0, 1.0, int(refine))
        cf = dgc.interp_ref(comp, xf, (d1, d2))
        Z = cf.transpose(p2, p1).values
        vmin, vmax, cmap = _color_range(Z, vmin, vmax, symmetric, cmap)
        if levels is None or np.isscalar(levels):
            levels = np.linspace(vmin, vmax, (levels or 20) + 1)
        m = ax.contourf(xf, xf, Z, levels=levels, cmap=cmap, extend="both", **kw)
    else:
        Z = comp.transpose(p2, p1).values
        vmin, vmax, cmap = _color_range(Z, vmin, vmax, symmetric, cmap)
        m = ax.pcolormesh(r1, r2, Z, shading="gouraud", cmap=cmap, vmin=vmin, vmax=vmax, **kw)
    if nodes:
        mk = dict(ls="none", marker=".", ms=3, color="k")
        mk.update(node_kw or {})
        X, Y = np.meshgrid(r1, r2)
        ax.plot(X.ravel(), Y.ravel(), **mk)
    ax.set_aspect("equal")
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.set_xlabel(dgc.REF_COORD[d1])
    ax.set_ylabel(dgc.REF_COORD[d2])
    if colorbar:
        ax.figure.colorbar(m, ax=ax, label=auto_label(comp, cbar_label))
    if title is not None:
        ax.set_title(title)
    return ax, m


def plot_composite1d(comp: xr.DataArray, dim: str = "x", ax=None, hue: Optional[str] = None,
                     refine: Optional[int] = 81, nodes: bool = True, vertical: Optional[bool] = None,
                     label: Optional[str] = None, value_label: Optional[str] = None,
                     title: Optional[str] = None, xlim="auto", ylim="auto",
                     node_kw: Optional[dict] = None, **plot_kw):
    """1-D element composite as a function of the reference coordinate of ``dim``.

    Node dims of the other directions present in ``comp`` are averaged over the
    element first (LGL quadrature, ``DG_common.element_mean``), e.g. the
    eta-average of a (py, px) composite for ``dim="x"``.

    Parameters
    ----------
    comp : DataArray with ``p<dim>`` (+ other ``p<d>``) (+ ``hue``).
    hue : optional extra dim (e.g. ``"ky"`` for per-row composites, or a case dim).
    refine : points for drawing the element polynomial (None: join the nodal values).
    xlim, ylim : axis ranges as in matplotlib.  ``"auto"`` (default) sets the
        reference-coordinate axis to (-1, 1) (x for horizontal, y for ``vertical``) and
        leaves the value axis to matplotlib; None leaves an axis untouched.
        E.g. ``ylim=(0, 0.5)`` fixes the value range of a horizontal plot.

    Returns
    -------
    ax, list of Line2D
    """
    pdim = dgc.ELEM_DIMS[dim][1]
    others = [d for d, (_, p) in dgc.ELEM_DIMS.items() if d != dim and p in comp.dims]
    if others:
        comp = dgc.element_mean(comp, others)
    comp = _check_dims(comp, [pdim], hue).compute()
    ax = _get_ax(ax)
    vertical = (dim == "z") if vertical is None else vertical
    r = dgc.ref_nodes_of(comp, dim)
    members = [(None, comp)] if hue is None else [(v, comp.sel({hue: v})) for v in comp[hue].values]
    lines = []
    for i, (hv, c) in enumerate(members):
        lab = label if hv is None else (f"{label} {hv}" if label else f"{hue}={hv}")
        kw = dict(plot_kw)
        if hue is not None:
            kw.setdefault("color", f"C{i % 10}")
        if refine:
            xf = np.linspace(-1.0, 1.0, int(refine))
            v = dgc.interp_ref(c, xf, dim).values
            ln, = (ax.plot(v, xf, label=lab, **kw) if vertical else ax.plot(xf, v, label=lab, **kw))
            lines.append(ln)
            if nodes:
                mk = dict(ls="none", marker="o", ms=4, color=ln.get_color())
                mk.update(node_kw or {})
                lines.append((ax.plot(c.values, r, **mk) if vertical else ax.plot(r, c.values, **mk))[0])
        else:
            kw.setdefault("marker", "o" if nodes else None)
            ln, = (ax.plot(c.values, r, label=lab, **kw) if vertical else ax.plot(r, c.values, label=lab, **kw))
            lines.append(ln)
    rl, vl = dgc.REF_COORD[dim], auto_label(comp, value_label)
    ax.set_xlabel(vl if vertical else rl)
    ax.set_ylabel(rl if vertical else vl)
    ref_lim, val_lim = (ylim, xlim) if vertical else (xlim, ylim)
    ref_lim = (-1, 1) if (isinstance(ref_lim, str) and ref_lim == "auto") else ref_lim
    val_lim = None if (isinstance(val_lim, str) and val_lim == "auto") else val_lim
    for setter, lim in (((ax.set_ylim, ref_lim), (ax.set_xlim, val_lim)) if vertical
                        else ((ax.set_xlim, ref_lim), (ax.set_ylim, val_lim))):
        if lim is not None:
            setter(lim)
    if label is not None or hue is not None:
        ax.legend()
    if title is not None:
        ax.set_title(title)
    return ax, lines


def plot_composite_summary(comp: xr.DataArray, dims: Tuple[str, str] = ("x", "y"),
                           refine2d: Optional[int] = None, refine1d: Optional[int] = 81,
                           vmin=None, vmax=None, ylim=None,
                           figsize=(11, 3.6), suptitle: Optional[str] = None,
                           kw1d: Optional[dict] = None, **kw2d):
    """Figure with the 2-D composite and the 1-D composites along each direction.

    Parameters
    ----------
    comp : composite with dims ``(p<d1>, p<d2>)``.
    vmin, vmax : colour range of the 2-D panel.
    ylim : value range of the 1-D panels: a tuple used for both, a dict per direction
        (``{"x": (0, 0.5), "y": (0, 0.6)}``), ``"shared"`` (same automatic range for
        both), or None (automatic for each).
    kw1d : other keyword arguments for :func:`plot_composite1d`.
    Other keyword arguments go to :func:`plot_composite2d` (``cmap``, ``levels``, ...).

    Returns
    -------
    fig, axes (3,)
    """
    plt = _plt()
    fig, axs = plt.subplots(1, 3, figsize=figsize, constrained_layout=True,
                            gridspec_kw={"width_ratios": [1.25, 1, 1]})
    kw2d.setdefault("cbar_label", "")          # the label goes to the title to save space
    kw2d.setdefault("title", auto_label(comp))
    plot_composite2d(comp, dims, ax=axs[0], refine=refine2d, vmin=vmin, vmax=vmax, **kw2d)
    for ax, d in zip(axs[1:], dims):
        lim = ylim.get(d) if isinstance(ylim, dict) else (None if isinstance(ylim, str) else ylim)
        plot_composite1d(comp, d, ax=ax, refine=refine1d, ylim=lim if lim is not None else "auto",
                         **(kw1d or {}))
        others = "".join(dd for dd in dims if dd != d)
        ax.set_title(f"mean over {dgc.REF_COORD[others]}" if others else "")
    if isinstance(ylim, str) and ylim == "shared":
        lo = min(ax.get_ylim()[0] for ax in axs[1:])
        hi = max(ax.get_ylim()[1] for ax in axs[1:])
        for ax in axs[1:]:
            ax.set_ylim(lo, hi)
    if suptitle:
        fig.suptitle(suptitle)
    return fig, axs


# ---------------------------------------------------------------------------
# Self test (draws with the Agg backend, checks only that everything runs)
# ---------------------------------------------------------------------------

def _selftest(verbose: bool = True) -> None:
    import matplotlib
    matplotlib.use("Agg")
    plt = _plt()
    f = lambda X, Y: np.sin(2 * np.pi * X / 4.0) * np.cos(2 * np.pi * Y / 3.0)
    da = dgc._make_test_array(func=f).rename("q")
    da.attrs.update(long_name="test field", units="1")
    da["x"].attrs["units"] = "m"
    da["y"].attrs["units"] = "m"
    for m in ("pcolormesh", "contourf"):
        plot_tone(da, ("x", "y"), method=m, scale=1.0)
    plot_line(da.isel(ky=0, py=0), "x", refine=16, scale=1.0)
    plot_line(da.isel(py=0).rename({"ky": "row"}).isel(row=[0, 1]), "x", hue="row", scale=1.0)
    plot_line(da.isel(kx=0, px=0).rename({"ky": "kz", "py": "pz", "y": "z", "eta": "zeta", "dy": "dz"}),
              "z", scale=1.0)
    ax, _ = plot_tone(dgc.composite_x(da), ("x", "y"), scale=1.0)     # (ky, py, px): xi vs y
    assert ax.get_xlabel() == "xi"
    comp = dgc.composite_xy(da)
    plot_composite2d(comp)
    plot_composite2d(comp, refine=21)
    plot_composite1d(comp, "x")
    plot_composite1d(dgc.composite_x(da), "x", hue="ky")
    plot_composite_summary(comp)
    fig, axs = plot_composite_summary(comp, vmin=-0.5, vmax=0.5, ylim=(-0.2, 0.3))
    assert axs[1].get_ylim() == (-0.2, 0.3) and axs[2].get_ylim() == (-0.2, 0.3)
    fig, axs = plot_composite_summary(comp, ylim={"x": (0, 1)})
    assert axs[1].get_ylim() == (0, 1)
    fig, axs = plot_composite_summary(comp, ylim="shared")
    assert axs[1].get_ylim() == axs[2].get_ylim()
    ax, _ = plot_composite1d(comp, "x", ylim=(0, 2))
    assert ax.get_ylim() == (0, 2) and ax.get_xlim() == (-1, 1)
    ax, _ = plot_composite1d(comp, "x", vertical=True, xlim=(0, 2))
    assert ax.get_xlim() == (0, 2) and ax.get_ylim() == (-1, 1)
    ax, m = plot_composite2d(comp, vmin=0, vmax=1, xlim=(-1, 0))
    assert m.get_clim() == (0, 1) and ax.get_xlim() == (-1, 0)
    # error for remaining dims
    try:
        plot_tone(da.expand_dims(time=2), ("x", "y"))
    except ValueError:
        pass
    else:
        raise AssertionError("extra dims must raise")
    u = _unique_nodes(io.to_node_form(da), "x")
    assert u.sizes["x"] == 4 * 7 + 1 and np.all(np.diff(u["x"].values) > 0)
    plt.close("all")
    if verbose:
        print("mkgraph_common: selftest passed")


if __name__ == "__main__":
    _selftest()
