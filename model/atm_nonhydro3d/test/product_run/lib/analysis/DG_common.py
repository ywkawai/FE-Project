"""Common DG utilities for analysing SCALE-DG output.

Reference-element operations (LGL/GL nodes, Lagrange/Legendre matrices) and
element-form operations on xarray objects (integrals, means, composites,
modal coefficients).

Element form
------------
Each horizontal/vertical direction ``d`` in ``("x", "y", "z")`` is represented
by a pair of dimensions ``(k<d>, p<d>)``: element index and node index in the
element.  Coordinates (set by :mod:`io_common`):

* ``x(kx, px)``  : physical node coordinate
* ``xi(px)``     : reference coordinate in [-1, 1] (``eta`` for y, ``zeta`` for z)
* ``dx(kx)``     : element width

Only numpy and xarray are required.  All matrix operations go through
``xr.dot`` so that dask-backed arrays stay lazy.
"""
from __future__ import annotations

from functools import lru_cache
from typing import Callable, Iterable, Optional, Sequence, Tuple, Union

import numpy as np
import xarray as xr
from numpy.polynomial import legendre as npleg

__all__ = [
    "ELEM_DIMS", "REF_COORD", "MODE_DIMS",
    "lgl_nodes", "gl_nodes", "ref_nodes", "legendre_vandermonde",
    "lagrange_matrix", "quad_weights", "mass_matrix",
    "ref_nodes_of", "element_width", "apply_along", "map_dataset",
    "interp_ref", "integrate", "mean", "element_mean",
    "composite", "composite_x", "composite_y", "composite_xy",
    "modal_coeffs", "composite_legendre", "lock_ratio",
]

#: element/node dimension names for each direction
ELEM_DIMS = {"x": ("kx", "px"), "y": ("ky", "py"), "z": ("kz", "pz")}
#: reference-coordinate names for each direction
REF_COORD = {"x": "xi", "y": "eta", "z": "zeta"}
#: modal (Legendre) dimension names for each direction
MODE_DIMS = {"x": "mx", "y": "my", "z": "mz"}

XrObj = Union[xr.DataArray, xr.Dataset]


# ---------------------------------------------------------------------------
# Reference element
# ---------------------------------------------------------------------------

def _readonly(a) -> np.ndarray:
    a = np.array(a, dtype=float)
    a.setflags(write=False)
    return a


@lru_cache(maxsize=None)
def lgl_nodes(Np: int) -> Tuple[np.ndarray, np.ndarray]:
    """Legendre-Gauss-Lobatto nodes and weights on [-1, 1] (ascending).

    Parameters
    ----------
    Np : int
        Number of nodes (= polynomial order + 1), Np >= 2.

    Returns
    -------
    xi, w : ndarray
        Nodes and quadrature weights (read-only).
    """
    if Np < 2:
        raise ValueError("LGL nodes require Np >= 2")
    N = Np - 1
    cN = np.zeros(N + 1)
    cN[N] = 1.0
    if N == 1:
        interior = np.empty(0)
    else:
        dc = npleg.legder(cN)
        d2c = npleg.legder(dc)
        x = np.sort(npleg.legroots(dc).real)
        for _ in range(3):  # Newton refinement
            x = x - npleg.legval(x, dc) / npleg.legval(x, d2c)
        interior = x
    xi = np.concatenate(([-1.0], interior, [1.0]))
    xi = 0.5 * (xi - xi[::-1])  # enforce exact symmetry
    w = 2.0 / (N * (N + 1) * npleg.legval(xi, cN) ** 2)
    return _readonly(xi), _readonly(w)


@lru_cache(maxsize=None)
def gl_nodes(Np: int) -> Tuple[np.ndarray, np.ndarray]:
    """Gauss-Legendre nodes and weights on [-1, 1] (ascending)."""
    if Np < 1:
        raise ValueError("GL nodes require Np >= 1")
    xi, w = npleg.leggauss(Np)
    return _readonly(xi), _readonly(w)


def ref_nodes(Np: int, node_type: str = "LGL") -> np.ndarray:
    """Reference nodes of the given type (``"LGL"`` or ``"GL"``)."""
    t = node_type.upper()
    if t == "LGL":
        return lgl_nodes(Np)[0]
    if t == "GL":
        return gl_nodes(Np)[0]
    raise ValueError(f"unknown node_type: {node_type}")


def legendre_vandermonde(xi: Sequence[float], nmode: int,
                         normalized: bool = True) -> np.ndarray:
    """Vandermonde matrix V[i, m] = P_m(xi_i) (orthonormal on [-1, 1] by default)."""
    xi = np.asarray(xi, dtype=float)
    V = npleg.legvander(xi, nmode - 1)
    if normalized:
        V = V * np.sqrt((2.0 * np.arange(nmode) + 1.0) / 2.0)
    return V


def lagrange_matrix(xi_src: Sequence[float], xi_dst: Sequence[float]) -> np.ndarray:
    """Interpolation matrix L (n_dst x n_src) of the Lagrange basis on ``xi_src``."""
    xi_src = np.asarray(xi_src, dtype=float)
    xi_dst = np.asarray(xi_dst, dtype=float)
    n = xi_src.size
    Vs = legendre_vandermonde(xi_src, n)
    Vd = legendre_vandermonde(xi_dst, n)
    return np.linalg.solve(Vs.T, Vd.T).T


def quad_weights(xi: Sequence[float]) -> np.ndarray:
    """Interpolatory quadrature weights on the nodes ``xi`` (= LGL/GL weights for those nodes).

    For LGL/GL nodes the tabulated weights are returned.  For other nodes (e.g. many
    equispaced points made by :func:`interp_ref` for plotting) the interpolatory rule
    can be unstable; a warning is issued when negative weights appear.  Integrate on
    the original DG nodes instead in that case.
    """
    xi = np.asarray(xi, dtype=float)
    n = xi.size
    for gen in (lgl_nodes, gl_nodes):
        if n >= (2 if gen is lgl_nodes else 1):
            x_ref, w_ref = gen(n)
            if np.allclose(xi, x_ref, atol=1e-10):
                return np.array(w_ref)
    invV = np.linalg.inv(legendre_vandermonde(xi, n))
    w = np.sqrt(2.0) * invV[0, :]
    if np.any(w < 0.0):
        import warnings
        warnings.warn(f"interpolatory quadrature on these {n} nodes has negative weights "
                      "(unstable, e.g. many equispaced points); integrate on the DG nodes instead",
                      RuntimeWarning, stacklevel=2)
    return w


def mass_matrix(xi: Sequence[float], lumped: bool = False) -> np.ndarray:
    """Mass matrix of the Lagrange basis on ``xi`` (exact, or lumped/diagonal)."""
    xi = np.asarray(xi, dtype=float)
    if lumped:
        return np.diag(quad_weights(xi))
    V = legendre_vandermonde(xi, xi.size)
    return np.linalg.inv(V @ V.T)


# ---------------------------------------------------------------------------
# Helpers for element-form xarray objects
# ---------------------------------------------------------------------------

def _parse_dirs(dims: Union[str, Iterable[str]]) -> Tuple[str, ...]:
    dirs = tuple(dims) if not isinstance(dims, str) else tuple(dims)
    for d in dirs:
        if d not in ELEM_DIMS:
            raise ValueError(f"unknown direction '{d}' (use 'x', 'y', 'z')")
    return dirs


def ref_nodes_of(obj: XrObj, d: str) -> np.ndarray:
    """Reference node coordinates along direction ``d`` (LGL is assumed if no coordinate)."""
    _, pdim = ELEM_DIMS[d]
    name = REF_COORD[d]
    if name in obj.coords and obj[name].dims == (pdim,):
        return np.asarray(obj[name].values, dtype=float)
    return lgl_nodes(obj.sizes[pdim])[0]


def element_width(obj: XrObj, d: str) -> xr.DataArray:
    """Element widths along ``d`` as a DataArray on ``k<d>``.

    Uses the coordinate ``d<d>`` if present, otherwise the physical node
    coordinate ``<d>(k<d>, p<d>)``.
    """
    kdim, pdim = ELEM_DIMS[d]
    wname = "d" + d
    if wname in obj.coords and obj[wname].dims == (kdim,):
        w = obj[wname]
        return xr.DataArray(w.values, dims=(kdim,),
                            coords={kdim: obj[kdim].values} if kdim in obj.coords else None)
    if d in obj.coords and set(obj[d].dims) == {kdim, pdim}:
        xi = ref_nodes_of(obj, d)
        c = obj[d].transpose(kdim, pdim).values
        w = (c[:, -1] - c[:, 0]) * 2.0 / (xi[-1] - xi[0])
        return xr.DataArray(w, dims=(kdim,),
                            coords={kdim: obj[kdim].values} if kdim in obj.coords else None)
    raise KeyError(f"cannot determine element width along '{d}': "
                   f"need coordinate '{wname}({kdim})' or '{d}({kdim},{pdim})'")


def apply_along(da: xr.DataArray, M: np.ndarray, d: str,
                out_pdim: Optional[str] = None,
                out_ref: Optional[Tuple[str, np.ndarray]] = None,
                transform_coord: bool = True) -> xr.DataArray:
    """Apply an element-local matrix ``M`` (n_out x n_in) along ``p<d>``.

    Parameters
    ----------
    da : DataArray in element form.
    M : ndarray (n_out, n_in) with n_in == da.sizes[p<d>].
    d : direction.
    out_pdim : name of the output node dimension (default: same as input).
    out_ref : (name, values) of the reference coordinate on the output dimension.
    transform_coord : apply ``M`` to the physical coordinate ``<d>(k<d>, p<d>)`` too
        (valid for interpolation and averaging operators, which map a linear
        function to its point values / cell means).
    """
    kdim, pdim = ELEM_DIMS[d]
    out_pdim = out_pdim or pdim
    M = np.asarray(M, dtype=float)
    if M.shape[1] != da.sizes[pdim]:
        raise ValueError(f"matrix shape {M.shape} does not match size of '{pdim}' ({da.sizes[pdim]})")

    p_coords = [c for c in da.coords if pdim in da[c].dims]
    phys = None
    if transform_coord and d in da.coords and set(da[d].dims) == {kdim, pdim}:
        phys = da[d].transpose(kdim, pdim).values
    base = da.drop_vars(p_coords)

    tmp = "__out_" + pdim
    Mda = xr.DataArray(M, dims=(tmp, pdim))
    out = xr.dot(base, Mda, dims=pdim).rename({tmp: out_pdim})
    out = out.transpose(*[out_pdim if dd == pdim else dd for dd in da.dims])
    if phys is not None:
        out = out.assign_coords({d: ((kdim, out_pdim), phys @ M.T)})
    if out_ref is not None:
        out = out.assign_coords({out_ref[0]: (out_pdim, np.asarray(out_ref[1], dtype=float))})
    out.attrs = dict(da.attrs)
    out.name = da.name
    return out


def map_dataset(func: Callable[..., xr.DataArray], obj: XrObj, dirs: Sequence[str],
                *args, **kwargs) -> XrObj:
    """Apply ``func(DataArray, ...)`` to a DataArray, or to every data variable of a Dataset
    that has all node dimensions of ``dirs``.

    Variables that have none of the element/node dimensions of ``dirs`` are kept
    unchanged; variables that have only some of them are dropped.
    """
    if isinstance(obj, xr.DataArray):
        return func(obj, *args, **kwargs)
    need = {ELEM_DIMS[d][1] for d in dirs}
    related = need | {ELEM_DIMS[d][0] for d in dirs}
    out = {}
    for name, var in obj.data_vars.items():
        dims = set(var.dims)
        if need <= dims:
            out[name] = func(var, *args, **kwargs)
        elif not (dims & related):
            out[name] = var
    ds = xr.Dataset(out)
    ds.attrs = dict(obj.attrs)
    return ds


# ---------------------------------------------------------------------------
# Interpolation, integrals and means
# ---------------------------------------------------------------------------

def interp_ref(obj: XrObj, xi_new: Sequence[float], dims: Union[str, Iterable[str]] = "x") -> XrObj:
    """Interpolate the element polynomials to new reference points ``xi_new`` (in [-1, 1])."""
    dirs = _parse_dirs(dims)
    xi_new = np.asarray(xi_new, dtype=float)

    def _f(da):
        for d in dirs:
            L = lagrange_matrix(ref_nodes_of(da, d), xi_new)
            if "d" + d not in da.coords:
                try:
                    da = da.assign_coords({"d" + d: element_width(da, d)})
                except KeyError:
                    pass
            da = apply_along(da, L, d, out_ref=(REF_COORD[d], xi_new))
        return da
    return map_dataset(_f, obj, dirs)


def integrate(obj: XrObj, dims: Union[str, Iterable[str]] = "x",
              weight: Optional[xr.DataArray] = None, exact: bool = False) -> XrObj:
    """Integral over the given directions (sum over elements of the element integrals).

    Parameters
    ----------
    obj : DataArray or Dataset in element form.
    dims : directions, e.g. ``"x"``, ``"xy"``, ``"xyz"``.
    weight : optional weight field (e.g. density) multiplied before integrating.
    exact : if True, interpolate to Np-point GL nodes before multiplying by ``weight``
        and integrating (exact for products of two degree-p polynomials).
        Without ``weight`` the nodal quadrature is already exact for a degree-p
        polynomial, so ``exact`` is only meaningful together with ``weight``.
    """
    dirs = _parse_dirs(dims)

    def _f(da):
        attrs = _integral_attrs(da, dirs)
        w_da = weight
        widths = {d: element_width(da, d) for d in dirs}
        if exact:
            for d in dirs:
                xg = gl_nodes(da.sizes[ELEM_DIMS[d][1]])[0]
                da = interp_ref(da, xg, d)
                if w_da is not None:
                    w_da = interp_ref(w_da, xg, d)
        if w_da is not None:
            da = da * w_da
        for d in dirs:
            kdim, pdim = ELEM_DIMS[d]
            w = xr.DataArray(quad_weights(ref_nodes_of(da, d)), dims=(pdim,))
            J = widths[d] * 0.5
            da = (da * w * J).sum((kdim, pdim), skipna=False)
        da.attrs = attrs
        return da
    return map_dataset(_f, obj, dirs)


def _integral_attrs(da: xr.DataArray, dirs: Sequence[str]) -> dict:
    """Attributes of the integral: long_name kept, units multiplied by the coordinate units."""
    attrs = dict(da.attrs)
    cu = [da[d].attrs.get("units", "") if d in da.coords else "" for d in dirs]
    cu = [u for u in cu if u]
    if cu:
        if all(u == cu[0] for u in cu) and len(cu) > 1:
            extra = f"{cu[0]}{len(cu)}"
        else:
            extra = " ".join(cu)
        attrs["units"] = f"{attrs['units']} {extra}".strip() if attrs.get("units") else extra
    return attrs


def _keep_attrs(out: XrObj, src: XrObj) -> XrObj:
    if isinstance(out, xr.DataArray):
        out.attrs = dict(src.attrs)
    else:
        for v in out.data_vars:
            if v in src.data_vars:
                out[v].attrs = dict(src[v].attrs)
        out.attrs = dict(src.attrs)
    return out


def mean(obj: XrObj, dims: Union[str, Iterable[str]] = "x",
         weight: Optional[xr.DataArray] = None, exact: bool = False) -> XrObj:
    """Mean over the given directions (integral divided by length/area/volume, or by the
    integral of ``weight`` when given, i.e. a weighted mean)."""
    dirs = _parse_dirs(dims)
    num = integrate(obj, dirs, weight=weight, exact=exact)
    if weight is None:
        ref = obj if isinstance(obj, xr.DataArray) else next(
            v for v in obj.data_vars.values() if {ELEM_DIMS[d][1] for d in dirs} <= set(v.dims))
        den = 1.0
        for d in dirs:
            den = den * float(element_width(ref, d).sum())
    else:
        den = integrate(weight, dirs, exact=exact)
    return _keep_attrs(num / den, obj)


def element_mean(obj: XrObj, dims: Union[str, Iterable[str]] = "x") -> XrObj:
    """Mean of the element polynomial over each element (``p<d>`` is reduced, ``k<d>`` kept)."""
    dirs = _parse_dirs(dims)

    def _f(da):
        for d in dirs:
            _, pdim = ELEM_DIMS[d]
            w = xr.DataArray(quad_weights(ref_nodes_of(da, d)) * 0.5, dims=(pdim,))
            attrs = dict(da.attrs)
            da = (da * w).sum(pdim, skipna=False)
            da.attrs = attrs
        return da
    return map_dataset(_f, obj, dirs)


# ---------------------------------------------------------------------------
# Element composites and Legendre analysis
# ---------------------------------------------------------------------------

def composite(obj: XrObj, dims: Union[str, Iterable[str]] = "xy",
              weights: Optional[xr.DataArray] = None, skipna: bool = False) -> XrObj:
    """Element composite: average over the element indices ``k<d>`` of the given directions.

    The result keeps the node dimensions ``p<d>`` (and the reference coordinates),
    i.e. it is the mean field as a function of the position inside the element.
    Take time / vertical means beforehand if needed.
    """
    dirs = _parse_dirs(dims)
    kdims = [ELEM_DIMS[d][0] for d in dirs]

    def _f(da):
        if weights is None:
            out = da.mean(kdims, skipna=skipna)
        else:
            out = da.weighted(weights).mean(kdims, skipna=skipna)
        out.attrs = dict(da.attrs)
        out.attrs["dg_op"] = "composite_" + "".join(dirs)
        return out
    return map_dataset(_f, obj, dirs)


def composite_x(obj: XrObj, **kw) -> XrObj:
    """Composite over ``kx`` only."""
    return composite(obj, "x", **kw)


def composite_y(obj: XrObj, **kw) -> XrObj:
    """Composite over ``ky`` only."""
    return composite(obj, "y", **kw)


def composite_xy(obj: XrObj, **kw) -> XrObj:
    """Horizontal composite over ``(kx, ky)``."""
    return composite(obj, "xy", **kw)


def modal_coeffs(obj: XrObj, dims: Union[str, Iterable[str]] = "xy") -> XrObj:
    """Orthonormal Legendre coefficients of the element polynomials.

    The node dimension ``p<d>`` is replaced by the mode dimension ``m<d>``
    (mode numbers 0..p).  The coefficients are obtained exactly with the
    inverse Vandermonde matrix (not by nodal quadrature).
    """
    dirs = _parse_dirs(dims)

    def _f(da):
        for d in dirs:
            _, pdim = ELEM_DIMS[d]
            n = da.sizes[pdim]
            invV = np.linalg.inv(legendre_vandermonde(ref_nodes_of(da, d), n))
            da = apply_along(da, invV, d, out_pdim=MODE_DIMS[d], transform_coord=False)
            da = da.assign_coords({MODE_DIMS[d]: np.arange(n)})
        return da
    return map_dataset(_f, obj, dirs)


def composite_legendre(comp: XrObj, dims: Union[str, Iterable[str]] = "xy",
                       normalize: bool = False) -> XrObj:
    """Legendre coefficients a_mn of an element composite.

    With ``normalize=True`` the coefficients are divided by a_00 (relative amplitude).
    """
    dirs = _parse_dirs(dims)
    a = modal_coeffs(comp, dirs)
    if normalize:
        a00 = a.isel({MODE_DIMS[d]: 0 for d in dirs})
        a = a / a00
    return a


def lock_ratio(a: XrObj, sample_dims: Sequence[str] = ("kx", "ky", "time"), ddof: int = 1) -> XrObj:
    """Element-locking ratio R = <a> / sigma(a) over ``sample_dims``.

    ``a`` are modal coefficients (see :func:`modal_coeffs`) per element and time.
    This is the 2-D extension of R_m of the toy model (simple_models/NOTES.md 4.5);
    R_e corresponds to the (m=2, n=0) and (m=0, n=2) components.  Subtracting the
    domain mean only changes the (0, 0) component.
    """
    dims = [d for d in sample_dims if d in a.dims]
    return a.mean(dims) / a.std(dims, ddof=ddof)


# ---------------------------------------------------------------------------
# Self test
# ---------------------------------------------------------------------------

def _make_test_array(NeX=4, NeY=3, Np=8, Lx=4.0, Ly=3.0, func=None) -> xr.DataArray:
    xi = lgl_nodes(Np)[0]
    dx, dy = Lx / NeX, Ly / NeY
    x = (np.arange(NeX)[:, None] + 0.5 * (xi[None, :] + 1.0)) * dx
    y = (np.arange(NeY)[:, None] + 0.5 * (xi[None, :] + 1.0)) * dy
    X = x[None, None, :, :]
    Y = y[:, :, None, None]
    f = func if func is not None else (lambda X, Y: np.ones(np.broadcast(X, Y).shape))
    data = f(X, Y) * np.ones((NeY, Np, NeX, Np))
    return xr.DataArray(
        data, dims=("ky", "py", "kx", "px"),
        coords={"kx": np.arange(NeX), "px": np.arange(Np), "ky": np.arange(NeY), "py": np.arange(Np),
                "xi": ("px", xi), "eta": ("py", xi),
                "x": (("kx", "px"), x), "y": (("ky", "py"), y),
                "dx": ("kx", np.full(NeX, dx)), "dy": ("ky", np.full(NeY, dy))})


def _selftest(verbose: bool = True) -> None:
    # LGL nodes/weights
    for Np in (2, 3, 5, 8, 12):
        xi, w = lgl_nodes(Np)
        assert abs(w.sum() - 2.0) < 1e-13
        assert np.allclose(quad_weights(xi), w, atol=1e-12)
        k = 2 * Np - 4  # even degree <= 2Np-3 is integrated exactly by LGL
        assert abs(np.dot(w, xi ** k) - 2.0 / (k + 1)) < 1e-12
    xi8 = lgl_nodes(8)[0]
    assert np.allclose(xi8[[0, -1]], [-1, 1]) and np.all(np.diff(xi8) > 0)
    # Lagrange interpolation exact for degree <= p
    xg = gl_nodes(8)[0]
    L = lagrange_matrix(xi8, xg)
    assert np.allclose(L @ xi8 ** 7, xg ** 7, atol=1e-12)
    # mass matrix
    V = legendre_vandermonde(xi8, 8)
    assert np.allclose(V.T @ mass_matrix(xi8) @ V, np.eye(8), atol=1e-10)

    # integrals on element-form arrays
    f = lambda X, Y: X ** 3 + 2.0 * Y
    da = _make_test_array(func=f)
    I = float(integrate(da, "xy"))
    I_exact = (4.0 ** 4 / 4) * 3.0 + 4.0 * 3.0 ** 2
    assert abs(I - I_exact) < 1e-10 * abs(I_exact), (I, I_exact)
    assert abs(float(mean(da, "xy")) - I_exact / 12.0) < 1e-10
    # weighted exact integral of a product of two degree-7 polynomials
    g = _make_test_array(func=lambda X, Y: X ** 7)
    Iw = float(integrate(g, "x", weight=g, exact=True).isel(ky=0, py=0))
    assert abs(Iw - 4.0 ** 15 / 15) < 1e-9 * 4.0 ** 15 / 15
    # element mean
    em = element_mean(_make_test_array(func=lambda X, Y: X), "x").isel(ky=0, py=0)
    assert np.allclose(em.values, np.arange(4) + 0.5)
    # interpolation keeps physical coordinates
    di = interp_ref(da, xg, "x")
    assert np.allclose(di["x"].values, (np.arange(4)[:, None] + 0.5 * (xg + 1)))
    assert np.allclose(di.values, f(di["x"].values[None, None], da["y"].values[:, :, None, None]))
    # composite and modal coefficients: field = P2(xi) in every element
    base = _make_test_array()
    p2 = xr.zeros_like(base) + 0.5 * (3 * base["xi"] ** 2 - 1)
    a = composite_legendre(composite_xy(p2))
    an = a.transpose("my", "mx").values
    expect = np.zeros((8, 8))
    expect[0, 2] = np.sqrt(2.0 / 5.0) * np.sqrt(2.0)  # P2 = sqrt(2/5)*P~2, times P~0 = 1/sqrt(2) in y
    assert np.allclose(an, expect, atol=1e-12), an[:3, :3]
    # attributes are kept by mean / element_mean; integrate multiplies the units
    da.attrs.update(long_name="q", units="kg")
    da["x"].attrs["units"] = "m"; da["y"].attrs["units"] = "m"
    assert mean(da, "xy").attrs["units"] == "kg" and element_mean(da, "x").attrs["long_name"] == "q"
    assert integrate(da, "xy").attrs["units"] == "kg m2" and integrate(da, "x").attrs["units"] == "kg m"
    # Dataset path
    ds = xr.Dataset({"a": da, "b": da * 2, "t": xr.DataArray([0.0, 1.0], dims="time")})
    m = mean(ds, "xy")
    assert abs(float(m["b"]) - 2 * I_exact / 12.0) < 1e-10 and "t" in m
    if verbose:
        print("DG_common: selftest passed")


if __name__ == "__main__":
    _selftest()
