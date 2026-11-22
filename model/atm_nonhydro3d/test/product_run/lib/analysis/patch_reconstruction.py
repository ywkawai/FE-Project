"""Patch reconstruction (FElib ``Reconstruction2``) in Python.

Port of ``MeshFieldFilterOperationBase_prepair_reconstruct2_matrix`` and
``..._apply_reconst1d_{x,y}_2`` in
``FElib/src/*/scale_meshfield_filter_operation_base.F90`` (centred patch
reconstruction, paper Sec. 4.1.1, Eq. 7).

Operator
--------
For each node x0 (reference coordinate) of an element, take the window of
width 2*Delta centred at x0 (Delta = element width).  It covers

* the left neighbour on [x0, 1]   (its own reference coordinate),
* the element itself on [-1, 1],
* the right neighbour on [-1, x0].

The DG polynomials on the window are L2-projected onto polynomials of degree
Nr-1 (exact mass matrix, as ``LineElement%Init(.., LumpedMassMatFlag=.false.)``)
and the projection is evaluated at the window centre.  The result is linear:

    q_new[i] = sum_j ( Ml[j,i] q_left[j] + Mc[j,i] q[j] + Mr[j,i] q_right[j] )

In 3-D SCALE-DG applies it horizontally only, x first and then y (with a halo
exchange in between).  Since the two 1-D operators act on different axes they
commute; the vertical is not filtered.
"""
from __future__ import annotations

from functools import lru_cache
from math import ceil
from typing import Iterable, Tuple, Union

import numpy as np
import xarray as xr
from numpy.polynomial import legendre as npleg

from . import DG_common as dgc

__all__ = ["reconst2_matrices", "apply_reconst2", "transfer_function_continuous",
           "transfer_function_discrete"]


@lru_cache(maxsize=None)
def reconst2_matrices(Np: int, Nr: int, node_type: str = "LGL") -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Matrices (Ml, Mc, Mr), each (Np x Np), of the Reconstruction2 operator.

    ``M[j, i]`` is the weight of node j of the left / centre / right element in the
    new value at node i of the centre element (same layout as ``Ml_tr`` etc. in FElib).

    Parameters
    ----------
    Np : nodes per element (= PolyOrder + 1).
    Nr : nodes of the reconstruction polynomial (``Nnode_h1D_reconst``; degree Nr-1).
    node_type : node type of the DG element ("LGL" in SCALE-DG).
    """
    xi = dgc.ref_nodes(Np, node_type)
    p, pr = Np - 1, Nr - 1
    nint = int(ceil(0.5 * (p + pr))) + 1           # same as FElib
    r, w = np.polynomial.legendre.leggauss(nint)
    # orthonormal Legendre of the reconstruction space evaluated at the centre
    P0 = dgc.legendre_vandermonde([0.0], Nr)[0]     # (Nr,)

    Ms = [np.zeros((Np, Np)) for _ in range(3)]
    for i, x0 in enumerate(xi):
        # (element coordinate range, reconstruction coordinate range) of each piece
        s_l1 = -1.0 + 0.5 * (1.0 - x0)
        pieces = (((x0, 1.0), (-1.0, s_l1)),          # left neighbour
                  ((-1.0, 1.0), (s_l1, s_l1 + 1.0)),  # centre element
                  ((-1.0, x0), (s_l1 + 1.0, 1.0)))    # right neighbour
        for M, ((xl, xr), (sl, sr)) in zip(Ms, pieces):
            if sr - sl <= 0.0:
                continue
            x_int = xl + 0.5 * (xr - xl) * (1.0 + r)
            s_int = sl + 0.5 * (sr - sl) * (1.0 + r)
            L = dgc.lagrange_matrix(xi, x_int)                  # (nint, Np)
            Ps = dgc.legendre_vandermonde(s_int, Nr)           # (nint, Nr)
            coef = 0.5 * (sr - sl)
            M[:, i] = L.T @ (coef * w * (Ps @ P0))
    for M in Ms:
        M.setflags(write=False)
    return Ms[0], Ms[1], Ms[2]


def _apply_1d(da: xr.DataArray, d: str, Nr: int, boundary: str, node_type: str) -> xr.DataArray:
    kdim, pdim = dgc.ELEM_DIMS[d]
    Np = da.sizes[pdim]
    Ml, Mc, Mr = reconst2_matrices(Np, Nr, node_type)

    p_coords = [c for c in da.coords if pdim in da[c].dims or kdim in da[c].dims]
    base = da.drop_vars(p_coords)
    qL = base.roll({kdim: 1}, roll_coords=False)    # qL[k] = q[k-1]
    qR = base.roll({kdim: -1}, roll_coords=False)   # qR[k] = q[k+1]
    tmp = "__out_" + pdim

    def _mat(M):  # new[i] = sum_j M[j, i] q[j]
        return xr.DataArray(M, dims=(pdim, tmp))
    out = (xr.dot(qL, _mat(Ml), dims=pdim) + xr.dot(base, _mat(Mc), dims=pdim)
           + xr.dot(qR, _mat(Mr), dims=pdim))
    out = out.rename({tmp: pdim}).transpose(*da.dims)
    if boundary == "nan":
        mask = xr.DataArray(np.r_[False, np.ones(da.sizes[kdim] - 2, bool), False], dims=(kdim,))
        out = out.where(mask)
    elif boundary != "periodic":
        raise ValueError("boundary must be 'periodic' or 'nan'")
    out = out.assign_coords({c: da.coords[c] for c in p_coords})
    out.attrs = dict(da.attrs)
    out.name = da.name
    return out


def apply_reconst2(obj, Nr: int, dims: Union[str, Iterable[str]] = "xy",
                   boundary: str = "periodic", node_type: str = "LGL"):
    """Apply Reconstruction2 to an element-form DataArray / Dataset.

    Parameters
    ----------
    obj : element-form DataArray or Dataset.
    Nr : nodes of the reconstruction polynomial (``Nnode_h1D_reconst`` /
        ``PostGLFilter_Nnodeh1D_reconst``).  E.g. "PR2p7" -> Nr = 8, "PR2p11" -> Nr = 12.
    dims : directions ("x", "y" or "xy"), applied in this order like SCALE-DG.
    boundary : ``"periodic"`` (the whole periodic domain must be loaded) or ``"nan"``
        (first/last elements become NaN; use for partial reads).
    """
    dirs = dgc._parse_dirs(dims)

    def _f(da):
        for d in dirs:
            da = _apply_1d(da, d, Nr, boundary, node_type)
        da.attrs["dg_op"] = (da.attrs.get("dg_op", "") + f" PR2_Nr{Nr}_{''.join(dirs)}").strip()
        return da
    return dgc.map_dataset(_f, obj, dirs)


# ---------------------------------------------------------------------------
# Response functions
# ---------------------------------------------------------------------------

def transfer_function_continuous(Nr: int, kdelta) -> np.ndarray:
    """Continuous-limit transfer function T_Nr(k) = sum_m (2m+1) i^m P_m(0) j_m(k Delta).

    Response of the projection-and-centre-evaluation on a window of half-width Delta
    to exp(i k x) (real because odd terms vanish).  kdelta = k * Delta = 2 pi Delta / lambda.
    """
    from scipy.special import spherical_jn
    kdelta = np.asarray(kdelta, dtype=float)
    T = np.zeros_like(kdelta)
    for m in range(0, Nr, 2):
        Pm0 = npleg.legval(0.0, np.eye(Nr)[m])
        T = T + (2 * m + 1) * (-1) ** (m // 2) * Pm0 * spherical_jn(m, kdelta)
    return T


def transfer_function_discrete(Np: int, Nr: int, kdelta, Ne: int = 64,
                               node_type: str = "LGL") -> np.ndarray:
    """Response of the discrete operator to exp(i k x) sampled at the DG nodes.

    The input is sampled on Ne periodic elements (k is rounded to a periodic
    wavenumber), filtered, and the response is the least-squares amplitude ratio.
    """
    xi = dgc.ref_nodes(Np, node_type)
    Ml, Mc, Mr = reconst2_matrices(Np, Nr, node_type)
    out = []
    for kd in np.atleast_1d(kdelta):
        n = max(1, int(round(kd * Ne / (2 * np.pi))))  # waves in the domain (Delta = 1)
        k = 2 * np.pi * n / Ne
        x = np.arange(Ne)[:, None] + 0.5 * (xi[None, :] + 1.0)
        q = np.exp(1j * k * x)
        qn = np.roll(q, 1, 0) @ Ml + q @ Mc + np.roll(q, -1, 0) @ Mr
        out.append(np.real(np.vdot(q, qn) / np.vdot(q, q)))
    return np.asarray(out)


# ---------------------------------------------------------------------------
# Self test
# ---------------------------------------------------------------------------

def _selftest(verbose: bool = True) -> None:
    for Np, Nr in ((8, 8), (8, 12), (8, 16), (5, 3)):
        Ml, Mc, Mr = reconst2_matrices(Np, Nr)
        # constants are preserved: sum of all weights = 1 for each output node
        assert np.allclose((Ml + Mc + Mr).sum(axis=0), 1.0, atol=1e-12)
        # polynomials of degree <= Nr-1 that are smooth across elements are reproduced
        xi = dgc.lgl_nodes(Np)[0]
        for deg in range(min(Nr, 6)):
            f = lambda x: (x / 3.0) ** deg
            qL, qC, qR = f(xi - 2.0), f(xi), f(xi + 2.0)
            assert np.allclose(qL @ Ml + qC @ Mc + qR @ Mr, qC, atol=1e-11), (Np, Nr, deg)

    # continuous-limit response at lambda = Delta (k Delta = 2 pi), table 1 of formulation
    T = transfer_function_continuous(8, [np.pi, 2 * np.pi, 4 * np.pi])
    assert abs(T[1] - 0.875) < 2e-3, T
    for Nr, ref in ((12, 0.999), (16, 1.000)):
        assert abs(transfer_function_continuous(Nr, 2 * np.pi) - ref) < 2e-3
    assert abs(transfer_function_continuous(8, 4 * np.pi) + 0.194) < 2e-3
    # the discrete operator approaches the continuous limit for well-resolved waves
    Td = transfer_function_discrete(16, 8, [2 * np.pi])[0]
    assert abs(Td - 0.875) < 1e-2, Td

    # xarray path: periodic sine, x then y
    da = dgc._make_test_array(NeX=8, NeY=8, Np=8, Lx=8.0, Ly=8.0,
                              func=lambda X, Y: np.sin(2 * np.pi * X / 8.0) * np.cos(2 * np.pi * Y / 8.0))
    out = apply_reconst2(da, 8)
    assert out.dims == da.dims and np.allclose(out["x"].values, da["x"].values)
    assert np.max(np.abs(out - da)) < 1e-3
    # domain integral is (approximately) conserved for a smooth field
    assert abs(float(dgc.integrate(out, "xy") - dgc.integrate(da, "xy"))) < 1e-6
    ds = xr.Dataset({"a": da})
    assert np.allclose(apply_reconst2(ds, 8, boundary="nan")["a"].isel(kx=slice(1, -1), ky=slice(1, -1)),
                       out.isel(kx=slice(1, -1), ky=slice(1, -1)))
    if verbose:
        print("patch_reconstruction: selftest passed "
              f"(T_8(2pi) = {T[1]:.3f}, T_8(4pi) = {T[2]:.3f}, discrete Np=16: {Td:.3f})")


if __name__ == "__main__":
    _selftest()
