#!/usr/bin/env python3
"""Uncertainty of the element-locking ratio R_e from time windows (e.g. 30 days).

Step 1 (slow, once per case): read the history of the runs and write the orthonormal
Legendre coefficients a of the modes (2,0) and (0,2) per element and time to
``coeffs_<case>.nc`` (PREC and the layer means of the 3-D variables of make_report.py).

Step 2 (fast): split the time axis into consecutive windows of ``--window`` days
(the incomplete last window is dropped) and compute, for each variable / layer:

  R_all     R_e = <a>/sigma(a) over (kx, ky, time) of the whole period (as make_report.py)
  R_w       R_e of each window (mean and sigma within the window)
  SE        std(R_w) / sqrt(N_w)
  r1        lag-1 autocorrelation of the series R_w
  SE_eff    SE with the effective number of windows N_eff = N_w (1 - r1) / (1 + r1)
            (only when r1 > 0)
  tau       integral time scale [day] of the element mean <a>_e(t) (sum of the
            autocorrelation of the 2.5-day series up to its first zero crossing)
  n_map, n_comp
            effective number of independent samples per day for the visibility of the
            locking in a T-day mean: the rms deviation of the window means of a from <a>
            is sigma(a)/sqrt(n T), per element (n_map: a time-mean map) or for the mean over
            the elements (n_comp: an element composite).  Used by the R_e criteria of
            summarize_Re.py (signal-to-noise |R_e| sqrt(n T)).

The R_e of each window is the mean of the (2,0) and (0,2) modes (they are taken in the
same window, so their correlation is included in SE); the modes are also given separately.

Usage (in Eh8Ez4P7/analysis/)::

    python3 window_Re.py ../asis ../PostPR2p7                  # -> window_Re.txt / .json
    python3 window_Re.py ../asis ../PostPR2p7 --window 30 60   # several window lengths
    python3 window_Re.py ../asis --recompute                   # re-read the history
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import sys
import time

import numpy as np
import xarray as xr

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_report as mr  # noqa: E402  (LIB path, VARIABLES, LAYERS, constants)
from analysis import io_common as io  # noqa: E402
from analysis import DG_common as dgc  # noqa: E402

MODES = [(2, 0), (0, 2)]


# ---------------------------------------------------------------------------
# Step 1: coefficients
# ---------------------------------------------------------------------------

def _mode_coeffs(da: xr.DataArray) -> xr.DataArray:
    a = dgc.modal_coeffs(da, "xy")
    return xr.concat([a.sel(mx=m, my=n, drop=True) for m, n in MODES],
                     dim=xr.DataArray([f"{m}{n}" for m, n in MODES], dims="mode", name="mode"))


def extract_coeffs(case_dir: str, tslice: slice, out: pathlib.Path) -> xr.Dataset:
    t0 = time.time()
    ds = io.open_dg_runs(case_dir, None, variables=mr.VARIABLES).isel(time=tslice)
    prec = (-ds["RAIN"] * 3600.0)
    pt = ds["T"] * (mr.PRE00 / ds["PRES"]) ** (mr.RDRY / mr.CPDRY)
    fields3d = {"DDENS": ds["DDENS"], "W": ds["W"], "PT": pt, "QV": ds["QV"], "MP_RHOH": ds["MP_RHOH"]}
    keep = ("mode", "kx", "ky", "time", "layer")
    res = {}

    def _clean(a):
        return a.drop_vars([c for c in a.coords if c not in keep])

    print(f"  [{time.time() - t0:7.1f} s] PREC", flush=True)
    res["PREC"] = _clean(_mode_coeffs(prec)).transpose("mode", "time", "ky", "kx").compute()
    lnames = list(mr.LAYERS)
    for v, f in fields3d.items():
        # all layers of a variable in one pass over the files
        print(f"  [{time.time() - t0:7.1f} s] {v}", flush=True)
        lms = xr.concat([dgc.mean_range(f, "z", *mr.LAYERS[ln]) for ln in lnames],
                        dim=xr.DataArray(lnames, dims="layer", name="layer"))
        a = _clean(_mode_coeffs(lms)).transpose("layer", "mode", "time", "ky", "kx").compute()
        for ln in lnames:
            res[f"{v}_{ln}"] = a.sel(layer=ln, drop=True)
    out_ds = xr.Dataset(res)
    out_ds.attrs = {"case_dir": os.path.abspath(case_dir), "runs": ds.attrs.get("dg_runs", ""),
                    "tslice": str(tslice), "layers": json.dumps(mr.LAYERS)}
    out_ds.to_netcdf(out)
    print(f"  written {out} ({time.time() - t0:.1f} s)", flush=True)
    return out_ds


# ---------------------------------------------------------------------------
# Step 2: window statistics
# ---------------------------------------------------------------------------

def _days(t: xr.DataArray) -> np.ndarray:
    v = t.values
    if np.issubdtype(v.dtype, np.datetime64):
        return (v - v[0]) / np.timedelta64(1, "D")
    # cftime objects
    return np.array([(x - v[0]).total_seconds() / 86400.0 for x in v])


def _R(a: np.ndarray) -> float:
    """R = mean/std over all samples (ddof=1)."""
    a = a.ravel()
    return float(a.mean() / a.std(ddof=1))


def _acf(x: np.ndarray, nlag: int) -> np.ndarray:
    x = x - x.mean()
    d = np.dot(x, x)
    return np.array([np.dot(x[: len(x) - k], x[k:]) / d for k in range(nlag + 1)])


def window_stats(a: xr.DataArray, days: np.ndarray, window: float) -> dict:
    """a: (mode, time, ky, kx)."""
    A = a.values
    dt = float(np.median(np.diff(days)))
    nper = int(round(window / dt))
    nw = len(days) // nper
    Rw = np.array([[_R(A[m, i * nper:(i + 1) * nper]) for i in range(nw)] for m in range(A.shape[0])])
    Rw_mean = Rw.mean(axis=0)                          # mean of the two modes per window
    out = {"window_days": window, "samples_per_window": nper, "n_windows": nw,
           "R_all": float(np.mean([_R(A[m]) for m in range(A.shape[0])])),
           "R_all_modes": [_R(A[m]) for m in range(A.shape[0])]}

    def _se(x):
        se = float(x.std(ddof=1) / np.sqrt(len(x)))
        r1 = float(_acf(x, 1)[1]) if len(x) > 2 else float("nan")
        neff = len(x) * (1 - r1) / (1 + r1) if r1 > 0 else float(len(x))
        se_eff = float(x.std(ddof=1) / np.sqrt(max(neff, 1.0)))
        return se, r1, neff, se_eff

    out["mean_Rw"] = float(Rw_mean.mean())
    out["SE"], out["r1"], out["N_eff"], out["SE_eff"] = _se(Rw_mean)
    out["modes"] = {}
    nmap, ncomp = [], []
    for m, lab in enumerate(a["mode"].values):
        se, r1, neff, se_eff = _se(Rw[m])
        # visibility noise of window means: per element (map) and element mean (composite)
        mu, sg = A[m].mean(), A[m].std(ddof=1)
        W = np.array([A[m, i * nper:(i + 1) * nper].mean(axis=0) for i in range(nw)])   # (nw, ky, kx)
        noise_map = float(np.sqrt(((W - mu) ** 2).mean()))
        noise_comp = float(np.sqrt(((W.mean(axis=(1, 2)) - mu) ** 2).mean()))
        n_map = float((sg / noise_map) ** 2 / window)       # effective samples per day
        n_comp = float((sg / noise_comp) ** 2 / window)
        nmap.append(n_map)
        ncomp.append(n_comp)
        out["modes"][str(lab)] = {"R_all": _R(A[m]), "mean_Rw": float(Rw[m].mean()), "SE": se,
                                  "r1": r1, "N_eff": neff, "SE_eff": se_eff,
                                  "noise_map_over_sigma": noise_map / float(sg),
                                  "noise_comp_over_sigma": noise_comp / float(sg),
                                  "n_map_per_day": n_map, "n_comp_per_day": n_comp}
    out["n_map_per_day"] = float(np.mean(nmap))
    out["n_comp_per_day"] = float(np.mean(ncomp))
    out["Rw"] = Rw_mean.round(4).tolist()
    return out


def time_scale(a: xr.DataArray, days: np.ndarray) -> dict:
    """Integral time scale of the element mean of a (mean over the two modes)."""
    s = a.mean(("kx", "ky")).values.mean(axis=0)     # (time,)
    dt = float(np.median(np.diff(days)))
    r = _acf(s, min(60, len(s) // 3))
    k0 = int(np.argmax(r <= 0)) if np.any(r <= 0) else len(r)
    tau = dt * (0.5 + float(np.sum(r[1:k0])))         # trapezoid-like: dt (1/2 + sum r_k)
    return {"tau_days": tau, "acf_lag": [round(float(x), 3) for x in r[:13]], "dt_days": dt}


# ---------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("cases", nargs="+", help="case directories, e.g. ../asis")
    p.add_argument("--window", type=float, nargs="+", default=[30.0], help="window lengths [day]")
    p.add_argument("--tslice", type=mr.parse_tslice, default=slice(1, None))
    p.add_argument("--recompute", action="store_true", help="re-read the history even if coeffs_<case>.nc exists")
    p.add_argument("-o", "--output", default="window_Re", help="basename of the .txt / .json output")
    args = p.parse_args(argv)

    results, lines = {}, []
    for cdir in args.cases:
        case = os.path.basename(os.path.normpath(cdir))
        nc = HERE / f"coeffs_{case}.nc"
        if nc.exists() and not args.recompute:
            print(f"{case}: reading {nc.name}", flush=True)
            ds = xr.open_dataset(nc).load()
        else:
            print(f"{case}: extracting coefficients from the history", flush=True)
            ds = extract_coeffs(cdir, args.tslice, nc)
        days = _days(ds["time"])
        results[case] = {"period": [str(ds["time"].values[0]), str(ds["time"].values[-1])],
                         "nsamples": int(ds.sizes["time"]), "vars": {}}
        lines.append(f"=== {case}: {results[case]['period'][0]} - {results[case]['period'][1]} "
                     f"({ds.sizes['time']} samples, dt = {np.median(np.diff(days)):.2f} day)")
        lines.append(f"{'variable':<12s}{'R_all':>8s}{'tau[d]':>8s}{'acf(30d)':>9s}  "
                     + "".join(f"| {'Nw':>3s}{'<R_w>':>8s}{'SE':>7s}{'r1':>7s}{'SE_eff':>8s}  " for _ in args.window))
        for v in ds.data_vars:
            ts = time_scale(ds[v], days)
            ws = [window_stats(ds[v], days, w) for w in args.window]
            results[case]["vars"][v] = {"time_scale": ts, "windows": ws}
            k30 = int(round(30.0 / ts["dt_days"]))
            acf30 = ts["acf_lag"][k30] if k30 < len(ts["acf_lag"]) else float("nan")
            lines.append(f"{v:<12s}{ws[0]['R_all']:+8.3f}{ts['tau_days']:8.1f}{acf30:+9.2f}  "
                         + "".join(f"| {w['n_windows']:3d}{w['mean_Rw']:+8.3f}{w['SE']:7.3f}"
                                   f"{w['r1']:+7.2f}{w['SE_eff']:8.3f}  " for w in ws))
        lines.append("")
    head = [f"windows: {', '.join(f'{w:g}' for w in args.window)} day",
            "R_all: R_e of the whole period (mean of (2,0) and (0,2));  <R_w>: mean of the window R_e",
            "tau: integral time scale of the element-mean coefficient;  acf(30d): its autocorrelation at lag 30 day",
            "SE = std(R_w)/sqrt(Nw);  r1: lag-1 autocorrelation of R_w;  SE_eff: SE with N_eff = Nw(1-r1)/(1+r1)",
            ""]
    txt = "\n".join(head + lines)
    print(txt)
    (HERE / f"{args.output}.txt").write_text(txt + "\n")
    with open(HERE / f"{args.output}.json", "w") as f:
        json.dump(results, f, indent=1)


if __name__ == "__main__":
    main()
