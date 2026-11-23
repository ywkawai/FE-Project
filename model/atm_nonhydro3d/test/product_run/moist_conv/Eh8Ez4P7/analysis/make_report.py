#!/usr/bin/env python3
"""Make a PDF report of one or more experiment cases of Eh8Ez4P7.

The figures and numbers are those of ``analysis_PrePR2p11_PostPR2p7.ipynb``:

  p.1  summary: case, runs, period, domain-mean precipitation time series,
       Legendre coefficients a_mn / a_00 and locking ratio R_e of the precipitation composite,
       and R_e of 3-D variables (density deviation, W, potential temperature, QV, MP_RHOH)
       averaged over representative height ranges (``LAYERS``: boundary layer,
       troposphere above it, layer above the troposphere)
  p.2  MP_RHOH: time and y mean (x-z)
  p.3  W: time and y mean (x-z)
  p.4  time-mean precipitation (x-y; pcolormesh and contourf)
  p.5  element composite of precipitation (2-D + 1-D summary)
  p.6  element polynomial of the composite, and composites per element row (ky)
  p.7  x-composites of the time and y mean MP_RHOH and W (xi-z)
  p.8  x-composites of the time and y mean DDENS, PT and QV (xi-z), as deviations from
       the horizontal mean at each height

The numbers of p.1 are also written to ``<output>.json`` for comparing cases.

Usage (in Eh8Ez4P7/analysis/)::

    python3 make_report.py ../PrePR2p11_PostPR2p7                 # -> report_PrePR2p11_PostPR2p7.pdf
    python3 make_report.py ../asis ../PrePR2p11_PostPR2p7          # one PDF per case
    python3 make_report.py ../asis --runs 1 2 --tstart 1           # select runs / skip the first step
    python3 make_report.py ../asis --tslice 20:                    # time index slice
    python3 make_report.py ../asis --range W=-0.01,0.01 --range PREC_COMP=auto
    python3 make_report.py ../asis -o reports/asis_run1.pdf
    python3 make_report.py ../asis --layer BL=0,1000 --layer UL=11000,14000   # height ranges [m]

Progress is printed while running, and the time spent in each part is summarised at the
end of each case (``--quiet`` prints only the summary line).  The data are read lazily,
so the reading of the files is included in the "compute" steps.

The plotting ranges are fixed by default (``RANGES`` below, the values used in the
notebook) so that reports of different cases can be compared directly; use
``--auto-range`` (or ``--range NAME=auto``) for the data range.
"""
from __future__ import annotations

import argparse
import contextlib
import datetime
import json
import os
import pathlib
import sys
import textwrap
import time
import warnings

# the analysis library must be imported before any netCDF file is opened
LIB = pathlib.Path(__file__).resolve().parent.joinpath("..", "..", "..", "lib").resolve()
sys.path.insert(0, str(LIB))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
import numpy as np  # noqa: E402
import xarray as xr  # noqa: E402

from analysis import io_common as io  # noqa: E402
from analysis import DG_common as dgc  # noqa: E402
from analysis import mkgraph_common as mkg  # noqa: E402

xr.set_options(keep_attrs=True)

# ---------------------------------------------------------------------------
# Settings (defaults = analysis_PrePR2p11_PostPR2p7.ipynb)
# ---------------------------------------------------------------------------

VARIABLES = ["W", "MP_RHOH", "RAIN", "DDENS", "QV", "T", "PRES"]

#: plotting ranges (vmin, vmax); None = data range.  Override with --range NAME=vmin,vmax
RANGES = {
    "MP_RHOH": (-0.01, 0.01),         # p.2  time & y mean
    "W": (-0.005, 0.005),             # p.3  time & y mean
    "PREC_MAP": (0.0, 1.0),           # p.4  time-mean precipitation (contourf)
    "PREC_COMP": (0.05, 0.25),        # p.5, p.6  precipitation composite (colour and 1-D value axis)
    "MP_RHOH_COMP": (-0.01, 0.01),    # p.7  x-composite of MP_RHOH
    "W_COMP": (-0.0025, 0.0025),      # p.7  x-composite of W
    "DDENS_COMP": (-5e-5,5e-5),       # p.8  x-composite of DDENS (deviation from the horizontal mean)
    "PT_COMP": (-0.1,0.1),            # p.8  x-composite of PT    (deviation from the horizontal mean)
    "QV_COMP": (-5e-4,5e-4),          # p.8  x-composite of QV    (deviation from the horizontal mean)
}
ZLIM = (0.0, 20.0)                    # height range of the x-z plots [km]
CMAP_XZ = "jet"
#: 2-D variables whose element locking (a_mn, R) is evaluated on p.1
LOCK_VARS = ["PREC"]

#: 3-D variables whose R_e is evaluated for the layer means of LAYERS (p.1).
#: DDENS is the deviation from the hydrostatic basic state; the basic state depends on z
#: only, so the horizontal modes (and R_e) are those of the full density.
LOCK3D_VARS = ["DDENS", "W", "PT", "QV", "MP_RHOH"]

#: 3-D variables whose x-composite (xi-z) is shown on p.8 as the deviation from the
#: horizontal mean at each height (MP_RHOH and W are on p.7 without subtraction)
DEV_COMP_VARS = ["DDENS", "PT", "QV"]

#: representative height ranges [m] for the layer means (override with --layer NAME=lo,hi)
#:   BL: boundary layer, up to ENV_PBL_Z of init.conf (1.2 km)
#:   FT: troposphere above the boundary layer, up to ENV_TRO_Z (11 km)
#:   UL: above the troposphere, up to the bottom of the sponge layer (SPONGE_HEIGHT = 15 km)
LAYERS = {
    "BL": (0.0, 1200.0),
    "FT": (1200.0, 11000.0),
    "UL": (11000.0, 15000.0),
}
#: R_e = R of the modes (mx, my) = (2, 0) and (0, 2)
RE_MODES = [(2, 0), (0, 2)]

# constants for the potential temperature (SCALE: CONST_Rdry, CONST_CPdry, CONST_PRE00)
RDRY, CPDRY, PRE00 = 287.04, 1004.64, 100000.0


# ---------------------------------------------------------------------------
# Progress and timing
# ---------------------------------------------------------------------------

class Progress:
    """Print the progress of the steps and record their elapsed times.

    Step names start with a category ("open:", "compute:", "plot:", "write:") used for
    the subtotals of :meth:`report`.
    """

    def __init__(self, label: str = "", verbose: bool = True):
        self.label = label
        self.verbose = verbose
        self.times = []                     # [(name, seconds)]
        self.t0 = time.perf_counter()

    @contextlib.contextmanager
    def step(self, name: str):
        if self.verbose:
            print(f"[{self.label}] {name} ...", flush=True)
        t = time.perf_counter()
        try:
            yield
        finally:
            dt = time.perf_counter() - t
            self.times.append((name, dt))
            if self.verbose:
                print(f"[{self.label}]   done in {dt:.1f} s", flush=True)

    def report(self, file=sys.stdout) -> None:
        total = time.perf_counter() - self.t0
        w = max([len(n) for n, _ in self.times] + [10])
        lines = [f"[{self.label}] elapsed time of each part",
                 f"  {'part':<{w}s} {'time [s]':>9s} {'share':>7s}"]
        for n, dt in self.times:
            lines.append(f"  {n:<{w}s} {dt:9.2f} {100 * dt / total:6.1f}%")
        cats = {}
        for n, dt in self.times:
            c = n.split(":")[0] if ":" in n else "other"
            cats[c] = cats.get(c, 0.0) + dt
        lines.append("  " + "-" * (w + 18))
        for c, dt in cats.items():
            lines.append(f"  {'subtotal ' + c:<{w}s} {dt:9.2f} {100 * dt / total:6.1f}%")
        other = total - sum(dt for _, dt in self.times)
        lines.append(f"  {'(not in any part)':<{w}s} {other:9.2f} {100 * other / total:6.1f}%")
        lines.append(f"  {'total':<{w}s} {total:9.2f}")
        print("\n".join(lines), file=file, flush=True)


#: progress of the case being processed (replaced for each case in make_report)
PROG = Progress(verbose=False)


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

def parse_tslice(text: str) -> slice:
    """'1:' -> slice(1, None), '20:80' -> slice(20, 80), ':-10' -> slice(None, -10)."""
    parts = [p.strip() for p in text.split(":")]
    if len(parts) not in (2, 3):
        raise argparse.ArgumentTypeError(f"invalid slice '{text}' (use start:stop[:step])")
    vals = [int(p) if p else None for p in parts]
    return slice(*vals)


def load_case(case_dir: str, runs, tslice: slice) -> dict:
    """Read the runs of a case and compute the fields used in the report."""
    with PROG.step("open: history files of the runs"), warnings.catch_warnings(record=True) as wlist:
        warnings.simplefilter("always")
        ds = io.open_dg_runs(case_dir, runs, variables=VARIABLES)
    notes = sorted({str(w.message) for w in wlist if issubclass(w.category, UserWarning)})
    ds_t = ds.isel(time=tslice)
    if ds_t.sizes["time"] == 0:
        raise ValueError(f"{case_dir}: no time step left after the time slice {tslice}")

    # precipitation rate: RAIN is negative downward; kg m-2 s-1 x 3600 s = mm/h
    prec = (-ds_t["RAIN"] * 3600.0).rename("PREC")
    prec.attrs.update(long_name="precipitation", units="mm/h")

    d = {"ds": ds, "ds_t": ds_t, "prec": prec, "notes": notes}
    with PROG.step("compute: MP_RHOH time & y mean"):
        d["rhoh_ty"] = dgc.mean(ds_t["MP_RHOH"].mean("time"), "y").compute()     # (kz, pz, kx, px)
    with PROG.step("compute: W time & y mean"):
        d["w_ty"] = dgc.mean(ds_t["W"].mean("time"), "y").compute()
    with PROG.step("compute: precipitation time mean & composites"):
        d["prec_tm"] = prec.mean("time").compute()                               # (ky, py, kx, px)
        d["prec_comp"] = dgc.composite_xy(d["prec_tm"])                          # (py, px)
        d["prec_by_row"] = dgc.composite_x(d["prec_tm"])                         # (ky, py, px)
        d["rhoh_comp_xz"] = dgc.composite_x(d["rhoh_ty"])                        # (kz, pz, px)
        d["w_comp_xz"] = dgc.composite_x(d["w_ty"])
    with PROG.step("compute: domain-mean precipitation series"):
        d["prec_dm_all"] = dgc.mean(-ds["RAIN"] * 3600.0, "xy").compute()        # (time,) incl. excluded steps

    # potential temperature (dry) from T and PRES at the nodes
    pt = (ds_t["T"] * (PRE00 / ds_t["PRES"]) ** (RDRY / CPDRY)).rename("PT")
    pt.attrs = {"long_name": "potential temperature", "units": "K"}
    fields3d = {"DDENS": ds_t["DDENS"], "W": ds_t["W"], "PT": pt, "QV": ds_t["QV"],
                "MP_RHOH": ds_t["MP_RHOH"]}

    # x-composites of the time & y mean, as deviations from the horizontal mean at each height
    d["dev_comp_xz"] = {}
    for v in DEV_COMP_VARS:
        with PROG.step(f"compute: {v} time & y mean and x-composite"):
            ty = dgc.mean(fields3d[v].mean("time"), "y").compute()       # (kz, pz, kx, px)
            comp = dgc.composite_x(ty)                                     # (kz, pz, px)
            dev = comp - dgc.element_mean(comp, "x")                       # minus the mean over xi
            dev.attrs = dict(fields3d[v].attrs)
            dev.attrs["long_name"] = f"{v} - horizontal mean"
            d["dev_comp_xz"][v] = dev

    # element locking of 3-D variables: R_e of the layer means (time, ky, py, kx, px).
    # All layer means of a variable are computed in one pass over the files.
    lock3d = {}
    lnames = list(LAYERS)
    for v in LOCK3D_VARS:
        with PROG.step(f"compute: R_e of {v} layer means"):
            lms = xr.concat([dgc.mean_range(fields3d[v], "z", *LAYERS[ln]) for ln in lnames],
                            dim=xr.DataArray(lnames, dims="layer", name="layer")).compute()
            lms.attrs = dict(fields3d[v].attrs)       # units for sigma(a)
            st = dgc.lock_stats_of(lms, "xy", sample_dims=("kx", "ky", "time"), modes=RE_MODES)
            lock3d[v] = {ln: st.sel(layer=ln) for ln in lnames}
    d["lock3d"] = lock3d

    # element locking of 2-D variables
    lock = {}
    for v in LOCK_VARS:
        with PROG.step(f"compute: a_mn and R_e of {v}"):
            var = {"PREC": prec}[v]
            comp = dgc.composite_xy(var.mean("time"))
            a = dgc.composite_legendre(comp, normalize=True).compute()
            R = dgc.lock_ratio(dgc.modal_coeffs(var), sample_dims=("kx", "ky", "time")).compute()
            st = dgc.lock_stats_of(var, "xy", sample_dims=("kx", "ky", "time"), modes=RE_MODES).compute()
            lock[v] = {"a": a, "R": R, "stats": st}
    d["lock"] = lock
    return d


def summary_numbers(d: dict, case_dir: str, tslice: slice) -> dict:
    ds, ds_t = d["ds"], d["ds_t"]
    out = {
        "case": ds.attrs.get("dg_case", os.path.basename(os.path.normpath(case_dir))),
        "case_dir": os.path.abspath(case_dir),
        "runs": ds.attrs.get("dg_runs", ""),
        "tslice": [tslice.start, tslice.stop, tslice.step],
        "period": [str(ds_t.time.values[0]), str(ds_t.time.values[-1])],
        "nsamples": int(ds_t.sizes["time"]),
        "prec_domain_mean_mm_per_h": float(dgc.mean(d["prec_tm"], "xy")),
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "notes": d["notes"],
        "layers": {k: list(v) for k, v in LAYERS.items()},
        "lock": {},
        "lock3d": {},
    }
    for v, per_layer in d["lock3d"].items():
        out["lock3d"][v] = {ln: {"R_e_20": float(R["R"].sel(mode="20")), "R_e_02": float(R["R"].sel(mode="02")),
                                 "sigma_20": float(R["std"].sel(mode="20")),
                                 "sigma_02": float(R["std"].sel(mode="02")),
                                 "units": R["std"].attrs.get("units", "")}
                            for ln, R in per_layer.items()}
    for v, L in d["lock"].items():
        a, R = L["a"], L["R"]
        out["lock"][v] = {
            "a_mn_over_a00": a.isel(mx=slice(0, 5), my=slice(0, 5)).transpose("my", "mx").values.round(5).tolist(),
            "R_e_20": float(R.sel(mx=2, my=0)),
            "R_e_02": float(R.sel(mx=0, my=2)),
            "sigma_20": float(L["stats"]["std"].sel(mode="20")),
            "sigma_02": float(L["stats"]["std"].sel(mode="02")),
            "units": L["stats"]["std"].attrs.get("units", ""),
        }
    return out


# ---------------------------------------------------------------------------
# Pages
# ---------------------------------------------------------------------------

def _rng(name):
    r = RANGES.get(name)
    return (None, None) if r is None else r


def page_summary(d: dict, info: dict):
    fig = plt.figure(figsize=(8.27, 11.69))          # A4 portrait
    fig.text(0.08, 0.95, f"Report: {info['case']}", fontsize=18, weight="bold")
    def _wrap(label, value):
        return textwrap.fill(str(value), width=88, initial_indent=f"{label:<15s}: ",
                             subsequent_indent=" " * 17, break_on_hyphens=False)
    lines = [
        _wrap("case directory", info["case_dir"]),
        f"runs           : {info['runs']}",
        f"period         : {info['period'][0]} - {info['period'][1]}  ({info['nsamples']} samples)",
        f"time slice     : {info['tslice'][0]}:{'' if info['tslice'][1] is None else info['tslice'][1]}"
        " (index along the concatenated time axis)",
        f"domain-mean precipitation (time mean) : {info['prec_domain_mean_mm_per_h']:.4f} mm/h",
        f"generated      : {info['generated']}",
    ]
    for n in info["notes"]:
        lines.append(_wrap("note", n))
    fig.text(0.08, 0.92, "\n".join(lines), fontsize=8.5, family="monospace", va="top")

    # time series of the domain-mean precipitation, coloured by run
    ax = fig.add_axes([0.10, 0.60, 0.82, 0.18])
    pdm = d["prec_dm_all"]
    t0 = pdm["time"].values[0]
    days = np.array([(t - t0).total_seconds() / 86400.0 for t in pdm["time"].values])
    for r in np.unique(pdm["run"].values):
        sel = (pdm["run"] == r).values
        ax.plot(days[sel], pdm.values[sel], ".-", ms=3, label=r)
    t_used = d["ds_t"]["time"].values
    ax.axvspan((t_used[0] - t0).total_seconds() / 86400.0, (t_used[-1] - t0).total_seconds() / 86400.0,
               color="0.9", zorder=0, label="analysed")
    ax.set_xlabel(f"days since {t0}")
    ax.set_ylabel("domain-mean precip. [mm/h]")
    ax.legend(fontsize=7, ncol=6)
    ax.set_title("domain-mean precipitation", fontsize=10)

    # element locking
    y = 0.53
    for v, L in info["lock"].items():
        a = np.array(L["a_mn_over_a00"])
        fig.text(0.08, y, f"element locking of {v}", fontsize=11, weight="bold")
        rows = ["      " + "".join(f"mx={m:<7d}" for m in range(a.shape[1]))]
        for n in range(a.shape[0]):
            rows.append(f"my={n:<2d} " + "".join(f"{a[n, m]:+8.3f} " for m in range(a.shape[1])))
        txt = ("a_mn / a_00 of the composite (orthonormal Legendre)\n" + "\n".join(rows) +
               f"\n\nR_e = <a>/sigma(a) over (kx, ky, time):  (2,0) = {L['R_e_20']:+.3f},  (0,2) = {L['R_e_02']:+.3f}"
               f"\nsigma(a) [{L['units']}]:                   (2,0) = {L['sigma_20']:.3e},  (0,2) = {L['sigma_02']:.3e}"
               "\n(positive R_e: edge-locked, negative: centre-locked)")
        fig.text(0.08, y - 0.02, txt, fontsize=8.5, family="monospace", va="top")
        y -= 0.20

    if info["lock3d"]:
        fig.text(0.08, y, "element locking of 3-D variables (layer means)", fontsize=11, weight="bold")
        layers = list(info["layers"])
        head = f"{'variable':<9s}" + "".join(f"{ln:>22s}" for ln in layers)
        rng = f"{'':<9s}" + "".join(
            f"{'(%g-%g km)' % (info['layers'][ln][0] / 1e3, info['layers'][ln][1] / 1e3):>22s}" for ln in layers)
        sub = f"{'':<9s}" + "".join(f"{'(2,0)':>11s}{'(0,2)':>11s}" for _ in layers)
        rows = [head, rng, sub]
        for v, per in info["lock3d"].items():
            rows.append(f"{v:<9s}" + "".join(f"{per[ln]['R_e_20']:+11.3f}{per[ln]['R_e_02']:+11.3f}"
                                             for ln in layers))
        srows = [f"{'variable':<9s}" + "".join(f"{ln:>22s}" for ln in layers) + "  units", sub]
        for v, per in info["lock3d"].items():
            u = per[layers[0]].get("units", "")
            srows.append(f"{v:<9s}" + "".join(f"{per[ln]['sigma_20']:11.3e}{per[ln]['sigma_02']:11.3e}"
                                              for ln in layers) + f"  {u}")
        txt = ("R_e = <a>/sigma(a) over (kx, ky, time) of the height mean over each layer\n"
               "DDENS: density deviation from the hydrostatic basic state (same horizontal modes as\n"
               "the full density); PT: potential temperature from T and PRES (dry constants)\n\n"
               + "\n".join(rows)
               + "\n\nsigma(a): standard deviation of the orthonormal Legendre coefficient a over\n"
               "(kx, ky, time), in the units of the variable\n\n"
               + "\n".join(srows))
        fig.text(0.08, y - 0.02, txt, fontsize=8.5, family="monospace", va="top")
    return fig


def page_xz(da, title, rng_name):
    fig, ax = plt.subplots(figsize=(10.0, 4.0), constrained_layout=True)
    vmin, vmax = _rng(rng_name)
    mkg.plot_tone(da, ("x", "z"), ax=ax, ylim=ZLIM, title=title, vmin=vmin, vmax=vmax, cmap=CMAP_XZ)
    return fig


def page_prec_map(d):
    fig, axs = plt.subplots(1, 2, figsize=(9, 4), constrained_layout=True)
    vmin, vmax = _rng("PREC_MAP")
    mkg.plot_tone(d["prec_tm"], ("x", "y"), ax=axs[0], edges=False, title="pcolormesh (nodal values)")
    mkg.plot_tone(d["prec_tm"], ("x", "y"), ax=axs[1], method="contourf", levels=20, edges=False,
                  title="contourf", vmin=vmin, vmax=vmax)
    fig.suptitle("time-mean precipitation")
    return fig


def page_comp_summary(d):
    vmin, vmax = _rng("PREC_COMP")
    ylim = (vmin, vmax) if vmin is not None else None
    fig, _ = mkg.plot_composite_summary(d["prec_comp"], suptitle="element composite of precipitation",
                                        vmin=vmin, vmax=vmax, ylim=ylim, figsize=(12, 4.2))
    return fig


def page_comp_rows(d):
    vmin, vmax = _rng("PREC_COMP")
    ylim = (vmin, vmax) if vmin is not None else "auto"
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    mkg.plot_composite2d(d["prec_comp"], ax=axs[0], refine=81, title="element polynomial",
                         vmin=vmin, vmax=vmax)
    mkg.plot_composite1d(d["prec_by_row"], "x", ax=axs[1], hue="ky", lw=0.8, nodes=False, ylim=ylim)
    mkg.plot_composite1d(d["prec_comp"], "x", ax=axs[1], color="k", lw=2, label="all elements", ylim=ylim)
    axs[1].set_title("mean over η (thin: composite per element row ky)")
    return fig


def _comp_xz_panels(entries, figsize):
    """One xi-z tone panel per (DataArray, range name, title)."""
    fig, axs = plt.subplots(1, len(entries), figsize=figsize, constrained_layout=True, squeeze=False)
    for ax, (da, name, title) in zip(axs[0], entries):
        vmin, vmax = _rng(name)
        mkg.plot_tone(da, ("x", "z"), ax=ax, ylim=ZLIM, vmin=vmin, vmax=vmax, cmap=CMAP_XZ, title=title)
    return fig


def page_comp_xz(d):
    return _comp_xz_panels([(d["rhoh_comp_xz"], "MP_RHOH_COMP", "MP_RHOH: x-composite of the time and y mean"),
                            (d["w_comp_xz"], "W_COMP", "W: x-composite of the time and y mean")],
                           figsize=(10, 4.6))


def page_dev_comp_xz(d):
    entries = [(d["dev_comp_xz"][v], f"{v}_COMP", f"{v} (deviation)")
               for v in DEV_COMP_VARS]
    fig = _comp_xz_panels(entries, figsize=(5.0 * len(entries), 4.6))
    fig.suptitle("x-composites of the time and y mean (deviation from the horizontal mean at each height)")
    return fig


def make_report(case_dir: str, runs, tslice: slice, output: str, write_json: bool = True,
                verbose: bool = True) -> dict:
    global PROG
    name = os.path.basename(os.path.normpath(case_dir))
    PROG = Progress(name, verbose)
    d = load_case(case_dir, runs, tslice)
    info = summary_numbers(d, case_dir, tslice)
    pages = [
        ("summary", lambda: page_summary(d, info)),
        ("MP_RHOH x-z", lambda: page_xz(d["rhoh_ty"], "MP_RHOH: time and y mean", "MP_RHOH")),
        ("W x-z", lambda: page_xz(d["w_ty"], "W: time and y mean", "W")),
        ("precipitation map", lambda: page_prec_map(d)),
        ("precipitation composite", lambda: page_comp_summary(d)),
        ("composite per row", lambda: page_comp_rows(d)),
        ("x-composites xi-z", lambda: page_comp_xz(d)),
        ("x-composites xi-z (deviations)", lambda: page_dev_comp_xz(d)),
    ]
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    pdf = PdfPages(output)
    try:
        for i, (title, mk) in enumerate(pages, 1):
            with PROG.step(f"plot: p.{i}/{len(pages)} {title}"):
                fig = mk()
                eng = fig.get_layout_engine()
                if eng is not None and hasattr(eng, "set"):      # leave room for the footer
                    try:
                        eng.set(rect=(0.0, 0.035, 1.0, 0.965))
                    except TypeError:
                        pass
                fig.text(0.99, 0.005, f"{info['case']}  ({info['runs']})  p.{i}", ha="right", va="bottom",
                         fontsize=7, color="0.4")
                pdf.savefig(fig)
                plt.close(fig)
        meta = pdf.infodict()
        meta["Title"] = f"Report: {info['case']}"
        meta["Subject"] = f"runs {info['runs']}, {info['period'][0]} - {info['period'][1]}"
    finally:
        with PROG.step("write: close PDF"):
            pdf.close()
    if write_json:
        with PROG.step("write: JSON"):
            with open(os.path.splitext(output)[0] + ".json", "w") as f:
                json.dump(info, f, indent=1, ensure_ascii=False)
    return info


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__.split("Usage (in Eh8Ez4P7/analysis/)::")[1])
    p.add_argument("cases", nargs="+", help="case directories containing run1, run2, ...")
    p.add_argument("--runs", nargs="+", default=None,
                   help="runs to read (numbers or names; default: all runs with history output)")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--tstart", type=int, default=None,
                   help="first time index used for the time means (default 1: skip step 0)")
    g.add_argument("--tslice", type=parse_tslice, default=None, help="time index slice, e.g. 20: or 10:80")
    p.add_argument("-o", "--output", default=None,
                   help="output PDF (only with one case). Default: report_<case>.pdf in this directory")
    p.add_argument("--range", action="append", default=[], metavar="NAME=VMIN,VMAX",
                   help=f"plotting range; NAME in {sorted(RANGES)}; VMIN,VMAX or 'auto'")
    p.add_argument("--auto-range", action="store_true", help="use the data range for all plots")
    p.add_argument("--layer", action="append", default=[], metavar="NAME=LO,HI",
                   help=f"height range [m] of a layer for the R_e of 3-D variables (default {LAYERS}); "
                        "a new NAME adds a layer")
    p.add_argument("--no-json", action="store_true", help="do not write the summary JSON")
    p.add_argument("-q", "--quiet", action="store_true",
                   help="do not print the progress (the timing summary is still printed)")
    args = p.parse_args(argv)

    if args.output and len(args.cases) > 1:
        p.error("--output can be used with one case only")
    if args.auto_range:
        for k in RANGES:
            RANGES[k] = None
    for item in args.range:
        name, _, val = item.partition("=")
        if name not in RANGES:
            p.error(f"unknown range name '{name}' (choose from {sorted(RANGES)})")
        if val.strip().lower() == "auto":
            RANGES[name] = None
        else:
            try:
                lo, hi = (float(v) for v in val.split(","))
            except ValueError:
                p.error(f"invalid range '{item}' (use NAME=VMIN,VMAX)")
            RANGES[name] = (lo, hi)

    for item in args.layer:
        name, _, val = item.partition("=")
        try:
            lo, hi = (float(v) for v in val.split(","))
        except ValueError:
            p.error(f"invalid layer '{item}' (use NAME=LO,HI in m)")
        if not hi > lo:
            p.error(f"invalid layer '{item}' (HI must be larger than LO)")
        LAYERS[name] = (lo, hi)

    runs = None
    if args.runs is not None:
        runs = [int(r) if r.isdigit() else r for r in args.runs]
    if args.tslice is not None:
        tslice = args.tslice
    else:
        tslice = slice(1 if args.tstart is None else args.tstart, None)

    here = pathlib.Path(__file__).resolve().parent
    status = 0
    for case in args.cases:
        name = os.path.basename(os.path.normpath(case))
        out = args.output or str(here / f"report_{name}.pdf")
        try:
            info = make_report(case, runs, tslice, out, write_json=not args.no_json,
                               verbose=not args.quiet)
        except Exception as e:  # keep going with the other cases
            print(f"[{name}] FAILED: {e}", file=sys.stderr)
            if PROG.label == name and PROG.times:
                PROG.report(file=sys.stderr)
            status = 1
            continue
        PROG.report()
        lk = info["lock"].get("PREC")
        extra = f", R_e(PREC) = {lk['R_e_20']:+.3f} / {lk['R_e_02']:+.3f}" if lk else ""
        print(f"[{name}] runs {info['runs']}, {info['nsamples']} samples "
              f"({info['period'][0]} - {info['period'][1]}){extra} -> {out}")
    return status


if __name__ == "__main__":
    sys.exit(main())
