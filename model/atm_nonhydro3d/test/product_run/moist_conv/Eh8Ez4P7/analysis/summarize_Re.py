#!/usr/bin/env python3
"""Collect the element-locking ratios R_e of all experiment cases into one PDF report.

Reads the ``report_<case>.json`` files written by ``make_report.py`` and the window
statistics ``window_Re.json`` written by ``window_Re.py`` (no model output is read, so it
runs in a few seconds), and typesets ``Re_summary.pdf`` (A4 portrait, Japanese text) with
XeLaTeX.

The document is ``Re_report.md``: the text (in Japanese) and the section structure are
written there by hand in a simple Markdown, and the tables and figures made from the data
are inserted at the markers ``<!-- auto:NAME -->``:

  meta            data sources, date and parameters
  cases           table of the cases (pre-/post-processing, runs, samples, period) and
                  the warnings about cases that are not directly comparable
  re-heatmap      heat map of R_e (mean of the (2,0) and (0,2) modes) +- SE
  criteria-def    definition of the R_e criteria and their table
  criteria-table  tables of the criteria for the --criteria-vars (default PREC)
  criteria-tvis   T_vis figure for the --criteria-vars
  criteria-map    status of the criteria for every variable / layer and case
  tvis-map        T_vis of the map and of the composite for every variable / layer and case,
                  coloured by the classes of the criterion periods
  tvis-summary    per case, the variable whose locking becomes visible first (min T_vis)
                  and the ratio T_vis(map) / T_vis(composite)
  variables       one figure per variable: R_e (+- SE) and sigma(a) of each case

The figures (English) are written to ``<output>_build/`` together with the LaTeX source.
Without XeLaTeX (or if it fails) only the figures are collected into the PDF.

Standard errors
---------------
SE of R_e comes from ``window_Re.json`` (time windows of ``--window-days`` days; SE_eff,
corrected for the lag-1 autocorrelation of the window values).  Cases without window
results are shown without SE.

R_e criteria
------------
The element locking cannot be removed completely at a finite resolution, so whether it is
acceptable depends on how the output is used.  A criterion is (kind, T): the locking must
not be visible in a T-day mean, either

  map        in a time-mean horizontal map (each element on its own), or
  composite  in the element composite (mean over all elements).

The locking is visible when the signal <a> exceeds kappa times the noise of the T-day mean,
sigma(a)/sqrt(n T), i.e. when |R_e| sqrt(n T) > kappa.  n is the effective number of
independent samples per day measured by window_Re.py (n_map per element, n_comp for the
element mean) of the window means, measured for the window lengths with at least
MIN_WINDOWS windows.  n depends on the averaging period, so n(T) is interpolated in log T
between the window lengths (constant outside their range):

  R_crit(T) = kappa / sqrt(n(T) T),   T_vis: |R_e| sqrt(n(T_vis) T_vis) = kappa

T_vis is the averaging period at which the locking essentially starts to be visible; a
criterion (kind, T) is met (best estimate) exactly when T_vis of that kind exceeds T.

The status uses |R_e| +- 2 SE: "pass" if even |R_e| + 2 SE is below R_crit, "fail" if
|R_e| - 2 SE is above it, otherwise "borderline".  Default criteria (--criterion):

  practical-30      map,         30 days
  practical-100     map,        100 days
  practical-1000    map,       1000 days
  intermediate-30   composite,   30 days
  intermediate-100  composite,  100 days
  strict            composite, 1000 days

i.e. the two kinds of figures times the averaging periods 30, 100 and 1000 days.

The values are also written to ``<output>.csv`` and the criteria to ``<output>_criteria.csv``.

Usage (in Eh8Ez4P7/analysis/)::

    python3 make_report.py ../asis ../PostPR2p7 ...        # make / update the JSON files first
    python3 window_Re.py ../asis ../PostPR2p7 ... --window 30 60 120   # SE and n
    python3 summarize_Re.py                                # -> Re_summary.pdf
    python3 summarize_Re.py --cases asis PostPR2p7 PrePR2p11_PostPR2p7   # subset and order
    python3 summarize_Re.py --kappa 1 --criterion practical=map:1000 --criterion strict=composite:1000
    python3 summarize_Re.py --criteria-vars PREC W_FT       # tables / T_vis for more variables
    python3 summarize_Re.py --document my_report.md -o reports/Re_summary.pdf
"""
from __future__ import annotations

import argparse
import csv
import datetime
import glob
import json
import os
import pathlib
import re
import shutil
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
import numpy as np  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent

#: reference case shown first by default
REFERENCE_CASE = "asis"
#: modes of R_e and their markers
MODES = (("20", "(2,0)", "o"), ("02", "(0,2)", "^"))
#: colour range of the heat map (|R_e|)
HEATMAP_VMAX = 0.6
#: |R_e| < SIG_FACTOR * SE is shown as not significant
SIG_FACTOR = 2.0
#: default R_e criteria: (name, kind, averaging period [day]); kind = "map" or "composite"
CRITERIA = [("practical-30", "map", 30.0),
            ("practical-100", "map", 100.0),
            ("practical-1000", "map", 1000.0),
            ("intermediate-30", "composite", 30.0),
            ("intermediate-100", "composite", 100.0),
            ("strict", "composite", 1000.0)]
#: Japanese names of the criteria (for the text and tables)
CRITERIA_JA = {"practical-30": "実用 30 日", "practical-100": "実用 100 日", "practical-1000": "実用 1000 日",
               "intermediate-30": "中間 30 日", "intermediate-100": "中間 100 日", "strict": "厳しい"}
KIND_JA = {"map": "時間平均の地図", "composite": "要素コンポジット"}
KIND_SHORT_JA = {"map": "地図", "composite": "コンポジット"}
#: visibility threshold of the signal-to-noise ratio |R_e| sqrt(n T)
KAPPA = 1.0
#: n (samples per day) is taken from window lengths with at least this many windows
MIN_WINDOWS = 8
#: status: symbol in the figures, symbol in the text, colour
STATUS = {"pass": ("o", "○", "#a6dba0"), "borderline": ("~", "△", "#fee08b"),
          "fail": ("x", "×", "#f4a582"), "n/a": ("-", "—", "#e0e0e0")}
#: width of the figures [inch] (A4 portrait text width)
FIGW = 6.6

plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
                     "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "legend.fontsize": 8,
                     "pdf.fonttype": 42})


# ---------------------------------------------------------------------------
# Collecting
# ---------------------------------------------------------------------------

def load_reports(paths):
    """{case: report dict} from the JSON files (later files win for duplicated cases)."""
    reports = {}
    for p in paths:
        with open(p) as f:
            d = json.load(f)
        d["_file"] = os.path.abspath(p)
        if d.get("case") in reports:
            print(f"warning: case '{d['case']}' appears twice; using {p}", file=sys.stderr)
        reports[d["case"]] = d
    return reports


def order_cases(reports, cases=None):
    if cases:
        missing = [c for c in cases if c not in reports]
        if missing:
            raise SystemExit(f"no JSON for case(s): {', '.join(missing)} (run make_report.py first)")
        return list(cases)
    names = sorted(reports)
    if REFERENCE_CASE in names:
        names.remove(REFERENCE_CASE)
        names.insert(0, REFERENCE_CASE)
    return names


def collect(reports, cases):
    """Rows of (variable, layer, case, mode, R_e, sigma, units); layer is '' for 2-D variables.

    Returns (rows, variables, layers) where variables keeps the order of appearance
    (2-D variables first) and layers the order of the first case that has them.
    """
    rows, variables, layers = [], [], []
    for c in cases:
        d = reports[c]
        for v, L in d.get("lock", {}).items():
            if v not in variables:
                variables.append(v)
            for m, _, _ in MODES:
                rows.append((v, "", c, m, L.get(f"R_e_{m}", np.nan), L.get(f"sigma_{m}", np.nan),
                             L.get("units", "")))
    for c in cases:
        d = reports[c]
        for ln in d.get("layers", {}):
            if ln not in layers:
                layers.append(ln)
        for v, per in d.get("lock3d", {}).items():
            if v not in variables:
                variables.append(v)
            for ln, L in per.items():
                for m, _, _ in MODES:
                    rows.append((v, ln, c, m, L.get(f"R_e_{m}", np.nan), L.get(f"sigma_{m}", np.nan),
                                 L.get("units", "")))
    return rows, variables, layers


def lookup(rows):
    return {(v, ln, c, m): (r, s, u) for v, ln, c, m, r, s, u in rows}


def row_keys(tab, variables, layers):
    """[(variable, layer)] in the display order."""
    keys = []
    for v in variables:
        lns = [""] if any(k[0] == v and k[1] == "" for k in tab) else layers
        keys += [(v, ln) for ln in lns]
    return keys


def _split_key(key):
    return (key, "") if "_" not in key else tuple(key.rsplit("_", 1))


def load_window_se(path, window_days):
    """SE of R_e from window_Re.json.

    Returns (se, info): se[(variable, layer, case, mode)] = SE_eff of that mode, with mode
    "mean" for the mean of the (2,0) and (0,2) modes; info[case] = dict(period, nsamples,
    n_windows).  Variable keys of window_Re.json are "PREC" or "<var>_<layer>".
    """
    se, info = {}, {}
    if not path or not os.path.exists(path):
        return se, info
    with open(path) as f:
        d = json.load(f)
    for c, dc in d.items():
        for key, x in dc.get("vars", {}).items():
            ws = [w for w in x.get("windows", []) if abs(w["window_days"] - window_days) < 1e-6]
            if not ws:
                continue
            w = ws[0]
            v, ln = _split_key(key)
            se[(v, ln, c, "mean")] = w.get("SE_eff", w["SE"])
            for m, mm in w.get("modes", {}).items():
                se[(v, ln, c, m)] = mm.get("SE_eff", mm["SE"])
            info[c] = {"period": dc.get("period"), "nsamples": dc.get("nsamples"),
                       "n_windows": w["n_windows"]}
    return se, info


def load_window_n(path, min_windows=MIN_WINDOWS):
    """Effective samples per day for the visibility: nv[(variable, layer, case)] =
    {window_days: (n_map, n_comp)} of the window lengths with >= min_windows windows."""
    nv = {}
    if not path or not os.path.exists(path):
        return nv
    with open(path) as f:
        d = json.load(f)
    for c, dc in d.items():
        for key, x in dc.get("vars", {}).items():
            ws = [w for w in x.get("windows", []) if w.get("n_windows", 0) >= min_windows
                  and "n_map_per_day" in w]
            if not ws:
                continue
            nv[(*_split_key(key), c)] = {w["window_days"]: (w["n_map_per_day"], w["n_comp_per_day"])
                                         for w in ws}
    return nv


def n_of_T(nwins, T, kind):
    """n(T) [per day] of ``kind`` ("map" or "composite"): linear interpolation in log T of the
    measured window lengths, constant outside their range."""
    ws = sorted(nwins)
    lw = np.log(ws)
    vals = [nwins[w][0 if kind == "map" else 1] for w in ws]
    return float(np.interp(np.log(T), lw, vals))


def solve_tvis(nwins, R, kind, kappa, tmin=1e-3, tmax=1e8):
    """T_vis: the averaging period with |R| sqrt(n(T) T) = kappa (bisection in log T)."""
    R = abs(R)
    if R == 0:
        return np.inf

    def f(T):
        return R * np.sqrt(n_of_T(nwins, T, kind) * T) - kappa
    if f(tmax) < 0:
        return np.inf
    if f(tmin) > 0:
        return tmin
    a, b = np.log(tmin), np.log(tmax)
    for _ in range(80):
        m = 0.5 * (a + b)
        if f(np.exp(m)) > 0:
            b = m
        else:
            a = m
    return float(np.exp(0.5 * (a + b)))


def evaluate_criteria(rows, se, nv, criteria, kappa):
    """crit[(variable, layer, case)] = dict(R, SE, n_map, n_comp, nwin (longest window, for the
    tables), T_vis_map, T_vis_comp (each (best, lo, hi), lo/hi from |R_e| -+ 2 SE),
    status {name: (status, R_crit, S/N)}).

    n depends on the averaging period: n(T) is interpolated in log T between the window lengths
    (constant outside), and is used both for the criteria and for T_vis, so that a criterion
    (kind, T) is met (best estimate) exactly when T_vis of that kind exceeds T."""
    tab = lookup(rows)
    keys = sorted({(v, ln, c) for v, ln, c, _ in tab})
    crit = {}
    for v, ln, c in keys:
        r = np.nanmean([tab.get((v, ln, c, m), (np.nan,))[0] for m, _, _ in MODES])
        e = se.get((v, ln, c, "mean"), np.nan)
        n = nv.get((v, ln, c))
        out = {"R": r, "SE": e, "status": {}}
        if n is None or not np.isfinite(r):
            out.update(n_map=np.nan, n_comp=np.nan, nwin=np.nan,
                       T_vis_map=(np.nan,) * 3, T_vis_comp=(np.nan,) * 3)
            for name, _, _ in criteria:
                out["status"][name] = ("n/a", np.nan, np.nan)
            crit[(v, ln, c)] = out
            continue
        wmax = max(n)
        ee = e if np.isfinite(e) else 0.0
        hi, lo = abs(r) + 2 * ee, max(abs(r) - 2 * ee, 0.0)
        out.update(n_map=n[wmax][0], n_comp=n[wmax][1], nwin=wmax)
        for kind, key in (("map", "T_vis_map"), ("composite", "T_vis_comp")):
            out[key] = tuple(solve_tvis(n, rr, kind, kappa) for rr in (abs(r), hi, lo))
        for name, kind, T in criteria:
            nn = n_of_T(n, T, kind)
            rc = kappa / np.sqrt(nn * T)
            st = "pass" if hi < rc else ("fail" if lo > rc else "borderline")
            out["status"][name] = (st, rc, abs(r) * np.sqrt(nn * T))
        crit[(v, ln, c)] = out
    return crit


def _translate_note(note):
    m = re.match(r"no history output in: (.*) \(skipped\)", note)
    if m:
        return f"history の出力がない run: {m.group(1)} (除外した)"
    return note


def warnings_for(reports, cases, min_samples, winfo=None, window_days=30.0):
    """Warnings (Japanese) about cases that are not directly comparable."""
    msgs = []
    for c in cases:
        n = reports[c].get("nsamples", 0)
        if n < min_samples:
            msgs.append(f"{c}: サンプル数が {n} しかない (< {min_samples})。$R_e$ はまだ信頼できない。")
        for note in reports[c].get("notes", []):
            msgs.append(f"{c}: {_translate_note(note)}")
    periods = {c: tuple(reports[c].get("period", [])) for c in cases}
    if len(set(periods.values())) > 1:
        msgs.append("解析期間がケースごとに違う (表 @tab:cases)。")
    layer_defs = {c: json.dumps(reports[c].get("layers", {}), sort_keys=True) for c in cases}
    if len(set(layer_defs.values())) > 1:
        msgs.append("層の定義がケースごとに違う: "
                    + "; ".join(f"{c}: {reports[c].get('layers')}" for c in cases))
    if winfo is not None:
        no_se = [c for c in cases if c not in winfo]
        if no_se:
            msgs.append(f"window_Re.json に結果がないケース: {', '.join(no_se)}。SE と criterion は示さない。")
        for c in cases:
            if c in winfo and list(winfo[c]["period"]) != list(reports[c].get("period", [])):
                msgs.append(f"{c}: SE ({window_days:g} 日窓) の期間 {winfo[c]['period'][0]} – "
                            f"{winfo[c]['period'][1]} が、$R_e$ の期間と違う。"
                            "make_report.py か window_Re.py を再実行すること。")
    return msgs


# ---------------------------------------------------------------------------
# Figures (English)
# ---------------------------------------------------------------------------

def fig_heatmap(rows, variables, layers, cases, se):
    tab = lookup(rows)
    keys = row_keys(tab, variables, layers)
    labels = [v if ln == "" else f"{v} {ln}" for v, ln in keys]
    vals = np.array([[np.nanmean([tab.get((v, ln, c, m), (np.nan,))[0] for m, _, _ in MODES])
                      for c in cases] for v, ln in keys], dtype=float)
    ses = np.array([[se.get((v, ln, c, "mean"), np.nan) for c in cases] for v, ln in keys], dtype=float)
    fig, ax = plt.subplots(figsize=(FIGW, 0.42 * len(labels) + 1.4), constrained_layout=True)
    im = ax.imshow(vals, cmap="RdBu_r", vmin=-HEATMAP_VMAX, vmax=HEATMAP_VMAX, aspect="auto")
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            if not np.isfinite(vals[i, j]):
                ax.text(j, i, "n/a", ha="center", va="center", fontsize=7, color="0.5")
                continue
            txt = f"{vals[i, j]:+.3f}"
            if np.isfinite(ses[i, j]):
                if abs(vals[i, j]) < SIG_FACTOR * ses[i, j]:
                    txt = f"({txt})"
                txt += f"\n\u00b1{ses[i, j]:.3f}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=7.5, linespacing=1.1,
                    color="white" if abs(vals[i, j]) > 0.6 * HEATMAP_VMAX else "black")
    ax.set_xticks(range(len(cases)))
    ax.set_xticklabels(cases, rotation=20, ha="right")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    i0 = 0
    for v in variables:
        n = sum(1 for k in keys if k[0] == v)
        if i0 > 0:
            ax.axhline(i0 - 0.5, color="k", lw=0.8)
        i0 += n
    ax.xaxis.tick_top()
    plt.setp(ax.get_xticklabels(), rotation=20, ha="left")
    fig.colorbar(im, ax=ax, label="R_e (mean of (2,0) and (0,2))", shrink=0.6, pad=0.01)
    return fig


def fig_variable(rows, v, layers, cases, reports, se):
    tab = lookup(rows)
    is2d = any(k[0] == v and k[1] == "" for k in tab)
    lns = [""] if is2d else layers
    units = next((u for (vv, _, _, _), (_, _, u) in tab.items() if vv == v and u), "")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(FIGW, 3.9), sharex=True, constrained_layout=True,
                                   gridspec_kw={"height_ratios": [2.2, 1.0]})
    x = np.arange(len(cases))
    nl = len(lns)
    width = 0.66
    for il, ln in enumerate(lns):
        off = (il - (nl - 1) / 2) * width / max(nl, 1)
        color = "k" if is2d else f"C{il}"
        for im_, (m, mlab, mk) in enumerate(MODES):
            dx = off + (im_ - 0.5) * 0.08
            r = [tab.get((v, ln, c, m), (np.nan, np.nan, ""))[0] for c in cases]
            sg = [tab.get((v, ln, c, m), (np.nan, np.nan, ""))[1] for c in cases]
            e = np.array([se.get((v, ln, c, m), np.nan) for c in cases], dtype=float)
            lab = f"{mlab}" if is2d else f"{ln} {mlab}"
            ax1.errorbar(x + dx, r, yerr=np.where(np.isfinite(e), e, 0.0), ls="none", marker=mk, ms=6,
                         color=color, label=lab, mfc=color if m == "20" else "white",
                         elinewidth=1.1, capsize=2.5)
            ax2.plot(x + dx, sg, ls="none", marker=mk, ms=5, color=color,
                     mfc=color if m == "20" else "white")
    ax1.axhline(0.0, color="0.4", lw=0.8)
    for xs in x[:-1] + 0.5:
        ax1.axvline(xs, color="0.88", lw=0.8)
        ax2.axvline(xs, color="0.88", lw=0.8)
    ax1.set_ylabel("R_e  (+: edge, -: centre)")
    ax1.legend(ncol=len(lns), loc="best", framealpha=0.9)
    ax1.set_title(f"{v}" + ("" if is2d else " (layer means)"), loc="left", weight="bold")
    ax2.set_ylabel(f"sigma(a) [{units}]" if units else "sigma(a)")
    sg_all = np.array([s for (vv, _, _, _), (_, s, _) in tab.items() if vv == v], dtype=float)
    if np.all(np.isfinite(sg_all)) and np.all(sg_all > 0) and sg_all.max() / sg_all.min() > 20:
        ax2.set_yscale("log")
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"{c}\n(n={reports[c].get('nsamples', 0)})" for c in cases], rotation=15, ha="right")
    return fig


def fig_tvis(crit, cases, criteria, var="PREC", layer=""):
    fig, ax = plt.subplots(figsize=(FIGW, 4.0), constrained_layout=True)
    x = np.arange(len(cases))
    tmax = 1e6
    shown = set()
    for kind, key, dx, mk, col in (("map", "T_vis_map", -0.12, "s", "C0"),
                                   ("composite", "T_vis_comp", 0.12, "o", "C1")):
        for i, c in enumerate(cases):
            xx = crit.get((var, layer, c))
            if xx is None or not np.isfinite(xx[key][0]):
                continue
            t, tlo, thi = (min(z, tmax) for z in xx[key])
            ax.errorbar(x[i] + dx, t, yerr=[[t - tlo], [thi - t]], marker=mk, color=col, ms=7,
                        capsize=3, label=None if kind in shown else f"T_vis ({kind})",
                        lolims=(thi >= tmax))
            shown.add(kind)
    groups = {}
    for name, kind, T in criteria:
        groups.setdefault(T, []).append((name, kind))
    for T, lst in groups.items():
        for name, kind in lst:
            ax.axhline(T, color="C0" if kind == "map" else "C1", ls="--" if kind == "map" else ":", lw=1.2)
        txt = " / ".join(f"{name} ({kind})" for name, kind in lst) + f": {T:g} d"
        ax.text(-0.45, T * 1.08, txt, color="0.25", va="bottom", ha="left", fontsize=7.5)
    ax.set_yscale("log")
    ax.set_ylim(0.1, tmax)
    ax.set_xticks(x)
    ax.set_xticklabels(cases, rotation=15, ha="right")
    ax.set_xlim(-0.5, len(cases) - 0.5)
    ax.set_ylabel("T_vis [day]")
    lab = var if not layer else f"{var} {layer}"
    ax.set_title(f"{lab}: averaging period beyond which the locking is visible", loc="left")
    ax.legend(loc="lower right")
    ax.grid(True, which="major", axis="y", color="0.9")
    return fig


def sample_dt(path, default=2.5):
    """Output interval [day] of the history (from window_Re.json)."""
    try:
        with open(path) as f:
            d = json.load(f)
        for dc in d.values():
            for x in dc.get("vars", {}).values():
                return float(x["time_scale"]["dt_days"])
    except (OSError, KeyError, ValueError, TypeError):
        pass
    return default


def _tvis_text(t, dt):
    if not np.isfinite(t):
        return "inf"
    if t < dt:
        return f"<{dt:g}"
    if t < 1e4:
        return f"{t:.0f}"
    mant, e = f"{t:.1e}".split("e")            # e.g. "1.2", "+04"
    return f"{mant.removesuffix('.0')}e{int(e)}"


def fig_tvis_map(crit, variables, layers, cases, criteria, tab, dt):
    """T_vis (map | composite) for every variable / layer and case, in classes of the
    criterion periods."""
    from matplotlib.colors import BoundaryNorm, ListedColormap
    keys = row_keys(tab, variables, layers)
    labels = [v if ln == "" else f"{v} {ln}" for v, ln in keys]
    periods = sorted({T for _, _, T in criteria})
    bounds = [dt] + periods + [periods[-1] * 10]
    edges = [1e-9] + bounds + [1e12]
    cols = ["#b2182b", "#ef8a62", "#fddbc7", "#d1e5f0", "#67a9cf", "#2166ac", "#053061"][:len(edges) - 1]
    cmap, norm = ListedColormap(cols), BoundaryNorm(edges, len(cols))
    fig, axes = plt.subplots(1, 2, figsize=(FIGW, 0.25 * len(labels) + 1.7), sharey=True,
                             constrained_layout=True)
    for ax, (kind, key) in zip(axes, (("map", "T_vis_map"), ("composite", "T_vis_comp"))):
        vals = np.full((len(keys), len(cases)), np.nan)
        for i, (v, ln) in enumerate(keys):
            for j, c in enumerate(cases):
                x = crit.get((v, ln, c))
                if x is None or not np.isfinite(x[key][0]) and not np.isinf(x[key][0]):
                    continue
                vals[i, j] = min(x[key][0], 1e11)
                insig = np.isfinite(x["SE"]) and abs(x["R"]) < SIG_FACTOR * x["SE"]
                txt = _tvis_text(x[key][0], dt)
                ax.text(j, i, f"({txt})" if insig else txt, ha="center", va="center", fontsize=7,
                        color="white" if vals[i, j] < periods[0] or vals[i, j] > periods[-1] * 10 else "black")
        im = ax.imshow(np.ma.masked_invalid(vals), cmap=cmap, norm=norm, aspect="auto")
        ax.set_facecolor("#e0e0e0")
        ax.set_xticks(range(len(cases)))
        ax.set_xticklabels(cases, rotation=35, ha="right", fontsize=7.5)
        ax.set_title(f"T_vis [day]: {kind}", fontsize=9.5)
        i0 = 0
        for v in variables:
            n = sum(1 for k in keys if k[0] == v)
            if i0 > 0:
                ax.axhline(i0 - 0.5, color="k", lw=0.6)
            i0 += n
    axes[0].set_yticks(range(len(labels)))
    axes[0].set_yticklabels(labels, fontsize=8)
    cb = fig.colorbar(im, ax=axes, shrink=0.7, pad=0.01, ticks=bounds)
    cb.ax.set_yticklabels([f"{b:g}" for b in bounds], fontsize=7.5)
    cb.set_label("T_vis [day]  (class boundaries: sampling interval and criterion periods)", fontsize=7.5)
    return fig


def auto_tvis_summary(crit, variables, layers, cases, tab, dt):
    """Per case: the variable / layer with the smallest T_vis (map, composite) and the ratio."""
    keys = row_keys(tab, variables, layers)
    tex = [r"\begin{table}[H]", r"\centering\footnotesize",
           r"\caption{各ケースで、要素ロックが最も早く見え出す変数・層と、その $T_\mathrm{vis}$。"
           r"比は全変数・層の $T_\mathrm{vis}$ (地図) $/$ $T_\mathrm{vis}$ (コンポジット) の中央値と範囲。}",
           r"\label{tab:tvis-summary}",
           r"\begin{tabular}{@{}llrllr@{}}", r"\toprule",
           r"ケース & \multicolumn{2}{c}{時間平均の地図} & \multicolumn{2}{c}{要素コンポジット} & "
           r"比 (地図/コンポ) \\", r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}",
           r" & 変数・層 & $T_\mathrm{vis}$ [日] & 変数・層 & $T_\mathrm{vis}$ [日] & \\", r"\midrule"]
    for c in cases:
        rows = [((v, ln), crit[(v, ln, c)]) for v, ln in keys if (v, ln, c) in crit
                and np.isfinite(crit[(v, ln, c)]["T_vis_map"][0])]
        if not rows:
            tex.append(f"{_tt(c)} & — & — & — & — & — \\\\")
            continue
        km, xm = min(rows, key=lambda r: r[1]["T_vis_map"][0])
        kc, xc = min(rows, key=lambda r: r[1]["T_vis_comp"][0])
        ratio = np.array([x["T_vis_map"][0] / x["T_vis_comp"][0] for _, x in rows
                          if np.isfinite(x["T_vis_comp"][0]) and x["T_vis_comp"][0] > 0])
        lab = lambda k: tex_escape(k[0] + (" " + k[1] if k[1] else ""))
        tex.append(f"{_tt(c)} & {lab(km)} & {tex_escape(_tvis_text(xm['T_vis_map'][0], dt))} & "
                   f"{lab(kc)} & {tex_escape(_tvis_text(xc['T_vis_comp'][0], dt))} & "
                   f"{np.median(ratio):.0f} ({ratio.min():.0f}–{ratio.max():.0f}) \\\\")
    tex += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return tex


def fig_criteria_map(crit, variables, layers, cases, criteria, tab):
    keys = row_keys(tab, variables, layers)
    labels = [v if ln == "" else f"{v} {ln}" for v, ln in keys]
    nc = len(criteria)
    fig, ax = plt.subplots(figsize=(FIGW, 0.3 * len(labels) + 1.6), constrained_layout=True)
    for i, (v, ln) in enumerate(keys):
        for j, c in enumerate(cases):
            x = crit.get((v, ln, c))
            for k, (name, _, _) in enumerate(criteria):
                st = x["status"][name][0] if x else "n/a"
                sym, _, col = STATUS[st]
                x0 = j - 0.45 + k * 0.9 / nc
                ax.add_patch(plt.Rectangle((x0, i - 0.42), 0.9 / nc, 0.84, facecolor=col, edgecolor="white"))
                ax.text(x0 + 0.45 / nc, i, sym, ha="center", va="center", fontsize=8)
    for j in range(len(cases)):
        for k, (name, _, _) in enumerate(criteria):
            ax.text(j - 0.45 + (k + 0.5) * 0.9 / nc, -0.75, str(k + 1), ha="center", va="bottom",
                    fontsize=7, color="0.3")
    ax.set_xlim(-0.5, len(cases) - 0.5)
    ax.set_ylim(len(labels) - 0.5, -1.0)
    ax.set_xticks(range(len(cases)))
    ax.set_xticklabels(cases, rotation=15, ha="right")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    byk = {}
    for k, (_, kind, T) in enumerate(criteria):
        byk.setdefault(kind, []).append(f"{k + 1}={T:g} d")
    legend = "    |    ".join(f"{kind}: " + ", ".join(v) for kind, v in byk.items())
    ax.set_title(f"{legend}\no pass   ~ borderline   x fail   - n/a", fontsize=8)
    return fig


# ---------------------------------------------------------------------------
# Markdown -> LaTeX
# ---------------------------------------------------------------------------

_TEX_ESC = {"\\": r"\textbackslash{}", "{": r"\{", "}": r"\}", "$": r"\$", "&": r"\&", "#": r"\#",
            "^": r"\^{}", "_": r"\_", "%": r"\%", "~": r"\textasciitilde{}"}


def tex_escape(s):
    return "".join(_TEX_ESC.get(ch, ch) for ch in s)


def inline(s):
    """Inline Markdown: $math$, `code`, **bold** (may contain math), @fig:/@tab:/@sec: references."""
    keep = []

    def _hold(m):
        tok = m.group(0)
        keep.append(tok if tok.startswith("$") else r"\texttt{" + tex_escape(tok[1:-1]) + "}")
        return f"\x00{len(keep) - 1}\x00"
    t = re.sub(r"\$[^$]+\$|`[^`]+`", _hold, s)
    t = tex_escape(t)
    t = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", t)
    t = re.sub(r"@((?:fig|tab|sec):[A-Za-z0-9\\_\-]+)",
               lambda m: r"\ref{" + m.group(1).replace(r"\_", "_") + "}", t)
    return re.sub("\x00(\\d+)\x00", lambda m: keep[int(m.group(1))], t)


def md_table(lines):
    cells = [[c.strip() for c in ln.strip().strip("|").split("|")] for ln in lines]
    body = [r for r in cells if not all(re.fullmatch(r":?-+:?", c) for c in r)]
    ncol = max(len(r) for r in body)
    tex = [r"\begin{center}\small", r"\begin{tabular}{" + "l" * ncol + "}", r"\toprule"]
    for k, r in enumerate(body):
        tex.append(" & ".join(inline(c) for c in r) + r" \\")
        if k == 0:
            tex.append(r"\midrule")
    tex += [r"\bottomrule", r"\end{tabular}", r"\end{center}"]
    return tex


def md_to_latex(text, auto):
    """Convert the simple Markdown of Re_report.md; ``auto`` maps marker names to LaTeX."""
    text = re.sub(r"<!--(?!\s*auto:).*?-->", "", text, flags=re.S)
    lines = text.split("\n")
    out, stack, para, i = [], [], [], 0     # stack: list of (indent, env)

    def flush_para():
        if para:
            out.append(inline(" ".join(para)))
            out.append("")
            para.clear()

    def close_lists(to_indent=-1):
        while stack and stack[-1][0] > to_indent:
            out.append(r"\end{" + stack.pop()[1] + "}")

    while i < len(lines):
        ln = lines[i]
        s = ln.strip()
        m_auto = re.fullmatch(r"<!--\s*auto:([\w-]+)\s*-->", s)
        m_head = re.match(r"(#{1,3})\s+(.*)", ln)
        m_item = re.match(r"(\s*)(-|\d+\.)\s+(.*)", ln)
        if not s:
            flush_para()
            i += 1
            continue
        if m_auto:
            flush_para()
            close_lists()
            name = m_auto.group(1)
            if name not in auto:
                print(f"warning: unknown auto block '{name}'", file=sys.stderr)
            out += auto.get(name, [])
            out.append("")
        elif m_head:
            flush_para()
            close_lists()
            lev, title = len(m_head.group(1)), m_head.group(2).strip()
            lab = re.search(r"\{#([\w:\-]+)\}\s*$", title)
            if lab:
                title = title[:lab.start()].strip()
            cmd = {1: r"\section", 2: r"\subsection", 3: r"\subsubsection*"}[lev]
            if title.startswith("付録") and lev == 1:
                out.append(r"\clearpage")
            out.append(cmd + "{" + inline(title) + "}" + (r"\label{" + lab.group(1) + "}" if lab else ""))
        elif s.startswith("|"):
            flush_para()
            close_lists()
            tbl = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                tbl.append(lines[i])
                i += 1
            out += md_table(tbl)
            continue
        elif s.startswith(">"):
            flush_para()
            close_lists()
            note = []
            while i < len(lines) and lines[i].strip().startswith(">"):
                note.append(lines[i].strip()[1:].strip())
                i += 1
            out.append(r"\notebox{" + inline(" ".join(note)) + "}")
            out.append("")
            continue
        elif m_item:
            flush_para()
            ind = len(m_item.group(1))
            env = "itemize" if m_item.group(2) == "-" else "enumerate"
            if not stack or ind > stack[-1][0]:
                stack.append((ind, env))
                out.append(r"\begin{" + env + "}")
            else:
                close_lists(ind)
                if not stack or stack[-1][0] < ind:
                    stack.append((ind, env))
                    out.append(r"\begin{" + env + "}")
            item = [m_item.group(3)]
            # continuation lines (indented more than the bullet, not a new item)
            while (i + 1 < len(lines) and lines[i + 1].strip()
                   and not re.match(r"\s*(-|\d+\.)\s+", lines[i + 1])
                   and len(lines[i + 1]) - len(lines[i + 1].lstrip()) > ind):
                i += 1
                item.append(lines[i].strip())
            out.append(r"\item " + inline(" ".join(item)))
        else:
            if stack and not ln.startswith(" "):
                close_lists()
            para.append(s)
        i += 1
    flush_para()
    close_lists()
    return "\n".join(out)


# ---------------------------------------------------------------------------
# Auto blocks (LaTeX)
# ---------------------------------------------------------------------------

def _tt(s):
    return r"\texttt{" + tex_escape(s) + "}"


def _num(x, fmt):
    return "—" if x is None or not np.isfinite(x) else format(x, fmt)


def _days(t):
    if t is None or not np.isfinite(t):
        return r"$\infty$"
    if t < 10:
        return f"{t:.1f}"
    if t < 1e5:
        return f"{t:,.0f}".replace(",", r"\,")
    e = int(np.floor(np.log10(t)))
    return f"${t / 10 ** e:.1f}\\times 10^{{{e}}}$"


def _processing(case):
    """(pre, post) of a case name like PrePR2p11_PostPR2p7."""
    def one(kind):
        m = re.search(kind + r"PR2p(\d+)", case)
        if m:
            return f"Reconstruction2 (Nr = {int(m.group(1)) + 1})"
        m = re.search(kind + r"([A-Za-z0-9]+)", case)
        return m.group(1) if m else "なし"
    return one("Pre"), one("Post")


def _figure(path, caption, label, width=r"\linewidth"):
    return [r"\begin{figure}[H]", r"\centering",
            r"\includegraphics[width=" + width + r",height=0.82\textheight,keepaspectratio]{" + path + "}",
            r"\caption{" + caption + "}", r"\label{" + label + "}", r"\end{figure}"]


def auto_meta(args, reports, cases, winfo):
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    gen = max((reports[c].get("generated", "") for c in cases), default="")
    rows = [("作成", now + " (summarize\\_Re.py)"),
            ("ケース数", f"{len(cases)} (" + ", ".join(_tt(c) for c in cases) + ")"),
            ("$R_e$ の出典", r"\texttt{report\_<case>.json} (make\_report.py、最新 " + tex_escape(gen) + ")"),
            ("SE と n の出典", _tt(os.path.basename(args.window_json or "なし"))
             + f" ({args.window_days:g} 日窓、window\\_Re.py)"),
            ("criterion", f"$\\kappa = {args.kappa:g}$、n は窓が {MIN_WINDOWS} 個以上とれる窓長から")]
    tex = [r"\begin{center}\small", r"\begin{tabular}{@{}p{0.2\linewidth}p{0.74\linewidth}@{}}", r"\toprule"]
    tex += [f"{a} & {b} \\\\" for a, b in rows]
    tex += [r"\bottomrule", r"\end{tabular}", r"\end{center}"]
    return tex


def auto_cases(reports, cases, msgs):
    tex = [r"\begin{table}[H]", r"\centering\footnotesize",
           r"\caption{比較したケース。サンプルは 2.5 日ごとの history 出力。}", r"\label{tab:cases}",
           r"\begin{tabular}{@{}lllrl@{}}", r"\toprule",
           r"ケース & 前処理 (MP 入力) & 後処理 (MP 加熱) & サンプル & 解析期間 \\", r"\midrule"]
    for c in cases:
        d = reports[c]
        pre, post = _processing(c)
        per = d.get("period", ["?", "?"])
        p0, p1 = (str(per[0])[:10], str(per[1])[:10])
        tex.append(f"{_tt(c)} & {tex_escape(pre)} & {tex_escape(post)} & {d.get('nsamples', 0)} & "
                   f"{p0} – {p1} \\\\")
    tex += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    if msgs:
        tex += [r"\begin{warnbox}", r"\textbf{比較の前に確認すること}", r"\begin{itemize}"]
        tex += [r"\item " + inline(m) for m in msgs]
        tex += [r"\end{itemize}", r"\end{warnbox}"]
    return tex


def auto_criteria_def(args, criteria):
    tex = [
        r"\paragraph{criterion の定義}",
        r"criterion は「種類」と「平均期間 $T$」の組で、$T$ 日平均のその種類の図で要素ロックが"
        r"\textbf{見えないこと}を要求する。種類は次の 2 つ。",
        r"\begin{itemize}",
        r"\item \textbf{地図}: 時間平均した水平分布。要素 1 つずつについて見えるかを判定する。",
        r"\item \textbf{コンポジット}: 要素コンポジット (全要素の平均)。",
        r"\end{itemize}",
        r"\paragraph{見えるかどうかの判定}",
        r"$T$ 日平均の係数のノイズは $\sigma(a)/\sqrt{nT}$ である。信号 $\langle a\rangle$ がその $\kappa$ 倍を超えると"
        r"「見える」とする。",
        r"\begin{equation}",
        r"\mathrm{S/N} = |R_e|\sqrt{nT} > \kappa"
        r"\quad\Longleftrightarrow\quad |R_e| > R_\mathrm{crit} = \frac{\kappa}{\sqrt{nT}}"
        r"\quad\Longleftrightarrow\quad T > T_\mathrm{vis} = \frac{\kappa^2}{R_e^2\, n}",
        r"\end{equation}",
        r"\begin{itemize}",
        r"\item $T_\mathrm{vis}$ は、ロックが見え出す平均期間である。",
        r"\item $n$ は 1 日あたりの実効的な独立サンプル数で、window\_Re.py が窓平均のばらつきから測る。"
        r"地図は要素ごと ($n_\mathrm{map}$)、コンポジットは要素平均 ($n_\mathrm{comp}$) について測る。",
        f"\\item $n$ には、窓が {MIN_WINDOWS} 個以上とれる窓長のうち、criterion の $T$ に (対数で) 最も近いものの値を使う"
        r"(30 日の criterion は 30 日窓の値で、30 日平均のノイズをそのまま測っている)。"
        r"$T_\mathrm{vis}$ と表の $n$ は、最長の窓長の値。",
        f"\\item 本資料では $\\kappa = {args.kappa:g}$ とした。",
        r"\item 判定は $|R_e| \pm 2\,\mathrm{SE}$ で行う。",
        r"\begin{itemize}",
        r"\item 合格 (○): $|R_e| + 2\,\mathrm{SE} < R_\mathrm{crit}$",
        r"\item 不合格 (×): $|R_e| - 2\,\mathrm{SE} > R_\mathrm{crit}$",
        r"\item 境界 (△): それ以外",
        r"\end{itemize}",
        r"\end{itemize}",
        r"\begin{table}[H]", r"\centering\footnotesize",
        r"\caption{用途別の criterion。図の種類 (地図・コンポジット) と平均期間 $T$ の組。"
        r"地図の criterion を実用、コンポジットの短・中期間を中間、長期間を厳しい criterion と呼ぶ。}",
        r"\label{tab:criteria}",
        r"\begin{tabular}{@{}llrl@{}}", r"\toprule",
        r"criterion & 種類 & $T$ [日] & 想定する用途 \\", r"\midrule"]
    use = {"practical-30": "月平均の水平分布の図で見えない",
           "practical-100": "季節平均程度の水平分布の図で見えない",
           "practical-1000": "長期平均の水平分布の図で見えない",
           "intermediate-30": "月平均の要素コンポジットで見えない",
           "intermediate-100": "季節平均程度の要素コンポジットで見えない",
           "strict": "長期平均の要素コンポジットで見えない"}
    for name, kind, T in criteria:
        tex.append(f"{CRITERIA_JA.get(name, tex_escape(name))} ({tex_escape(name)}) & {KIND_JA[kind]} & "
                   f"{T:g} & {use.get(name, '')} \\\\")
    tex += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return tex


def auto_criteria_table(crit, cases, criteria, var, layer):
    lab = var if not layer else f"{var} {layer}"
    key = (var + ("_" + layer if layer else "")).replace("_", "-")
    tex = [r"\begin{table}[H]", r"\centering\footnotesize",
           r"\caption{" + tex_escape(lab) + r" の $R_e$ と、ロックが見え出す平均期間 $T_\mathrm{vis}$。"
           r"括弧内は $R_e \mp 2\,\mathrm{SE}$ による範囲。$n$ は 1 日あたりの実効サンプル数、"
           r"「窓」は $n$ を求めた窓長 [日]。}",
           r"\label{tab:tvis-" + key + "}",
           r"\begin{tabular}{@{}lrrrrrll@{}}", r"\toprule",
           r"ケース & $R_e$ & SE & $n_\mathrm{map}$ & $n_\mathrm{comp}$ & 窓 & "
           r"$T_\mathrm{vis}$ 地図 [日] & $T_\mathrm{vis}$ コンポジット [日] \\", r"\midrule"]
    for c in cases:
        x = crit.get((var, layer, c))
        if x is None:
            continue
        tm, tcm = x["T_vis_map"], x["T_vis_comp"]

        def rng(t):
            if not np.isfinite(t[0]):
                return "—"
            return f"{_days(t[0])} \\,{{\\footnotesize [{_days(t[1])}–{_days(t[2])}]}}"
        tex.append(f"{_tt(c)} & ${x['R']:+.3f}$ & {_num(x['SE'], '.3f')} & {_num(x['n_map'], '.3f')} & "
                   f"{_num(x['n_comp'], '.1f')} & {_num(x['nwin'], '.0f')} & {rng(tm)} & {rng(tcm)} \\\\")
    tex += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]

    groups = []                               # [(kind, [criteria of that kind])] in order
    for cr in criteria:
        if groups and groups[-1][0] == cr[1]:
            groups[-1][1].append(cr)
        else:
            groups.append((cr[1], [cr]))
    grp_name = {"map": "時間平均の地図 (実用)", "composite": "要素コンポジット (中間・厳しい)"}

    def header():
        h1, rules, col = ["ケース"], [], 2
        for kind, lst in groups:
            h1.append(r"\multicolumn{" + str(len(lst)) + "}{c}{" + grp_name.get(kind, kind) + "}")
            rules.append(r"\cmidrule(lr){" + f"{col}-{col + len(lst) - 1}" + "}")
            col += len(lst)
        h2 = [""] + [f"$T$ = {T:g} 日" for _, _, T in criteria]
        return [" & ".join(h1) + r" \\", "".join(rules), " & ".join(h2) + r" \\", r"\midrule"]

    def table(caption, label, cell):
        t = [r"\begin{table}[H]", r"\centering\footnotesize", r"\caption{" + caption + "}",
             r"\label{" + label + "}", r"\begin{tabular}{@{}l" + "c" * len(criteria) + "@{}}", r"\toprule"]
        t += header()
        for c in cases:
            x = crit.get((var, layer, c))
            if x is not None:
                t.append(f"{_tt(c)} & " + " & ".join(cell(x["status"][name]) for name, _, _ in criteria) + r" \\")
        return t + [r"\bottomrule", r"\end{tabular}", r"\end{table}"]

    def cell_status(stv):
        st, rc, sn = stv
        body = STATUS[st][1] if not np.isfinite(sn) else f"{STATUS[st][1]} {sn:.2f}"
        return f"\\cellcolor[HTML]{{{STATUS[st][2].lstrip('#')}}}{body}"

    def cell_rcrit(stv):
        return "—" if not np.isfinite(stv[1]) else f"{stv[1]:.3f}"

    tex += table(tex_escape(lab) + r" の criterion の判定。各欄は判定と S/N ($= |R_e|\sqrt{nT}$)。"
                 r"○ 合格、△ 境界、× 不合格、— 判定なし。", "tab:status-" + key, cell_status)
    tex += table(tex_escape(lab) + r" で、各 criterion を満たすのに必要な $|R_e|$ の上限 "
                 r"$R_\mathrm{crit} = \kappa/\sqrt{nT}$。$n$ がケースごとに違うので、値もケースごとに違う。",
                 "tab:rcrit-" + key, cell_rcrit)
    return tex


# ---------------------------------------------------------------------------
# Document
# ---------------------------------------------------------------------------

PREAMBLE = r"""\documentclass[10pt,a4paper]{article}
\usepackage{fontspec}
\setmainfont{VL PGothic}[AutoFakeBold=2.0]
\setsansfont{VL PGothic}[AutoFakeBold=2.0]
\setmonofont{VL Gothic}[Scale=0.92]
\XeTeXlinebreaklocale "ja"
\XeTeXlinebreakskip=0pt plus 0.1em
\usepackage[a4paper,top=22mm,bottom=24mm,left=21mm,right=21mm,headheight=14pt]{geometry}
\usepackage{amsmath,amssymb}
\usepackage[table]{xcolor}
\usepackage{graphicx,booktabs,array,float}
\usepackage[font=small,labelfont=bf,labelsep=quad,skip=4pt]{caption}
\usepackage{enumitem}
\usepackage{fancyhdr}
\usepackage[hidelinks,unicode]{hyperref}
\linespread{1.3}
\setlength{\parindent}{0pt}
\setlength{\parskip}{0.55em}
\setlist{itemsep=0.15em,topsep=0.3em,parsep=0pt,leftmargin=1.6em}
\setlist[itemize,1]{label=\textbullet}
\setlist[itemize,2]{label=--}
\definecolor{headc}{HTML}{1F4E79}
\definecolor{notec}{HTML}{EEF3F8}
\definecolor{warnc}{HTML}{FDF1E7}
\renewcommand{\figurename}{図}
\renewcommand{\tablename}{表}
\renewcommand{\contentsname}{目次}
\makeatletter
\renewcommand\section{\@startsection{section}{1}{\z@}{-2.2ex \@plus -.5ex \@minus -.2ex}{1.2ex}%
  {\normalfont\Large\bfseries\color{headc}}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}{-1.8ex \@plus -.4ex \@minus -.2ex}{0.8ex}%
  {\normalfont\large\bfseries\color{headc}}}
\renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}{-1.4ex \@plus -.3ex}{0.5ex}%
  {\normalfont\normalsize\bfseries}}
\renewcommand\paragraph{\@startsection{paragraph}{4}{\z@}{1.0ex \@plus .3ex}{-1em}%
  {\normalfont\normalsize\bfseries\color{headc}}}
\makeatother
\newcommand{\notebox}[1]{\par\noindent\colorbox{notec}{\parbox{\dimexpr\linewidth-2\fboxsep}{\small #1}}\par}
\newenvironment{warnbox}{\par\noindent\begin{lrbox}{0}\begin{minipage}{\dimexpr\linewidth-2\fboxsep}\small}%
  {\end{minipage}\end{lrbox}\colorbox{warnc}{\usebox{0}}\par}
\renewcommand{\arraystretch}{1.25}
\setlength{\tabcolsep}{5pt}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small\color{gray} 要素ロック比 $R_e$ のまとめ (Eh8Ez4P7)}
\fancyhead[R]{\small\color{gray} \leftmark}
\fancyfoot[C]{\small \thepage}
\renewcommand{\headrulewidth}{0.4pt}
\renewcommand{\sectionmark}[1]{\markboth{\thesection\quad #1}{}}
\setcounter{tocdepth}{2}
\begin{document}
\begin{center}
{\LARGE\bfseries\color{headc} 要素ロック構造に対する前処理・後処理の効果}\par\vspace{0.6em}
{\large 要素ロック比 $R_e$・不確かさ・用途別 criterion のまとめ — Eh8Ez4P7 理想化湿潤対流実験}\par\vspace{0.4em}
{\small DATE}
\end{center}
\vspace{0.5em}
{\small\tableofcontents}
\vspace{0.5em}
"""


def build_latex(body, date):
    return PREAMBLE.replace("DATE", date) + body + "\n\\end{document}\n"


def run_xelatex(tex_path):
    cwd = os.path.dirname(tex_path)
    name = os.path.basename(tex_path)
    for _ in range(2):        # twice for the table of contents and references
        r = subprocess.run(["xelatex", "-interaction=nonstopmode", "-halt-on-error", name],
                           cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if r.returncode != 0:
            return False, r.stdout[-3000:]
    log = open(os.path.join(cwd, os.path.splitext(name)[0] + ".log"), errors="replace").read()
    missing = sorted(set(re.findall(r"Missing character: There is no (.) in font", log)))
    return True, ("missing characters: " + " ".join(missing)) if missing else ""


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_csv(rows, path, se=None):
    se = se or {}
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["variable", "layer", "case", "mode", "R_e", "SE", "sigma", "units"])
        for v, ln, c, m, r, s, u in rows:
            e = se.get((v, ln, c, m), np.nan)
            w.writerow([v, ln, c, f"({m[0]},{m[1]})", f"{r:.6g}", f"{e:.6g}", f"{s:.6g}", u])


def write_criteria_csv(crit, criteria, path):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["variable", "layer", "case", "R_e", "SE", "n_map_per_day", "n_comp_per_day",
                    "window_days", "T_vis_map", "T_vis_comp"]
                   + [f"{name}_{x}" for name, _, _ in criteria for x in ("status", "R_crit", "SN")])
        for (v, ln, c), x in sorted(crit.items()):
            row = [v, ln, c, f"{x['R']:.4g}", f"{x['SE']:.3g}", f"{x['n_map']:.4g}", f"{x['n_comp']:.4g}",
                   f"{x['nwin']:g}", f"{x['T_vis_map'][0]:.4g}", f"{x['T_vis_comp'][0]:.4g}"]
            for name, _, _ in criteria:
                st, rc, sn = x["status"][name]
                row += [st, f"{rc:.4g}", f"{sn:.3g}"]
            w.writerow(row)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__.split("Usage (in Eh8Ez4P7/analysis/)::")[1])
    p.add_argument("--json", nargs="+", default=None,
                   help="report JSON files (default: report_*.json in this directory)")
    p.add_argument("--cases", nargs="+", default=None,
                   help=f"cases to show, in this order (default: all, '{REFERENCE_CASE}' first)")
    p.add_argument("-o", "--output", default=str(HERE / "Re_summary.pdf"), help="output PDF")
    p.add_argument("--document", default=str(HERE / "Re_report.md"),
                   help="Markdown text of the report (default Re_report.md)")
    p.add_argument("--min-samples", type=int, default=100,
                   help="warn about cases with fewer samples (default 100)")
    p.add_argument("--window-json", default=str(HERE / "window_Re.json"),
                   help="SE and n from window_Re.py (default window_Re.json; '' to skip)")
    p.add_argument("--window-days", type=float, default=30.0, help="window length of the SE [day]")
    p.add_argument("--kappa", type=float, default=KAPPA,
                   help=f"visibility threshold of |R_e| sqrt(n T) (default {KAPPA:g})")
    p.add_argument("--criterion", action="append", default=None, metavar="NAME=KIND:DAYS",
                   help="R_e criterion (repeatable; replaces the defaults), KIND = map or composite, "
                        "e.g. practical=map:1000")
    p.add_argument("--criteria-vars", nargs="+", default=["PREC"],
                   help="variables (VAR or VAR_LAYER) with criterion tables and a T_vis figure")
    args = p.parse_args(argv)

    paths = args.json or sorted(glob.glob(str(HERE / "report_*.json")))
    if not paths:
        p.error("no report_*.json found; run make_report.py first")
    reports = load_reports(paths)
    cases = order_cases(reports, args.cases)
    rows, variables, layers = collect(reports, cases)
    tab = lookup(rows)
    se, winfo = load_window_se(args.window_json, args.window_days)
    msgs = warnings_for(reports, cases, args.min_samples, winfo if args.window_json else None,
                        args.window_days)
    criteria = CRITERIA
    if args.criterion:
        criteria = []
        for cdef in args.criterion:
            try:
                name, rest = cdef.split("=")
                kind, days = rest.split(":")
                assert kind in ("map", "composite")
                criteria.append((name, kind, float(days)))
            except (ValueError, AssertionError):
                p.error(f"invalid --criterion '{cdef}' (use NAME=map:DAYS or NAME=composite:DAYS)")
    crit = evaluate_criteria(rows, se, load_window_n(args.window_json), criteria, args.kappa)

    # figures
    stem = os.path.splitext(os.path.abspath(args.output))[0]
    build = stem + "_build"
    os.makedirs(build, exist_ok=True)
    figs = {}

    def save(name, fig):
        path = os.path.join(build, f"fig_{name}.pdf")
        fig.savefig(path)
        plt.close(fig)
        figs[name] = path
        return os.path.basename(path)

    save("re_heatmap", fig_heatmap(rows, variables, layers, cases, se))
    save("criteria_map", fig_criteria_map(crit, variables, layers, cases, criteria, tab))
    dt = sample_dt(args.window_json)
    save("tvis_map", fig_tvis_map(crit, variables, layers, cases, criteria, tab, dt))
    cvars = [_split_key(cv) for cv in args.criteria_vars]
    for v, ln in cvars:
        save(f"tvis_{v}{'_' + ln if ln else ''}", fig_tvis(crit, cases, criteria, v, ln))
    for v in variables:
        save(f"var_{v}", fig_variable(rows, v, layers, cases, reports, se))

    # auto blocks
    auto = {
        "meta": auto_meta(args, reports, cases, winfo),
        "cases": auto_cases(reports, cases, msgs),
        "re-heatmap": _figure("fig_re_heatmap.pdf",
                              r"全変数・全ケースの $R_e$ (モード (2,0) と (0,2) の平均)。赤は端寄り、青は中央寄り。"
                              r"各欄に $R_e$ と SE を示す。括弧付きの値は $|R_e| < 2\,\mathrm{SE}$ で、0 と区別できない。"
                              r"PREC は降水、ほかは層平均 (BL 0–1.2 km、FT 1.2–11 km、UL 11–15 km)。",
                              "fig:re-heatmap"),
        "criteria-def": auto_criteria_def(args, criteria),
        "criteria-table": sum((auto_criteria_table(crit, cases, criteria, v, ln) for v, ln in cvars), []),
        "criteria-tvis": sum((_figure(f"fig_tvis_{v}{'_' + ln if ln else ''}.pdf",
                                      tex_escape(v + (' ' + ln if ln else ''))
                                      + r" の $T_\mathrm{vis}$ (■ 地図、● コンポジット)。誤差棒は "
                                      r"$R_e \mp 2\,\mathrm{SE}$ による範囲。横線は各 criterion の平均期間 "
                                      r"(破線: 地図、点線: コンポジット)。マーカーが、同じ種類の線より上にあれば合格。",
                                      f"fig:tvis-{v}{'-' + ln if ln else ''}", r"0.92\linewidth")
                              for v, ln in cvars), []),
        "criteria-map": _figure("fig_criteria_map.pdf",
                                r"全変数・全ケースの criterion の判定。各ケースの " + str(len(criteria)) + r" つの小列は、左から "
                                + "、".join(CRITERIA_JA.get(n, tex_escape(n)) for n, _, _ in criteria)
                                + r"。緑 (o) 合格、黄 (\textasciitilde) 境界、赤 (x) 不合格、灰 (-) 判定なし。",
                                "fig:criteria-map"),
        "tvis-map": _figure("fig_tvis_map.pdf",
                            r"全変数・全ケースの、要素ロックが見え出す平均期間 $T_\mathrm{vis}$ [日] "
                            r"(左: 時間平均の地図、右: 要素コンポジット)。色の境界は、出力間隔 ("
                            + f"{dt:g}" + r" 日) と criterion の平均期間 ("
                            + "・".join(f"{T:g}" for T in sorted({T for _, _, T in criteria}))
                            + r" 日)。赤いほど早く見える。括弧付きは $|R_e| < 2\,\mathrm{SE}$ で、値の不確かさが大きい。"
                            r"解析期間 (約 1000 日) より長い $T_\mathrm{vis}$ は外挿。",
                            "fig:tvis-map"),
        "tvis-summary": auto_tvis_summary(crit, variables, layers, cases, tab, dt),
        "variables": sum((_figure(f"fig_var_{v}.pdf",
                                  tex_escape(v) + r" の $R_e$ (上、誤差棒 ±1 SE) と $\sigma(a)$ (下)。"
                                  r"横軸の n はサンプル数。", f"fig:var-{v}")
                          for v in variables), []),
    }
    with open(args.document, encoding="utf-8") as f:
        body = md_to_latex(f.read(), auto)
    date = datetime.datetime.now().strftime("%Y-%m-%d")
    tex_path = os.path.join(build, "Re_summary.tex")
    with open(tex_path, "w", encoding="utf-8") as f:
        f.write(build_latex(body, date))

    # typeset
    ok, info = (False, "xelatex not found")
    if shutil.which("xelatex"):
        ok, info = run_xelatex(tex_path)
    if ok:
        shutil.copyfile(os.path.splitext(tex_path)[0] + ".pdf", args.output)
        if info:
            print(f"  warning: {info}")
    else:
        print(f"  LaTeX failed ({tex_path}); writing the figures only:\n{info}", file=sys.stderr)
        with PdfPages(args.output) as pdf:
            for fig in [fig_heatmap(rows, variables, layers, cases, se),
                        fig_criteria_map(crit, variables, layers, cases, criteria, tab)] + \
                       [fig_tvis(crit, cases, criteria, v, ln) for v, ln in cvars] + \
                       [fig_variable(rows, v, layers, cases, reports, se) for v in variables]:
                pdf.savefig(fig)
                plt.close(fig)

    csv_path = stem + ".csv"
    write_csv(rows, csv_path, se)
    crit_csv = stem + "_criteria.csv"
    write_criteria_csv(crit, criteria, crit_csv)
    print(f"{len(cases)} cases ({', '.join(cases)}), {len(variables)} variables"
          + ("" if ok else "  [figures only]"))
    for m in msgs:
        print("  check: " + re.sub(r"@(tab|fig|sec):\S+?(?=\))", "1", m).replace("$", ""))
    print(f"-> {args.output}\n-> {csv_path}\n-> {crit_csv}\n   (LaTeX source and figures: {build}/)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
