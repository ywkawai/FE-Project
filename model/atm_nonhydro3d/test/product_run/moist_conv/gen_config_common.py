"""Common part of gen_config.py for the moist_conv experiments (independent of resolution / polynomial order).

Each resolution directory (e.g. Eh8Ez4P7/) has a small gen_config.py that sets the
resolution-dependent values in a Config object and calls main(cfg):

    import sys
    from pathlib import Path
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from gen_config_common import Config, main

    cfg = Config(here=Path(__file__).resolve().parent, NPRC=8, ...)
    if __name__ == "__main__":
        main(cfg, __doc__)

Files in the resolution directory (cfg.here):

    run_common_template.conf     common part of run.conf
    run_<CASE>_template.conf     case-specific part (namelist groups)
    job_run_template.sh          job script template (Odyssey)
    phy_preoperation_*NN.conf    configuration of variable containers (copied when needed)
    run_spinup2/                 spin-up; its last restart is read by run1

    <CASE>/run<N>/run.conf       generated
    <CASE>/run<N>/run.sh         generated (job script)
    <CASE>/run<N>/phy_preoperation_*NN.conf   copied when the case uses container NN >= 2
    <CASE>/submit_chain.sh       generated: submits run<first>..run<last> as a step job

How run.conf is built:
  1. A namelist group in run_<CASE>_template.conf replaces the group of the same name in the
     common template. The whole group is replaced, so write all variables you need in it.
     Groups that appear more than once (HISTORY_ITEM, MONITOR_ITEM) are never replaced.
  2. Other groups and comments in the case template are inserted at the line containing
     @@CASE_SPECIFIC@@ in the common template (appended at the end if there is no such line).
  3. Placeholders written as a KEY in curly braces are replaced with the values from build_params().
     A brace preceded by '$' (bash ${VAR}) is left alone. An unknown KEY is an error.

Usage (in the resolution directory):
    python3 gen_config.py asis                 # generate one case (asis/run1..run4)
    python3 gen_config.py asis case1           # generate several cases at once
    python3 gen_config.py --dry-run PostPR2p7  # show what would be done
    python3 gen_config.py --force case1        # overwrite run.conf / run.sh that differ
    python3 gen_config.py --list               # list the available case templates
"""

import argparse
import os
import re
import sys
from dataclasses import dataclass, field
from datetime import date, timedelta
from pathlib import Path

CASE_MARKER = "@@CASE_SPECIFIC@@"


# ============================================================================
# Settings (defaults = values common to all resolutions)
# ============================================================================
@dataclass
class Config:
    here: Path  # resolution directory (where gen_config.py and the templates are)

    # --- Resolution-dependent: must be given by each gen_config.py
    START_DATE: date  # start of the main run = end of run_spinup2
    NPRC: int  # number of MPI processes
    ELAPSE: str  # elapse time limit of one job
    JOB_PREFIX: str  # prefix of the job name

    NPRC_PER_NODE: int = 4
    NTHREAD: int = 12

    # --- Templates (in `here`)
    COMMON_TEMPLATE: str = "run_common_template.conf"
    CASE_TEMPLATE: str = "run_{case}_template.conf"
    JOB_TEMPLATE: str = "job_run_template.sh"

    # --- Time integration
    TINTEG_DAYS: int = 1000  # total length of the main run
    TINTEG_PER_RUN: int = 250  # length of one job
    TIME_DT_RESTART_DAY: float = 125  # intermediate restarts within a job (for recovery)
    HISTORY_TINTERVAL_DAY: float = 2.5

    # --- Spin-up whose last restart is read by run1 (relative to `here`)
    SPINUP_DIR: str = "run_spinup2"

    # --- T_ref of the Newtonian relaxation (relative to `here`)
    #     Same as run_spinup2 (horizontal mean of T in the initial state = step 1 of run_spinup1/history),
    #     so that T_ref does not change between the spin-up and the main run.
    #     The Fugaku template used '../run_spinup2/history' instead.
    TREF_FILEBASE: str = "run_spinup1/history"
    TREF_STEP: int = 1

    # --- Job (Odyssey)
    GROUP: str = "gz06"
    RSCGRP: str = "regular-o"
    BIN: str = "../scale-dg"  # relative to `here`

    # --- Additional files copied into every run directory of a case (relative to `here`).
    #     phy_preoperation_*NN.conf files are detected automatically from run.conf.
    CASE_EXTRA_FILES: dict = field(default_factory=dict)

    @property
    def nrun(self):
        return (self.TINTEG_DAYS + self.TINTEG_PER_RUN - 1) // self.TINTEG_PER_RUN


# ============================================================================
# Namelist helpers
# ============================================================================
def f90(x):
    """Convert a Python value to a Fortran namelist literal."""
    if isinstance(x, bool):
        return ".true." if x else ".false."
    if isinstance(x, float):
        return f"{x:.8e}".replace("e", "D")
    return str(x)


def strip_comment(line):
    """Remove a Fortran '!' comment (quote-aware)."""
    quote = None
    for i, c in enumerate(line):
        if quote:
            if c == quote:
                quote = None
        elif c in "'\"":
            quote = c
        elif c == "!":
            return line[:i]
    return line


NML_BEGIN = re.compile(r"^\s*&(\w+)")


def split_items(text, fname):
    """Split a conf text into items: {'kind': 'nml'|'text', 'name': str|None, 'lines': [...]}."""
    lines = [l if l.endswith("\n") else l + "\n" for l in text.splitlines(keepends=True)]
    items = []
    i = 0
    while i < len(lines):
        m = NML_BEGIN.match(lines[i])
        if not m:
            items.append({"kind": "text", "name": None, "lines": [lines[i]]})
            i += 1
            continue
        name = m.group(1).upper()
        blk = [lines[i]]
        first = strip_comment(lines[i]).rstrip()
        if not first.endswith("/"):  # multi-line group
            i += 1
            while True:
                if i >= len(lines):
                    sys.exit(f"ERROR: {fname}: namelist group &{name} is not terminated by '/'")
                blk.append(lines[i])
                if strip_comment(lines[i]).strip() == "/":
                    break
                i += 1
        items.append({"kind": "nml", "name": name, "lines": blk})
        i += 1
    return items


def attach_comments(items):
    """Attach comment lines directly above a namelist group to that group (case template)."""
    out = []
    pending = []
    for it in items:
        if it["kind"] == "text" and it["lines"][0].lstrip().startswith("#"):
            pending.append(it)
            continue
        if it["kind"] == "nml":
            it = dict(it, lines=sum((p["lines"] for p in pending), []) + it["lines"])
        else:
            out.extend(pending)
        pending = []
        out.append(it)
    out.extend(pending)
    return out


def merge(common_text, case_text, common_fname, case_fname):
    """Merge the case template into the common template (see the module docstring)."""
    common = split_items(common_text, common_fname)
    case = attach_comments(split_items(case_text, case_fname))

    count = {}
    for it in common:
        if it["kind"] == "nml":
            count[it["name"]] = count.get(it["name"], 0) + 1

    replaced = {}
    inserted = []
    for it in case:
        if it["kind"] == "nml" and count.get(it["name"]) == 1:
            if it["name"] in replaced:
                sys.exit(f"ERROR: {case_fname}: &{it['name']} appears twice")
            replaced[it["name"]] = it["lines"]
        else:
            inserted.extend(it["lines"])

    out = []
    has_marker = False
    for it in common:
        if it["kind"] == "nml" and it["name"] in replaced:
            out.extend(replaced[it["name"]])
        elif it["kind"] == "text" and CASE_MARKER in it["lines"][0]:
            has_marker = True
            out.extend(inserted)
        else:
            out.extend(it["lines"])
    if not has_marker:
        out.append("\n")
        out.extend(inserted)
    return "".join(out), sorted(replaced)


PLACEHOLDER = re.compile(r"(?<![$])\{([A-Za-z_]\w*)\}")


def render(text, params, fname):
    missing = sorted({k for k in PLACEHOLDER.findall(text) if k not in params})
    if missing:
        sys.exit(f"ERROR: {fname}: no value for placeholder(s): {', '.join(missing)}")
    return PLACEHOLDER.sub(lambda m: params[m.group(1)], text)


def nml_value(conf, group, var):
    """Return the value of `var` in namelist `group` of conf (str without quotes) or None."""
    for it in split_items(conf, "run.conf"):
        if it["kind"] != "nml" or it["name"] != group.upper():
            continue
        for line in it["lines"]:
            m = re.match(rf"^\s*{var}\s*=\s*(.+)$", strip_comment(line), re.IGNORECASE)
            if m:
                return m.group(1).strip().rstrip(",").strip().strip("'\"")
    return None


def required_container_files(conf):
    """phy_preoperation files that scale-dg reads for the variable containers used in conf."""
    base = nml_value(conf, "PARAM_ATMOS_VARS", "PHY_PREPROC_FILE_BASENAME") or "phy_preoperation"
    ids = set()
    for group, var in [("PARAM_ATMOS_PHY_MP", "atm_var_container_typeid"),
                       ("PARAM_USER", "BL_VARS_CONTAINER_ID")]:
        v = nml_value(conf, group, var)
        if v is not None and int(v) >= 2:
            ids.add(int(v))
    return [f"{base}{i:02d}.conf" for i in sorted(ids)]


# ============================================================================
# Run parameters
# ============================================================================
def restart_name(d):
    return f"restart_{d.year:04d}{d.month:02d}{d.day:02d}-000000.000"


def rel(cfg, target, run_dir):
    return os.path.relpath(cfg.here / target, run_dir)


def build_params(cfg, case, irun, run_dir):
    start = cfg.START_DATE + timedelta(days=(irun - 1) * cfg.TINTEG_PER_RUN)
    duration = min(cfg.TINTEG_PER_RUN, cfg.TINTEG_DAYS - (irun - 1) * cfg.TINTEG_PER_RUN)
    prev_dir = cfg.SPINUP_DIR if irun == 1 else f"{case}/run{irun - 1}"
    p = {
        "CASE": case,
        "IRUN": irun,
        "START_YEAR": start.year,
        "START_MONTH": start.month,
        "START_DAY": start.day,
        "TIME_DURATION_DAY": float(duration),
        "TIME_DT_RESTART_DAY": float(cfg.TIME_DT_RESTART_DAY),
        "HISTORY_TINTERVAL_DAY": float(cfg.HISTORY_TINTERVAL_DAY),
        "RESTART_IN_BASENAME": f"{rel(cfg, prev_dir, run_dir)}/{restart_name(start)}",
        "TREF_FILEBASE": rel(cfg, cfg.TREF_FILEBASE, run_dir),
        "TREF_STEP": cfg.TREF_STEP,
        # job script
        "GROUP": cfg.GROUP,
        "RSCGRP": cfg.RSCGRP,
        "NPRC": cfg.NPRC,
        "NPRC_PER_NODE": cfg.NPRC_PER_NODE,
        "NODE": (cfg.NPRC + cfg.NPRC_PER_NODE - 1) // cfg.NPRC_PER_NODE,
        "NTHREAD": cfg.NTHREAD,
        "ELAPSE": cfg.ELAPSE,
        "JOB_NAME": f"{cfg.JOB_PREFIX}_{case}_r{irun}",
        "RUN_DIR": str(run_dir),
        "BIN": rel(cfg, cfg.BIN, run_dir),
        "CONF": "run.conf",
    }
    return {k: f90(v) for k, v in p.items()}, start, duration, prev_dir


# ============================================================================
# Output
# ============================================================================
class Writer:
    """Writes files, honoring --force / --dry-run, and remembers files that were kept."""

    def __init__(self, force, dry_run):
        self.force = force
        self.dry_run = dry_run
        self.kept = []  # existing files that differ from the generated text and were not overwritten

    def write(self, path, text, mode=None):
        """Write text to path. Returns a status string."""
        if path.exists():
            if path.read_text() == text:
                return "unchanged"
            if not self.force:
                self.kept.append(path)
                return "DIFFERS (kept; use --force to overwrite)"
            status = "overwritten"
        else:
            status = "created"
        if not self.dry_run:
            path.write_text(text)
            if mode is not None:
                path.chmod(mode)
        return status


def submit_chain_script(case, nrun):
    return f"""#! /bin/bash
# Submit {case}/run<first>..run<last> as one step job (each step starts after the previous one ends;
# the remaining steps are cancelled if a step ends with a non-zero exit code).
# Generated by gen_config.py.
#   usage: ./submit_chain.sh [first_run (default 1)] [last_run (default {nrun})]
set -e
cd "$(dirname "$0")"
first=${{1:-1}}
last=${{2:-{nrun}}}
jid=""
for irun in $(seq ${{first}} ${{last}}); do
  if [ -z "${{jid}}" ]; then
    out=$(pjsub -z jid --step run${{irun}}/run.sh)
  else
    out=$(pjsub -z jid --step --sparam "jid=${{jid}},sd=ec!=0:after" run${{irun}}/run.sh)
  fi
  echo "run${{irun}}: ${{out}}"
  jid=${{out%%_*}}
done
"""


def generate_case(cfg, case, common_text, job_text, writer, templates):
    case_fname = cfg.CASE_TEMPLATE.format(case=case)
    if not (cfg.here / case_fname).exists():
        sys.exit(f"ERROR: {case_fname} not found. available cases: {' '.join(templates)}")
    merged, replaced = merge(common_text, (cfg.here / case_fname).read_text(),
                             cfg.COMMON_TEMPLATE, case_fname)
    print(f"=== case {case}: {case_fname} replaces "
          f"{', '.join('&' + r for r in replaced) or '(nothing)'}")

    case_dir = cfg.here / case
    for irun in range(1, cfg.nrun + 1):
        run_dir = case_dir / f"run{irun}"
        params, start, duration, prev_dir = build_params(cfg, case, irun, run_dir)
        header = (f"#--- Generated by gen_config.py from {cfg.COMMON_TEMPLATE} + {case_fname}.\n"
                  f"#    Edit the templates and regenerate instead of editing this file.\n")
        conf = header + render(merged, params, case_fname)
        job = render(job_text, params, cfg.JOB_TEMPLATE)

        if not writer.dry_run:
            run_dir.mkdir(parents=True, exist_ok=True)
        st_conf = writer.write(run_dir / "run.conf", conf)
        st_job = writer.write(run_dir / "run.sh", job, mode=0o750)
        print(f"  run{irun}: {start.isoformat()} + {duration} days, "
              f"restart from {prev_dir}/ | run.conf {st_conf}, run.sh {st_job}")

        extra = required_container_files(conf) + cfg.CASE_EXTRA_FILES.get(case, [])
        for f in extra:
            src = cfg.here / f
            if not src.exists():
                sys.exit(f"ERROR: {f} is required by case {case} but not found in {cfg.here}")
            st = writer.write(run_dir / Path(f).name, src.read_text())
            print(f"        {Path(f).name} {st}")

    st = writer.write(case_dir / "submit_chain.sh", submit_chain_script(case, cfg.nrun), mode=0o750)
    print(f"  submit_chain.sh {st}")


# ============================================================================
# Main
# ============================================================================
def main(cfg, description=None):
    ap = argparse.ArgumentParser(description=description or __doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cases", nargs="*", help="case names (run_<CASE>_template.conf)")
    ap.add_argument("--force", action="store_true", help="overwrite run.conf/run.sh that differ")
    ap.add_argument("--dry-run", action="store_true", help="show what would be done")
    ap.add_argument("--list", action="store_true", help="list available case templates")
    args = ap.parse_args()

    prefix, suffix = cfg.CASE_TEMPLATE.split("{case}")
    templates = sorted(p.name[len(prefix):-len(suffix)]
                       for p in cfg.here.glob(f"{prefix}*{suffix}") if p.name != cfg.COMMON_TEMPLATE)
    if args.list or not args.cases:
        print("available cases:", " ".join(templates))
        return

    common_text = (cfg.here / cfg.COMMON_TEMPLATE).read_text()
    job_text = (cfg.here / cfg.JOB_TEMPLATE).read_text()

    first_restart = cfg.here / cfg.SPINUP_DIR / f"{restart_name(cfg.START_DATE)}.pe000000.nc"
    if not first_restart.exists():
        print(f"WARNING: {first_restart.relative_to(cfg.here)} does not exist (needed by run1)")
    if not (cfg.here / cfg.BIN).exists():
        print(f"WARNING: binary {cfg.BIN} does not exist")

    writer = Writer(args.force, args.dry_run)
    for case in args.cases:
        generate_case(cfg, case, common_text, job_text, writer, templates)

    if writer.kept:
        print(f"\nWARNING: {len(writer.kept)} existing file(s) differ from the templates and were NOT updated:")
        for f in writer.kept:
            print(f"  {f.relative_to(cfg.here)}")
        print("  -> rerun with --force to overwrite them.")
    if args.dry_run:
        print("(dry run: nothing was written)")
