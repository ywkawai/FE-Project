#!/usr/bin/env python3
"""Generate run directories, run.conf and job scripts for the main runs of Eh8Ez4P7.

Only the resolution-dependent settings are here. The templates (run_common_template.conf,
run_<CASE>_template.conf, job_run_template.sh, phy_preoperation_*NN.conf) are in this directory.
The generator itself and the default settings are in ../gen_config_common.py.

Usage:
    python3 gen_config.py asis                 # generate one case (asis/run1..run4)
    python3 gen_config.py asis case1           # generate several cases at once
    python3 gen_config.py --dry-run PostPR2p7  # show what would be done
    python3 gen_config.py --force case1        # overwrite run.conf / run.sh that differ
    python3 gen_config.py --list               # list the available case templates
"""

import sys
from datetime import date
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from gen_config_common import Config, main  # noqa: E402

cfg = Config(
    here=Path(__file__).resolve().parent,
    # Start of the main run = end of run_spinup2 (run_spinup2/restart_00010625-000000.000)
    START_DATE=date(1, 6, 25),
    # Job: 8 MPI processes (2 nodes x 4 processes x 12 threads)
    NPRC=8,
    NPRC_PER_NODE=4,
    NTHREAD=12,
    ELAPSE="04:00:00",  # run_spinup2: 150 days took about 50 min (asis-like setting)
    JOB_PREFIX="MC_E8P7",
    # Other defaults (TINTEG_DAYS, TINTEG_PER_RUN, TIME_DT_RESTART_DAY, HISTORY_TINTERVAL_DAY,
    # SPINUP_DIR, TREF_FILEBASE, GROUP, RSCGRP, BIN, CASE_EXTRA_FILES, ...) are in
    # ../gen_config_common.py and can be overridden here.
)

if __name__ == "__main__":
    main(cfg, __doc__)
