#!/bin/bash -x

echo "+make directory"
mkdir -p analysis

### check error norm ###
# echo "+mkgraph monitor"
# for var in DDENS ENGT ENGP ENGK ENGI; do
#    python ../common/cmd_analysis_monitor.py monitor.peall ${var} 0.25 analysis/monitor_${var}.png
# done

### make figures ###
echo "+mkgraph"
python ../common/cmd_mkgraph.py history.pe00\*.nc@T,x=0e0,y=0e0 analysis/T.png --prc_num_xy 1 1 --figsize 8 4
python ../common/cmd_mkgraph.py history.pe00\*.nc@BL_RHOU_t,x=0e0,y=0e0 analysis/BL_RHOU_t.png --prc_num_xy 1 1 --figsize 8 4 --range -4e-4 4e-4
python ../common/cmd_mkgraph.py history.pe00\*.nc@BL_RHOT_t,x=0e0,y=0e0 analysis/BL_RHOT_t.png --prc_num_xy 1 1 --figsize 8 4
python ../common/cmd_mkgraph.py history.pe00\*.nc@NU,x=0e0,y=0e0 analysis/NU.png --prc_num_xy 1 1 --figsize 8 4 --range 0 30
python ../common/cmd_mkgraph.py history.pe00\*.nc@KH,x=0e0,y=0e0 analysis/KH.png --prc_num_xy 1 1 --figsize 8 4 --range 0 30
python ../common/cmd_mkgraph.py history.pe00\*.nc@U,x=0e0,y=0e0 analysis/U.png --prc_num_xy 1 1 --figsize 8 4

### make animation ###
# echo "+make animation"
# python visualize/mkanim.py analysis/advdiff1d.mp4 -1.1 1.1

### check error norm ###
# echo "+check numerical errors"
# python visualize/mkgraph_numerror.py LOG_NUMERROR.peall analysis/numerror_
