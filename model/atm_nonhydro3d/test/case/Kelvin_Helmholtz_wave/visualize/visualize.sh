#! /bin/bash -x

echo "+make directory"
mkdir -p analysis

### check error norm ###
echo "+mkgraph monitor"
for var in DDENS ENGT ENGP ENGK ENGI; do
   python ../common/cmd_analysis_monitor.py monitor.peall ${var} 0.25 analysis/monitor_${var}.png
done

### make figures ###
echo "+mkgraph"
for time in 0 900 1800 2400 3600 5400; do
   python ../common/cmd_mkgraph.py history.pe00\*.nc@PT,y=0e0,time=${time} analysis/PT_y0_t${time}s.png --prc_num_xy 4 2 --figsize 8 3 --range 300.0 301.0 --int 0.1
   python ../common/cmd_mkgraph.py history.pe00\*.nc@U,y=0e0,time=${time} analysis/U_y0_t${time}s.png --prc_num_xy 4 2 --figsize 8 3 --range -10 60 --int 2.5
   python ../common/cmd_mkgraph.py history.pe00\*.nc@W,y=0e0,time=${time} analysis/W_y0_t${time}s.png --prc_num_xy 4 2 --figsize 8 3 --range -10 10 --int 1
done

### make animation ###
# echo "+make animation"
# python visualize/mkanim.py analysis/advdiff1d.mp4 -1.1 1.1

### check error norm ###
# echo "+check numerical errors"
# python visualize/mkgraph_numerror.py LOG_NUMERROR.peall analysis/numerror_
