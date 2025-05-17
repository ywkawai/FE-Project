#! /bin/bash -x

echo "+make directory"
mkdir -p analysis

### check error norm ###
for var in DDENS ENGT ENGP ENGK ENGI; do
   python ../common/cmd_analysis_monitor.py monitor.peall ${var} 180.0 analysis/monitor_${var}.png
done

### make figures ###
echo "+mkgraph"
for time in 0 432000 864000; do
   python ../common/cmd_mkgraph.py history.pe00\*.nc@PRES,z=500e0,time=${time} analysis/PRES_t${time}s.png --prc_num_xy 4 2 --figsize 11 4 --interp --range 929e2 949e2 --int 2e2
   python ../common/cmd_mkgraph.py history.pe00\*.nc@T,z=500e0,time=${time} analysis/T_t${time}s.png --prc_num_xy 4 2 --figsize 11 4 --interp --range 268 304 --int 2
done

### make animation ###
# echo "+make animation"
# python visualize/mkanim.py analysis/advdiff1d.mp4 -1.1 1.1

### check error norm ###
# echo "+check numerical errors"
# python visualize/mkgraph_numerror.py LOG_NUMERROR.peall analysis/numerror_
