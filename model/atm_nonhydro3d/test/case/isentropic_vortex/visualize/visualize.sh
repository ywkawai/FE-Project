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
for time in 0 1 2 3 4 5 6 7 8; do
   python ../common/cmd_mkgraph.py history.pe00\*.nc@DDENS,z=0e0,time=${time} analysis/DDENS_t${time}s.png --prc_num_xy 2 2 --figsize 8 8 --range " -0.7" 0.7
   python ../common/cmd_mkgraph.py history.pe00\*.nc@U,z=0e0,time=${time} analysis/U_t${time}s.png --prc_num_xy 2 2 --figsize 8 8 --range 4 6
   python ../common/cmd_mkgraph.py history.pe00\*.nc@V,z=0e0,time=${time} analysis/V_t${time}s.png --prc_num_xy 2 2 --figsize 8 8 --range 4 6
   python ../common/cmd_mkgraph.py history.pe00\*.nc@PT,z=0e0,time=${time} analysis/PT_t${time}s.png --prc_num_xy 2 2 --figsize 8 8
done

### make animation ###
# echo "+make animation"
# python visualize/mkanim.py analysis/advdiff1d.mp4 -1.1 1.1

### check error norm ###
# echo "+check numerical errors"
# python visualize/mkgraph_numerror.py LOG_NUMERROR.peall analysis/numerror_
