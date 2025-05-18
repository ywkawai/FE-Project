#! /bin/bash -x

echo "+make directory"
mkdir -p analysis

### check error norm ###
# for var in DDENS ENGT ENGP ENGK ENGI; do
#    python ../common/cmd_analysis_monitor.py monitor.peall ${var} 180.0 analysis/monitor_${var}.png
# done

### make figures ###
echo "+mkgraph"
for time in 0 300 600 900; do
   python ../common/cmd_mkgraph.py history.pe00\*.nc@PT_diff,y=0e0,time=${time} analysis/PT_diff_t${time}s.png --prc_num_xy 4 1 --figsize 11 4 --range " -9.5" 0.5 --int 1.0 --xlim 0e0 19e3 --ylim 0e0 4e3
   python ../common/cmd_mkgraph.py history.pe00\*.nc@U,y=0e0,time=${time} analysis/U_t${time}s.png --prc_num_xy 4 1 --figsize 11 4 --range " -24" 24 --int 2.0 --xlim 0e0 19e3 --ylim 0e0 4e3
   python ../common/cmd_mkgraph.py history.pe00\*.nc@W,y=0e0,time=${time} analysis/W_t${time}s.png --prc_num_xy 4 1 --figsize 11 4 --range " -24" 24 --int 2.0 --xlim 0e0 19e3 --ylim 0e0 4e3
done

### make animation ###
# echo "+make animation"
# python visualize/mkanim.py analysis/advdiff1d.mp4 -1.1 1.1

### check error norm ###
# echo "+check numerical errors"
# python visualize/mkgraph_numerror.py LOG_NUMERROR.peall analysis/numerror_
