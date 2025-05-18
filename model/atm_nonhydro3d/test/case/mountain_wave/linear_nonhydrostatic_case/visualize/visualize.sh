#! /bin/bash -x

echo "+make directory"
mkdir -p analysis

### check error norm ###
for var in DDENS ENGT ENGP ENGK ENGI; do
   python ../../common/cmd_analysis_monitor.py monitor.peall ${var} 9.6 analysis/monitor_${var}.png
done

### make figures ###
echo "+mkgraph"
for time in 0 3600 7200 18000; do
   python ../../common/cmd_mkgraph.py reg_outdata/history.pe00\*.nc@U,y=0e0,time=${time} analysis/U_t${time}s.png --prc_num_xy 4 1 --figsize 10 5 --range 9.992 10.008 --int 1e-3 --xlim 60e3 100e3 --ylim 0e3 16e3
   python ../../common/cmd_mkgraph.py reg_outdata/history.pe00\*.nc@W,y=0e0,time=${time} analysis/W_t${time}s.png --prc_num_xy 4 1 --figsize 10 5 --range " -6e-3" 6e-3 --int 1e-3 --xlim 60e3 100e3 --ylim 0e0 16e3
done

### make animation ###
# echo "+make animation"
# python visualize/mkanim.py analysis/advdiff1d.mp4 -1.1 1.1

### check error norm ###
# echo "+check numerical errors"
# python visualize/mkgraph_numerror.py LOG_NUMERROR.peall analysis/numerror_
