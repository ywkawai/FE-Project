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
python ../common/cmd_mkgraph.py history.pe00\*.nc@W,x=0e0,y=0e0 analysis/W.png --prc_num_xy 1 1 --exch --figsize 8 6 --range " -1.2e-4" "1.2e-4" 
python ../common/cmd_mkgraph.py history.pe00\*.nc@PRES_diff,x=0e0,y=0e0 analysis/PRES_diff.png --prc_num_xy 1 1 --exch --figsize 8 6 --range " -1.0e-1" "1.0e-1" 

### make animation ###
# echo "+make animation"
# python visualize/mkanim.py analysis/advdiff1d.mp4 -1.1 1.1

### check error norm ###
# echo "+check numerical errors"
# python visualize/mkgraph_numerror.py LOG_NUMERROR.peall analysis/numerror_
