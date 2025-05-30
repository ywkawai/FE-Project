#! /bin/bash -x

echo "+make directory"
mkdir -p analysis

### check error norm ###
echo "+mkgraph monitor"
for var in DDENS QV QC QR ENGT ENGP ENGK ENGI; do
   python ../../common/cmd_analysis_monitor.py monitor.peall ${var} 3.0 analysis/monitor_${var}.png
done

### make figures ###
echo "+mkgraph"
for time in 300 600 900 1200 1500 1800; do
   python ../../common/cmd_mkgraph.py history.pe00*.nc@U,y=16e3,time=${time} analysis/U_y16km_t${time}s.png --prc_num_xy 5 2 --figsize 12 4 --range -20 20 --int 2
   python ../../common/cmd_mkgraph.py history.pe00*.nc@W,y=16e3,time=${time} analysis/W_y16km_t${time}s.png --prc_num_xy 5 2 --figsize 12 4 --range -40 40 --int 5
   python ../../common/cmd_mkgraph.py history.pe00*.nc@QC,y=16e3,time=${time} analysis/QC_y16km_t${time}s.png --prc_num_xy 5 2 --figsize 12 4
   python ../../common/cmd_mkgraph.py history.pe00*.nc@QR,y=16e3,time=${time} analysis/QR_y16km_t${time}s.png --prc_num_xy 5 2 --figsize 12 4
done

for time in 2700 3600 5400; do
   python ../../common/cmd_mkgraph.py history.pe00*.nc@U,y=16e3,time=${time} analysis/U_y16km_t${time}s.png --prc_num_xy 5 2 --figsize 12 4 --range -20 20 --int 2
   python ../../common/cmd_mkgraph.py history.pe00*.nc@W,y=16e3,time=${time} analysis/W_y16km_t${time}s.png --prc_num_xy 5 2 --figsize 12 4 --range -20 20 --int 2.5
   python ../../common/cmd_mkgraph.py history.pe00*.nc@QC,y=16e3,time=${time} analysis/QC_y16km_t${time}s.png --prc_num_xy 5 2 --figsize 12 4
   python ../../common/cmd_mkgraph.py history.pe00*.nc@QR,y=16e3,time=${time} analysis/QR_y16km_t${time}s.png --prc_num_xy 5 2 --figsize 12 4
done

python ../../common/cmd_mkgraph.py history.pe00*.nc@PREC,y=16e3 analysis/PREC_y16km.png --exch --prc_num_xy 5 2 --figsize 12 6

### make animation ###
# echo "+make animation"
# python visualize/mkanim.py analysis/advdiff1d.mp4 -1.1 1.1

### check error norm ###
# echo "+check numerical errors"
# python visualize/mkgraph_numerror.py LOG_NUMERROR.peall analysis/numerror_
