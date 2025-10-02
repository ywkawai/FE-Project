#! /bin/bash -x

echo "+make directory"
mkdir -p analysis

### check error norm ###
echo "+mkgraph monitor"
for var in DDENS ENGT ENGP ENGK ENGI; do
   python ../common/cmd_analysis_monitor.py monitor.peall ${var} 0.15 analysis/monitor_${var}.png
done

### make figures ###
echo "+mkgraph"
for time in 300 900 1800 2400 3600; do
   python ../common/cmd_mkgraph.py history.pe00\*.nc@U,z=50e0,time=${time} analysis/U_z50m_t${time}s.png --prc_num_xy 4 2 --figsize 8 8 --range 2 8 --int 0.3
   python ../common/cmd_mkgraph.py history.pe00\*.nc@V,z=50e0,time=${time} analysis/V_z50m_t${time}s.png --prc_num_xy 4 2 --figsize 8 8 --range -3 3 --int 0.3
   python ../common/cmd_mkgraph.py history.pe00\*.nc@W,z=50e0,time=${time} analysis/W_z50m_t${time}s.png --prc_num_xy 4 2 --figsize 8 8 --range -3 3 --int 0.3
   python ../common/cmd_mkgraph.py history.pe00\*.nc@U,z=250e0,time=${time} analysis/U_z250m_t${time}s.png --prc_num_xy 4 2 --figsize 8 8 --range 2 8 --int 0.3
   python ../common/cmd_mkgraph.py history.pe00\*.nc@V,z=250e0,time=${time} analysis/V_z250m_t${time}s.png --prc_num_xy 4 2 --figsize 8 8 --range -3 3 --int 0.3
   python ../common/cmd_mkgraph.py history.pe00\*.nc@W,z=250e0,time=${time} analysis/W_z250m_t${time}s.png --prc_num_xy 4 2 --figsize 8 8 --range -3 3 --int 0.3
   python ../common/cmd_mkgraph.py history.pe00\*.nc@W,y=0e0,time=${time} analysis/W_y0m_t${time}s.png --prc_num_xy 4 2 --figsize 6 8 --range -3 3 --int 0.3
done

### make animation ###
# echo "+make animation"
# python visualize/mkanim.py analysis/advdiff1d.mp4 -1.1 1.1

### check error norm ###
# echo "+check numerical errors"
# python visualize/mkgraph_numerror.py LOG_NUMERROR.peall analysis/numerror_
