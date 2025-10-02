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
for time in 600 1800; do
   python ../common/cmd_mkgraph.py history.pe00\*.nc@U,z=800e0,time=${time} analysis/U_z800m_t${time}s.png --prc_num_xy 4 4 --figsize 8 8 --range -0.8 0.8 --int 0.08
   python ../common/cmd_mkgraph.py history.pe00\*.nc@V,z=800e0,time=${time} analysis/V_z800m_t${time}s.png --prc_num_xy 4 4 --figsize 8 8 --range -0.8 0.8 --int 0.08
   python ../common/cmd_mkgraph.py history.pe00\*.nc@W,z=800e0,time=${time} analysis/W_z800m_t${time}s.png --prc_num_xy 4 4 --figsize 8 8 --range -0.8 0.8 --int 0.08
   python ../common/cmd_mkgraph.py history.pe00\*.nc@PT,z=800e0,time=${time} analysis/PT_z800m_t${time}s.png --prc_num_xy 4 4 --figsize 8 8 --range 299.6 300.4 --int 0.04
   python ../common/cmd_mkgraph.py history.pe00\*.nc@U,y=0e0,time=${time} analysis/U_y0m_t${time}s.png --prc_num_xy 4 4 --figsize 8 4 --range -0.8 0.8 --int 0.08
   python ../common/cmd_mkgraph.py history.pe00\*.nc@V,y=0e0,time=${time} analysis/V_y0m_t${time}s.png --prc_num_xy 4 4 --figsize 8 4 --range -0.8 0.8 --int 0.08
   python ../common/cmd_mkgraph.py history.pe00\*.nc@W,y=0e0,time=${time} analysis/W_y0m_t${time}s.png --prc_num_xy 4 4 --figsize 8 4 --range -0.8 0.8 --int 0.08
   python ../common/cmd_mkgraph.py history.pe00\*.nc@PT,y=0e0,time=${time} analysis/PT_y0m_t${time}s.png --prc_num_xy 4 4 --figsize 8 4 --range 299.6 300.4 --int 0.04
done
for time in 3600 5400; do
   python ../common/cmd_mkgraph.py history.pe00\*.nc@U,z=800e0,time=${time} analysis/U_z800m_t${time}s.png --prc_num_xy 4 4 --figsize 8 8 --range -1.5 1.5 --int 0.15
   python ../common/cmd_mkgraph.py history.pe00\*.nc@V,z=800e0,time=${time} analysis/V_z800m_t${time}s.png --prc_num_xy 4 4 --figsize 8 8 --range -1.5 1.5 --int 0.15
   python ../common/cmd_mkgraph.py history.pe00\*.nc@W,z=800e0,time=${time} analysis/W_z800m_t${time}s.png --prc_num_xy 4 4 --figsize 8 8 --range -1.5 1.5 --int 0.15
   python ../common/cmd_mkgraph.py history.pe00\*.nc@PT,z=800e0,time=${time} analysis/PT_z800m_t${time}s.png --prc_num_xy 4 4 --figsize 8 8 --range 299.6 300.4 --int 0.04
   python ../common/cmd_mkgraph.py history.pe00\*.nc@U,y=0e0,time=${time} analysis/U_y0m_t${time}s.png --prc_num_xy 4 4 --figsize 8 4 --range -1.5 1.5 --int 0.15
   python ../common/cmd_mkgraph.py history.pe00\*.nc@V,y=0e0,time=${time} analysis/V_y0m_t${time}s.png --prc_num_xy 4 4 --figsize 8 4 --range -1.5 1.5 --int 0.15
   python ../common/cmd_mkgraph.py history.pe00\*.nc@W,y=0e0,time=${time} analysis/W_y0m_t${time}s.png --prc_num_xy 4 4 --figsize 8 4 --range -1.5 1.5 --int 0.15
   python ../common/cmd_mkgraph.py history.pe00\*.nc@PT,y=0e0,time=${time} analysis/PT_y0m_t${time}s.png --prc_num_xy 4 4 --figsize 8 4 --range 299.6 300.4 --int 0.04
done

### make animation ###
# echo "+make animation"
# python visualize/mkanim.py analysis/advdiff1d.mp4 -1.1 1.1

### check error norm ###
# echo "+check numerical errors"
# python visualize/mkgraph_numerror.py LOG_NUMERROR.peall analysis/numerror_
