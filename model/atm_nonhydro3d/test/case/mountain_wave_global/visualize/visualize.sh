#! /bin/bash -x

echo "+make directory"
mkdir -p analysis

### check error norm ###
echo "+mkgraph monitor"
# for var in DDENS ENGT ENGP ENGK ENGI; do
#    python ../common/cmd_analysis_monitor.py monitor.peall ${var} 300 analysis/monitor_${var}.png
# done

### make figures ###
echo "+mkgraph"

for time in 0 43200 86400 129600 172800; do
    python ../common/cmd_mkgraph.py outdata/history.pe00\*.nc@Umet,time=${time},z=5000 analysis/Umet_z5000m_t${time}s.png --interp --range 0 50 --int 2 --prc_num_xy 2 2 --merge_coords lat lon --figsize 11 4
    python ../common/cmd_mkgraph.py outdata/history.pe00\*.nc@W,time=${time},z=5000 analysis/W_z5000m_t${time}s.png --interp --range -2.5e-2 2.5e-2 --int 2.5e-3 --prc_num_xy 2 2 --merge_coords lat lon --figsize 11 4
    python ../common/cmd_mkgraph.py outdata/history.pe00\*.nc@Umet,time=${time},lat=0 analysis/Umet_lat0_t${time}s.png --range 30 50 --int 2 --prc_num_xy 2 2 --merge_coords lat lon --figsize 11 4
    python ../common/cmd_mkgraph.py outdata/history.pe00\*.nc@W,time=${time},lat=0 analysis/W_lat0_t${time}s.png --range -5e-2 5e-2 --int 5e-3 --prc_num_xy 2 2 --merge_coords lat lon --figsize 11 4
done

### make animation ###
# echo "+make animation"
# python visualize/mkanim.py analysis/advdiff1d.mp4 -1.1 1.1

### check error norm ###
# echo "+check numerical errors"
# python visualize/mkgraph_numerror.py LOG_NUMERROR.peall analysis/numerror_
