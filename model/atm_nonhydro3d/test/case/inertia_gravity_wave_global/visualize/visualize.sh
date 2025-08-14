#! /bin/bash -x

echo "+make directory"
mkdir -p analysis

### check error norm ###
echo "+mkgraph monitor"
for var in DDENS ENGT ENGP ENGK ENGI; do
   python ../common/cmd_analysis_monitor.py monitor.peall ${var} 360 analysis/monitor_${var}.png
done

### make figures ###
echo "+mkgraph"
for time in 0 43200 86400; do
    python ../common/cmd_mkgraph.py outdata/history.pe00\*.nc@Umet,time=${time},lat=0 analysis/Umet_lat0_t${time}s.png --range -15 15 --int 2.5 --prc_num_xy 2 2 --merge_coords lat lon --figsize 11 4 
    python ../common/cmd_mkgraph.py outdata/history.pe00\*.nc@W,time=${time},lat=0 analysis/W_lat0_t${time}s.png --range -0.03 0.03 --int 0.003 --prc_num_xy 2 2 --merge_coords lat lon --figsize 11 4 
    python ../common/cmd_mkgraph.py outdata/history.pe00\*.nc@PT_diff,time=${time},lat=0 analysis/PT_diff_lat0_t${time}s.png --range -4 4 --int .4 --prc_num_xy 2 2 --merge_coords lat lon --figsize 11 4 
    python ../common/cmd_mkgraph.py outdata/history.pe00\*.nc@PT_diff,time=${time},z=5000 analysis/PT_diff_z5000_t${time}s.png --range -4 4 --int .4 --prc_num_xy 2 2 --merge_coords lat lon
done
for time in 172800 259200 345600; do
    python ../common/cmd_mkgraph.py outdata/history.pe00\*.nc@Umet,time=${time},lat=0 analysis/Umet_lat0_t${time}s.png --range -8 8 --int 1 --prc_num_xy 2 2 --merge_coords lat lon --figsize 11 4 
    python ../common/cmd_mkgraph.py outdata/history.pe00\*.nc@W,time=${time},lat=0 analysis/W_lat0_t${time}s.png --range -0.025 0.025 --int 0.0025 --prc_num_xy 2 2 --merge_coords lat lon --figsize 11 4 
    python ../common/cmd_mkgraph.py outdata/history.pe00\*.nc@PT_diff,time=${time},lat=0 analysis/PT_diff_lat0_t${time}s.png --range -2 2 --int .2 --prc_num_xy 2 2 --merge_coords lat lon --figsize 11 4 
    python ../common/cmd_mkgraph.py outdata/history.pe00\*.nc@PT_diff,time=${time},z=5000 analysis/PT_diff_z5000_t${time}s.png --range -2 2 --int .2 --prc_num_xy 2 2 --merge_coords lat lon
done

### make animation ###
# echo "+make animation"
# python visualize/mkanim.py analysis/advdiff1d.mp4 -1.1 1.1

### check error norm ###
# echo "+check numerical errors"
# python visualize/mkgraph_numerror.py LOG_NUMERROR.peall analysis/numerror_
