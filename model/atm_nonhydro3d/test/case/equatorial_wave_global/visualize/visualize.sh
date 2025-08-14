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

for z in 250 9750; do
    python ../common/cmd_mkgraph.py outdata/history.pe00\*.nc@Umet,lat=0,z=${z} analysis/Umet_lat0_z${z}m.png --interp --range -11 11 --int 1 --prc_num_xy 2 2 --merge_coords lat lon --figsize 5 8 
done

for time in 0 432000 864000 1728000 2592000; do
    python ../common/cmd_mkgraph.py outdata/history.pe00\*.nc@Umet,time=${time},z=1250 analysis/Umet_z1250m_t${time}s.png --interp --range -11 11 --int 1 --prc_num_xy 2 2 --merge_coords lat lon --figsize 10 5
    python ../common/cmd_mkgraph.py outdata/history.pe00\*.nc@Vmet,time=${time},z=1250 analysis/Vmet_z1250m_t${time}s.png --interp --range -2.4 2.4 --int 0.2 --prc_num_xy 2 2 --merge_coords lat lon --figsize 10 5 
    python ../common/cmd_mkgraph.py outdata/history.pe00\*.nc@Umet,time=${time},z=8750 analysis/Umet_z8750m_t${time}s.png --interp --range -11 11 --int 1 --prc_num_xy 2 2 --merge_coords lat lon --figsize 10 5 
    python ../common/cmd_mkgraph.py outdata/history.pe00\*.nc@Vmet,time=${time},z=8750 analysis/Vmet_z8750m_t${time}s.png --interp --range -2.4 2.4 --int .2 --prc_num_xy 2 2 --merge_coords lat lon --figsize 10 5 
    python ../common/cmd_mkgraph.py outdata/history.pe00\*.nc@Umet,time=${time},lat=0 analysis/Umet_lat0_t${time}s.png --interp --range -11 11 --int 1 --prc_num_xy 2 2 --merge_coords lat lon --figsize 11 4 
    python ../common/cmd_mkgraph.py outdata/history.pe00\*.nc@W,time=${time},lat=0 analysis/W_lat0_t${time}s.png --interp --range -1e-2 1e-2 --int 1e-3 --prc_num_xy 2 2 --merge_coords lat lon --figsize 11 4
done

### make animation ###
# echo "+make animation"
# python visualize/mkanim.py analysis/advdiff1d.mp4 -1.1 1.1

### check error norm ###
# echo "+check numerical errors"
# python visualize/mkgraph_numerror.py LOG_NUMERROR.peall analysis/numerror_
