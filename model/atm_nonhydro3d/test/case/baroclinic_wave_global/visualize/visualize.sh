#! /bin/bash -x

echo "+make directory"
mkdir -p analysis

### check error norm ###
echo "+mkgraph monitor"
for var in DDENS ENGT ENGP ENGK ENGI; do
   python ../common/cmd_analysis_monitor.py monitor.peall ${var} 600 analysis/monitor_${var}.png
done

### make figures ###
for time in 0 345600 518400; do
    python ../common/cmd_mkgraph.py outdata/history.pe00\*.nc@PRES,time=${time},z=0,lat=0:90,lon=45:360 analysis/Ps_t${time}s.png  --range 992e2 1006e2 --int 1e2 --prc_num_xy 4 2 --merge_coords lat lon --figsize 11 4
done 

time=691200
python ../common/cmd_mkgraph.py outdata/history.pe00\*.nc@PRES,time=${time},z=0,lat=0:90,lon=45:360 analysis/Ps_t${time}s.png  --range 980e2 1012e2 --int 4e2 --prc_num_xy 4 2 --merge_coords lat lon --figsize 11 4

for time in 777600 8664000; do
    python ../common/cmd_mkgraph.py outdata/history.pe00\*.nc@PRES,time=${time},z=0,lat=0:90,lon=45:360 analysis/Ps_t${time}s.png  --range 930e2 1030e2 --int 10e2 --prc_num_xy 4 2 --merge_coords lat lon --figsize 11 4
done

for time in 0 345600 518400 691200 777600 8664000; do
    python ../common/cmd_mkgraph.py outdata_p/history.pe00\*.nc@T,time=${time},p=850e2,lat=0:90,lon=45:360 analysis/T_p850hPa_t${time}s.png --range 220 310 --int 10 --prc_num_xy 4 2 --merge_coords lat lon --figsize 11 4 
done 

### make animation ###
# echo "+make animation"
# python visualize/mkanim.py analysis/advdiff1d.mp4 -1.1 1.1

### check error norm ###
# echo "+check numerical errors"
# python visualize/mkgraph_numerror.py LOG_NUMERROR.peall analysis/numerror_
