#! /bin/bash -x

echo "+make directory"
mkdir -p analysis

### check error norm ###
echo "+mkgraph monitor"
python ../common/cmd_analysis_monitor.py monitor.peall "PTracer" 2.0 analysis/monitor_PTracer.png

### make figures ###
echo "+mkgraph"
for time in 0 300 600 900; do
    python ../common/cmd_mkgraph.py history.pe00\*.nc@PTracer,time=${time},z=5e3 analysis/PTracer_t${time}s.png --range 0.0 1.0 --prc_num_xy 2 2
done

### make animation ###
# echo "+make animation"
# python visualize/mkanim.py analysis/advdiff1d.mp4 -1.1 1.1

### check error norm ###
# echo "+check numerical errors"
# python visualize/mkgraph_numerror.py LOG_NUMERROR.peall analysis/numerror_
