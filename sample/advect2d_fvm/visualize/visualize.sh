#! /bin/bash -x

echo "+make directory"
mkdir -p analysis

### make animation ###
echo "+make animation"
python visualize/mkanim.py analysis/advect2d.mp4 -1.1 1.1

### check error norm ###
echo "+check numerical errors"
python visualize/mkgraph_numerror.py LOG_NUMERROR.peall analysis/numerror_
