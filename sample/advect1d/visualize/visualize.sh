#! /bin/bash -x

echo "+make directory"
mkdir -p analysis

### make animation ###
echo "+make animation"
python visualize/mkanim.py analysis/advect1d.mp4 -0.1 1.1

### check error norm ###
echo "+check numerical errors"
python visualize/mkgraph_numerror.py LOG_NUMERROR.peall analysis/numerror_
