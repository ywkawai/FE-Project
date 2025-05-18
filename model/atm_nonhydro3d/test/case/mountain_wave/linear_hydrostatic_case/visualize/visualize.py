import numpy as np
import os
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
import visualize_common

#-------------------------------------
exp_case="../linear_hydrostatic_case/"; 
topo_params = {"name": "bell_shape", "h0": 1.0, "a": 10e3, "xc": 120e3}
U0 = 20.0; PTEMP0=250.0; BruntFreq = 1.9567954884127995E-2; 

# Parameters for linear solution
Nx=4096*4; Nz=300; Lx=480e3; Lx_lin=Lx*8; Lz=30e3; 

# Paramters for visualization
umin=-3e-2; umax=3e-2; wmin=-5e-3; wmax=5e-3; 
u_lev = np.linspace(-3e-2,3e-2,13); w_lev = np.linspace(-5e-3,5e-3,21); 

#-------------------------------------
DATA_DIR=f"{exp_case}/outdata/"
VIS_OUT_DIR=f"{exp_case}/vis"

os.makedirs(VIS_OUT_DIR, exist_ok=True)

u = visualize_common.open_tmp_nc(f"{DATA_DIR}", "U")
w = visualize_common.open_tmp_nc(f"{DATA_DIR}", "W")

u_lin, w_lin = visualize_common.gen_linsol(topo_params, U0, PTEMP0, BruntFreq, Lx_lin, Nx, Lz, Nz)

xc = topo_params["xc"]
xlim = [xc-40e3,xc+40e3]
zlim = [0.0, 12e3]

u_ = u.isel(y=0).sel(time=3600*10)
visualize_common.plot_var_xz("U'", f"{VIS_OUT_DIR}/U_dash.png", u_-U0, umin, umax, u_lev, xlim=xlim, zlim=zlim, v_lin=u_lin)
w_ = w.isel(y=0).sel(time=3600*10)
visualize_common.plot_var_xz("W", f"{VIS_OUT_DIR}/W.png", w_, wmin, wmax, w_lev, xlim=xlim, zlim=zlim, v_lin=w_lin)
