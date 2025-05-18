import numpy as np
import os
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
import visualize_common

#-------------------------------------
exp_case="../Schaer_case/"; 
topo_params = {"name": "Schaer", "h0": 250.0, "a": 5e3, "lam":4e3, "xc": 25e3}
U0 = 10.0; PTEMP0=280.0; BruntFreq = 1e-2; 

# Parameters for linear solution
Nx=4096*4; Nz=300; Lx=200e3; Lx_lin=Lx*8; Lz=30e3; 

# Paramters for visualization
umin=-2e-0; umax=2e-0; wmin=-2.2e-0; wmax=2.2e-0; 
u_lev = np.linspace(-2.2e-0,2.2e-0,23); w_lev = np.linspace(-2.2e-0,2.2e-0,45); 

#-------------------------------------
DATA_DIR=f"{exp_case}/outdata/"
VIS_OUT_DIR=f"{exp_case}/vis"

os.makedirs(VIS_OUT_DIR, exist_ok=True)

u = visualize_common.open_tmp_nc(f"{DATA_DIR}", "U")
w = visualize_common.open_tmp_nc(f"{DATA_DIR}", "W")

u_lin, w_lin = visualize_common.gen_linsol(topo_params, U0, PTEMP0, BruntFreq, Lx_lin, Nx, Lz, Nz)

xc = topo_params["xc"]
xlim = [xc-10e3,xc+10e3]
zlim = [0.0, 10e3]

u_ = u.isel(y=0).sel(time=3600*10)
visualize_common.plot_var_xz("U'", f"{VIS_OUT_DIR}/U_dash.png", u_-U0, umin, umax, u_lev, xlim=xlim, zlim=zlim, v_lin=u_lin)
w_ = w.isel(y=0).sel(time=3600*10)
visualize_common.plot_var_xz("W", f"{VIS_OUT_DIR}/W.png", w_, wmin, wmax, w_lev, xlim=xlim, zlim=zlim, v_lin=w_lin)
