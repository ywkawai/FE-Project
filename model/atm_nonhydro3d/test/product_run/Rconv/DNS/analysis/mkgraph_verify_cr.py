import os
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
sys.path.append(os.path.join(os.path.dirname('__file__'), '../../analysis_lib/'))
import energy_spectra_common as common

OUTDIR="fig_verify_cr/"
outfigext = "svg"

exp_name_list = [
  'Dx3.1m_P7',           
]
dir_ind_list = { 
  'Dx3.1m_P7': np.arange(13,101, 1), 
}
exp_ncut = {
  'Dx3.1m_P7': 512, 
}

ZLEVEL_list = [800]
L = 3.2e3

#---------------------------------------------------------------------------------------
def read_spectra_data(exp_name, tmp_dir, dir_ind, ke_spectra_3dvel_list, 
                      ncut, start_ind=1, OUTNC_suffix=""):
  di = dir_ind[0]
  fac = 1.0 / float(len(dir_ind))
  ke_spectra_hvel = xr.open_mfdataset(f'{tmp_dir}_{di}/KE_spectra_momh{OUTNC_suffix}.nc', decode_times=False, combine='by_coords').KE_spectra_momh
  ke_spectra_wvel = xr.open_mfdataset(f'{tmp_dir}_{di}/KE_spectra_momz{OUTNC_suffix}.nc', decode_times=False, combine='by_coords').KE_spectra_momz
  ke_spectra_3dvel_tmp = ke_spectra_hvel[start_ind:,:,1:ncut] + ke_spectra_wvel[start_ind:,:,1:ncut]
  time = ke_spectra_3dvel_tmp.time
  dtime = time[1] - time[0]
  time_s = time[-1] + dtime
  for i in range(1,len(dir_ind)):
    di = dir_ind[i]
    ke_spectra_hvel = xr.open_mfdataset(f'{tmp_dir}_{di}/KE_spectra_momh{OUTNC_suffix}.nc', decode_times=False, combine='by_coords').KE_spectra_momh
    ke_spectra_wvel = xr.open_mfdataset(f'{tmp_dir}_{di}/KE_spectra_momz{OUTNC_suffix}.nc', decode_times=False, combine='by_coords').KE_spectra_momz
    tmp = ke_spectra_hvel[start_ind:,:,1:ncut] + ke_spectra_wvel[start_ind:,:,1:ncut]
    tmp = tmp.assign_coords(time=tmp.time+time_s)
    ke_spectra_3dvel_tmp = xr.concat([ke_spectra_3dvel_tmp, tmp ], "time")
    time_s = tmp.time[-1] + dtime
    
  ke_spectra_3dvel_list[exp_name] = ke_spectra_3dvel_tmp
  return

def eval_Neff( ke_spectra_3dvel ):
  k = ke_spectra_3dvel.k
  dk = k[1] - k[0]
  ke_tot = (ke_spectra_3dvel*dk).sum(["k"])
  kc = ke_tot - ke_tot.mean(["time"])
  conv1 = ( kc * kc.shift(time=1)).mean("time", skipna=True)
  var0 = (kc * kc).mean("time", skipna=True)
  rho1 = (conv1/var0).compute().item()

  N = int(len(ke_tot.time.values))
  Neff = N * (1-rho1) / (1+rho1)
  return Neff, rho1

#-
ke_spectra_3dvel_list = {}

for exp_name in exp_name_list:
  tmp_dir = f"../tmp_data_energy_spectra/tmp_{exp_name}"
  print(f"tmp_dir={tmp_dir}")

  read_spectra_data( exp_name, tmp_dir, dir_ind_list[exp_name], ke_spectra_3dvel_list, exp_ncut[exp_name], start_ind=0 )
  print(ke_spectra_3dvel_list[exp_name].time)
  
ke_spectra_3dvel_ref = ke_spectra_3dvel_list[exp_name_list[0]].sel(z=ZLEVEL_list[0])

#--

time = ke_spectra_3dvel_ref.time
Nt = len(time.values)
Neff, rho1 = eval_Neff(ke_spectra_3dvel_ref)
print(f"Neff={Neff}, rho1={rho1}")
print(Nt)

vari = ( (ke_spectra_3dvel_ref - ke_spectra_3dvel_ref.mean(["time"]))**2 ).sum(["time"])  / float(Nt-1)
stddev = ( vari )**0.5
stderr = ( vari / float(Neff)  )**0.5

#--

ke_tmean = ke_spectra_3dvel_ref.mean(["time"])
ke_ratio_stddev_0p8 = (ke_tmean - 0.67*stddev)/ke_tmean
ke_ratio_stddev_1 = (ke_tmean - stddev)/ke_tmean
ke_ratio_stddev_2 = (ke_tmean - 2*stddev)/ke_tmean
ke_ratio_stderr_5 = (ke_tmean - 5*stderr)/ke_tmean
ke_ratio_stderr_8 = (ke_tmean - 8*stderr)/ke_tmean
ke_ratio_stderr_10 = (ke_tmean - 10*stderr)/ke_tmean
fig, ax = plt.subplots(1,1,figsize=(8,6))
ax.set_xscale("log")
ax.set_xlabel("inverse of wavelength [m$^{-1}$]", fontsize=15)
ax.xaxis.set_major_formatter(FuncFormatter(common.xaxis_txt))
ax.xaxis.set_minor_formatter(FuncFormatter(common.xaxis_txt_minor))

# ax.set_yscale("log")
kw = ke_ratio_stddev_1.k/(2.0*np.pi)
ax.plot(kw, ke_ratio_stddev_0p8, label="stddev_0.67sig", color="red", linestyle="--")
ax.plot(kw, ke_ratio_stddev_1, label="stddev_1sig", color="red")
ax.plot(kw, ke_ratio_stddev_2, label="stddev_2sig", color="red", linestyle="--")
ax.plot(kw, ke_ratio_stderr_5, label="stderr_5sig", color="blue", alpha=0.2)
# ax.plot(kw, ke_ratio_stderr_8, label="stderr_8sig", color="blue", alpha=0.2)
ax.plot(kw, ke_ratio_stderr_10, label="stderr_10sig", color="blue", alpha=0.2)
ax.set_ylim(0.65,0.93)
# ax.legend()
ax.grid(axis="y")

os.makedirs(OUTDIR, exist_ok=True)
plt.savefig(f"{OUTDIR}/verify_cr.{outfigext}", bbox_inches="tight")