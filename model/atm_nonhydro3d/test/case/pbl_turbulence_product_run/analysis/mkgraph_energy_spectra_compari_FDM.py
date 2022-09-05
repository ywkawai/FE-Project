import matplotlib.pyplot as plt
import numpy as np
from numpy.fft import fftn
import xarray as xr
import glob
import os
from joblib import Parallel, delayed
from matplotlib.ticker import FuncFormatter

ZLEVEL_list=[500]
N=960

exp_name_list = [
  "E80P11",   
  "E120P7",   
  "E160P5",   
  "E192P4",     
  "E240P3", 
  "E480P1",   
  "E480P1_MFoff", 
]
exp_name_list_FDM = [
  "RK3CD2_ND2Gam2e-4",     
  "RK4UD3",   
  "RK6UD5",     
  "RK8UD7",   
]
OUT_DIR="./compari_FDM/"

exp_label_list = {
  "E80P11": "P11",  
  "E120P7": "P7", 
  "E160P5": "P5", 
  "E192P4": "P4",   
  "E240P3": "P3", 
  "E480P1": "P1", 
  "E480P1_MFoff": "P1_MFoff", 
  "reference": "reference (CD8ND8)", 
  "RK3CD2_ND2Gam2e-4": "CD2ND2",   
  "RK4UD3": "UD3",     
  "RK6UD5": "UD5",   
  "RK8UD7": "UD7",       
}
exp_color_list = {
  "E80P11": "black",   
  "E120P7": "goldenrod", 
  "E160P5": "green",
  "E192P4": "red", #"forestgreen",  
  "E240P3": "blue", 
  "E240P3_RK3": "violet", 
  "E480P1": "cyan",   
  "E480P1_MFoff": "lightskyblue",   
  "reference": "black", 
  "RK3CD2_ND2Gam2e-4": "cyan",    
  "RK4UD3": "lightskyblue",    
  "RK6UD5": "blue",    
  "RK8UD7": "red",  
}
exp_ltype_list = {
  "E80P11": "-",
  "E120P7": "-",
  "E160P5": "-",
  "E192P4": "-",  
  "E240P3": "-",
  "E480P1": "-",
  "E480P1_MFoff": "-.",
  "reference": "-",  
  "RK3CD2_ND2Gam2e-4": ":",   
  "RK4UD3": ":",   
  "RK6UD5": ":",     
  "RK8UD7": ":",     
}
exp_ltype_width = {
  "E80P11": 1,  
  "E120P7": 1,
  "E160P5": 1,
  "E192P4": 4,    
  "E240P3": 4,
  "E480P1": 1,
  "E480P1_MFoff": 1,
  "reference": 1, 
  "RK3CD2_ND2Gam2e-4": 3,  
  "RK4UD3": 3,     
  "RK6UD5": 3,       
  "RK8UD7": 3, 
}

OUTNC_suffix=""

#---

Nproc = 16
CpDry=1004.64
L = 9.6e3

def xaxis_txt(inv_lam, pos=None):
  l = str(int(1.0/inv_lam))
  return "  "+l+'$^{-1}$'

def xaxis_txt_minor(inv_lam, pos=None):
  l = str(int(1.0/inv_lam))
  if l=="20" or l=="50" or l=="200" or l=="500" or l=="2000" or l=="5000":
    return "  "+l+'$^{-1}$'
  else:
    return ""


def create_fig_energy_spectra_compari_slopem35(ke_spectra_list, ke_spectra_list_FDM, zlev, figname, slope_m35_ampl = 1.7e-5):
  print(f"create_fig_energy_spectra: {figname}")

  fig, ax = plt.subplots(figsize=(18,8))
  ax.set_xlim(1.0/5000.0,1.0/18.0)
  ax.set_ylim(4e-1, 1.4e0)
  ax.set_yscale('log')
  ax.set_xscale('log')

  ke_spectra_ref = ke_spectra_list["reference"].sel(z=zlev)
  for exp_name in ke_spectra_list.keys(): 
    ke_spectra = ke_spectra_list[exp_name].sel(z=zlev)
    slope_m35 = slope_m35_ampl*ke_spectra.k**(-5.0/3.0) 
    ax.plot(ke_spectra.k/(2.0*np.pi), ke_spectra/slope_m35, 
      linestyle=exp_ltype_list[exp_name], color=exp_color_list[exp_name], label=exp_label_list[exp_name])

  for exp_name in ke_spectra_list_FDM.keys(): 
    ke_spectra = ke_spectra_list_FDM[exp_name].sel(z=zlev)
    slope_m35 = slope_m35_ampl*ke_spectra.k**(-5.0/3.0) 
    ax.plot(ke_spectra.k/(2.0*np.pi), ke_spectra/slope_m35, 
      linestyle=exp_ltype_list[exp_name], color=exp_color_list[exp_name], label=exp_label_list[exp_name])

  # kx = ke_spectra_list[exp_name_list[0]].k
  # ax.plot(kx/(2.0*np.pi), ke_spectra_ref/slope_m35, 
  #         linestyle=exp_ltype_list["reference"], color=exp_color_list["reference"], label=exp_label_list["reference"])

  kx = ke_spectra_list[exp_name_list[0]].k
  ax.plot(kx/(2.0*np.pi), slope_m35/slope_m35, color="grey", linestyle="-.", label="-5/3")

  ax.set_xlabel("inverse of wavelength [m$^{-1}$]", fontsize=20)
  ax.xaxis.set_major_formatter(FuncFormatter(xaxis_txt))
  ax.xaxis.set_minor_formatter(FuncFormatter(xaxis_txt_minor))
  ax.tick_params(which="both", labelsize=18, length=3)

  ax.legend(loc='lower left')
  plt.savefig(f"{OUT_DIR}/{figname}")

def create_fig_energy_spectra_compari_ref(ke_spectra_list, ke_spectra_list_FDM, zlev, figname, slope_m35_ampl = 1.7e-5):
  print(f"create_fig_energy_spectra: {figname}")

  fig, ax = plt.subplots(figsize=(10,8))
  ax.set_xlim(1.0/5000.0,1.0/18.0)
  ax.set_ylim(.2e-1, 1.8e0)
  ax.set_yscale('log')
  ax.set_xscale('log')

  ke_spectra_ref = ke_spectra_list["reference"].sel(z=zlev)
  for exp_name in ke_spectra_list.keys(): 
    ke_spectra = ke_spectra_list[exp_name].sel(z=zlev)
    slope_m35 = slope_m35_ampl*ke_spectra.k**(-5.0/3.0) 
    if exp_name=="reference":
      alpha=1.0
    else:
      alpha=0.4
    ax.plot(ke_spectra.k/(2.0*np.pi), ke_spectra/ke_spectra_ref, 
      linestyle=exp_ltype_list[exp_name], color=exp_color_list[exp_name], label=exp_label_list[exp_name], 
      linewidth=exp_ltype_width[exp_name], alpha=alpha )

  for exp_name in ke_spectra_list_FDM.keys(): 
    ke_spectra = ke_spectra_list_FDM[exp_name].sel(z=zlev)
    slope_m35 = slope_m35_ampl*ke_spectra.k**(-5.0/3.0) 
    ax.plot(ke_spectra.k/(2.0*np.pi), ke_spectra/ke_spectra_ref, 
      linestyle=exp_ltype_list[exp_name], color=exp_color_list[exp_name], label=exp_label_list[exp_name], 
      linewidth=exp_ltype_width[exp_name] )
    
  # kx = ke_spectra_list[exp_name_list[0]].k
  # ax.plot(kx/(2.0*np.pi), ke_spectra_ref/ke_spectra_ref, 
  #         linestyle=exp_ltype_list["reference"], color=exp_color_list["reference"], label=exp_label_list["reference"])

  kx = ke_spectra_list[exp_name_list[0]].k
  ax.plot(kx/(2.0*np.pi), slope_m35/ke_spectra_ref, color="grey", linestyle="-.", label="-5/3")

  ax.set_xlabel("inverse of wavelength [m$^{-1}$]", fontsize=20)
  ax.xaxis.set_major_formatter(FuncFormatter(xaxis_txt))
  ax.xaxis.set_minor_formatter(FuncFormatter(xaxis_txt_minor))
  ax.tick_params(which="both", labelsize=18, length=3)

  ax.legend(loc='lower left')
  plt.savefig(f"{OUT_DIR}/{figname}")

def create_fig_energy_spectra_compari_ref_zoom(ke_spectra_list, ke_spectra_list_FDM, zlev, figname, slope_m35_ampl = 1.7e-5):
  print(f"create_fig_energy_spectra: {figname}")

  fig, ax = plt.subplots(figsize=(10,8))
  ax.set_xlim(1.0/400.0,1.0/18.0)
  ax.set_ylim(2e-1, 1.4e0)
  ax.set_yscale('log')
  ax.set_xscale('log')

  ke_spectra_ref = ke_spectra_list["reference"].sel(z=zlev)
  for exp_name in ke_spectra_list.keys(): 
    ke_spectra = ke_spectra_list[exp_name].sel(z=zlev)
    slope_m35 = slope_m35_ampl*ke_spectra.k**(-5.0/3.0) 
    if exp_name=="reference":
      alpha=1.0
    else:
      alpha=0.4
    ax.plot(ke_spectra.k/(2.0*np.pi), ke_spectra/ke_spectra_ref, 
      linestyle=exp_ltype_list[exp_name], color=exp_color_list[exp_name], label=exp_label_list[exp_name], 
      linewidth=exp_ltype_width[exp_name], alpha=alpha )

  for exp_name in ke_spectra_list_FDM.keys(): 
    ke_spectra = ke_spectra_list_FDM[exp_name].sel(z=zlev)
    slope_m35 = slope_m35_ampl*ke_spectra.k**(-5.0/3.0) 
    ax.plot(ke_spectra.k/(2.0*np.pi), ke_spectra/ke_spectra_ref, 
      linestyle=exp_ltype_list[exp_name], color=exp_color_list[exp_name], label=exp_label_list[exp_name], 
      linewidth=exp_ltype_width[exp_name] )
    
  # kx = ke_spectra_list[exp_name_list[0]].k
  # ax.plot(kx/(2.0*np.pi), ke_spectra_ref/ke_spectra_ref, 
  #         linestyle=exp_ltype_list["reference"], color=exp_color_list["reference"], label=exp_label_list["reference"])

  kx = ke_spectra_list[exp_name_list[0]].k
  ax.plot(kx/(2.0*np.pi), slope_m35/ke_spectra_ref, color="grey", linestyle="-.", label="-5/3")

  ax.set_xlabel("inverse of wavelength [m$^{-1}$]", fontsize=20)
  ax.xaxis.set_major_formatter(FuncFormatter(xaxis_txt))
  ax.xaxis.set_minor_formatter(FuncFormatter(xaxis_txt_minor))
  ax.tick_params(which="both", labelsize=14, length=3)

  ax.legend(loc='lower left')
  plt.savefig(f"{OUT_DIR}/{figname}")
  
def read_spectra_data(exp_name, tmp_dir, ke_spectra_hvel_list, ke_spectra_wvel_list, ke_spectra_3dvel_list):
  ke_spectra_hvel = xr.open_mfdataset(f'{tmp_dir}/KE_spectra_momh{OUTNC_suffix}.nc', decode_times=False, combine='by_coords').KE_spectra_momh
  ke_spectra_hvel_list[exp_name] = np.mean(ke_spectra_hvel[:,:,1:Ncut], axis=0)

  ke_spectra_wvel = xr.open_mfdataset(f'{tmp_dir}/KE_spectra_momz{OUTNC_suffix}.nc', decode_times=False, combine='by_coords').KE_spectra_momz
  ke_spectra_wvel_list[exp_name] = np.mean(ke_spectra_wvel[:,:,1:Ncut], axis=0)

  ke_spectra_3dvel_list[exp_name] = ke_spectra_hvel_list[exp_name] + ke_spectra_wvel_list[exp_name]
  return

#-----------------------------
os.makedirs(OUT_DIR, exist_ok=True)


Ncut = int(N)

ke_spectra_hvel_list = {}
ke_spectra_wvel_list = {}
ke_spectra_3dvel_list = {}

ke_spectra_hvel_list_FDM = {}
ke_spectra_wvel_list_FDM = {}
ke_spectra_3dvel_list_FDM = {}

for exp_name in exp_name_list:
  tmp_dir = f"tmp_data/tmp_{exp_name}"
  print(f"tmp_dir={tmp_dir}")

  read_spectra_data( exp_name, tmp_dir, ke_spectra_hvel_list, ke_spectra_wvel_list, ke_spectra_3dvel_list )


for exp_name in exp_name_list_FDM:
  tmp_dir = f"KT2021_SCALERM/tmp3_{exp_name}"
  print(f"tmp_dir={tmp_dir}")

  read_spectra_data( exp_name, tmp_dir, ke_spectra_hvel_list_FDM, ke_spectra_wvel_list_FDM, ke_spectra_3dvel_list_FDM )

tmp_dir_FDM = f"KT2021_SCALERM/tmp3_RK8CD8_ND8Gam2e-4"
read_spectra_data( "reference", tmp_dir_FDM, ke_spectra_hvel_list, ke_spectra_wvel_list, ke_spectra_3dvel_list )

for zind, zlev in enumerate(ZLEVEL_list):
  create_fig_energy_spectra_compari_slopem35(ke_spectra_3dvel_list, ke_spectra_3dvel_list_FDM, zlev, f"3dvel_spectra_z{zlev}m_FDM_compari_slopem35.png")  
  create_fig_energy_spectra_compari_ref(ke_spectra_3dvel_list, ke_spectra_3dvel_list_FDM, zlev, f"3dvel_spectra_z{zlev}m_FDM_compari_CD8ND8.png")  
  create_fig_energy_spectra_compari_ref_zoom(ke_spectra_3dvel_list, ke_spectra_3dvel_list_FDM, zlev, f"3dvel_spectra_z{zlev}m_FDM_compari_CD8ND8_zoom.png")  
  
  plt.close()
