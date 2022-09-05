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

#--
exp_name_list = [
  "E120P7",   
  "E120P7_MF32OrdAlphx10",   
  "E192P4", 
  "E192P4_MF32OrdAlphx10",     
  "E240P3", 
  "E240P3_MF32OrdAlphx10",   
  "E480P1",     
  "E480P1_MF32OrdAlphx10",   
  "E480P1_MFoff", 
]
OUT_DIR="./compari_modalfilter_coefdep/"


exp_label_list = {
  "E120P7": "P7", 
  "E120P7_MF32OrdAlphx10": r"P7_$p_m$32_$\alpha_m$x10", 
  "E192P4": "P4",   
  "E192P4_MF32OrdAlphx10": r"P4_$p_m$32_$\alpha_m$x10", 
  "E240P3": "P3", 
  "E240P3_MF32OrdAlphx10": r"P3_$p_m$32_$\alpha_m$x10", 
  "E480P1": "P1", 
  "E480P1_MF32OrdAlphx10": r"P1_$p_m$32_$\alpha_m$x10", 
  "E480P1_MFoff": "P1_MFoff", 
  "reference": "reference (CD8ND8)", 
}
exp_color_list = {
  "E120P7": "goldenrod", 
  "E120P7_MF32OrdAlphx10": "goldenrod", 
  "E192P4": "red", #"forestgreen",  
  "E192P4_MF32OrdAlphx10": "red", 
  "E240P3": "blue", 
  "E240P3_MF32OrdAlphx10": "blue", 
#  "E480P1": "lightcyan", 
  "E480P1": "cyan",   
  "E480P1_MF32OrdAlphx10": "cyan", 
  "E480P1_MFoff": "lightskyblue", 
  "reference": "black",   
}
exp_ltype_list = {
  "E120P7": "-",
  "E120P7_MF32OrdAlphx10": "--",
  "E192P4": "-",  
  "E192P4_MF32OrdAlphx10": "--",
  "E240P3": "-",
  "E240P3_MF32OrdAlphx10": "--",
  "E480P1": "-",
  "E480P1_MF32OrdAlphx10": "--",
  "E480P1_MFoff": "-.",
  "reference": "-",   
}
exp_ltype_width = {
  "E120P7": 2,
  "E120P7_MF32OrdAlphx10": 1,
  "E192P4": 2,  
  "E192P4_MF32OrdAlphx10": 1,
  "E240P3": 2,
  "E240P3_MF32OrdAlphx10": 1,
  "E480P1": 2,
  "E480P1_MF32OrdAlphx10": 1,
  "E480P1_MFoff": 2,
  "reference": 1,   
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

def create_fig_energy_spectra_zoom(ke_spectra_list, zlev, figname, slope_m35_ampl = 1.7e-5):
  print(f"create_fig_energy_spectra: {figname}")

  fig, ax = plt.subplots(figsize=(5,10))
#pcm = plt.pcolormesh(kx[1:Ncut], ky[1:Ncut], np.log10(s_f_abs[1:Ncut,1:Ncut]), vmin=-6, vmax=6)
#plt.colorbar(pcm)
  ax.set_xlim(1.0/100.0,1.0/18.0)
  ax.set_ylim(4e-7, 2e-3)
  #ax.set_ylim(1e-3, 1.5e-2)
  ax.set_yscale('log')
  ax.set_xscale('log')

  for exp_name in ke_spectra_list.keys(): 
    ke_spectra = ke_spectra_list[exp_name].sel(z=zlev)
    ax.plot(ke_spectra.k/(2.0*np.pi), ke_spectra, 
            linestyle=exp_ltype_list[exp_name], 
            linewidth=exp_ltype_width[exp_name], 
            color=exp_color_list[exp_name], 
            label=exp_label_list[exp_name] )

  kx = ke_spectra_list["E120P7"].k
  slope_m35 = slope_m35_ampl * kx**(-5.0/3.0) 
  ax.plot(kx/(2.0*np.pi), slope_m35, color="grey", linestyle="-.", label="-5/3", linewidth=2)
  ax.set_xlabel("inverse of wavelength [m$^{-1}$]", fontsize=15)
  ax.set_ylabel("E(k)", fontsize=15)
  ax.xaxis.set_major_formatter(FuncFormatter(xaxis_txt))
  ax.xaxis.set_minor_formatter(FuncFormatter(xaxis_txt_minor))
  ax.tick_params(which="both", labelsize=15)

  ax.legend(loc='lower left')
  plt.savefig(f"{OUT_DIR}/{figname}")

def create_fig_energy_spectra(ke_spectra_list, zlev, figname, slope_m35_ampl = 1.7e-5):
  print(f"create_fig_energy_spectra: {figname}")

  fig, ax = plt.subplots(figsize=(8,6))
#pcm = plt.pcolormesh(kx[1:Ncut], ky[1:Ncut], np.log10(s_f_abs[1:Ncut,1:Ncut]), vmin=-6, vmax=6)
#plt.colorbar(pcm)
  ax.set_xlim(1.0/5000.0,1.0/18.0)
  #ax.set_ylim(.7e-6, 5e-1)
  ax.set_ylim(3e-8, 2e-1)
  ax.set_yscale('log')
  ax.set_xscale('log')

  for exp_name in ke_spectra_list.keys(): 
    ke_spectra = ke_spectra_list[exp_name].sel(z=zlev)
    ax.plot(ke_spectra.k/(2.0*np.pi), ke_spectra, 
            linestyle=exp_ltype_list[exp_name], 
            linewidth=exp_ltype_width[exp_name], 
            color=exp_color_list[exp_name], 
            label=exp_label_list[exp_name])

#  kx = ke_spectra_list["reference"].k
  kx = ke_spectra_list[exp_name_list[0]].k
  slope_m35 = slope_m35_ampl * kx**(-5.0/3.0) 
  ax.plot(kx/(2.0*np.pi), slope_m35, color="grey", linestyle="-.", label="-5/3", linewidth=2)
  ax.set_xlabel("inverse of wavelength [m$^{-1}$]", fontsize=15)
#  ax.set_ylabel("E(k)", fontsize=14)
  ax.xaxis.set_major_formatter(FuncFormatter(xaxis_txt))
  ax.xaxis.set_minor_formatter(FuncFormatter(xaxis_txt_minor))
  ax.tick_params(which="both", labelsize=13)

  ax.legend(loc='lower left')
  plt.savefig(f"{OUT_DIR}/{figname}")


def create_fig_energy_spectra_diffm53(ke_spectra_list, zlev, figname, slope_m35_ampl = 1.7e-5):
  print(f"create_fig_energy_spectra: {figname}")

  fig, ax = plt.subplots(figsize=(10,7))
  ax.set_xlim(1.0/5000.0,1.0/18.0)
  ax.set_ylim(4e-1, 1.4e0)
  ax.set_yscale('log')
  ax.set_xscale('log')

  for exp_name in ke_spectra_list.keys(): 
    ke_spectra = ke_spectra_list[exp_name].sel(z=zlev)
    slope_m35 = slope_m35_ampl*ke_spectra.k**(-5.0/3.0) 
    ax.plot(ke_spectra.k/(2.0*np.pi), ke_spectra/slope_m35, 
      linestyle=exp_ltype_list[exp_name], color=exp_color_list[exp_name], label=exp_label_list[exp_name])

  kx = ke_spectra_list[exp_name_list[0]].k
  ax.plot(kx/(2.0*np.pi), slope_m35/slope_m35, color="grey", linestyle="-.", label="-5/3")

  ax.set_xlabel("inverse of wavelength [m$^{-1}$]", fontsize=20)
  ax.xaxis.set_major_formatter(FuncFormatter(xaxis_txt))
  ax.xaxis.set_minor_formatter(FuncFormatter(xaxis_txt_minor))
  ax.tick_params(which="both", labelsize=18, length=3)

  plt.savefig(f"{OUT_DIR}/{figname}")
  
def create_fig_energy_spectra_compari_ref_zoom(ke_spectra_list, zlev, figname, slope_m35_ampl = 1.7e-5):
  print(f"create_fig_energy_spectra: {figname}")

  fig, ax = plt.subplots(figsize=(10,7))
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
    # if exp_name=="E480P1_MFoff" or exp_name=="E240P3" or exp_name=="E192P4" or exp_name=="E120P7":
    #   alpha=0.4 
    else:
      alpha=1.0
    ax.plot(ke_spectra.k/(2.0*np.pi), ke_spectra/ke_spectra_ref, 
      linestyle=exp_ltype_list[exp_name], color=exp_color_list[exp_name], label=exp_label_list[exp_name], 
      linewidth=exp_ltype_width[exp_name], alpha=alpha )

  kx = ke_spectra_list[exp_name_list[0]].k
  ax.plot(kx/(2.0*np.pi), slope_m35/ke_spectra_ref, color="grey", linestyle="-.", label="-5/3")

  ax.set_xlabel("inverse of wavelength [m$^{-1}$]", fontsize=20)
  ax.xaxis.set_major_formatter(FuncFormatter(xaxis_txt))
  ax.xaxis.set_minor_formatter(FuncFormatter(xaxis_txt_minor))
  ax.tick_params(which="both", labelsize=14, length=3)

#  ax.legend(loc='lower left')
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



for exp_name in exp_name_list:
  tmp_dir = f"tmp_data/tmp_{exp_name}"
  print(f"tmp_dir={tmp_dir}")

  read_spectra_data( exp_name, tmp_dir, ke_spectra_hvel_list, ke_spectra_wvel_list, ke_spectra_3dvel_list )

for zind, zlev in enumerate(ZLEVEL_list):
  create_fig_energy_spectra(ke_spectra_3dvel_list, zlev, f"3dvel_spectra_z{zlev}m.png", 1.7e-5) # 1.7e-5  
  create_fig_energy_spectra_zoom(ke_spectra_3dvel_list, zlev, f"3dvel_spectra_z{zlev}m_zoom.png")  
  create_fig_energy_spectra_diffm53(ke_spectra_3dvel_list, zlev, f"3dvel_spectra_z{zlev}m_diff.png")  
  
  plt.close()

#----
tmp_dir_FDM = f"KT2021_SCALERM/tmp3_RK8CD8_ND8Gam2e-4"
read_spectra_data( "reference", tmp_dir_FDM, ke_spectra_hvel_list, ke_spectra_wvel_list, ke_spectra_3dvel_list )

for zind, zlev in enumerate(ZLEVEL_list):
  create_fig_energy_spectra_compari_ref_zoom(ke_spectra_3dvel_list, zlev, f"3dvel_spectra_z{zlev}m_compari_CD8ND8_zoom.png")  
  plt.close()
