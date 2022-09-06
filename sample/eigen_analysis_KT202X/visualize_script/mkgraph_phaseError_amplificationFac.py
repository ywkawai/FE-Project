import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os

TARGET_DIR = "stdupwind_semidiscrete/"
DATA_DIR   = "data/"
VIS_DIR    = "vis/phaseError_amplfac/"
BASIS_TYPE = "nodal"
FORM_TYPE  = 'strong'
TSCHEME_NAME = ""
MF_SUFFIX    = ""
PMAX         = 8
Cr_DG_dxeff  = 5.0 * 0.0125 / 10.0
TLEV_list = [10, 100, 1000]

#---
file_prefix=f"{TARGET_DIR}/{DATA_DIR}/eigenval_{BASIS_TYPE}_{FORM_TYPE}_"
OUT_DIR=f"{TARGET_DIR}/{VIS_DIR}"

def set_xax(ax):
  ax.set_xlim(0.0, np.pi)
  ax.set_xticks([0, 0.25*np.pi, 0.5*np.pi, 0.75*np.pi, np.pi])
  ax.set_xticklabels(["0", "$\pi/4$", "$\pi/2$", "$3\pi/4$", "$\pi$"])
  ax.grid(color='k', linestyle='dotted', linewidth=1)


def create_fig_amplification_phase_error(cr_dg_dxeff, k_nondim,  
  km_real, km_aimg, amplifiFac_list, phaseError_list, 
  porder, fname):
  fig = plt.figure(figsize=(12,4))

  color_list = ["black", "blue", "cyan", "red"]

  ax1 = fig.add_subplot(1,2,1, xlabel='K', ylabel='$|\psi|$')
  ax1.set_yscale("log")
  ax1.set_title("Dispersion error")
  ax1.set_ylim(1e-10, 1e2)
  set_xax(ax1)
  for i in range(0,len(TLEV_list)):
    n = TLEV_list[i]; color = color_list[i]
    ax1.plot(k_nondim, n * cr_dg_dxeff * (np.abs(km_real-k_nondim)), color=color, linestyle=":", label=f"n={n} (physical mode)")
    ax1.plot(k_nondim, np.abs(phaseError_list[:,i]), color=color, linestyle="-", label=f"n={n} (combined mode)")    

  ax2 = fig.add_subplot(1,2, 2, xlabel='K', ylabel='G')
  ax2.set_title("Amplification factor")
  ax2.set_ylim(1e-12, 1e0)  
  #ax2.set_yscale("log")  
  set_xax(ax2)
  for i in range(0,len(TLEV_list)):
    n = TLEV_list[i]; color = color_list[i]
    ax2.plot(k_nondim, np.exp(n*cr_dg_dxeff*km_aimg), color=color, linestyle=":", label=f"n={n} (physical mode)")
    ax2.plot(k_nondim, amplifiFac_list[:,i], color=color,  linestyle="-", label=f"n={n} (combined mode)")    

  ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)

  plt.subplots_adjust(wspace=0.4)
  plt.tight_layout()  
  plt.savefig(fname)


#-------------

print(f"graph porderdep: {file_prefix} ..")
os.makedirs(OUT_DIR, exist_ok=True)

for porder in range(1, PMAX):
 
  dat_suffix = ""
  if len(TSCHEME_NAME) > 0: 
    dat_suffix = f"{dat_suffix}_{TSCHEME_NAME}"
  if len(MF_SUFFIX) > 0: 
    dat_suffix = f"{dat_suffix}_{MF_SUFFIX}"
  

  print(f"{file_prefix}P{porder:02}_combinedmode{dat_suffix}.dat")
  k_nondim = np.loadtxt(f"{file_prefix}P{porder:02}_Mode01.dat", skiprows=1, usecols=0)
  amplifiFac_list = np.loadtxt(f"{file_prefix}P{porder:02}_combinedmode{dat_suffix}.dat", skiprows=1, usecols=[3,7,13])
  phaseError_list = np.loadtxt(f"{file_prefix}P{porder:02}_combinedmode{dat_suffix}.dat", skiprows=1, usecols=[4,8,14])

  for m in range(1,porder+2):
    km_real = np.loadtxt(f"{file_prefix}P{porder:02}_Mode{m:02}{dat_suffix}.dat", skiprows=1, usecols=1)
    km_aimg = np.loadtxt(f"{file_prefix}P{porder:02}_Mode{m:02}{dat_suffix}.dat", skiprows=1, usecols=2)

    if ( abs(km_real[0]) < 1.0E-10 and km_real[1]-km_real[0] > 0.0 ):  
      k_nondim = np.loadtxt(f"{file_prefix}P01_Mode01{dat_suffix}.dat", skiprows=1, usecols=0)
      create_fig_amplification_phase_error( Cr_DG_dxeff, 
        k_nondim, km_real, km_aimg, amplifiFac_list, phaseError_list, 
        porder, f"{OUT_DIR}/P{porder:02}_{BASIS_TYPE}_{FORM_TYPE}_G_phi{dat_suffix}.png")
