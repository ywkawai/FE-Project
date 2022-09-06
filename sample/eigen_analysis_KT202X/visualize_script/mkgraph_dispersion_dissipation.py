import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import os

TARGET_DIR = "stdupwind_semidiscrete/"
DATA_DIR   = "data/"
VIS_DIR    = "vis/dispersion_dissipation/"
BASIS_TYPE = "nodal"
FORM_TYPE  = 'strong'
TSCHEME_NAME = ""
MF_SUFFIX    = ""
PMAX         = 8


#---
file_prefix=f"{TARGET_DIR}/{DATA_DIR}/eigenval_{BASIS_TYPE}_{FORM_TYPE}_"
OUT_DIR=f"{TARGET_DIR}/{VIS_DIR}"

secondary_mode_color_list = ["red", "orange", "green", "magenta", "blue", "cyan", "gold"]


def set_xax(ax):
  ax.set_xlim(0.0, np.pi)
  ax.set_xticks([0, 0.25*np.pi, 0.5*np.pi, 0.75*np.pi, np.pi])
  ax.set_xticklabels(["0", "$\pi/4$", "$\pi/2$", "$3\pi/4$", "$\pi$"])
  ax.grid(color='k', linestyle='dotted', linewidth=1)


def create_fig_disp_diff(k_nondim, km_real_list, km_aimg_list, fname):
  fig = plt.figure(figsize=(10,4))
  primary_mode_flag = [False]*len(km_real_list)

  ax1 = fig.add_subplot(1,2, 1, xlabel='K', ylabel='Re(Km)')
  ax1.set_title("Numerical dispersion")
  set_xax(ax1)
  ax1.set_ylim(-3.2, 3.2)  
  ax1.plot(k_nondim, k_nondim, 'k-.', label="exact", linewidth=2)
  i=0
  for m in range(0,len(km_real_list)):
    km_real = km_real_list[m]
    if km_real[1] <  0.01 and km_real[1] > 0.0: 
      ax1.plot(k_nondim, km_real, label=f"mode{m}", color="black", linewidth=2)
      primary_mode_flag[m] = True 
    else:
      ax1.plot(k_nondim, km_real, label=f"mode{m}", color=secondary_mode_color_list[i], alpha=0.4)
      i = i + 1

  ax2 = fig.add_subplot(1,2, 2, xlabel='K', ylabel='Im(Km)')
  ax2.set_title("Numerical dissipation")
  ax2.plot(k_nondim, 0.0*k_nondim, 'k-.', label="exact", linewidth=2)
  set_xax(ax2)
  ax2.set_ylim(-5.1, 0.2)
#  ax2.set_ylim(-400, 0.2)
#  ax2.set_ylim(-0.3, 0.05)

  i=0
  for m in range(0,len(km_aimg_list)):
    if primary_mode_flag[m]:
      ax2.plot(k_nondim, km_aimg_list[m], label=f"mode{m}", color="black", linewidth=2)
    else:
      ax2.plot(k_nondim, km_aimg_list[m], label=f"mode{m}", color=secondary_mode_color_list[i], alpha=0.4)
      i = i + 1

  ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)

  plt.subplots_adjust(wspace=0.4)
  plt.tight_layout()  
  plt.savefig(fname)


# Create figure for each order of polynomial

os.makedirs(OUT_DIR, exist_ok=True)
for porder in range(1,PMAX):
  dat_suffix = ""
  if len(TSCHEME_NAME) > 0: 
    dat_suffix = f"{dat_suffix}_{TSCHEME_NAME}"
  if len(MF_SUFFIX) > 0: 
    dat_suffix = f"{dat_suffix}_{MF_SUFFIX}"

  print(f"graph for each p: {file_prefix}P{porder:02}_Mode0*{dat_suffix} ..")

  k_nondim = np.loadtxt(f"{file_prefix}P{porder:02}_Mode01{dat_suffix}.dat", skiprows=1, usecols=0)
  km_real_list = {}; km_aimg_list = {}
  for m in range(1,porder+2):
    km_real_list[m-1] = np.loadtxt(f"{file_prefix}P{porder:02}_Mode{m:02}{dat_suffix}.dat", skiprows=1, usecols=1)
    km_aimg_list[m-1] = np.loadtxt(f"{file_prefix}P{porder:02}_Mode{m:02}{dat_suffix}.dat", skiprows=1, usecols=2)

  create_fig_disp_diff(k_nondim, km_real_list, km_aimg_list, f"{OUT_DIR}/P{porder:02}_{BASIS_TYPE}_{FORM_TYPE}{dat_suffix}.png")
