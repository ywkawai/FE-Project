import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

DATA_DIR = "data/"
OUT_DIR  = "visualize/"

def set_xax(ax):
  ax.set_xlim(0.0, np.pi)
  ax.set_xticks([0, 0.25*np.pi, 0.5*np.pi, 0.75*np.pi, np.pi])
  ax.set_xticklabels(["0", "$\pi/4$", "$\pi/2$", "$3\pi/4$", "$\pi$"])

def create_fig_disp_diff(k_nondim, km_real_list, km_aimg_list, fname):
  fig = plt.figure(figsize=(10,4))

  ax1 = fig.add_subplot(1,2, 1, xlabel='K', ylabel='Re(Km)')
  ax1.set_title("dispersion")
  set_xax(ax1)
  ax1.plot(k_nondim, k_nondim, 'k--', label="exact")
  for m in range(0,len(km_real_list)):
    ax1.plot(k_nondim, km_real_list[m], label=f"mode{m}")

  ax2 = fig.add_subplot(1,2, 2, xlabel='K', ylabel='Im(Km)')
  ax2.set_title("dissipation")
  ax2.plot(k_nondim, 0.0*k_nondim, 'k--', label="exact")
  set_xax(ax2)
  for m in range(0,len(km_aimg_list)):
    ax2.plot(k_nondim, km_aimg_list[m], label=f"mode{m}")
  ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)

  plt.subplots_adjust(wspace=0.4)
  plt.tight_layout()  
  plt.savefig(fname)

def create_fig_disp_diff_porder_dep(porder_list, k_nondim, km_real_list, km_aimg_list, fname):
  fig = plt.figure(figsize=(10,4))

  ax1 = fig.add_subplot(1,2, 1, xlabel='K', ylabel='Re(Km)')
  ax1.set_title("dispersion")
  set_xax(ax1)
  ax1.plot(k_nondim, k_nondim, 'k--', label="exact")
  for m in porder_list:
    ax1.plot(k_nondim, km_real_list[m], label=f"p={m}")

  ax2 = fig.add_subplot(1,2, 2, xlabel='K', ylabel='Im(Km)')
  ax2.plot(k_nondim, 0.0*k_nondim, 'k--', label="exact")
  ax2.set_title("dissipation")
  set_xax(ax2)
  for m in porder_list:
    ax2.plot(k_nondim, km_aimg_list[m], label=f"p={m}")
  ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)

  plt.subplots_adjust(wspace=0.4)
  plt.tight_layout()  
  plt.savefig(fname)

# Create figure for each order of polynomial
for basis_type in ['modal', 'nodal']:
  form_type='weak'
  file_prefix=f"{DATA_DIR}/eigenval_{basis_type}_{form_type}_"
  print(f"graph for each p: {file_prefix} ..")

  for porder in range(1,8):
    k_nondim = np.loadtxt(f"{file_prefix}P{porder:02}_Mode01.dat", skiprows=1, usecols=0)
    km_real_list = {}; km_aimg_list = {}
    for m in range(1,porder+2):
      km_real_list[m-1] = np.loadtxt(f"{file_prefix}P{porder:02}_Mode{m:02}.dat", skiprows=1, usecols=1)
      km_aimg_list[m-1] = np.loadtxt(f"{file_prefix}P{porder:02}_Mode{m:02}.dat", skiprows=1, usecols=2)

    create_fig_disp_diff(k_nondim, km_real_list, km_aimg_list, f"{OUT_DIR}/P{porder:02}_{basis_type}_{form_type}.png")

# Create figure to show the dependence of polynomial order on the numerical dispersion/dissipation
for basis_type in ['modal', 'nodal']:
  form_type='weak'
  file_prefix=f"{DATA_DIR}/eigenval_{basis_type}_{form_type}_"
  print(f"graph porderdep: {file_prefix} ..")

  porder_list = [2, 4, 3, 4, 5]
  km_real_list_pysmode = {}; km_aimg_list_pysmode = {}
  for porder in porder_list:
    for m in range(1,porder+2):
      km_real = np.loadtxt(f"{file_prefix}P{porder:02}_Mode{m:02}.dat", skiprows=1, usecols=1)
      km_aimg = np.loadtxt(f"{file_prefix}P{porder:02}_Mode{m:02}.dat", skiprows=1, usecols=2)
      if ( abs(km_real[0]) < 1.0E-10 and km_real[1]-km_real[0] > 0.0 ):
        km_real_list_pysmode[porder] = km_real 
        km_aimg_list_pysmode[porder] = km_aimg

  k_nondim = np.loadtxt(f"{file_prefix}P{porder_list[0]:02}_Mode01.dat", skiprows=1, usecols=0)
  create_fig_disp_diff_porder_dep( porder_list, k_nondim, km_real_list_pysmode, km_aimg_list_pysmode, 
    f"{OUT_DIR}/PolyOrderDep_{basis_type}_{form_type}.png")