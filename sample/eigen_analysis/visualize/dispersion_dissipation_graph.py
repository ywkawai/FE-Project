import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

DATA_DIR = "data/"
OUT_DIR  = "visualize/"

# Stability limit based on Table 3 of Alhawwary and Wang (2018). 
CFL_max = {
  "RK3": [0.409, 0.209, 0.13, 0.089, 0.066],
  "RK4": [0.464, 0.235, 0.145, 0.1, 0.073]
}
CourantNumberFac = 0.5

def set_xax(ax):
  ax.set_xlim(0.0, np.pi)
  ax.set_xticks([0, 0.25*np.pi, 0.5*np.pi, 0.75*np.pi, np.pi])
  ax.set_xticklabels(["0", "$\pi/4$", "$\pi/2$", "$3\pi/4$", "$\pi$"])
  ax.grid(color='k', linestyle='dotted', linewidth=1)

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

def create_fig_sgscompari_porder_dep(porder_list, k_nondim, km_real_list, km_aimg_list, fname):
  fig = plt.figure(figsize=(10,4))

  ax1 = fig.add_subplot(1,2, 1, xlabel='K', ylabel='Re(Km)')
  ax1.set_title("dispersion")
  set_xax(ax1)
  ax1.plot(k_nondim, k_nondim, 'k--', label="exact")
  for m in porder_list:
    ax1.plot(k_nondim, km_real_list[m], label=f"p={m}")

  ax2 = fig.add_subplot(1,2, 2, xlabel='K', ylabel='he/(a*(p+1)*Im(Km))')
  ax2.set_title("decay time scale [s]")
  set_xax(ax2)
  ax2.set_yscale("log")
  ax2.set_xlim(np.pi/16.0, np.pi)
  ax2.set_ylim(5e-1, 1e5)

  he_ov_pplus1 = 200/8.0#10.0
  a = 5.0
  Cs = 0.13
  Dsgs=2.0*he_ov_pplus1
  sabs = 0.15 * (np.pi/Dsgs)**(2.0/3.0)
  for m in porder_list:
    #print(-he_ov_pplus1/(a * (1e-15 + km_aimg_list[m])))
    ax2.plot(k_nondim, -he_ov_pplus1/(a * (1e-15 + km_aimg_list[m])), label=f"p={m}")
  ax2.plot(k_nondim, 1.0/( (k_nondim/he_ov_pplus1)**2 * (Cs*Dsgs)**2*sabs ), label=f"SGS time scale")
  ax2.plot(k_nondim, 10.0/( (k_nondim/he_ov_pplus1)**2 * (Cs*Dsgs)**2*sabs ), label=f"10xSGS time scale")
  ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)

  plt.subplots_adjust(wspace=0.4)
  plt.tight_layout()  
  plt.savefig(fname)

def create_fig_amplification_phase_error(k_nondim,  
  km_real, km_aimg, amplifiFac_list, phaseError_list, 
  porder, tscheme_name, fname):
  fig = plt.figure(figsize=(10,4))

  tlev_list = [1, 10, 100]
  color_list = ["black", "blue", "cyan"]

  courant_num = CFL_max[tscheme_name][porder-1] * CourantNumberFac 
  ax1 = fig.add_subplot(1,2, 1, xlabel='K', ylabel='Re(Km)')
  ax1.set_yscale("log")
  ax1.set_title("Dispersion error")
  ax1.set_ylim(1e-10, 1e2)
  set_xax(ax1)
  for i in range(0,len(tlev_list)):
    n = tlev_list[i]; color = color_list[i]
    ax1.plot(k_nondim, n*courant_num*(np.abs(km_real-k_nondim)), color=color, linestyle=":", label=f"n={n} (physical mode)")
    ax1.plot(k_nondim, np.abs(phaseError_list[:,i]), color=color, linestyle="-", label=f"n={n} (combined mode)")    

  ax2 = fig.add_subplot(1,2, 2, xlabel='K', ylabel='G')
  #ax2.plot(k_nondim, np.ones(len(k_nondim)), 'k', label="exact")
  ax2.set_title("Amplification factor")
  ax2.set_ylim(1e-12, 1e0)  
  ax2.set_yscale("log")  
  set_xax(ax2)
  for i in range(0,len(tlev_list)):
    n = tlev_list[i]; color = color_list[i]
    ax2.plot(k_nondim, np.exp(n*(porder+1.0)*courant_num*km_aimg), color=color, linestyle=":", label=f"n={n} (physical mode)")
    ax2.plot(k_nondim, amplifiFac_list[:,i], color=color,  linestyle="-", label=f"n={n} (combined mode)")    

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

    # create_fig_disp_diff(k_nondim, km_real_list, km_aimg_list, f"{OUT_DIR}/P{porder:02}_{basis_type}_{form_type}.png")

# Create figure to show the dependence of polynomial order on the numerical dispersion/dissipation
for basis_type in ['modal', 'nodal']:
  form_type='weak'
  file_prefix=f"{DATA_DIR}/eigenval_{basis_type}_{form_type}_"
  print(f"graph porderdep: {file_prefix} ..")

  porder_list = [2, 3, 4, 5, 6, 7]
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

  # create_fig_sgscompari_porder_dep( porder_list, k_nondim, km_real_list_pysmode, km_aimg_list_pysmode, 
  #  f"{OUT_DIR}/PolyOrderDep_{basis_type}_{form_type}_sgscompari.png")

# Create figure of full-discrete analysis for p=2
for tscheme_name in ['RK3']:
  for basis_type in ['modal']:
    form_type='weak'
    file_prefix=f"{DATA_DIR}/eigenval_{basis_type}_{form_type}_"

    print(f"graph porderdep: {file_prefix} ..")

    porder_list = [2]

    for porder in porder_list:
      k_nondim = np.loadtxt(f"{file_prefix}P{porder:02}_Mode01_{tscheme_name}.dat", skiprows=1, usecols=0)
      amplifiFac_list = np.loadtxt(f"{file_prefix}P{porder:02}_combinedmode_{tscheme_name}.dat", skiprows=1, usecols=[1,3,5])
      phaseError_list = np.loadtxt(f"{file_prefix}P{porder:02}_combinedmode_{tscheme_name}.dat", skiprows=1, usecols=[2,4,6])

      for m in range(1,porder+2):
        km_real = np.loadtxt(f"{file_prefix}P{porder:02}_Mode{m:02}_{tscheme_name}.dat", skiprows=1, usecols=1)
        km_aimg = np.loadtxt(f"{file_prefix}P{porder:02}_Mode{m:02}_{tscheme_name}.dat", skiprows=1, usecols=2)

        if ( abs(km_real[0]) < 1.0E-10 and km_real[1]-km_real[0] > 0.0 ):  
          k_nondim = np.loadtxt(f"{file_prefix}P{porder_list[0]:02}_Mode01_{tscheme_name}.dat", skiprows=1, usecols=0)
          create_fig_amplification_phase_error( 
            k_nondim, km_real, km_aimg, amplifiFac_list, phaseError_list, 
            porder, tscheme_name, f"{OUT_DIR}/P{porder:02}_{basis_type}_{form_type}_{tscheme_name}_G_phi.png")
