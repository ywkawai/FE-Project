import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors
import os


U_VEL        = 5.0

# TARGET_DIR = "stdupwind_semidiscrete/"; TSCHEME_NAME = ""; MF_SUFFIX = ""
# PORDER_LIST = [1, 3, 4, 5, 7, 11 ]
# Cr_DG_dxeff  = U_VEL * 0.0125 / 10.0

# TARGET_DIR = "stdupwind_MF32Alph1E-3_RK4/"; TSCHEME_NAME = "RK4"; MF_SUFFIX = "MFMF32Alph1E-3"
# TARGET_DIR = "stdupwind_MF32Alph1E-2_RK4/"; TSCHEME_NAME = "RK4"; MF_SUFFIX = "MFMF32Alph1E-2"
# TARGET_DIR = "stdupwind_MF16Alph1E-2_RK4/"; TSCHEME_NAME = "RK4"; MF_SUFFIX = "MFMF16Alph1E-2"
# TARGET_DIR = "stdupwind_MF08Alph1E-2_RK4/"; TSCHEME_NAME = "RK4"; MF_SUFFIX = "MFMF08Alph1E-2"
# TARGET_DIR = "stdupwind_MF32Alph1E-1_RK4/"; TSCHEME_NAME = "RK4"; MF_SUFFIX = "MFMF32Alph1E-1"
# PORDER_LIST = [1, 3, 4, 5, 7 ]
# Cr_DG_dxeff  = U_VEL * 0.0125 / 10.0

#TARGET_DIR = "stdupwind_RK3/"; TSCHEME_NAME = "RK3"; MF_SUFFIX = ""
# TARGET_DIR = "stdupwind_RK4/"; TSCHEME_NAME = "RK4"; MF_SUFFIX = ""
# TARGET_DIR = "stdupwind_RKo4s10/"; TSCHEME_NAME = "RKo4s10"; MF_SUFFIX = ""
# PORDER_LIST = [1, 3, 4, 7]
# PORDER_LIST = [9, 11]
# Cr_DG_dxeff  = U_VEL * 0.0125 / 10.0

# TARGET_DIR = "stdupwind_Crx50_semidiscrete/"; TSCHEME_NAME = ""; MF_SUFFIX = ""
#TARGET_DIR = "stdupwind_Crx50_RK3/"; TSCHEME_NAME = "RK3"; MF_SUFFIX = ""
#TARGET_DIR = "stdupwind_Crx50_RK4/"; TSCHEME_NAME = "RK4"; MF_SUFFIX = ""
# TARGET_DIR = "stdupwind_Crx50_RKo4s10/"; TSCHEME_NAME = "RKo4s10"; MF_SUFFIX = ""
# PORDER_LIST = [3, 4, 7]
# Cr_DG_dxeff  = 50.0 * U_VEL * 0.0125 / 10.0

# TARGET_DIR = "stdupwind_Crx190_RKo4s10/"; TSCHEME_NAME = "RKo4s10"; MF_SUFFIX = ""
# PORDER_LIST = [3, 4, 7]
# Cr_DG_dxeff  = 190.0 * U_VEL * 0.0125 / 10.0

# TARGET_DIR = "ovupwind_RKo4s10/"; TSCHEME_NAME = "RKo4s10"; MF_SUFFIX = ""
# TARGET_DIR = "ovupwind_MF32Alph1E-3_RKo4s10/"; TSCHEME_NAME = "RKo4s10"; MF_SUFFIX = "MFMF32Alph1E-3"
# PORDER_LIST = [1, 3, 4, 5, 7 ]
# Cr_DG_dxeff  = U_VEL * 0.0125 / 10.0

TARGET_DIR = "stdupwind_MF32Alph1E-3_RKo4s10/"; TSCHEME_NAME = "RKo4s10"; MF_SUFFIX = "MFMF32Alph1E-3"
PORDER_LIST = [2] #[1, 3, 4, 5, 7, 11 ]
Cr_DG_dxeff  = U_VEL * 0.0125 / 10.0

#---------
DATA_DIR   = "data/"
VIS_DIR    = "vis/Rdiff_Rdisp/"
BASIS_TYPE = "nodal"
FORM_TYPE  = 'strong'


# 1, 10, 50, 100, 200, 500, 1000, 3000 

TSLOT_list = [ 1, 10, 50, 100, 1000, 3000 ]
TSLOT_ary_ind = [0, 1, 2, 3, 6, 7]

Cs = 0.13
Eta = 0.15
Eta_dash = (2.0/15.0)**0.5 * Eta
m = 2.0

mask_min_value = 8e-9

#--------------------------------------

file_prefix=f"{TARGET_DIR}/{DATA_DIR}/eigenval_{BASIS_TYPE}_{FORM_TYPE}_"
OUT_DIR=f"{TARGET_DIR}/{VIS_DIR}"

dx_ = np.array([0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 
       12, 15, 17.5, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000,  
       1200, 1500, 1750, 2000, 2500,  3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]) #np.arange(0.5, 1000.1, 0.1)

#l_ = np.arange(2.0, 100.2, 0.2)

def set_xax(ax):
  ax.set_xlim(1.0, 1000.0)
  ax.set_xticks([0, 0.25*np.pi, 0.5*np.pi, 0.75*np.pi, np.pi])
  ax.set_xticklabels(["0", "$\pi/4$", "$\pi/2$", "$3\pi/4$", "$\pi$"])
  ax.grid(color='k', linestyle='dotted', linewidth=1)

def create_fig_rdiff_rdisp(dx, l, rdiff, rdisp, fname):
  
  fig = plt.figure(figsize=(12, 10))
  ax1 = fig.add_subplot(111)

  ax1.set_xlabel('$\Delta x_{eff}$ [m]', fontsize=30)
  ax1.set_ylabel('$l$ [grid]', fontsize=30)
  ax1.set_xscale('log')
  ax1.set_yscale('log')#, basey=2)
  ax1.tick_params(which="major", length=8, width=2)
  ax1.tick_params(which="minor", length=4, width=2)
  ax1.set_ylim(2.0, 100.0)

  rdisp_mask = rdisp.copy()
  rdisp_mask[rdisp < mask_min_value] = np.nan
  cnt_lv = [0.00000001, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000]
  cnt_opts = {'levels':cnt_lv, 'colors': ['black'], 'linewidths': 4, 'alpha':0.8, "linestyles":"dotted"}
  cnt = ax1.contour(dx, l, rdisp_mask, **cnt_opts)
  zc = cnt.collections[7]
  plt.setp(zc, linewidth=8)
  fmt = ticker.LogFormatterSciNotation()
  fmt.create_dummy_axis()
  lbl = cnt.clabel(fmt=fmt, fontsize=28)
  #lbl = cnt.clabel(fmt=fmt, fontsize=28, inline=True, manual=clabel_loc_Rdisp[scheme_type])

  rdiff_mask = rdiff.copy()
  rdiff_mask[rdiff < mask_min_value] = np.nan
  pcm = ax1.pcolor(dx, l, rdiff_mask, norm=colors.LogNorm(vmin=1e-2, vmax=2e2), cmap='YlGnBu')
  cnt_lv = [0.00000001, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000]
  cnt_opts = {'levels':cnt_lv, 'colors': ['red'], 'linewidths': 4, 'alpha':1.0}
  cnt = ax1.contour(dx, l, rdiff_mask, **cnt_opts)
  zc = cnt.collections[7]
  plt.setp(zc, linewidth=7)
  fmt = ticker.LogFormatterSciNotation()
  fmt.create_dummy_axis()
  lbl = cnt.clabel(fmt=fmt, fontsize=32) 
  #lbl = cnt.clabel(fmt=fmt, fontsize=32, inline=True, manual=clabel_loc_Rdiff[scheme_type])

  cbar = plt.colorbar(pcm, extend='both', orientation='vertical', shrink=1.0, aspect=30.0)
  cbar.ax.tick_params(labelsize=30)
  ax1.tick_params(labelsize=25)

  plt.savefig(fname)

print(f"graph Rdiff & Rdisp: {file_prefix} ..")

for porder in PORDER_LIST:
  
  dat_suffix = ""
  if len(TSCHEME_NAME) > 0: 
    dat_suffix = f"{dat_suffix}_{TSCHEME_NAME}"
  if len(MF_SUFFIX) > 0: 
    dat_suffix = f"{dat_suffix}_{MF_SUFFIX}"
  
  print(f"{file_prefix}P{porder:02}_combinedmode{dat_suffix}.dat")
  k_nondim = np.loadtxt(f"{file_prefix}P{porder:02}_Mode01{dat_suffix}.dat", skiprows=1, usecols=0)

  amplifiFac_list = np.loadtxt(f"{file_prefix}P{porder:02}_combinedmode{dat_suffix}.dat", skiprows=1, usecols=[2*i+1 for i in TSLOT_ary_ind])
  phaseError_list = np.loadtxt(f"{file_prefix}P{porder:02}_combinedmode{dat_suffix}.dat", skiprows=1, usecols=[2*(i+1) for i in TSLOT_ary_ind])

  k_nondim = np.loadtxt(f"{file_prefix}P{porder:02}_Mode01{dat_suffix}.dat", skiprows=1, usecols=0)
  km_real_list_combmode = {}
  km_aimg_list_combmode = {}
  
  
  l_ = (2.0*np.pi) / k_nondim[1:]
  dx, l = np.meshgrid(dx_, l_)
  print(l[:,0])
  
  cr_DG_elem = Cr_DG_dxeff /  ( porder + 1.0 )
  

  for tind in range(0,len(TSLOT_list)):
    n = TSLOT_list[tind]
    rdiff = 1.0 / Eta * (np.pi/m)**(4.0/3.0) * U_VEL * (1.0 / cr_DG_elem) / (2.0*np.pi*Cs)**2 * (l/np.pi)**2 \
            * np.tile( -np.log( amplifiFac_list[1:,tind] ).reshape(-1,1), (1,len(dx_)) ) / ( n * (porder+1.0) ) * dx**(-1.0/3.0)
            
    rdisp = 1.0 / Eta_dash  * (np.pi/m)**(1.0/3.0) * U_VEL * (1.0 / cr_DG_elem) / (np.pi*Cs)**2  \
            * np.tile( ( np.abs(phaseError_list[1:,tind]) ).reshape(-1,1), (1,len(dx_)) ) / ( n * (porder+1.0) ) * dx**(-1.0/3.0)
            
    os.makedirs(f"{OUT_DIR}/n{int(n)}", exist_ok=True)
    create_fig_rdiff_rdisp( dx, l, rdiff, rdisp,  
      f"{OUT_DIR}/n{int(n)}/Rdiff_Rdisp_P{porder:02}_{BASIS_TYPE}_{FORM_TYPE}{dat_suffix}_m{int(m)}.png")
            
  