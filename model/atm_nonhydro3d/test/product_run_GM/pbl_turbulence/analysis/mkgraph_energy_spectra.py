import pyshtools as pysh
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import xarray as xr
import numpy as np
import os

TOP_DIR="rp3.4km"

# EXP_list= [
#     "Eh64Ez34P7",   
#     "Eh100Ez52P4",                  
#     "Eh128Ez64P3",
# ]
EXP_list= [
    "Eh64Ez34P7_deepatm",     
    "Eh100Ez52P4_deepatm",  
    "Eh128Ez64P3_deepatm",      
]

TARGET_RUNDIR_NO_LIST = {
  "Eh128Ez64P3": list(range(8,9)),     
  "Eh100Ez52P4": list(range(8,9)),       
  "Eh64Ez34P7": list(range(8,9)), 
  "Eh128Ez64P3_deepatm": list(range(8,9)),       
  "Eh100Ez52P4_deepatm": list(range(8,9)),         
  "Eh64Ez34P7_deepatm": list(range(8,9))  
}

exp_label_list = {
  "Eh128Ez64P3": "$p=3$", 
  "Eh100Ez52P4": "$p=4$", 
  "Eh64Ez34P7": "$p=7$",   
  "Eh128Ez64P3_deepatm": "$p=3$", 
  "Eh100Ez52P4_deepatm": "$p=4$",   
  "Eh64Ez34P7_deepatm": "$p=7$",     
}

exp_color_list = {
  "Eh128Ez64P3": "blue", 
  "Eh100Ez52P4": "red", 
  "Eh64Ez34P7": "goldenrod",   
  "Eh128Ez64P3_deepatm": "blue", 
  "Eh100Ez52P4_deepatm": "red",   
  "Eh64Ez34P7_deepatm": "goldenrod",     
}

LMAX_list = {
    "Eh128Ez64P3": 1024,
    "Eh100Ez52P4": 1024,  
    "Eh64Ez34P7": 1024, 
    "Eh64Ez34P7_deepatm": 1024, 
    "Eh100Ez52P4_deepatm": 1024,      
    "Eh128Ez64P3_deepatm": 1024,      
}
exp_ltype_width = {
  "Eh128Ez64P3": 2,  
  "Eh100Ez52P4": 2,   
  "Eh64Ez34P7": 2,
  "Eh128Ez64P3_deepatm": 2,      
  "Eh100Ez52P4_deepatm": 2,  
  "Eh64Ez34P7_deepatm": 2,    
}

ANALYSIS_OUT_DIR="analysis_out"

# fig_suffix="_rp3.4km_shallow_atm_approx"
# Z_INTERP=500; Y_lim = [5e-6,1.1e0];  Y_lim_zoom = [1e-5,.5e-1]; 
# SLOPE_m3_ampl=2.4e5; SLOPE_m5div3_ampl=.8e2

fig_suffix="_rp3.4km_no_shallow_atm_approx"
Z_INTERP=500; Y_lim = [2e-6,0.8e0]; Y_lim_zoom = [1e-5,.3e-1]; 
SLOPE_m3_ampl=1.5e5; SLOPE_m5div3_ampl=.53e2

#---------
def mkgraph( hke_spectra_list, wke_spectra_list, slope_m3_ampl, slope_m5div3_ampl, pngname,):
    fig, ax = plt.subplots(1, 1, figsize=(9,6))
    
    for expname, hspectra in hke_spectra_list.items():
        spectra = hspectra + wke_spectra_list[expname]
        spectra_tmp = spectra.sel(n=slice(0,LMAX_list[expname])).copy()
        spectra_tmp[:,1:-1] = 2.0 * spectra_tmp[:,1:-1]        
        spectra_ = spectra_tmp.sum(["m"])
        print(spectra_.values)
        ax.plot(spectra_.n[0:], spectra_.values, label=exp_label_list[expname], color=exp_color_list[expname], linewidth=exp_ltype_width[expname])

#    ax.plot(spectra_.n, slope_m3_ampl*spectra_.n**(-3.0),label="-3", linestyle="-.", color="black")
    ax.plot(spectra_.n, slope_m5div3_ampl*spectra_.n**(-5.0/3.0),label="-5/3", linestyle="-.", color="gray")
    # ax.plot([1024.0]*3,[0.0,0.01,10.0], linestyle="--", color="lightgray", alpha=0.5)
    # ax.plot([512.0]*3,[0.0,0.01,10.0], linestyle="--", color="lightgray", alpha=0.5)    
    # ax.plot([256.0]*3,[0.0,0.01,10.0], linestyle="--", color="lightgray", alpha=0.5)    
    ax.tick_params(axis="both", which="major", length=6, labelsize=14)
    ax.tick_params(axis="both", which="minor", length=4)
    ax.set(yscale='log', xscale='log')
    ax.set_xlabel('Spherical harmonic degree', fontsize=18)
    ax.set_ylabel("$E(k)$", fontsize=18)
#    ax.grid()
    ax.set_xlim([4e0,1.024e3])
    ax.set_ylim(Y_lim)
#    ax.legend(fontsize=15)
    plt.savefig(pngname)

def mkgraph_diffm53(hke_spectra_list, wke_spectra_list, slope_m3_ampl, slope_m5div3_ampl, pngname,):
  fig, ax = plt.subplots(1, 1, figsize=(10,5.8))
  
  for expname, hspectra in hke_spectra_list.items():
      spectra = hspectra + wke_spectra_list[expname]
      spectra_tmp = spectra.sel(n=slice(0,LMAX_list[expname])).copy()
      spectra_tmp[:,1:-1] = 2.0 * spectra_tmp[:,1:-1]        
      spectra_ = spectra_tmp.sum(["m"])
      slope_m35 = slope_m5div3_ampl*spectra_.n**(-5.0/3.0)
      print(spectra_.values)
      ax.plot(spectra_.n[0:], spectra_.values/slope_m35, label=expname, color=exp_color_list[expname], linewidth=exp_ltype_width[expname])

#    ax.plot(spectra_.n, slope_m3_ampl*spectra_.n**(-3.0),label="-3", linestyle="-.", color="black")
  ax.plot(spectra_.n, slope_m35/slope_m35, label="-5/3", linestyle="-.", color="gray")    
  ax.tick_params(axis="both", which="major", length=6, labelsize=14)
  ax.tick_params(axis="both", which="minor", length=4)  
  ax.set(yscale='log', xscale='log')
  ax.set_xlabel('Spherical harmonic degree', fontsize=18)  
  ax.set_ylabel("$E(k)*k^{5/3}$", fontsize=18)
#  ax.grid()
  ax.set_xlim([4e0,1.024e3])
  ax.set_ylim(4e-1, 1.4e0)
#  ax.legend(fontsize=15)
  plt.savefig(pngname)

def mkgraph_zoom(hke_spectra_list, wke_spectra_list, slope_m5div3_ampl, pngname):
  fig, ax = plt.subplots(figsize=(5,12))

  for expname, hspectra in hke_spectra_list.items():
      spectra = hspectra + wke_spectra_list[expname]
      spectra_tmp = spectra.sel(n=slice(0,LMAX_list[expname])).copy()
      spectra_tmp[:,1:-1] = 2.0 * spectra_tmp[:,1:-1]        
      spectra_ = spectra_tmp.sum(["m"])
      print(spectra_.values)
      ax.plot(spectra_.n[0:], spectra_.values, label=exp_label_list[expname], color=exp_color_list[expname], linewidth=exp_ltype_width[expname])

  ax.plot(spectra_.n, slope_m5div3_ampl*spectra_.n**(-5.0/3.0),label="-5/3", linestyle="-.", color="gray")
  ax.plot([1024.0]*3,[0.0,0.01,10.0], linestyle="--", color="lightgray", alpha=0.5)
  ax.plot([512.0]*3,[0.0,0.01,10.0], linestyle="--", color="lightgray", alpha=0.5)    
  ax.plot([256.0]*3,[0.0,0.01,10.0], linestyle="--", color="lightgray", alpha=0.5)    
  ax.tick_params(axis="both", which="major", length=6, labelsize=14)
  ax.tick_params(axis="both", which="minor", length=4)
  ax.set(yscale='log', xscale='log')
  ax.set_xlabel('Spherical harmonic degree', fontsize=16)
  ax.set_ylabel("$E(k)$", fontsize=18)
#  ax.grid()
  ax.set_xlim([1e2,1.024e3])
  ax.set_ylim(Y_lim_zoom)
  ax.legend(fontsize=15)
  
  plt.savefig(pngname)
  
def get_index(var, coord_name, target_pos):
  print(var)
  coord = var.coords[coord_name].values
  return np.argmin((coord-target_pos)**2)

def cal_kinetic_energy(dir, runno_s, runno_e):
    run_num = runno_e - runno_s + 1
    ncpath = f'{dir}/run{runno_s}/analysis/spectral_data.nc'
    print(ncpath)
    ds = xr.open_mfdataset(ncpath, decode_times=False, combine='by_coords')
    level_ipos = get_index(ds["U_r"], "level", Z_INTERP)
    ds = ds.isel(level=level_ipos)

    # u_r = ds["U_r"]; u_i = ds["U_i"]; 
    # v_r = ds["V_r"]; v_i = ds["V_i"];
    # w_r = ds["W_r"]; w_i = ds["W_i"];
    u_r = ds["RsqrtUmet_r"]; u_i = ds["RsqrtUmet_i"]; 
    v_r = ds["RsqrtVmet_r"]; v_i = ds["RsqrtVmet_i"];
    w_r = ds["RsqrtW_r"]; w_i = ds["RsqrtW_i"];                    
     
    print(f"runno: {runno_s}", u_r.time.values)
    hke = 0.5 * ( u_r**2+u_i**2 + v_r**2+v_i**2 ).mean(["time"]) / float(run_num)
    vke = 0.5 * ( w_r**2+w_i**2 ).mean(["time"]) / float(run_num)    
    for runno in range(runno_s+1,runno_e+1):
        ds = xr.open_mfdataset(f'{dir}/run{runno}/analysis/spectral_data.nc', decode_times=False, combine='by_coords').isel(time=slice(1,-1))
        level_ipos = get_index(ds["U_r"], "level", Z_INTERP)
        ds = ds.isel(level=level_ipos)
        
        # u_r = ds["U_r"]; u_i = ds["U_i"]; 
        # v_r = ds["V_r"]; v_i = ds["V_i"];
        # w_r = ds["W_r"]; w_i = ds["W_i"];            
        u_r = ds["RsqrtUmet_r"]; u_i = ds["RsqrtUmet_i"]; 
        v_r = ds["RsqrtVmet_r"]; v_i = ds["RsqrtVmet_i"];
        w_r = ds["RsqrtW_r"]; w_i = ds["RsqrtW_i"];                    
#        print(f"runno: {runno}", u_r.time.values)
        hke = hke + 0.5 * ( u_r**2+u_i**2 + v_r**2+v_i**2 ).mean(["time"]) / float(run_num)
        vke = vke + 0.5 * ( w_r**2+w_i**2 ).mean(["time"]) / float(run_num)
        
    return hke.transpose().rename("s_hke"), vke.transpose().rename("s_vke")
#---------

hke_spectra_tavg_list = {}
vke_spectra_tavg_list = {}

for exp in EXP_list:
    runno_list = TARGET_RUNDIR_NO_LIST[exp]
    hke_spectra_tavg_list[exp], vke_spectra_tavg_list[exp] = cal_kinetic_energy(f"{TOP_DIR}/{exp}", runno_list[0], runno_list[-1])

os.makedirs(f"{ANALYSIS_OUT_DIR}/energy_spectra", exist_ok=True)
mkgraph( hke_spectra_tavg_list, vke_spectra_tavg_list, 
        SLOPE_m3_ampl, SLOPE_m5div3_ampl, 
        f"{ANALYSIS_OUT_DIR}/energy_spectra/KE_spectra_z{int(Z_INTERP)}m{fig_suffix}.pdf")
mkgraph_diffm53( hke_spectra_tavg_list, vke_spectra_tavg_list, 
        SLOPE_m3_ampl, SLOPE_m5div3_ampl, 
        f"{ANALYSIS_OUT_DIR}/energy_spectra/KE_comp_spectra_z{int(Z_INTERP)}m{fig_suffix}.pdf")
mkgraph_zoom( hke_spectra_tavg_list, vke_spectra_tavg_list, 
        SLOPE_m5div3_ampl, 
        f"{ANALYSIS_OUT_DIR}/energy_spectra/KE_spectra_zoom_z{int(Z_INTERP)}m{fig_suffix}.pdf")
