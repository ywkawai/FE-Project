import pyshtools as pysh
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import os

# day_per_run: 
# Eh12Ez8P3: 250 day, Eh24Ez16P3: 125 day, Eh48Ez24P3: 50 day
# Eh6Ez4P7: 250 day, Eh12Ez8P7: 125 day, Eh24Ez16P7: 50 day, Eh48Ez32P7: 10 day
# Eh4Ez3P11: 250 day, Eh8Ez6P11: 125 day,  Eh16Ez12P11: 50 day

EXP_listlist= {
  "P3": ["Eh12Ez8P3", "Eh24Ez16P3", "Eh48Ez32P3"],
  "P7": ["Eh6Ez4P7", "Eh12Ez8P7", "Eh24Ez16P7", "Eh48Ez32P7"], 
  "P11": ["Eh4Ez3P11", "Eh8Ez6P11", "Eh16Ez12P11"], 
}
TARGET_RUNDIR_NO_LIST = {
  "Eh12Ez8P3": list(range(1,5)), "Eh24Ez16P3": list(range(1,9)), "Eh48Ez32P3": list(range(1,21)), 
  "Eh6Ez4P7": list(range(1,5)), "Eh12Ez8P7": list(range(1,9)), "Eh24Ez16P7": list(range(1,21)), "Eh48Ez32P7": list(range(1,31)), 
  "Eh4Ez3P11": list(range(1,5)), "Eh8Ez6P11": list(range(1,9)), "Eh16Ez12P11": list(range(1,21)), 
}
EXP_ref="Eh48Ez32P7"

PRES_INTERP=250e2; Y_lim = [1e-4,1e3]; Y_comps_lim = [8e-2,.14e1]; SLOPE_m3_ampl=2.5e5; SLOPE_m5div3_ampl=1e2
#PRES_INTERP=850e2; Y_lim = [.5e-4,2e2]; Y_comps_lim = [1e-1,.2e1]; SLOPE_m3_ampl=4.5e4; SLOPE_m5div3_ampl=2.5e1

LMAX_list = {
#    "Eh12Ez8P3": 64, "Eh24Ez16P3": 128,
    "Eh12Ez8P3": 56, "Eh24Ez16P3": 107, "Eh48Ez32P3": 214, 
    "Eh6Ez4P7": 64, "Eh12Ez8P7": 128, "Eh24Ez16P7":256, "Eh48Ez32P7":512, 
    "Eh4Ez3P11": 64, "Eh8Ez6P11": 128, "Eh16Ez12P11": 256}

EXP_color_list = {
    "Eh12Ez8P3": "blue", "Eh24Ez16P3": "blue", "Eh48Ez32P3": "blue",
    "Eh6Ez4P7": "red", "Eh12Ez8P7": "red", "Eh24Ez16P7":"red", "Eh48Ez32P7": "black", 
    "Eh4Ez3P11": "green", "Eh8Ez6P11": "green", "Eh16Ez12P11": "green" }

EXP_ltype_list = {
    "Eh12Ez8P3": "--", "Eh24Ez16P3": "-.", "Eh48Ez32P3": ":", 
    "Eh6Ez4P7": "--", "Eh12Ez8P7": "-.", "Eh24Ez16P7":":", "Eh48Ez32P7": "-", 
    "Eh4Ez3P11": "--", "Eh8Ez6P11": "-.", "Eh16Ez12P11": ":" }

EXP_lbl_list = {
    "Eh12Ez8P3": "$p=3$, $\Delta_{h,eq}\sim 208$ km, VDOF=32", 
    "Eh24Ez16P3": "$p=3$, $\Delta_{h,eq}\sim 104$ km, VDOF=64", 
    "Eh48Ez32P3": "$p=3$, $\Delta_{h,eq}\sim 52$ km, VDOF=128",     
    "Eh6Ez4P7": "$p=7$, $\Delta_{h,eq}\sim 208$ km, VDOF=32", 
    "Eh12Ez8P7": "$p=7$, $\Delta_{h,eq}\sim 104$ km, VDOF=64", 
    "Eh24Ez16P7": "$p=7$, $\Delta_{h,eq}\sim 52$ km, VDOF=128", 
    "Eh48Ez32P7": "$p=7$, $\Delta_{h,eq}\sim 26$ km, VDOF=256", 
    "Eh4Ez3P11": "$p=11$, $\Delta_{h,eq}\sim 208$ km, VDOF=36", 
    "Eh8Ez6P11": "$p=11$, $\Delta_{h,eq}\sim 104$ km, VDOF=72", 
    "Eh16Ez12P11": "$p=11$, $\Delta_{h,eq}\sim 52$ km, VDOF=144",     
    }

ANALYSIS_OUT_DIR="analysis_out"

#---------
def set_axis(ax, ylabel, ylim):
    ax.tick_params(labelsize=18, length=8)    
    ax.set(yscale='log', xscale='log', ylabel='Power')
    ax.set_xlabel('Spherical harmonic degree', fontsize=24)
    ax.set_ylabel(ylabel, fontsize=24)
    ax.grid()
    ax.set_xlim([1e0,6e2])
    ax.set_ylim(ylim)

def mkgraph( hke_spectra_listlist, wke_spectra_listlist, slope_m3_ampl, slope_m5div3_ampl, pngname,):
    fig, ax = plt.subplots(1, 1, figsize=(12,11))
    set_axis(ax, "E(k)", Y_lim)
    
    i = 0
    for explist_key, EXP_list in EXP_listlist.items():
        hke_spectra_list = hke_spectra_listlist[i]
        for exp_name, spectra in hke_spectra_list.items():
            spectra_tmp = spectra.sel(n=slice(0,LMAX_list[exp_name])).copy()
            spectra_tmp[:,1:-1] = 2.0 * spectra_tmp[:,1:-1]        
            spectra_ = spectra_tmp.sum(["m"])
            print(spectra_)
            ax.plot( spectra_.n[0:], spectra_.values, 
                     color=EXP_color_list[exp_name], 
                     linestyle=EXP_ltype_list[exp_name], linewidth=3,
                     label=EXP_lbl_list[exp_name] )
            if EXP_ref == exp_name:
                n_ref = spectra_.n
        i = i + 1
    
        
    ax.plot(n_ref, slope_m3_ampl*n_ref**(-3.0), linestyle="-", color="gray")        
#    ax.plot(n_ref, slope_m5div3_ampl*n_ref**(-5.0/3.0),label="-5/3", linestyle="-", color="lightgray")        
    ax.legend(fontsize=18)
    
    plt.savefig(pngname)

def mkgraph_compensate( hke_spectra_listlist, wke_spectra_listlist, slope_m3_ampl, slope_m5div3_ampl, pngname,):
    fig, ax = plt.subplots(1, 1, figsize=(12,6))
    fig.subplots_adjust(bottom=0.15)    
    set_axis(ax, "$E(k)*k^3$", Y_comps_lim)
    
    i = 0
    for explist_key, EXP_list in EXP_listlist.items():
        hke_spectra_list = hke_spectra_listlist[i]
        for exp_name, spectra in hke_spectra_list.items():
            spectra_tmp = spectra.sel(n=slice(0,LMAX_list[exp_name])).copy()
            spectra_tmp[:,1:-1] = 2.0 * spectra_tmp[:,1:-1]        
            spectra_ = spectra_tmp.sum(["m"]) 
            comp_spectra_ = spectra_ / (slope_m3_ampl*spectra_.n[0:]**(-3.0))
            print(spectra_)
            ax.plot(comp_spectra_.n[0:], comp_spectra_.values, 
                    label=EXP_lbl_list[exp_name], color=EXP_color_list[exp_name], 
                    linestyle=EXP_ltype_list[exp_name], 
                    linewidth=3 )
            if EXP_ref == exp_name:
                n_ref = comp_spectra_.n
        i = i + 1
    
    ax.plot(n_ref, n_ref**(0), linestyle="-", color="gray")        
#    ax.plot(n_ref, slope_m5div3_ampl*n_ref**(-5.0/3.0+3.0),label="-5/3", linestyle="-", color="lightgray")        
#    ax.legend(fontsize=18)
    
    plt.savefig(pngname)

def get_index(var, coord_name, target_pos):
  print(var)
  coord = var.coords[coord_name].values
  return np.argmin((coord-target_pos)**2)

def cal_kinetic_energy(dir, runno_s, runno_e):
    run_num = runno_e - runno_s + 1
    ds = xr.open_mfdataset(f'{dir}/run{runno_s}/outdata_p_spectra/spectral_data.nc', decode_times=False, combine='by_coords')
    level_ipos = get_index(ds["U_r"], "level", PRES_INTERP)
    ds = ds.isel(level=level_ipos)

    u_r = ds["U_r"]; u_i = ds["U_i"]; 
    v_r = ds["V_r"]; v_i = ds["V_i"];
    w_r = ds["W_r"]; w_i = ds["W_i"];
     
    print(f"runno: {runno_s}", u_r.time.values)
    hke = 0.5 * ( u_r**2+u_i**2 + v_r**2+v_i**2 ).mean(["time"]) / float(run_num)
    vke = 0.5 * ( w_r**2+w_i**2 ).mean(["time"]) / float(run_num)    
    for runno in range(runno_s+1,runno_e+1):
        ds = xr.open_mfdataset(f'{dir}/run{runno}/outdata_p_spectra/spectral_data.nc', decode_times=False, combine='by_coords').isel(time=slice(1,-1))
        level_ipos = get_index(ds["U_r"], "level", PRES_INTERP)
        ds = ds.isel(level=level_ipos)
        
        u_r = ds["U_r"]; u_i = ds["U_i"]; 
        v_r = ds["V_r"]; v_i = ds["V_i"];
        w_r = ds["W_r"]; w_i = ds["W_i"];            
#        print(f"runno: {runno}", u_r.time.values)
        hke = hke + 0.5 * ( u_r**2+u_i**2 + v_r**2+v_i**2 ).mean(["time"]) / float(run_num)
        vke = vke + 0.5 * ( w_r**2+w_i**2 ).mean(["time"]) / float(run_num)
    
    return hke.transpose().rename("s_hke"), vke.transpose().rename("s_vke")

#---------

hke_spectra_tavg_listlist = []
vke_spectra_tavg_listlist = []

for exp_list_key, exp_list in EXP_listlist.items():
    hke_spectra_tavg_list = {}
    vke_spectra_tavg_list = {}

    for exp in exp_list:
        runno_list = TARGET_RUNDIR_NO_LIST[exp]
        hke_spectra_tavg_list[exp], vke_spectra_tavg_list[exp] = cal_kinetic_energy(exp, runno_list[0], runno_list[-1])

    hke_spectra_tavg_listlist.append(hke_spectra_tavg_list)
    vke_spectra_tavg_listlist.append(vke_spectra_tavg_list)


os.makedirs(f"{ANALYSIS_OUT_DIR}/energy_spectra", exist_ok=True)
mkgraph( hke_spectra_tavg_listlist, vke_spectra_tavg_listlist, 
        SLOPE_m3_ampl, SLOPE_m5div3_ampl, 
        f"{ANALYSIS_OUT_DIR}/energy_spectra/KE_spectra_p{int(PRES_INTERP/1e2)}hPa_compari.png")

mkgraph_compensate( hke_spectra_tavg_listlist, vke_spectra_tavg_listlist, 
        SLOPE_m3_ampl, SLOPE_m5div3_ampl, 
        f"{ANALYSIS_OUT_DIR}/energy_spectra/KE_compensated_spectra_p{int(PRES_INTERP/1e2)}hPa_compari.png")
