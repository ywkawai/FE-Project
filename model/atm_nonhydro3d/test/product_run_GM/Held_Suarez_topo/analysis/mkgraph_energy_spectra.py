import pyshtools as pysh
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import os

EXP_list= [
    "Eh12Ez8P3", "Eh24Ez16P3",
#    "Eh6Ez4P7", "Eh12Ez8P7", "Eh24Ez16P7", "Eh48Ez32P7", 
#    "Eh4Ez3P11", "Eh8Ez6P11", 
]
TARGET_RUNDIR_NO_LIST = {
  "Eh12Ez8P3": list(range(1,5)), "Eh24Ez16P3": list(range(1,5)),
  "Eh6Ez4P7": list(range(1,5)), "Eh12Ez8P7": list(range(1,9)), "Eh24Ez16P7": list(range(1,9)), "Eh48Ez32P7": list(range(1,7)), 
  "Eh4Ez3P11": list(range(1,5)), "Eh8Ez6P11": list(range(1,9))
}

PRES_INTERP=250e2; Y_lim = [1e-4,1e3]; SLOPE_m3_ampl=2.5e5; SLOPE_m5div3_ampl=1e2
#PRES_INTERP=850e2; Y_lim = [.5e-4,2e2]; SLOPE_m3_ampl=4.5e4; SLOPE_m5div3_ampl=2.5e1

LMAX_list = {
#    "Eh12Ez8P3": 64, "Eh24Ez16P3": 128,
    "Eh12Ez8P3": 42, "Eh24Ez16P3": 85,    
    "Eh6Ez4P7": 64, "Eh12Ez8P7": 128, "Eh24Ez16P7":256, "Eh48Ez32P7":512, 
    "Eh4Ez3P11": 64, "Eh8Ez6P11": 128, }

ANALYSIS_OUT_DIR="analysis_out"
fig_suffix="_p3_resoldep"
#fig_suffix="_p7_resoldep"
#fig_suffix="_p11_resoldep"

#---------
def mkgraph( hke_spectra_list, wke_spectra_list, slope_m3_ampl, slope_m5div3_ampl, pngname,):
    fig, ax = plt.subplots(1, 1, figsize=(8,6))
    
    for resol, spectra in hke_spectra_list.items():
        spectra_tmp = spectra.sel(n=slice(0,LMAX_list[resol])).copy()
        spectra_tmp[:,1:-1] = 2.0 * spectra_tmp[:,1:-1]        
        spectra_ = spectra_tmp.sum(["m"])
        print(spectra_)
        ax.plot(spectra_.n[0:], spectra_.values, label=resol)

    ax.plot(spectra_.n, slope_m3_ampl*spectra_.n**(-3.0),label="-3", linestyle="-.", color="black")
    ax.plot(spectra_.n, slope_m5div3_ampl*spectra_.n**(-5.0/3.0),label="-5/3", linestyle="--", color="lightgray")    
    ax.set(yscale='log', xscale='log', xlabel='Spherical harmonic degree', ylabel='Power')
    ax.set_ylabel("E(k)", fontsize=15)
    ax.grid()
    ax.set_xlim([1e0,6e2])
    ax.set_ylim(Y_lim)
    ax.legend()
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

hke_spectra_tavg_list = {}
vke_spectra_tavg_list = {}

for exp in EXP_list:
    runno_list = TARGET_RUNDIR_NO_LIST[exp]
    hke_spectra_tavg_list[exp], vke_spectra_tavg_list[exp] = cal_kinetic_energy(exp, runno_list[0], runno_list[-1])

os.makedirs(f"{ANALYSIS_OUT_DIR}/energy_spectra", exist_ok=True)
mkgraph( hke_spectra_tavg_list, vke_spectra_tavg_list, 
        SLOPE_m3_ampl, SLOPE_m5div3_ampl, 
        f"{ANALYSIS_OUT_DIR}/energy_spectra/KE_spectra_p{int(PRES_INTERP/1e2)}hPa{fig_suffix}.png")
