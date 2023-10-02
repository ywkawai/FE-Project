import pyshtools as pysh
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import os

EXP_DIR_LIST = {
#  "Eh6Ez4P7", 
# "Eh12Ez8P7", 
 "Eh24Ez16P7", 
#  "Eh48Ez32P7", 
#  "Eh4Ez3P11", 
#  "Eh8Ez6P11"
}
TARGET_RUNDIR_NO_LIST = {
  "Eh6Ez4P7": list(range(1,5)), "Eh12Ez8P7": list(range(1,6)), "Eh24Ez16P7": list(range(1,4)), "Eh48Ez32P7": list(range(1,3)), 
  "Eh4Ez3P11": list(range(1,5)), "Eh8Ez6P11": list(range(1,8))
}
PRES_INTERP=850e2

ANALYSIS_OUT_DIR="analysis_out"

#---------
def interp_lat(var):
  lat = var.lat
  dlat = 180.0/len(lat)
  lat_interp = np.arange(-90.0+dlat,90.0+dlat/2.0,dlat)
  
  print(f"Interp lat .. {var.name}")
  var_ip = var.interp(lat=lat_interp, kwargs={"fill_value": "extrapolate"})
  return var_ip.rename(var.name)

def get_index(var, coord_name, target_pos):
  print(var)
  coord = var.coords[coord_name].values
  return np.argmin((coord-target_pos)**2)

def get_power_spectra_tavg(umet, vmet, w, time_list):
    power_per_l_list = []
    
    for t in  time_list:
        print(f"time={t/3600.0} hour")    
        u_ = umet.sel(time=t).values
        v_ = vmet.sel(time=t).values        
        w_ = w.sel(time=t).values        
        u_coeffs = pysh.expand.SHExpandDH(u_, sampling=2)
        v_coeffs = pysh.expand.SHExpandDH(v_, sampling=2)        
        w_coeffs = pysh.expand.SHExpandDH(w_, sampling=2)        
        power_per_l_list.append( 
                                0.5*(  pysh.spectralanalysis.spectrum(u_coeffs) 
                                     + pysh.spectralanalysis.spectrum(v_coeffs)
                                     + pysh.spectralanalysis.spectrum(w_coeffs) ) )

    degrees = np.arange(u_coeffs.shape[1])
    
    tlen = len(time_list)
    power_per_l_tavg = power_per_l_list[0] / float(tlen)
    for i in range(0,tlen):
        power_per_l_tavg = power_per_l_tavg + power_per_l_list[i] / float(tlen)
    
    return degrees, power_per_l_tavg

def analyze_spectra(exp_dir, run_no_list, out_dir):
    ke_spectra_tavg = None
    RUN_NUM = len(run_no_list)
    
    out_dir_tmp = f"{out_dir}/tmp_data"
    os.makedirs(out_dir_tmp, exist_ok=True)
    
    for run_no in run_no_list:
        dir = f"{exp_dir}/run{run_no}/outdata_p_uniform"
        ds = xr.open_mfdataset(f'{dir}/history.pe*.nc', decode_times=False, combine='by_coords')
        plev_ind = get_index( ds["Umet"], "p", PRES_INTERP)

        print("Interp & Output tmp data")
        umet_ip_ = interp_lat( ds["Umet"].isel(p=plev_ind) ).to_netcdf(f"{out_dir_tmp}/tmp_umet_ip_{run_no}.nc")
        vmet_ip_ = interp_lat( ds["Vmet"].isel(p=plev_ind) ).to_netcdf(f"{out_dir_tmp}/tmp_vmet_ip_{run_no}.nc")
        w_ip_ = interp_lat( ds["W"].isel(p=plev_ind) ).to_netcdf(f"{out_dir_tmp}/tmp_w_ip_{run_no}.nc")

        print("Spectral analysis")
        umet_ip_ = xr.open_mfdataset(f"{out_dir_tmp}/tmp_umet_ip_{run_no}.nc", decode_times=False, combine='by_coords')["Umet"]
        vmet_ip_ = xr.open_mfdataset(f"{out_dir_tmp}/tmp_vmet_ip_{run_no}.nc", decode_times=False, combine='by_coords')["Vmet"]
        w_ip_ = xr.open_mfdataset(f"{out_dir_tmp}/tmp_w_ip_{run_no}.nc", decode_times=False, combine='by_coords')["W"]
        time_list = umet_ip_.time.values

        degrees, ke_spectra_ = get_power_spectra_tavg(umet_ip_, vmet_ip_, w_ip_, time_list)
        if run_no > 1:
            ke_spectra_tavg = ke_spectra_tavg + ke_spectra_ / float(RUN_NUM)        
        else:
            ke_spectra_tavg = ke_spectra_ / float(RUN_NUM)

    xr_ke_spectra_tavg = xr.DataArray(ke_spectra_tavg, 
                                    dims=["l"], coords={'l': degrees}).rename("KE_spectra")
    xr_ke_spectra_tavg.to_netcdf(f"{out_dir_tmp}/KE_spectra_p{int(PRES_INTERP/1e2)}hPa.nc")
    
for exp_dir in EXP_DIR_LIST:
  rundir_no_list = TARGET_RUNDIR_NO_LIST[exp_dir]
  analyze_spectra(exp_dir, rundir_no_list, f"{ANALYSIS_OUT_DIR}/{exp_dir}")    