import numpy as np
import pandas as pd
import xarray as xr
import os
import joblib

# Memo: day_per_run: ------
# Eh12Ez8P3: 250 day, Eh24Ez16P3: 125 day, Eh48Ez24P3: 50 day
# Eh6Ez4P7: 250 day, Eh12Ez8P7: 125 day, Eh24Ez16P7: 50 day, Eh48Ez32P7: 10 day
# Eh4Ez3P11: 250 day, Eh8Ez6P11: 125 day,  Eh16Ez12P11: 50 day
#-----------------------------------------------------------------


EXP_TOP_DIR = "./rhot_hevi"
EXP_DIR_LIST = {
#  "Eh12Ez8P3", 
# "Eh24Ez16P3",   
#  "Eh48Ez32P3"
#  "Eh6Ez4P7", "Eh12Ez8P7",  
# "Eh24Ez16P7", 
 "Eh48Ez32P7", 
#  "Eh4Ez3P11", 
#  "Eh8Ez6P11",
#"Eh16Ez12P11",
}
# TARGET_RUNDIR_NO_LIST = {
#   "Eh12Ez8P3": list(range(1,5)), "Eh24Ez16P3": list(range(1,16)), "Eh48Ez32P3": list(range(1,21)), 
#   "Eh6Ez4P7": list(range(1,5)), "Eh12Ez8P7": list(range(1,9)), "Eh24Ez16P7": list(range(1,21)), "Eh48Ez32P7": list(range(1,24)), 
#   "Eh4Ez3P11": list(range(1,5)), "Eh8Ez6P11": list(range(1,16)), "Eh16Ez12P11": list(range(1,13)),
# }
TARGET_RUNDIR_NO_LIST = {
  "Eh12Ez8P3": list(range(1,5)), "Eh24Ez16P3": list(range(1,9)), "Eh48Ez32P3": list(range(1,21)), 
  "Eh6Ez4P7": list(range(1,5)), "Eh12Ez8P7": list(range(1,9)), "Eh24Ez16P7": list(range(1,21)), "Eh48Ez32P7": list(range(1,31)), 
  "Eh4Ez3P11": list(range(1,5)), "Eh8Ez6P11": list(range(1,9)), "Eh16Ez12P11": list(range(1,21)),
}


ANALYSIS_OUT_DIR="analysis_out"

#--------------------------------------

def mean_lon( var ):
  lon = var["lon"]
  lon_weight = [ 0.166666666666667, 0.833333333333333, 0.833333333333333, 0.166666666666667 ]
  nelem_lon = int(len(lon)/4)

  df = pd.DataFrame({ 'lon': lon.values,  'lon_intweight': lon_weight  *  nelem_lon } ).set_index("lon")
  ds = df.to_xarray()
  var_mean_lon = ( var * ds["lon_intweight"] ).sum("lon") / (2.0 * nelem_lon )
  var_mean_lon.name = var.name
  return var_mean_lon

def analysis_xtmean_sub( vars_list, exp_dir, n, runno ):
  dir = f"{exp_dir}/run{runno}/outdata_p"
  print(dir)
  nc = xr.open_mfdataset(f"{dir}/history.pe00*.nc", decode_times=False, use_cftime=False, combine='by_coords')
  time = nc["time"]
  #print(time)

  if n > 0: 
      nc = nc.isel(time=slice(1,len(time.values)))
      # print(  nc[vname].time )    
  
  vars_list = {}
  
  for vname in ["Umet", "Vmet", "T"]:
    v = nc[vname].mean("lon")
    v.attrs["units"] = nc[vname].attrs["units"]
    vars_list[vname] = v
  
  up = nc["Umet"] - mean_lon(nc["Umet"])
  vp = nc["Vmet"] - mean_lon(nc["Vmet"])
  tp = nc["T"] - mean_lon(nc["T"])

  up_vp = mean_lon( vp * up )
  up_vp.name = "merid_eddy_momflx"
  up_vp.attrs["units"] = "m2/s2"
  vars_list["merid_eddy_momflx"] = up_vp

  vp_tp = mean_lon( vp * tp )
  vp_tp.name = "merid_eddy_hflx"
  vp_tp.attrs["units"] = "K.m/s"    
  vars_list["merid_eddy_hflx"] = vp_tp

  eke = mean_lon( 0.5*(up**2 + vp**2) )
  eke.name = "eddy_kinetic_energy"
  eke.attrs["units"] = "m2/s2"        
  vars_list["eddy_kinetic_energy"] = eke

  tp_tp = mean_lon( tp * tp )
  tp_tp.name = "eddy_temp_variance"
  tp_tp.attrs["units"] = "K2"        
  vars_list["eddy_temp_variance"] = tp_tp
  
  return vars_list
  
def merge_var(vlist):
  time_list = []
  t_offset = 0.0
  for v in vlist:
    time_list.extend(v["time"].values+t_offset)
    t_offset = time_list[-1]
  #print(time_list)
  
  v0 = list(vlist)[0]
  lat_axis = v0.coords["lat"]
  p_axis = v0.coords["p"]
  coords = {
    "lat": (["latitude"], lat_axis.values, lat_axis.attrs), 
    "p": (["pressure"], p_axis.values, p_axis.attrs), 
    "time": (["time"], time_list, v0.coords["time"].attrs),
  }
  np_data = np.zeros((len(time_list),len(p_axis.values),len(lat_axis.values)))
  i_toffset = 0
  for v in vlist:
    tlen = len(v["time"].values)
    np_data[i_toffset:i_toffset+tlen,:,:] = v.values[:,:,:]
    i_toffset = i_toffset + tlen
  
  #print(v0.attrs)
  return xr.DataArray( np_data, dims=["time", "p", "lat"], coords=coords, attrs=v0.attrs, name=v0.name)
  
def merge_var(key, vlist_ori, out_dir_tmp):
  print(f"{key}: concat & output")
  vlist = dict(sorted(vlist_ori.items())).values()

  time_list = []
  t_offset = 0.0
  for v in vlist:
    time_list.extend(v["time"].values+t_offset)
    t_offset = time_list[-1]
  #print(time_list)
  
  v0 = list(vlist)[0]
  lat_axis = v0.coords["lat"]
  p_axis = v0.coords["p"]
  coords = {
    "lat": (["lat"], lat_axis.values, lat_axis.attrs), 
    "p": (["p"], p_axis.values, p_axis.attrs), 
    "time": (["time"], time_list, v0.coords["time"].attrs),
  }
  np_data = np.zeros((len(time_list),len(p_axis.values),len(lat_axis.values)))
  i_toffset = 0
  for v in vlist:
    tlen = len(v["time"].values)
    np_data[i_toffset:i_toffset+tlen,:,:] = v.values[:,:,:]
    i_toffset = i_toffset + tlen
  
  #print(v0.attrs)
  da = xr.DataArray( np_data, dims=["time", "p", "lat"], coords=coords, attrs=v0.attrs, name=v0.name)
  da.to_netcdf(f"{out_dir_tmp}/{key}.nc")
  
def analysis_xtmean(exp_dir, runno_list, out_dir):  
  vars_list = {
    "Umet":{}, "Vmet":{}, "T":{}, 
    "merid_eddy_momflx":{}, "merid_eddy_hflx":{}, 
    "eddy_kinetic_energy":{}, "eddy_temp_variance":{}
    }

  vars_listlist = joblib.Parallel(n_jobs=4, verbose=2)(
    joblib.delayed(analysis_xtmean_sub)( vars_list, exp_dir, n, runno_list[n] ) for n in range(0,len(runno_list)) )

  for n, vars_list0 in enumerate(vars_listlist):
    for varname, var in vars_list0.items():
      vars_list[varname][n] = var      
  # print(vars_list)
  
  out_dir_tmp = f"{out_dir}/tmp_data/"
  os.makedirs(out_dir_tmp, exist_ok=True)
  keys = list(vars_list.keys())
  print(keys)
  result = joblib.Parallel(n_jobs=4, verbose=2)(
    joblib.delayed(merge_var)( keys[n], vars_list[keys[n]], out_dir_tmp ) for n in range(0,len(keys) ) )
  
  # for key, vlist_ori in vars_list.items():
  #   print(f"{key}: concat & output")
  #   vlist = dict(sorted(vlist_ori.items()))
  #   #print(vlist)
  #   merge_var(vlist.values()).to_netcdf(f"{out_dir_tmp}/{key}.nc")

for exp_dir in EXP_DIR_LIST:
  rundir_no_list = TARGET_RUNDIR_NO_LIST[exp_dir]
  analysis_xtmean(f"{EXP_TOP_DIR}/{exp_dir}", rundir_no_list, f"{ANALYSIS_OUT_DIR}/{exp_dir}")