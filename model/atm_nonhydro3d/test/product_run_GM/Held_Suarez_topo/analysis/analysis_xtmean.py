import numpy as np
import pandas as pd
import xarray as xr
import os

EXP_DIR_LIST = {
#  "Eh12Ez8P3",  "Eh24Ez16P3",   
# "Eh6Ez4P7_topo_egn64", "Eh12Ez8P7_topo_egn64",  
  "Eh24Ez16P7_topo_egn64", 
#  "Eh48Ez32P7", 
#  "Eh4Ez3P11", "Eh8Ez6P11"
}
TARGET_RUNDIR_NO_LIST = {
  "Eh12Ez8P3_topo_egn64": list(range(1,5)), "Eh24Ez16P3_topo_egn64": list(range(1,9)), 
  "Eh6Ez4P7_topo_egn64": list(range(1,5)), "Eh12Ez8P7_topo_egn64": list(range(1,6)), "Eh24Ez16P7_topo_egn64": list(range(1,10)), "Eh48Ez32P7_topo_egn64": list(range(1,15)), 
  "Eh4Ez3P11_topo_egn64": list(range(1,5)), "Eh8Ez6P11_topo_egn64": list(range(1,9))
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

def analysis_xtmean(exp_dir, runno_list, out_dir):  
  vars_list = {
    "Umet":[], "Vmet":[], "T":[], "merid_eddy_momflx":[], "merid_eddy_hflx":[], 
    "eddy_kinetic_energy":[], "eddy_temp_variance":[]}

  for n in range(0,len(runno_list)):
    runno = runno_list[n]
    dir = f"{exp_dir}/run{runno}/outdata_p"
    print(dir)
    nc = xr.open_mfdataset(f"{dir}/history.pe000*.nc", decode_times=True, use_cftime=True, combine='by_coords')
    time = nc["time"]
    # print(time)

    if n > 0: 
        nc = nc.isel(time=slice(1,len(time.values)))
        # print(  nc[vname].time )    
    
    for vname in ["Umet", "Vmet", "T"]:
        vars_list[vname].append( nc[vname].mean("lon") )

    up = nc["Umet"] - mean_lon(nc["Umet"])
    vp = nc["Vmet"] - mean_lon(nc["Vmet"])
    tp = nc["T"] - mean_lon(nc["T"])

    up_vp = mean_lon( vp * up )
    up_vp.name = "merid_eddy_momflx"
    vars_list["merid_eddy_momflx"].append( up_vp )

    vp_tp = mean_lon( vp * tp )
    vp_tp.name = "merid_eddy_hflx"
    vars_list["merid_eddy_hflx"].append( vp_tp )

    eke = mean_lon( 0.5*(up**2 + vp**2) )
    eke.name = "eddy_kinetic_energy"
    vars_list["eddy_kinetic_energy"].append( eke )  

    tp_tp = mean_lon( tp * tp )
    tp_tp.name = "eddy_temp_variance"
    vars_list["eddy_temp_variance"].append( tp_tp )  

  out_dir_tmp = f"{out_dir}/tmp_data/"
  os.makedirs(out_dir_tmp, exist_ok=True)
  for key, vlist in vars_list.items():
    print(f"{key}: concat & output")
    var = xr.concat(vlist, "time")
    var.to_netcdf(f"{out_dir_tmp}/{key}.nc")

for exp_dir in EXP_DIR_LIST:
  rundir_no_list = TARGET_RUNDIR_NO_LIST[exp_dir]
  analysis_xtmean(exp_dir, rundir_no_list, f"{ANALYSIS_OUT_DIR}/{exp_dir}")