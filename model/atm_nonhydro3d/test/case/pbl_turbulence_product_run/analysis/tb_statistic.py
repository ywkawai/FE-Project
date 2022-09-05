#import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import glob
import os
from joblib import Parallel, delayed

RESOL="E120P7"
#RESOL="E240P3"
#RESOL="E480P1"
exp_name = f"{RESOL}" 
#exp_name = f"{RESOL}_MF16Ord" 
#exp_name = f"{RESOL}_MF32OrdAlphx10" 
#exp_name = f"{RESOL}_hypupw" 
#exp_name = f"{RESOL}_MFoff" 

dir=f"../run8_{exp_name}/outdata/"  # the directory which has simulation results
tmp_dir=f"tmp3_data/tmp_{exp_name}/"
out_dir=f"{exp_name}/horidist"
Nproc = 8
PEnum = 144

#time_list = [12600, 13200, 13800, 14400]
time_list = np.arange(3.5*3600, 14440, 60)

def create_tmp_nc(data_dir, tmp_dir):  
  Parallel(n_jobs=Nproc)( [delayed(create_tmp_nc_sub1)(pe,data_dir,tmp_dir) for pe in range(0,PEnum) ] )
  for varname in ["DENS", "MOMZ"]:
    merge_nc(varname)

def create_tmp_nc_sub1(pe, dir, tmp_dir):
  print(f'{dir}history.pe{pe:06}.nc')
  ds_basic = xr.open_mfdataset(f'{dir}basic_state.pe{pe:06}.nc', decode_times=False, combine='by_coords')
  ds = xr.open_mfdataset(f'{dir}history.pe{pe:06}.nc', decode_times=False, combine='by_coords')
  print(f"pe={pe}..")

  vars = {}
  for varname in ["U", "V", "W", "DDENS"]:
    vars[varname] = ds[varname].sel(time=time_list, method='nearest')

  vars["DENS"] = ds_basic["DENS_hyd"] + vars["DDENS"]
  vars["MOMZ"] = vars["DENS"] * vars["W"]

  for varname in ["DENS", "MOMZ"]:
    v = vars[varname].mean(["x", "y"]).rename(varname)
    v.to_netcdf(f"{tmp_dir}/tmp_{varname}.pe{pe:06}.nc")

def merge_nc(varname):
  print(f'Merge nc: {varname}')

  v = xr.open_mfdataset(f"{tmp_dir}/tmp_{varname}.pe000000.nc", decode_times=False, combine='by_coords')[varname]
  v_new = v*1.0
  for pe in range(1,PEnum):
    v = xr.open_mfdataset(f"{tmp_dir}/tmp_{varname}.pe{pe:06}.nc", decode_times=False, combine='by_coords')[varname]
    v_new = v_new + v

  ( v_new / float(PEnum) ).to_netcdf(f"{tmp_dir}/{varname}_XYmean.nc")

#-------------------

def cal_tb_statics(data_dir, tmp_dir):  
  Parallel(n_jobs=Nproc)( [delayed(cal_tb_statics_sub)(pe,data_dir,tmp_dir) for pe in range(0,PEnum) ] )
  for varname in ["W_PRIM2", "W_PRIM3"]:
    merge_nc(varname)

def cal_tb_statics_sub(pe, dir, tmp_dir):
  print(f'{dir}history.pe{pe:06}.nc')

  ds_basic = xr.open_mfdataset(f'{dir}basic_state.pe{pe:06}.nc', decode_times=False, combine='by_coords')
  ds = xr.open_mfdataset(f'{dir}history.pe{pe:06}.nc', decode_times=False, combine='by_coords')
  dens_xymean = xr.open_mfdataset(f"{tmp_dir}/DENS_XYmean.nc", decode_times=False, combine='by_coords')['DENS']
  momz_xymean = xr.open_mfdataset(f"{tmp_dir}/MOMZ_XYmean.nc", decode_times=False, combine='by_coords')['MOMZ']

  print(f"pe={pe}..")
  vars = {}
  for varname in ["W", "DDENS"]:
    vars[varname] = ds[varname].sel(time=time_list, method='nearest')

  vars["DENS"] = ds_basic["DENS_hyd"] + vars["DDENS"]
  w_xymean = momz_xymean / dens_xymean
  w_prim = vars["W"] - w_xymean
  vars["W_PRIM2"] = ( vars["DENS"] * w_prim**2 ).mean(["x", "y"]) / dens_xymean
  vars["W_PRIM3"] = ( vars["DENS"] * w_prim**3 ).mean(["x", "y"]) / dens_xymean

  for varname in ["W_PRIM2", "W_PRIM3"]:
    vars[varname].rename(varname).to_netcdf(f"{tmp_dir}/tmp_{varname}.pe{pe:06}.nc")

# make directory
os.makedirs(tmp_dir, exist_ok=True)
os.makedirs(out_dir, exist_ok=True)

create_tmp_nc(dir, tmp_dir)
cal_tb_statics(dir, tmp_dir)
