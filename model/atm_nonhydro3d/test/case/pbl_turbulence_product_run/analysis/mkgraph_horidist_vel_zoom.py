import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from joblib import Parallel, delayed
from matplotlib.ticker import FuncFormatter
import matplotlib.colors as colors
import scipy.stats as st
import os

exp_name_list = [
#  "E80P11", 
#  "E120P7", 
#  "E160P5"
  "E192P4"  
#   "E240P3", 
#   "E480P1", 
#   "E480P1_MFoff", 
]

Nproc = 16
PEnum = 16

time_list = [14400] 
zlev_list = [500]

def create_tmp_nc(data_dir, tmp_dir):  
  Parallel(n_jobs=Nproc)( [delayed(create_tmp_nc_sub1)(pe,data_dir,tmp_dir) for pe in range(0,PEnum) ] )
  Parallel(n_jobs=Nproc)( [delayed(create_tmp_nc_sub2)(pe,data_dir,tmp_dir) for pe in range(0,PEnum) ] )

def create_tmp_nc_sub1(pe, dir, tmp_dir):
  print(f'{dir}history.pe{pe:06}.nc')
  ds = xr.open_mfdataset(f'{dir}history.pe{pe:06}.nc', decode_times=False, combine='by_coords')

  for varname in ["U", "V"]:
    for z in zlev_list:
      var = ds[varname].sel(z=z, method='nearest').sel(time=time_list, method='nearest')
      print(f"Output {varname} z={z} pe={pe}..")
      var.to_netcdf(f"{tmp_dir}/{varname}_z{z}.pe{pe:06}.nc")
  
def create_tmp_nc_sub2(pe, dir, tmp_dir):
  print(f'{dir}history.pe{pe:06}.nc')
  ds = xr.open_mfdataset(f'{dir}history.pe{pe:06}.nc', decode_times=False, combine='by_coords')

  # W full level
  momz = ds.MOMZ.sel(time=time_list, method='nearest')
  dens = ds.DENS.sel(time=time_list, method='nearest')
  z = dens.z
  for zlev in zlev_list: 
    w = (momz.interp(zh=zlev)/dens.interp(z=zlev)).rename("W")
    print(f"Output W z={zlev} pe={pe}..")
    w.to_netcdf(f"{tmp_dir}/W_z{zlev}.pe{pe:06}.nc")

    
def m2km_txt(m, pos=None):
  km = int(m/1000)
  return f"{km}"

def m2km_txt2(m, pos=None):
  km = m/1000.0
  return f"{km:.1f}"


def create_fig_UVW(lev, vmin, vmax, out_dir):
  print(f"{tmp_dir}/tmp_W_z{lev}_interp.pe*.nc")  
  nc_w = xr.open_mfdataset(f"{tmp_dir}/tmp_W_z{lev}_interp.pe*.nc", combine='by_coords', decode_times=False)

  for tsec in time_list:
    fig, ax = plt.subplots()
#    w = nc_w.W.sel(time=tsec, method='nearest').sel(x=slice(2500,5500), y=slice(3000, 6000))
    w = nc_w.W.sel(time=tsec, method='nearest').sel(x=slice(4500,7500), y=slice(3000, 6000))

    w.plot(vmin=-4.0, vmax=4.0, cmap='seismic', 
      cbar_kwargs={"shrink": 0.6, "extend": "both"})
    ax.set_xlabel("x [km]")
    ax.xaxis.set_major_formatter(FuncFormatter(m2km_txt2))
    ax.set_ylabel("y [km]")
    ax.yaxis.set_major_formatter(FuncFormatter(m2km_txt2))
    plt.title("w [m/s]")
    plt.savefig(f"{out_dir}/UVW_zoom_z{lev}_t{tsec}.png")

#-------------------------------------------------------------

for exp_name in exp_name_list:

  data_dir=f"./data/{exp_name}/"  # the directory which has simulation results
  tmp_dir=f"tmp_data/tmp_{exp_name}/"
  out_dir=f"{exp_name}/horidist"

  print("create merged nc ..")

  os.makedirs(tmp_dir, exist_ok=True)
#  create_tmp_nc(data_dir, tmp_dir)

  print("create figures..")
  os.makedirs(out_dir, exist_ok=True)

  # Fig.7
  create_fig_UVW(500, -4.0, 4.0, out_dir)
