import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from joblib import Parallel, delayed
from matplotlib.ticker import FuncFormatter
import matplotlib.colors as colors
import scipy.stats as st
import os

#RESOL="E80P11"
#RESOL="E120P7"
#RESOL="E160P5"
RESOL="E192P4"
#RESOL="E240P3"
#RESOL="E480P1"
exp_name = f"{RESOL}" 
#exp_name = f"{RESOL}_MF16Ord" 
#exp_name = f"{RESOL}_MF32OrdAlphx10" 
#exp_name = f"{RESOL}_hypupw" 
#exp_name = f"{RESOL}_MFoff" 

dir=f"../run8_{exp_name}/outdata/"  # the directory which has simulation results
tmp_dir=f"tmp_data/tmp_{exp_name}/"
out_dir=f"{exp_name}/horidist"
Nproc = 8
PEnum = 144
if exp_name == 'E80P11':
  PeNum = 288

time_list = [12600,12900, 13200, 13500, 13800, 14100, 14400]
zlev_list = [1200]

def create_tmp_nc():  
  Parallel(n_jobs=Nproc)( [delayed(create_tmp_nc_sub1)(pe) for pe in range(0,PEnum) ] )

def create_tmp_nc_sub1(pe):
  print(f'{dir}history.pe{pe:06}.nc')
  ds = xr.open_mfdataset(f'{dir}history.pe{pe:06}.nc', decode_times=False, combine='by_coords')
  print(f"pe={pe}..")

  for varname in ["U", "V", "W"]:
    for z in zlev_list:
      var = ds[varname].sel(z=z, method='nearest').sel(time=time_list, method='nearest')
      print(f"Output {varname} z={z} pe={pe}..")
      var.to_netcdf(f"{tmp_dir}/tmp_{varname}_z{z}_interp.pe{pe:06}.nc")
  
    
def m2km_txt(m, pos=None):
  km = int(m/1000)
  return f"{km}"

def m2km_txt2(m, pos=None):
  km = m/1000.0
  return f"{km:.1f}"


def create_fig_UVW(lev, vmin, vmax):
  print(f"{tmp_dir}/tmp_W_z{lev}_interp.pe*.nc")
  nc_w = xr.open_mfdataset(f"{tmp_dir}/tmp_W_z{lev}_interp.pe*.nc", combine='by_coords', decode_times=False)

  for tsec in time_list:
    fig, ax = plt.subplots(figsize=(15,12))

    w = nc_w["W"].sel(time=int(tsec))
    im = ax.pcolormesh(w.coords['x'], w.coords['y'], w, vmin=vmin, vmax=vmax, cmap='seismic')
    cbar = fig.colorbar(im, ax=ax, aspect=30, extend='both', shrink=0.7, pad=0.015)
    cbar.ax.tick_params(labelsize=26, length=6, width=2)

    ax.set_xlabel("x [km]", fontsize=26)
    ax.xaxis.set_major_formatter(FuncFormatter(m2km_txt))
    ax.set_ylabel("y [km]", fontsize=26)
    ax.yaxis.set_major_formatter(FuncFormatter(m2km_txt))
    ax.tick_params(labelsize=26, length=6, width=2)

    spines = 2
    ax.spines["top"].set_linewidth(spines)
    ax.spines["left"].set_linewidth(spines)
    ax.spines["bottom"].set_linewidth(spines)
    ax.spines["right"].set_linewidth(spines)

    plt.title("w [m/s]", fontsize=30)
    plt.savefig(f"{out_dir}/W_z{lev}_t{int(tsec)}.png", dpi=600)

def create_fig_W_PRIM2(lev, vmin, vmax):
  print(f"{tmp_dir}/tmp_W_PRIM3_z{lev}_interp.pe*.nc")
  nc_w = xr.open_mfdataset(f"{tmp_dir}/tmp_W_z{lev}_interp.pe*.nc", combine='by_coords', decode_times=False)

  for tsec in time_list:
    fig, ax = plt.subplots(figsize=(15,12))

    w = nc_w["W"].sel(time=int(tsec))
    w_prim = w - w.mean(['x','y'])
    im = ax.pcolormesh(w.coords['x'], w.coords['y'], 
      w_prim**2, vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(im, ax=ax, aspect=30, extend='both', shrink=0.7, pad=0.015)
    cbar.ax.tick_params(labelsize=26, length=6, width=2)

    ax.set_xlabel("x [km]", fontsize=26)
    ax.xaxis.set_major_formatter(FuncFormatter(m2km_txt))
    ax.set_ylabel("y [km]", fontsize=26)
    ax.yaxis.set_major_formatter(FuncFormatter(m2km_txt))
    ax.tick_params(labelsize=26, length=6, width=2)

    spines = 2
    ax.spines["top"].set_linewidth(spines)
    ax.spines["left"].set_linewidth(spines)
    ax.spines["bottom"].set_linewidth(spines)
    ax.spines["right"].set_linewidth(spines)

    plt.title("w_prim2 [m/s]", fontsize=30)
    plt.savefig(f"{out_dir}/W_PRIM2_z{lev}_t{int(tsec)}.png", dpi=600)

def create_fig_W_PRIM3(lev, vmin, vmax):
  print(f"{tmp_dir}/tmp_W_PRIM3_z{lev}_interp.pe*.nc")
  nc_w = xr.open_mfdataset(f"{tmp_dir}/tmp_W_z{lev}_interp.pe*.nc", combine='by_coords', decode_times=False)

  for tsec in time_list:
    fig, ax = plt.subplots(figsize=(15,12))

    w = nc_w["W"].sel(time=int(tsec))
    w_prim = w - w.mean(['x','y'])
    im = ax.pcolormesh(w.coords['x'], w.coords['y'], 
      w_prim**3, vmin=vmin, vmax=vmax, cmap='seismic')
    cbar = fig.colorbar(im, ax=ax, aspect=30, extend='both', shrink=0.7, pad=0.015)
    cbar.ax.tick_params(labelsize=26, length=6, width=2)

    ax.set_xlabel("x [km]", fontsize=26)
    ax.xaxis.set_major_formatter(FuncFormatter(m2km_txt))
    ax.set_ylabel("y [km]", fontsize=26)
    ax.yaxis.set_major_formatter(FuncFormatter(m2km_txt))
    ax.tick_params(labelsize=26, length=6, width=2)

    spines = 2
    ax.spines["top"].set_linewidth(spines)
    ax.spines["left"].set_linewidth(spines)
    ax.spines["bottom"].set_linewidth(spines)
    ax.spines["right"].set_linewidth(spines)

    plt.title("w_prim3 [m/s]", fontsize=30)
    plt.savefig(f"{out_dir}/W_PRIM3_z{lev}_t{int(tsec)}.png", dpi=600)

def create_fig_Sabs(lev, vmin, vmax):
  nc_sabs = xr.open_mfdataset(f"{tmp_dir}/tmp_Sabs_hCD10_cX7_z500_t007.pe0000*.nc", combine='by_coords', decode_times=False)
  #print(nc_u.time)

  fig, ax = plt.subplots()
  sabs = nc_sabs.Sabs_hCD10_cX7.sel(time=14400, method='nearest')
  sabs.plot( 
    norm=colors.LogNorm(vmin=vmin, vmax=vmax), 
    cmap='hot_r', 
    cbar_kwargs={"shrink": 0.6, "extend": "min"})

  ax.set_xlabel("x [km]")
  ax.xaxis.set_major_formatter(FuncFormatter(m2km_txt))
  ax.set_ylabel("y [km]")
  ax.yaxis.set_major_formatter(FuncFormatter(m2km_txt))
  plt.title("|S| [s$^{-1}$]")
  plt.savefig(f"{out_dir}/Sabs_z500_t14400.png")


##---------------------------------
print("create merged nc ..")

os.makedirs(tmp_dir, exist_ok=True)
create_tmp_nc()

##---------------------------------
print("create figures..")

os.makedirs(out_dir, exist_ok=True)

# Fig.1(a)
create_fig_UVW(500, -5.0, 5.0)
create_fig_UVW(1200, -5.0, 5.0)
# create_fig_W_PRIM2(500, 0.0, 20.0)
# create_fig_W_PRIM2(1200, 0.0, 20.0)
# create_fig_W_PRIM3(500, -60.0, 60.0)
# create_fig_W_PRIM3(1200, -60.0, 60.0)

# Fig.1(b)
#create_fig_Sabs(500,  3e-3, 1e-1)
