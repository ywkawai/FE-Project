import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import os
import cartopy
import cartopy.crs as ccrs
import joblib

#TOP_EXP_DIR="./rp3.4km/Eh64Ez34P7"
TOP_EXP_DIR="./rp3.4km/Eh64Ez34P7_deepatm"
OUT_DIR=f"analysis_out/anim/{TOP_EXP_DIR}/"

RUN_NO_S=1
RUN_NO_E=8
StartTimeSec=0
TimePerRun = 1800
TimeInterval=300
RPlanet=3.4e3

def open_bs_tmp_nc(dir, varname):
  print(f'{dir}/bs.pe*.nc')
  return xr.open_mfdataset(f'{dir}/bs.pe*.nc', decode_times=False, combine='by_coords')[varname]

def open_tmp_nc(dir, varname):
  print(f'{dir}/history.pe*.nc')
  return xr.open_mfdataset(f'{dir}/history.pe*.nc', decode_times=False, combine='by_coords')[varname]


def gen_graph_sub(fig, X_list, Z_list, vmin, vmax, var_list, time, time_ind, png_name):
  print("time="+str(time[time_ind]))

  fig = plt.figure(figsize=(12,12))
  ax1 = fig.add_subplot(1,1,1, projection="polar")
  ax_list = [ax1]

  for i, ax in enumerate(ax_list):
    print(var_list[i].shape)
    ax.set_rlim(0, RPlanet)
    ax.set_ylim(0, RPlanet+3e3)
    cax = ax.pcolormesh(X_list[i]/180.0*np.pi, Z_list[i], var_list[i],
        vmin=vmin, vmax=vmax, 
        cmap = "jet" ) 
  
    fig.colorbar(cax, ax=ax, shrink=0.62)
        
    # title = ax.text(150.0, 3100, "time="+str(time.values[time_ind]) + " [s]", backgroundcolor="white", size="large")
    # artists.append([cax, title])
  plt.savefig(png_name)

  
def gen_graph(var_lonz_list, vmin, vmax, anim_file):
  X1, R1 = np.meshgrid(var_lonz_list[0].lon.values, RPlanet+var_lonz_list[0].z.values)
#  X2, R2 = np.meshgrid(var_latz_list[0].lat.values, RPlanet+var_latz_list[0].z.values)
  
  plt.subplots_adjust(wspace=0.1, hspace=0.1)

  istart=0
  itime_offset = 0

  for var in var_lonz_list:
    time_ = var.time
    time = np.arange(time_[0],time_[0]+1830,300)
    print(time)
    print(itime_offset)
    
    for time_ind in range(istart,len(time)):
      if time[time_ind] >= StartTimeSec:
        fig = plt.figure(figsize=(12,12))

        out_png = f"{anim_file}_t{itime_offset+time_ind:06}.png"
        print(out_png)
        gen_graph_sub( fig, [X1], [R1], vmin, vmax, 
                      [var.sel(time=time[time_ind]).values], time, time_ind, 
                      out_png )
    
    itime_offset = itime_offset + len(time) - 1
    istart=1
    


def get_index(var, coord_name, target_pos):
  print(var)
  coord = var.coords[coord_name].values
  return np.argmin((coord-target_pos)**2)
  
def get_var_lonz_list(varname, dir_list):
  var_list = []
  for dir in dir_list:
    var = open_tmp_nc(dir, varname)
    lat_pos = get_index(var, "lat", 0.0)
    var_list.append(var.isel(lat=lat_pos))
  return var_list


dir_list = []
for runno in range(RUN_NO_S,RUN_NO_E+1):
  dir_list.append(f"{TOP_EXP_DIR}/run{runno}/outdata")

os.makedirs(OUT_DIR, exist_ok=True)
gen_graph(get_var_lonz_list("W", dir_list), -1.5, 1.5, f'{OUT_DIR}/W_lonz')
