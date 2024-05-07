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

def open_bs_tmp_nc(dir, varname):
  print(f'{dir}/bs.pe*.nc')
  return xr.open_mfdataset(f'{dir}/bs.pe*.nc', decode_times=False, combine='by_coords')[varname]

def open_tmp_nc(dir, varname):
  print(f'{dir}/history.pe*.nc')
  return xr.open_mfdataset(f'{dir}/history.pe*.nc', decode_times=False, combine='by_coords')[varname]


def gen_graph_sub(fig, ax_list, X, Y, vmin, vmax, var, time, time_ind, png_name):
  print("time="+str(time[time_ind]))

  for ax in ax_list:
    cax = ax.pcolormesh(X, Y, var,
        vmin=vmin, vmax=vmax, 
        cmap = "jet", transform=ccrs.PlateCarree(), 
    )
    fig.colorbar(cax, ax=ax, shrink=0.62)
        
    # title = ax.text(150.0, 3100, "time="+str(time.values[time_ind]) + " [s]", backgroundcolor="white", size="large")
    # artists.append([cax, title])
    ax.relim()
    ax.autoscale_view()
  plt.savefig(png_name)
  
def gen_graph(var_list, vmin, vmax, anim_file):
  X, Y = np.meshgrid(var_list[0].lon.values, var_list[0].lat.values)
  plt.subplots_adjust(wspace=0.1, hspace=0.1)

  istart=0
  itime_offset = 0

  for var in var_list:
    time_ = var.time
    time = np.arange(time_[0],time_[0]+1830,300)
    print(time)
    print(itime_offset)
    
    for time_ind in range(istart,len(time)):
      if time[time_ind] >= StartTimeSec:
        fig = plt.figure(figsize=(12,12))

        inclination = 60
        plotcrs = ccrs.Orthographic(0, 90 - inclination)
        ax1 = fig.add_subplot(2,2,1, projection=plotcrs)

        inclination = 60
        plotcrs = ccrs.Orthographic(180, 90 - inclination)
        ax2 = fig.add_subplot(2,2,2, projection=plotcrs)

        plotcrs = ccrs.Orthographic(0, 90)
        ax3 = fig.add_subplot(2,2,3, projection=plotcrs)

        plotcrs = ccrs.Orthographic(0, -90)
        ax4 = fig.add_subplot(2,2,4, projection=plotcrs)

        for ax in [ax1, ax2, ax3, ax4]:
          gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, linestyle='--', color='black')
          gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 30))
          gl.ylocator = mticker.FixedLocator(np.arange(-90, 90, 30))

        out_png = f"{anim_file}_t{itime_offset+time_ind:06}.png"
        print(out_png)
        gen_graph_sub( fig, [ax1, ax2, ax3, ax4], X, Y, vmin, vmax, 
                      var.sel(time=time[time_ind]).values, time, time_ind, 
                      out_png )
      # result = joblib.Parallel(n_jobs=1)(joblib.delayed(gen_graph_sub) 
      #   ( fig, [ax1, ax2, ax3, ax4], X, Y, vmin, vmax, 
      #                 var.isel(time=time_ind).values, time, time_ind, 
      #                 f"{anim_file}_t{itime_offset+time_ind:06}.png" ) for time_ind in range(istart,len(time.values)))
    
    itime_offset = itime_offset + len(time) - 1
    istart=1
    


def get_index(var, coord_name, target_pos):
  print(var)
  coord = var.coords[coord_name].values
  return np.argmin((coord-target_pos)**2)
  
def get_var_list(varname, dir_list):
  var_list = []
  for dir in dir_list:
    var = open_tmp_nc(dir, varname)
    z_pos = get_index(var, "z", 500.0)
    var_list.append(var.isel(z=z_pos))
  return var_list


dir_list = []
for runno in range(RUN_NO_S,RUN_NO_E+1):
  dir_list.append(f"{TOP_EXP_DIR}/run{runno}/outdata")

os.makedirs(OUT_DIR, exist_ok=True)
gen_graph(get_var_list("W", dir_list), -1.5, 1.5, f'{OUT_DIR}/Wsphere')
