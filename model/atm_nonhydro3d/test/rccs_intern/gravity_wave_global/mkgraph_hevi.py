import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import os

#---
run_dir = "Eh8Ez6P7_hevi"
dist_dir = f"figs/{run_dir}"
#---
# run_dir = "Eh16Ez12P3_hevi"
# dist_dir = f"figs/{run_dir}"

#------------------
def open_nc(dir, varname):
    return  xr.open_mfdataset(f"{dir}/history.pe000*.nc", decode_times=False, combine='by_coords')[varname]

def v_1Daxis_fmt(tick_val, pos):
  val = int(tick_val/1000)
  return f'{val}'

def set_fig_Xz_axis(ax):
  ax.tick_params(labelsize=20, length=8)
  ax.set_xlabel('longitude [deg]', fontsize=22)  
  ax.xaxis.set_major_locator(tick.MultipleLocator(30.0))
  ax.xaxis.set_minor_locator(tick.MultipleLocator(10.0))
  ax.set_ylabel('height [km]', fontsize=22)
  ax.yaxis.set_major_locator(tick.MultipleLocator(2000.0))
  ax.yaxis.set_minor_locator(tick.MultipleLocator(500.0))
  ax.yaxis.set_major_formatter(tick.FuncFormatter(v_1Daxis_fmt))

def plot_var_xz(v, fig_title, vmin, vmax, cnt_levels, png_name):
  x = v.coords["lon"]
  z = v.coords["z"]

  X, Z = np.meshgrid(x,z)
  fig = plt.figure(figsize=(15,10)) 
  ax = fig.add_subplot(1,1,1)
  set_fig_Xz_axis(ax)
  ax.set_title(fig_title, fontsize=26)

  print(f"output: {png_name} ..")
  pcm = ax.pcolormesh(X, Z, v.values, vmin=vmin, vmax=vmax, cmap='jet')
  levels = cnt_levels[cnt_levels != 0.0]
  cont = ax.contour(X, Z, v.values, levels=levels, colors=['black'])  
  cont = ax.contour(X, Z, v.values, levels=[0.0], colors=['lightgreen'], linewidths=2.0)  

  cbar = plt.colorbar(pcm, aspect=50.0, 
                      extend='both', shrink=0.8, orientation='horizontal', pad=0.12)
  cbar.ax.ticklabel_format(style='sci', scilimits=(-3,3))
  cbar.formatter.set_useMathText(True)
  cbar.ax.xaxis.get_offset_text().set(size=16) 
  cbar.ax.tick_params(labelsize=20)

  plt.savefig(png_name)

  #----------------------------------------------------

os.makedirs(dist_dir, exist_ok=True)    

vars = {}
for varname in ["Umet", "W", "DDENS"]:
    vars[varname] = open_nc(f"{run_dir}/outdata", varname)

levels = np.arange(-8e-4, 8e-4+1e-4, 1e-4)    
plot_var_xz(vars["Umet"].sel(lat=0).isel(time=-1)[:,0,:], 'Zonal wind [m/s]', -8e-4, 8e-4, levels,  f"{dist_dir}/Umet.png")

levels = np.arange(-3e-6, 3e-6+2e-7, 2e-7)    
plot_var_xz(vars["W"].sel(lat=0).isel(time=-1)[:,0,:], 'Vertical wind [m/s]', -3e-6, 3e-6, levels,  f"{dist_dir}/W.png")

