import numpy as np
import xarray as xr
import matplotlib
#from IPython.display import set_matplotlib_formats
#set_matplotlib_formats('retina')
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import matplotlib.animation as animation
import netCDF4

def hori_1Daxis_fmt(tick_val, pos):
  val = int(tick_val/1000)
  return f'{val}'

def set_fig_XY_axis(ax):
  ax.set_xlabel('X (km)')  
  ax.xaxis.set_major_formatter(tick.FuncFormatter(hori_1Daxis_fmt))
  ax.xaxis.set_major_locator(tick.MultipleLocator(1000e3))
  ax.xaxis.set_minor_locator(tick.MultipleLocator(250e3))
  ax.set_ylabel('Y (km)')
  ax.yaxis.set_major_formatter(tick.FuncFormatter(hori_1Daxis_fmt))
  ax.yaxis.set_major_locator(tick.MultipleLocator(1000e3))
  ax.yaxis.set_minor_locator(tick.MultipleLocator(250e3))

def set_fig_YZ_axis(ax):
  ax.set_xlabel('Y (km)')  
  ax.xaxis.set_major_formatter(tick.FuncFormatter(hori_1Daxis_fmt))
  ax.xaxis.set_major_locator(tick.MultipleLocator(1000e3))
  ax.xaxis.set_minor_locator(tick.MultipleLocator(250e3))
  ax.set_ylabel('Z (km)')
  ax.yaxis.set_major_formatter(tick.FuncFormatter(hori_1Daxis_fmt))
  ax.yaxis.set_major_locator(tick.MultipleLocator(1e3))
  ax.yaxis.set_minor_locator(tick.MultipleLocator(250.0))

def plot_var(TIME_cut, varname, vmin, vmax, zisel, nc, ini_diff=False):
  #x_isel=6*24+1
  x_isel=30*8-1

  v = nc.sel(time=TIME_cut).isel(x=slice(0,x_isel),z=zisel)[varname]
  if ini_diff:
    v = v - nc.sel(time=0).isel(x=x_isel,z=zisel)[varname]
  
  x = v.coords["x"]
  y = v.coords["y"]
  z = v.coords["z"]
  time = nc.time

  zlev = z.values

  X, Y = np.meshgrid(x, y)
  fig = plt.figure(figsize=(60,9)) 
  ax = fig.add_subplot(1,1,1)
  set_fig_XY_axis(ax)
  ax.set_title(f"baroclinic wave (z={int(zlev)}m)")

  pcm = ax.pcolormesh(X, Y, v, vmin=vmin, vmax=vmax, cmap='seismic')
  fmt = tick.ScalarFormatter(useMathText=True)
  fmt.set_powerlimits((0,0))
  cbar = plt.colorbar(pcm, aspect=40.0, extend='both', orientation='horizontal', shrink=0.5, format=fmt)
  cbar.ax.xaxis.get_offset_text().set_fontsize(6)

  plt.savefig(f"{v.name}_t{TIME_cut}_z{int(zlev)}.png")

def plot_var_yz(TIME_cut, varname, vmin, vmax, x_isel, nc, ini_diff=False):

  v = nc.sel(time=TIME_cut).isel(x=x_isel)[varname]
  if ini_diff:
    v = v - nc.sel(time=0).isel(x=x_isel)[varname]

  x = v.coords["x"]
  y = v.coords["y"]
  z = v.coords["z"]
  time = nc.time
  x_ = x.values

  Y, Z = np.meshgrid(y,z)
  fig = plt.figure(figsize=(10,8)) 
  ax = fig.add_subplot(1,1,1)
  set_fig_YZ_axis(ax)
  ax.set_title(f"baroclinic wave (x={int(x_/1000e3)}km)")

  pcm = ax.pcolormesh(Y, Z, v, vmin=vmin, vmax=vmax, cmap='jet')
  fmt = tick.ScalarFormatter(useMathText=True)
  fmt.set_powerlimits((0,0))
  cbar = plt.colorbar(pcm, aspect=40.0, extend='both', orientation='horizontal', shrink=0.5, format=fmt)
  cbar.ax.xaxis.get_offset_text().set_fontsize(6)

  plt.savefig(f"{v.name}_t{TIME_cut}_x{int(x_/1000e3)}km.png")

nc = xr.open_mfdataset('./history.pe000*.nc', decode_times=False, combine='by_coords')
#'''
for n_t in range(0,17):
  for ke_z in range(0,3):
    nz = 0#+5*ke_z
    TIME_cut = 900*n_t
    #plot_var(TIME_cut, "MOMX", -0.3, 0.3, nz, nc, True)
    #plot_var(TIME_cut, "MOMY", -0.3, 0.3, nz, nc)
    plot_var(TIME_cut, "MOMY", -0.03, 0.03, nz, nc)
    #plot_var(TIME_cut, "DDENS", -0.0005, 0.0005, nz, nc)    
    #plot_var(TIME_cut, "MOMZ", -0.001, 0.001, nz, nc)    
    #plot_var(TIME_cut, "DRHOT", -0.06, 0.06, nz, nc)
#'''

'''
for n_t in range(0,17):
  TIME_cut = 900*n_t
  plot_var_yz(TIME_cut, "MOMX", -3e-2, 3e-2, 6*8+1, nc, True)
  plot_var_yz(TIME_cut, "MOMY", -3e-2, 3e-2, 6*8+1, nc) 
  plot_var_yz(TIME_cut, "MOMZ", -0.00024, 0.00024, 6*8+1, nc) 
  plot_var_yz(TIME_cut, "W", -0.005, 0.005, 6*8+1, nc) 
  plot_var_yz(TIME_cut, "V", -0.3, 0.3, 6*8+1, nc) 
  plot_var_yz(TIME_cut, "DDENS", -0.0005, 0.0005, 6*8+1, nc) 
  plot_var_yz(TIME_cut, "DRHOT", -0.06, 0.06, 6*8+1, nc) 
'''