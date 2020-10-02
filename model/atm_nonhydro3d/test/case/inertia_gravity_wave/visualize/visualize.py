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
  ax.xaxis.set_major_locator(tick.MultipleLocator(20e3))
  ax.xaxis.set_minor_locator(tick.MultipleLocator(10e3))
  ax.set_ylabel('Y (km)')
  ax.yaxis.set_major_formatter(tick.FuncFormatter(hori_1Daxis_fmt))
  ax.yaxis.set_major_locator(tick.MultipleLocator(20e3))
  ax.yaxis.set_minor_locator(tick.MultipleLocator(10e3))

def plot_var_xy(TIME_cut, varname, vmin, vmax, zisel, nc, ini_diff=False):

  PRES_hyd = nc.sel(time=0).isel(z=zisel)["PRES_hyd"]
  DENS_hyd = nc.sel(time=0).isel(z=zisel)["DENS_hyd"]

  RHOT_hyd = 1e5/287.0*(PRES_hyd/1e5)**(1.0/1.41)
  DRHOT = nc.sel(time=TIME_cut).isel(z=zisel)["DRHOT"]
  DDENS = nc.sel(time=TIME_cut).isel(z=zisel)["DDENS"]


  DTHETA = (RHOT_hyd + DRHOT)/(DENS_hyd + DDENS).rename("DTHETA")
  DTHETA = DTHETA - RHOT_hyd/DENS_hyd

  print(DTHETA)
  lv = list(map(lambda i: (i-10)*2e-4, range(0,21)))
  print(lv)
  
  x = DRHOT.coords["x"]
  y = DRHOT.coords["y"]
  z = DRHOT.coords["z"]
  time = nc.time

  zlev = z.values

  X, Y = np.meshgrid(x, y)
  fig = plt.figure(figsize=(10,10)) 
  ax = fig.add_subplot(1,1,1)
  set_fig_XY_axis(ax)
  ax.set_title(f"inertia gravity wave (z={int(zlev)}m)")

  pcm = ax.pcolormesh(X, Y, DTHETA, vmin=vmin, vmax=vmax, cmap='seismic')
  fmt = tick.ScalarFormatter(useMathText=True)
  fmt.set_powerlimits((0,0))
  cbar = plt.colorbar(pcm, aspect=40.0, extend='both', orientation='horizontal', shrink=0.5, format=fmt)
  cbar.ax.xaxis.get_offset_text().set_fontsize(6)

  ax.contour(X,Y,DTHETA, levels=lv, colors=['black'])  

  plt.savefig(f"DTHETA_t{TIME_cut}_z{int(zlev)}.png")

def plot_var_xz(TIME_cut, varname, vmin, vmax, yisel, nc, ini_diff=False):

  PRES_hyd = nc.sel(time=0).isel(y=yisel)["PRES_hyd"]
  DENS_hyd = nc.sel(time=0).isel(y=yisel)["DENS_hyd"]

  RHOT_hyd = 1e5/287.0*(PRES_hyd/1e5)**(1.0/1.41)
  DRHOT = nc.sel(time=TIME_cut).isel(y=yisel)["DRHOT"]
  DDENS = nc.sel(time=TIME_cut).isel(y=yisel)["DDENS"]


  DTHETA = (RHOT_hyd + DRHOT)/(DENS_hyd + DDENS).rename("DTHETA")
  DTHETA = DTHETA - RHOT_hyd/DENS_hyd

  print(DTHETA)
  lv = list(map(lambda i: (i-10)*2e-4, range(0,21)))
  print(lv)
  
  x = DRHOT.coords["x"]
  y = DRHOT.coords["y"]
  z = DRHOT.coords["z"]
  time = nc.time

  zlev = z.values

  X, Z = np.meshgrid(x, z)
  fig = plt.figure(figsize=(10,10)) 
  ax = fig.add_subplot(1,1,1)
  set_fig_XY_axis(ax)
  ax.set_title(f"inertia gravity wave (y=50km)")

  pcm = ax.pcolormesh(X, Z, DTHETA, vmin=vmin, vmax=vmax, cmap='seismic')
  fmt = tick.ScalarFormatter(useMathText=True)
  fmt.set_powerlimits((0,0))
  cbar = plt.colorbar(pcm, aspect=40.0, extend='both', orientation='horizontal', shrink=0.5, format=fmt)
  cbar.ax.xaxis.get_offset_text().set_fontsize(6)

  ax.contour(X,Z,DTHETA, levels=lv, colors=['black'])  

  plt.savefig(f"DTHETA_t{TIME_cut}_y50km.png")

nc = xr.open_mfdataset('./history.pe000*.nc', decode_times=False, combine='by_coords')
for TIME_cut in [700]: #, 12600, 14400, 18000, 21600, 25200, 28800]:
  for ke_z in [5]:
    nz = 1+5*(ke_z-1)
    # plot_var(TIME_cut, "MOMX", -0.3, 0.3, nz, nc, True)
    # plot_var(TIME_cut, "MOMY", -0.3, 0.3, nz, nc)
    # plot_var(TIME_cut, "DDENS", -0.0005, 0.0005, nz, nc)    
    # plot_var(TIME_cut, "MOMZ", -0.001, 0.001, nz, nc) 
    plot_var_xy(TIME_cut, "DTHETA", -2e-3, 2e-3, nz, nc)    
    plot_var_xz(TIME_cut, "DTHETA", -2e-3, 2e-3, 75, nc)    