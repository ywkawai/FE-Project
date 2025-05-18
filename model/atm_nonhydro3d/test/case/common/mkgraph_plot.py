import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from mpl_toolkits.axes_grid1 import make_axes_locatable
#import matplotlib.animation as animation

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

def plot_var_xy(var, vmin, vmax, zisel, nc):
  x = var.coords["x"]
  y = var.coords["y"]

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

def plot(da, exch=False, vmin=None, vmax=None, xlim=None, ylim=None, vint=None, cmap="jet", title=None, figsize=None):
  cut_dim_info_str = []
  for k, v in da.coords.items(): 
    if k not in da.dims:
      cut_dim_info_str.append(f"{k}={v.values:.2f}")
  if (len(cut_dim_info_str) > 0 ):
    cut_dim_info_str = ", ".join(cut_dim_info_str)
    title_pad = 18
  else:
    cut_dim_info_str = ""
    title_pad = 5
#  print(cut_dim_info_str)

  if exch:
    coord1 = da.coords[da.dims[0]]
    coord2 = da.coords[da.dims[1]]
  else:
    coord1 = da.coords[da.dims[1]]
    coord2 = da.coords[da.dims[0]]

  yaxis_name = coord2.name

  if figsize:
    fig_size = figsize
  else:
    w = coord1.values[-1]-coord1.values[0]
    h = coord2.values[-1]-coord2.values[0]
    h_ov_w = h/w
    print(f"h_ov_w: {h_ov_w}")
    width = 8
    fig_size = (width,width*h_ov_w)
    
  print(f"fig_size: {fig_size}")
  print(f"yaxis_name: {yaxis_name}")
  plt.rcParams["font.size"] = 14
  fig, ax = plt.subplots(1,1, figsize=fig_size) 

  if vint:
    lv = np.arange(vmin, vmax + vint, vint)
    cs = da.plot.contourf(y=yaxis_name, ax=ax, levels=lv, cmap=cmap, add_colorbar=False)
  else:
    cs = da.plot.contourf(y=yaxis_name, ax=ax, vmin=vmin, vmax=vmax, cmap=cmap, add_colorbar=False)
  
  cs2 = da.plot.contour(y=yaxis_name, ax=ax, levels=cs.levels, colors='black', add_colorbar=False)
  ax.text(0.5, 1.0, cut_dim_info_str, transform=ax.transAxes, ha='center', va='bottom', fontsize=12)

  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="3%", pad=0.2)  
  fig.colorbar(cs, cax=cax)

  if title:
    print(title)
    ax.set_title(title, fontsize=18, pad=title_pad)

  if xlim:
    ax.set_xlim( [xlim[0], xlim[1]] )
  if ylim:
    ax.set_ylim( [ylim[0], ylim[1]] )

  ax.set_xlabel(f"{coord1.name} [{coord1.units}]", fontsize=16)
  ax.set_ylabel(f"{coord2.name} [{coord2.units}]", fontsize=16)
  # ax.set_aspect('equal')

  return fig