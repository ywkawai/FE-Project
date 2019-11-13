import numpy as np
import xarray as xr
import matplotlib
#from IPython.display import set_matplotlib_formats
#set_matplotlib_formats('retina')
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import matplotlib.animation as animation
import netCDF4
import glob

def hori_1Daxis_fmt(tick_val, pos):
  val = int(tick_val)
  return f'{val}'

def set_fig_XY_axis(ax):
  ax.set_xlabel('X (m)')  
  ax.set_xlim(0.0, 1.0e3)
  ax.xaxis.set_major_formatter(tick.FuncFormatter(hori_1Daxis_fmt))
  ax.xaxis.set_major_locator(tick.MultipleLocator(0.5e3))
  ax.xaxis.set_minor_locator(tick.MultipleLocator(0.1e3))
  ax.set_ylabel('Z (m)')
  ax.set_ylim(0.0, 1.5e3)
  ax.yaxis.set_major_formatter(tick.FuncFormatter(hori_1Daxis_fmt))
  ax.yaxis.set_major_locator(tick.MultipleLocator(0.5e3))
  ax.yaxis.set_minor_locator(tick.MultipleLocator(0.1e3))

#nc = xr.open_dataset('history.pe000000.nc', decode_times=False)
nc = xr.open_mfdataset('history.pe000*.nc', decode_times=False, combine='by_coords')
dtheta = nc.DTHETA
u = nc.U
w = nc.W
x = nc.x
z = nc.y
time = nc.time

X, Z = np.meshgrid(x, z)
fig = plt.figure(figsize=(9,10)) 
ax = fig.add_subplot(1,1,1)
set_fig_XY_axis(ax)
ax.set_title("Rising thermal bubble ($\Theta', U, W$)")

pcm = ax.pcolormesh(X, Z, dtheta.isel(time=0), vmin=-0.2, vmax=+0.6, cmap='jet')#cmap='YlGnBu_r')
fmt = tick.ScalarFormatter(useMathText=True)
fmt.set_powerlimits((0,0))
cbar = plt.colorbar(pcm, extend='both', orientation='vertical', shrink=0.8, aspect=30.0)#, format=fmt)
cbar.ax.xaxis.get_offset_text().set_fontsize(6)

cnt_lv = np.arange(-0.1, 0.5, 0.05)
cnt_lv = np.delete(cnt_lv, 2)
cnt_opts = {'levels': cnt_lv, 'colors': ['black'], 'linewidths': 1, 'alpha':0.5}
matplotlib.rcParams['contour.negative_linestyle']= 'solid'
cnt = ax.contour(X, Z, dtheta.isel(time=0), **cnt_opts)

u0 = 0.0
v_xslice_int = 5
v_zslice_int = 5
q = ax.quiver(X[::v_zslice_int,::v_xslice_int], Z[::v_zslice_int,::v_xslice_int], 
              (u.isel(time=0)[::v_zslice_int,::v_xslice_int]-u0), 
              w.isel(time=0)[::v_zslice_int,::v_xslice_int], 
              pivot='mid', angles='xy', scale_units='xy', scale=0.03, alpha=0.6, linewidth=0.2 )

def update(i):
  global cnt, pcm
  global v_xslice_int, v_zslice_int, u0

  print(f"time={i}[s]")
  dtheta_ = dtheta.sel(time=i)
  u_ = u.sel(time=i)[::v_zslice_int,::v_xslice_int]
  w_ = w.sel(time=i)[::v_zslice_int,::v_xslice_int]
  
  pcm.set_array(dtheta_.values[:-1,:-1].ravel())
  q.set_UVC((u_-u0), w_)

  for c in cnt.collections:
    c.remove()
  cnt = ax.contour(X, Z, dtheta_, **cnt_opts)

  ax.set_title(f"Rising thermal bubble ($\Theta', U, W$) t={i}[s]")

  

ani = animation.FuncAnimation(fig, update, frames=time.values, interval=1000*10)
ani.save('rising_thermal_bubble.mp4', writer="ffmpeg", fps=8)
#ani.save('rising_thermal_bubble.gif', writer="imagemagick", fps=8)
#plt.show()
