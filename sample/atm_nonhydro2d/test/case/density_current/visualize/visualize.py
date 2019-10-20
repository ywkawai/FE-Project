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
  ax.set_xlim(0.0, 19.2e3)
  ax.xaxis.set_major_formatter(tick.FuncFormatter(hori_1Daxis_fmt))
  ax.xaxis.set_major_locator(tick.MultipleLocator(10e3))
  ax.xaxis.set_minor_locator(tick.MultipleLocator(2.5e3))
  ax.set_ylabel('Z (km)')
  ax.set_ylim(0.0, 4.8e3)
  ax.yaxis.set_major_formatter(tick.FuncFormatter(hori_1Daxis_fmt))
  ax.yaxis.set_major_locator(tick.MultipleLocator(1e3))
  ax.yaxis.set_minor_locator(tick.MultipleLocator(0.25e3))

nc = xr.open_dataset('history.pe000000.nc', decode_times=False)
dtheta = nc.DTHETA
u = nc.U
w = nc.W
x = nc.x
z = nc.y
time = nc.time

X, Z = np.meshgrid(x, z)
fig = plt.figure(figsize=(10,5)) 
ax = fig.add_subplot(1,1,1)
set_fig_XY_axis(ax)
ax.set_title("gravity wave ($\Theta', U', W$)")

pcm = ax.pcolormesh(X, Z, dtheta.isel(time=0), vmin=-12.0, vmax=+2.0, cmap='jet')#cmap='YlGnBu_r')
fmt = tick.ScalarFormatter(useMathText=True)
fmt.set_powerlimits((0,0))
cbar = plt.colorbar(pcm, aspect=40.0, extend='both', orientation='horizontal', shrink=0.5)#, format=fmt)
cbar.ax.xaxis.get_offset_text().set_fontsize(6)

cnt_lv = np.arange(-9, -0.25, 0.25)
#cnt_lv = np.delete(cnt_lv, 10)
cnt_opts = {'levels': cnt_lv, 'colors': ['black'], 'linewidths': 1, 'alpha':0.5}
matplotlib.rcParams['contour.negative_linestyle']= 'solid'
cnt = ax.contour(X, Z, dtheta.isel(time=0), **cnt_opts)

u0 = 0.0
v_xslice_int = 5
v_zslice_int = 5
q = ax.quiver(X[::v_zslice_int,::v_xslice_int], Z[::v_zslice_int,::v_xslice_int], 
              2.0*(u.isel(time=0)[::v_zslice_int,::v_xslice_int]-u0), 
              w.isel(time=0)[::v_zslice_int,::v_xslice_int], 
              pivot='mid', angles='xy', scale_units='xy', scale=0.075, alpha=0.6, linewidth=0.2 )

def update(i):
  global cnt, pcm
  global v_xslice_int, v_zslice_int, u0

  print(f"time={i}[s]")
  dtheta_ = dtheta.sel(time=i)
  u_ = u.sel(time=i)[::v_zslice_int,::v_xslice_int]
  w_ = w.sel(time=i)[::v_zslice_int,::v_xslice_int]
  
  pcm.set_array(dtheta_.values[:-1,:-1].ravel())
  q.set_UVC(2.0*(u_-u0), w_)

  for c in cnt.collections:
    c.remove()
  cnt = ax.contour(X, Z, dtheta_, **cnt_opts)

  ax.set_title(f"density current ($\Theta', U', W$) t={i}[s]")
  

ani = animation.FuncAnimation(fig, update, frames=time.values, interval=1000*10)
ani.save('density_current.mp4', writer="ffmpeg", fps=5)
#plt.show()
