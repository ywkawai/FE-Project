import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as animation

args = sys.argv
out_filename = args[1]
vmin = float(args[2])
vmax = float(args[3])

levels_ = np.arange(-1.2,1.2+0.2,0.2)
levels = np.delete(levels_, np.argmin(levels_**2))

u = xr.open_mfdataset("history.pe000000.nc", decode_times=False, combine='by_coords')["q"]
x = u.coords["x"]
y = u.coords["y"]
time = u.coords["time"]

fig, ax = plt.subplots(figsize=(5,4.5))

ax.set_xlim(0.0, 1.0)
ax.set_ylim(0.0, 1.0)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("2D linear advection") 

def update(n):
  print(f"n={n}")
  ax.cla()
  ax.contour(x, y, u.isel(time=n), levels=levels, colors="black", alpha=0.5)
  ax.contourf(x, y, u.isel(time=n), vmin=vmin, vmax=vmax, cmap="jet")
  
  time_txt = "{:.2f}".format(time.values[n])
  ax.text(0.01, 1.01, "time="+time_txt)

ani = animation.FuncAnimation(fig, update, frames=len(time), interval=200)
print( f'generate: {out_filename}' )

fig.tight_layout()
ani.save( out_filename, writer="ffmpeg" )