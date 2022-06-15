import sys
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as animation

args = sys.argv
out_filename = args[1]
ylim_min = float(args[2])
ylim_max = float(args[3])

u = xr.open_mfdataset("history.pe000000.nc", decode_times=False, combine='by_coords')["q"]
x = u.coords["x"]
time = u.coords["time"]

fig, ax = plt.subplots(figsize=(6,4))

ax.set_xlim(0.0, 1.0)
ax.set_ylim(ylim_min , ylim_max)
ax.set_xlabel("x")
ax.set_ylabel("u")
ax.set_title("1D linear advection") 

ims = []
for n in range(0,len(time)):
  print(f"n={n}")
  im = plt.plot(x, u.isel(time=n), color="black")
  time_txt = "{:.2f}".format(time.values[n])
  title = ax.text(0.8, ylim_max*0.95, "time="+time_txt)
  ims.append(im + [title])

ani = animation.ArtistAnimation(fig, ims, interval=200)
print( f'generate: {out_filename}' )
ani.save( out_filename, writer="ffmpeg" )