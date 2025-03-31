import sys
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as animation

args = sys.argv
out_filename = args[1]
ylim_min = float(args[2])
ylim_max = float(args[3])

u = xr.open_mfdataset("history.pe000000.nc", decode_times=False, combine='by_coords')["q"]
uex = xr.open_mfdataset("history.pe000000.nc", decode_times=False, combine='by_coords')["qexact"]
x = u.coords["x"]
time = u.coords["time"]

fig, ax = plt.subplots(figsize=(7,4))

ax.set_xlim(0.0, 1.0)
ax.set_ylim(ylim_min , ylim_max)
ax.set_xlabel("x")
ax.set_ylabel("u")
ax.set_title("1D linear advection") 

ims = []
flag_legend = True

for n in range(1,len(time)):
  print(f"n={n}")
  im1 = plt.plot(x, uex.isel(time=n), color="red", linestyle = "dashed", label="exact", alpha=0.5, linewidth = 3.0) 
  im2 = plt.plot(x, u.isel(time=n), color="black", label="numsol")
  
  time_txt = "{:.2f}".format(time.values[n])
  title = ax.text(0.01, ylim_max*0.95, "time="+time_txt)
  if flag_legend:
    plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
    flag_legend = False
  ims.append(im2 + im1 + [title])

ani = animation.ArtistAnimation(fig, ims, interval=200)
print( f'generate: {out_filename}' )

fig.tight_layout()
ani.save( out_filename, writer="ffmpeg" )