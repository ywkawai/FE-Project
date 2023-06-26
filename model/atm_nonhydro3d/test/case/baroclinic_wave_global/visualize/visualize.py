import xarray as xr
import numpy as np
import matplotlib.ticker as tick
import matplotlib.pyplot as plt
import os

REGRID_DATA_DIR="../outdata"
REGRID_PCOORD_DATA__DIR="../outdata_p"
OUT_FIG_dir=f"./"

#---------------------------------------------
def open_bs_tmp_nc(dir, varname):
  print(f'{dir}/bs.pe*.nc')
  return xr.open_mfdataset(f'{dir}/bs.pe*.nc', decode_times=False, combine='by_coords')[varname]

def open_tmp_nc(dir, varname):
  print(f'{dir}/history.pe*.nc')
  return xr.open_mfdataset(f'{dir}/history.pe*.nc', decode_times=False, combine='by_coords')[varname]


def hori_1Daxis_fmt(tick_val, pos):
  val = int(tick_val)
  return f'{val}'

def set_fig_XY_axis(ax):
  ax.set_xlabel('longitude [degrees]')  
  ax.xaxis.set_major_formatter(tick.FuncFormatter(hori_1Daxis_fmt))
  ax.xaxis.set_major_locator(tick.MultipleLocator(90.0))
  ax.xaxis.set_minor_locator(tick.MultipleLocator(5.0))
  ax.set_ylabel('latitude [degrees]')
  ax.yaxis.set_major_formatter(tick.FuncFormatter(hori_1Daxis_fmt))
  ax.yaxis.set_major_locator(tick.MultipleLocator(30.0))
  ax.yaxis.set_minor_locator(tick.MultipleLocator(5.0))


def create_graph_zoom(data, vmin, vmax, levels, out_pngname):
  print(f"mkgraph {out_pngname}")
  fig = plt.figure(figsize=(15,4))
  ax = fig.add_subplot(111)

  data.plot(ax=ax, cmap="jet", vmin=vmin,vmax=vmax)
  data.plot.contour(ax=ax,levels=levels,colors='black')

  ax.set_xlim(45,360)
  ax.set_ylim(0,90)
  set_fig_XY_axis(ax)
  
  plt.savefig(out_pngname)

def create_graph_global(data, vmin, vmax, levels, out_pngname):
  print(f"mkgraph {out_pngname}")
  fig = plt.figure(figsize=(15,8))
  ax = fig.add_subplot(111)

  data.plot(ax=ax, cmap="jet", vmin=vmin,vmax=vmax)
  data.plot.contour(ax=ax,levels=levels,colors='black')

  ax.set_xlim(0,360)
  ax.set_ylim(-90,90)
  set_fig_XY_axis(ax)
    
  plt.savefig(out_pngname)

def create_fig_sfcpre(dir, out_fig_dir):
  sfc_pres = open_tmp_nc(dir, "PRES").isel(z=0)

  for day in [4]:
      create_graph_zoom( sfc_pres.sel(time=day*86400.0, method="nearest") / 1e2, 
                  992.0, 1006.0, np.linspace(992,1006,15),
                  f"{out_fig_dir}/ps_day{day}_zoom.png")
      
  for day in [6]:
      create_graph_zoom( sfc_pres.sel(time=day*86400.0, method="nearest") / 1e2, 
                  992.0, 1006.0, np.linspace(992,1006,15),
                  f"{out_fig_dir}/ps_day{day}_zoom.png")
  for day in [7]:        
      create_graph_global( sfc_pres.sel(time=day*86400.0, method="nearest") / 1e2, 
                  980.0, 1012.0, np.linspace(980,1012,9),
                  f"{out_fig_dir}/ps_day{day}_global.png")
      
  for day in [8, 9, 10]:
      create_graph_zoom( sfc_pres.sel(time=day*86400.0, method="nearest") / 1e2, 
                  930.0, 1030.0, np.linspace(930,1030,11),
                  f"{out_fig_dir}/ps_day{day}_zoom.png")
      create_graph_global( sfc_pres.sel(time=day*86400.0, method="nearest") / 1e2, 
                  930.0, 1030.0, np.linspace(930,1030,11),
                  f"{out_fig_dir}/ps_day{day}_global.png")

def search_index_plev(p, target_plev):
  return np.argmin( (p-target_plev)**2 )
  
def create_fig_temp_850hPa(dir, out_fig_dir):
  temp = open_tmp_nc(dir, "T")
  temp_ = temp.isel(p=search_index_plev(temp.coords["p"].values, 850e2))

  for day in [0, 4]:
      create_graph_zoom( temp_.sel(time=day*86400.0, method="nearest"),
                  220.0, 310.0, np.linspace(220,310,10),
                  f"{out_fig_dir}/temp_p850hPa_day{day}_zoom.png")
      
  for day in [6]:
      create_graph_zoom( temp_.sel(time=day*86400.0, method="nearest"), 
                  220.0, 310.0, np.linspace(220,310,10),
                  f"{out_fig_dir}/temp_p850hPa_day{day}_zoom.png")
      
  for day in [8, 9, 10]:
      create_graph_zoom( temp_.sel(time=day*86400.0, method="nearest"), 
                  220.0, 310.0, np.linspace(220,310,10),
                  f"{out_fig_dir}/temp_p850hPa_day{day}_zoom.png")

#-------------

os.makedirs(OUT_FIG_dir, exist_ok=True)

create_fig_sfcpre(REGRID_DATA_DIR, OUT_FIG_dir)
create_fig_temp_850hPa(REGRID_PCOORD_DATA__DIR, OUT_FIG_dir)