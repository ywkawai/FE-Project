#import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.ticker as tick
import matplotlib.pyplot as plt
import os

EXP_DIR="./rhot_hevi"
#EXP_DIR="./rhot_hevi"

RESOL_LIST = [
# p=3
    # "Eh10Ez8P3", 
    # "Eh20Ez8P3",     
    # "Eh40Ez8P3",     
    # "Eh80Ez8P3",     
    # "Eh160Ez8P3",         
# p=7
    # "Eh5Ez4P7", 
    # "Eh10Ez4P7", 
    # "Eh20Ez4P7", 
    # "Eh40Ez4P7",
    #"Eh80Ez4P7",     
# p=11
    "Eh3Ez3P11",     
    # "Eh6Ez3P11",     
    # "Eh12Ez3P11",         
]

#---------------------------------------------
def open_bs_tmp_nc(dir, varname):
  print(f'{dir}/bs.pe*.nc')
  return xr.open_mfdataset(f'{dir}/bs.pe*.nc', decode_times=False, combine='by_coords')[varname]

def open_tmp_nc(dir, varname):
  print(f'{dir}/history.pe*.nc')
  return xr.open_mfdataset(f'{dir}/history.pe*.nc', decode_times=False, combine='by_coords')[varname]


def hori_1Daxis_fmt(tick_val, pos):
  # val = int(tick_val/1000)
  val = int(tick_val)
  return f'{val}'

def set_fig_XY_axis(ax):
  ax.tick_params(labelsize=18, length=10)  
  ax.set_xlabel('longitude [deg]', fontsize=18)
  ax.xaxis.set_major_formatter(tick.FuncFormatter(hori_1Daxis_fmt))
  ax.xaxis.set_major_locator(tick.MultipleLocator(30.0))
  ax.xaxis.set_minor_locator(tick.MultipleLocator(10.0))
  ax.set_ylabel('latitude [deg]', fontsize=18)
  ax.yaxis.set_major_formatter(tick.FuncFormatter(hori_1Daxis_fmt))
  ax.yaxis.set_major_locator(tick.MultipleLocator(30.0))
  ax.yaxis.set_minor_locator(tick.MultipleLocator(5.0))


def create_graph_zoom(data, vmin, vmax, cnt_levels, cnt_zero_val, out_pngname):
    print(f"mkgraph {out_pngname}")
    
    lon = data.coords["lon"]
    lat = data.coords["lat"]

    Lon, Lat = np.meshgrid(lon,lat)    
    fig = plt.figure(figsize=(13,5))
    fig.subplots_adjust(bottom=0.2)
    
    ax = fig.add_subplot(111)
    ax.set_xlim(45,360)
    ax.set_ylim(0,89)
    set_fig_XY_axis(ax)

    pcm = ax.pcolormesh(Lon, Lat, data.values, vmin=vmin, vmax=vmax, cmap='jet')
    levels = cnt_levels[cnt_levels != cnt_zero_val]    
    cont = ax.contour(Lon, Lat, data.values, levels=levels, colors=['black'])  
#    cont.clabel(fmt='%2.0f', fontsize=12)
    cont = ax.contour(Lon, Lat, data.values, levels=[cnt_zero_val], colors=['darkgray'], linewidths=2.0)  
    
    fmt = tick.ScalarFormatter(useMathText=True)
    cbar = plt.colorbar(pcm, aspect=30.0, extend='both', orientation='vertical', shrink=0.8, format=fmt, pad=0.02)
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.yaxis.get_offset_text().set_fontsize(14)
    
    plt.savefig(out_pngname)

def create_graph_global(data, vmin, vmax, levels, out_pngname):
    print(f"mkgraph {out_pngname}")
    fig = plt.figure(figsize=(15,8))
    ax = fig.add_subplot(111)

    data.plot(ax=ax, cmap="jet", vmin=vmin,vmax=vmax)
    data.plot.contour(ax=ax,levels=levels,colors='black')

    ax.set_xlim(0,360)
    ax.set_ylim(-89,89)
    set_fig_XY_axis(ax)
    
    plt.savefig(out_pngname)

def create_fig_sfcpre(dir1, dir2, out_fig_dir):
    sfc_pres1 = open_tmp_nc(dir1, "PRES").isel(z=0)
    sfc_pres2 = open_tmp_nc(dir2, "PRES").isel(z=0)

    for day in [4]:
        create_graph_zoom( sfc_pres1.sel(time=day*86400.0, method="nearest") / 1e2, 
                    992.0, 1006.0, np.linspace(992,1006,15), 1000.0, 
                    f"{out_fig_dir}/ps_day{day}_zoom.png")
        
    for day in [6]:
        create_graph_zoom( sfc_pres2.sel(time=day*86400.0, method="nearest") / 1e2, 
                    992.0, 1006.0, np.linspace(992,1006,15), 1000.0, 
                    f"{out_fig_dir}/ps_day{day}_zoom.png")
    for day in [7]:        
        create_graph_global( sfc_pres2.sel(time=day*86400.0, method="nearest") / 1e2, 
                    980.0, 1012.0, np.linspace(980,1012,9), 
                    f"{out_fig_dir}/ps_day{day}_global.png")
        
    for day in [8, 9, 10]:
        create_graph_zoom( sfc_pres2.sel(time=day*86400.0, method="nearest") / 1e2, 
                    930.0, 1030.0, np.linspace(930,1030,11), 1000.0, 
                    f"{out_fig_dir}/ps_day{day}_zoom.png")
        create_graph_global( sfc_pres2.sel(time=day*86400.0, method="nearest") / 1e2, 
                    930.0, 1030.0, np.linspace(930,1030,11), 
                    f"{out_fig_dir}/ps_day{day}_global.png")
        

for resol in RESOL_LIST:
    dir1=f"{EXP_DIR}/{resol}_1/outdata"
    dir2=f"{EXP_DIR}/{resol}_2/outdata"
    dir_ini=f"{EXP_DIR}/{resol}_1/outdata"
    out_fig_dir=f"analysis_out/{EXP_DIR}/{resol}"
    
    os.makedirs(out_fig_dir, exist_ok=True)    
    create_fig_sfcpre(dir1, dir2, out_fig_dir)