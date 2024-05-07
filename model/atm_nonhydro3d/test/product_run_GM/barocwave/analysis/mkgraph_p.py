#import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.ticker as tick
import matplotlib.pyplot as plt
import os

EXP_DIR="./rhot_hevi"

RESOL_LIST = [
# p=3
    # "Eh10Ez8P3", 
    # "Eh20Ez8P3",     
    # "Eh40Ez8P3",     
    # "Eh80Ez8P3",         
# p=7
    #"Eh5Ez4P7", 
    # "Eh10Ez4P7", 
    # "Eh20Ez4P7", 
    "Eh40Ez4P7",
    #"Eh80Ez4P7",     
# p=11
    #"Neh4_Nez3_P11",     
    #"Neh8_Nez3_P11", 
  #  "Neh16_Nez3_P11",              
    #  "Neh32_Nez3_P11",                
]

#---------------------------------------------
def open_bs_tmp_nc(dir, varname):
  print(f'{dir}/bs.pe*.nc')
  return xr.open_mfdataset(f'{dir}/bs.pe*.nc', decode_times=False, combine='by_coords')[varname]

def open_tmp_nc(dir, varname, ncfname="history"):
  print(f'{dir}/history.pe*.nc')
  return xr.open_mfdataset(f'{dir}/{ncfname}.pe*.nc', decode_times=False, combine='by_coords')[varname]

def cal_temp(dens_hyd, pres_hyd_, ddens, drhot):
  gam = 1004.64/717.6
  rdry = 287.04
  pre00 = 1e5
  
  pres_hyd = pres_hyd_.p
  rhot_hyd = (pres_hyd / pre00)**(1.0/gam) * pre00 / rdry

  pres = pre00 * (rdry/pre00 * (rhot_hyd + drhot))**gam
  dens = dens_hyd + ddens
  return pres / (dens * rdry)

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

def create_graph_zoom(data, vmin, vmax, levels, out_pngname, contour_flag=True):
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
    if contour_flag:
      cont = ax.contour(Lon, Lat, data.values, levels=levels, colors=['black'])  
  
    fmt = tick.ScalarFormatter(useMathText=True)
    cbar = plt.colorbar(pcm, aspect=30.0, extend='both', orientation='vertical', shrink=0.8, format=fmt, pad=0.02)
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.yaxis.get_offset_text().set_fontsize(14)  
    
    plt.savefig(out_pngname)


def create_fig_temp(dir1, dir2, out_fig_dir):
    # ddens1 = open_tmp_nc(dir1, "DDENS")
    # drhot1 = open_tmp_nc(dir1, "THERM")
    # ddens2 = open_tmp_nc(dir2, "DDENS")
    # drhot2 = open_tmp_nc(dir2, "THERM")

    # dens_hyd = open_bs_tmp_nc(dir_ini, "DENS_hyd")
    # pres_hyd = open_bs_tmp_nc(dir_ini, "PRES_hyd")

    # temp1 = cal_temp(dens_hyd, pres_hyd, ddens1, drhot1).isel(p=7)
    # temp2 = cal_temp(dens_hyd, pres_hyd, ddens2, drhot2).isel(p=7)

    temp1 = open_tmp_nc(dir1, "T").isel(p=7)
    temp2 = open_tmp_nc(dir2, "T").isel(p=7)

    for day in [0, 4]:
        create_graph_zoom( temp1.sel(time=day*86400.0, method="nearest"),
                    220.0, 310.0, np.linspace(220,310,10),
                    f"{out_fig_dir}/temp_p850hPa_day{day}_zoom.png")
        
    for day in [6]:
        create_graph_zoom( temp2.sel(time=day*86400.0, method="nearest"), 
                    220.0, 310.0, np.linspace(220,310,10),
                    f"{out_fig_dir}/temp_p850hPa_day{day}_zoom.png")
        
    for day in [8, 9, 10]:
        create_graph_zoom( temp2.sel(time=day*86400.0, method="nearest"), 
                    220.0, 310.0, np.linspace(220,310,10),
                    f"{out_fig_dir}/temp_p850hPa_day{day}_zoom.png")
        
def create_fig_vor(dir1, dir2, out_fig_dir):
    #vor1 = open_tmp_nc(dir1, "VOR", "diag").isel(p=7)
    vor2 = open_tmp_nc(dir2, "VOR", "diag").isel(p=7)

    for day in [0, 4]:
        create_graph_zoom( vor1.sel(time=day*86400.0, method="nearest"),
                    220.0, 310.0, np.linspace(220,310,10),
                    f"{out_fig_dir}/temp_p850hPa_day{day}_zoom.png")
        
    for day in [7]:
        create_graph_zoom( vor2.sel(time=day*86400.0, method="nearest"), 
                    -3e-5, 6e-5, np.linspace(-3e-5,6e-5,10),
                    f"{out_fig_dir}/vor_p850hPa_day{day}_zoom.png")
        
    for day in [8, 9, 10]:
        create_graph_zoom( vor2.sel(time=day*86400.0, method="nearest"), 
                    -1e-4, 3.5e-4, np.linspace(-1e-4,3.5e-4,10),
                    f"{out_fig_dir}/ vor_p850hPa_day{day}_zoom.png", 
                    False )

for resol in RESOL_LIST:
    dir1=f"{EXP_DIR}/{resol}_1/outdata_p"
    dir2=f"{EXP_DIR}/{resol}_2/outdata_p"
    dir_ini=f"{EXP_DIR}/{resol}_1/outdata_p"
    out_fig_dir=f"analysis_out/{EXP_DIR}/{resol}"
    
    os.makedirs(out_fig_dir, exist_ok=True)    
    create_fig_temp(dir1, dir2, out_fig_dir)
#    create_fig_vor(dir1, dir2, out_fig_dir)    