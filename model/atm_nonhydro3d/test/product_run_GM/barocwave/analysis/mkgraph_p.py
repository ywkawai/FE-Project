#import pandas as pd
import xarray as xr
import numpy as np
import os
import mkgraph_sub as mg

EXP_DIR="./rhot_hevi"

RESOL_LIST = [
# p=3
    # "Eh10Ez8P3", 
    # "Eh20Ez8P3",     
    # "Eh40Ez8P3",     
    # "Eh80Ez8P3",         
# p=7
    #"Eh5Ez4P7", 
    #"Eh10Ez4P7", 
    #"Eh20Ez4P7", 
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

def create_fig_temp(dir1, dir2, out_fig_dir):
    temp1 = open_tmp_nc(dir1, "T").isel(p=7)
    temp2 = open_tmp_nc(dir2, "T").isel(p=7)

    for day in [0, 4]:
        mg.create_graph_zoom( temp1.sel(time=day*86400.0, method="nearest"),
                    220.0, 310.0, 2.5, np.linspace(220,310,10), 20, 0.0, 
                    f"{out_fig_dir}/temp_p850hPa_day{day}_zoom.pdf")
        
    for day in [6]:
        mg.create_graph_zoom( temp2.sel(time=day*86400.0, method="nearest"), 
                    220.0, 310.0, 2.5, np.linspace(220,310,10), 20, 0.0, 
                    f"{out_fig_dir}/temp_p850hPa_day{day}_zoom.pdf")
        
    for day in [8, 9, 10]:
        mg.create_graph_zoom( temp2.sel(time=day*86400.0, method="nearest"), 
                    220.0, 310.0, 2.5, np.linspace(220,310,10), 20, 0.0,
                    f"{out_fig_dir}/temp_p850hPa_day{day}_zoom.pdf")
        
# def create_fig_vor(dir1, dir2, out_fig_dir):
#     #vor1 = open_tmp_nc(dir1, "VOR", "diag").isel(p=7)
#     vor2 = open_tmp_nc(dir2, "VOR", "diag").isel(p=7)

#     for day in [0, 4]:
#         create_graph_zoom( vor1.sel(time=day*86400.0, method="nearest"),
#                     220.0, 310.0, np.linspace(220,310,10),
#                     f"{out_fig_dir}/temp_p850hPa_day{day}_zoom.pdf")
        
#     for day in [7]:
#         create_graph_zoom( vor2.sel(time=day*86400.0, method="nearest"), 
#                     -3e-5, 6e-5, np.linspace(-3e-5,6e-5,10),
#                     f"{out_fig_dir}/vor_p850hPa_day{day}_zoom.pdf")
        
#     for day in [8, 9, 10]:
#         create_graph_zoom( vor2.sel(time=day*86400.0, method="nearest"), 
#                     -1e-4, 3.5e-4, np.linspace(-1e-4,3.5e-4,10),
#                     f"{out_fig_dir}/ vor_p850hPa_day{day}_zoom.pdf", 
#                     False )

for resol in RESOL_LIST:
    dir1=f"{EXP_DIR}/{resol}_1/outdata_p"
    dir2=f"{EXP_DIR}/{resol}_2/outdata_p"
    dir_ini=f"{EXP_DIR}/{resol}_1/outdata_p"
    out_fig_dir=f"analysis_out/{EXP_DIR}/{resol}"
    
    os.makedirs(out_fig_dir, exist_ok=True)    
    create_fig_temp(dir1, dir2, out_fig_dir)
#    create_fig_vor(dir1, dir2, out_fig_dir)    