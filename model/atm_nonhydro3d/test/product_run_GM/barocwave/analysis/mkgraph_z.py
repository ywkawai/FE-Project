#import pandas as pd
import xarray as xr
import numpy as np
import os
import mkgraph_sub as mg

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
    #"Eh5Ez4P7", 
    #"Eh10Ez4P7", 
    #"Eh20Ez4P7", 
    "Eh40Ez4P7",
    #"Eh80Ez4P7",     
# p=11
    #"Eh3Ez3P11",     
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


def create_fig_sfcpre(dir1, dir2, out_fig_dir):
    sfc_pres1 = open_tmp_nc(dir1, "PRES").isel(z=0)
    sfc_pres2 = open_tmp_nc(dir2, "PRES").isel(z=0)

    for day in [4]:
        mg.create_graph_zoom( sfc_pres1.sel(time=day*86400.0, method="nearest") / 1e2, 
                    992.0, 1006.0, 0.5, np.linspace(992,1006,15), 2, 1000.0, 
                    f"{out_fig_dir}/ps_day{day}_zoom.pdf")
        
    for day in [6]:
        mg.create_graph_zoom( sfc_pres2.sel(time=day*86400.0, method="nearest") / 1e2, 
                    992.0, 1006.0, 0.5, np.linspace(992,1006,15), 2, 1000.0, 
                    f"{out_fig_dir}/ps_day{day}_zoom.pdf")
    for day in [7]:        
        mg.create_graph_global( sfc_pres2.sel(time=day*86400.0, method="nearest") / 1e2, 
                    980.0, 1012.0, np.linspace(980,1012,9), 
                    f"{out_fig_dir}/ps_day{day}_global.pdf")
        
    for day in [8, 9, 10]:
        mg.create_graph_zoom( sfc_pres2.sel(time=day*86400.0, method="nearest") / 1e2, 
                    930.0, 1030.0, 2.0, np.linspace(930,1030,11), 10, 1000.0, 
                    f"{out_fig_dir}/ps_day{day}_zoom.pdf")
        mg.create_graph_global( sfc_pres2.sel(time=day*86400.0, method="nearest") / 1e2, 
                    930.0, 1030.0, np.linspace(930,1030,11), 
                    f"{out_fig_dir}/ps_day{day}_global.pdf")
        

for resol in RESOL_LIST:
    dir1=f"{EXP_DIR}/{resol}_1/outdata"
    dir2=f"{EXP_DIR}/{resol}_2/outdata"
    dir_ini=f"{EXP_DIR}/{resol}_1/outdata"
    out_fig_dir=f"analysis_out/{EXP_DIR}/{resol}"
    
    os.makedirs(out_fig_dir, exist_ok=True)    
    create_fig_sfcpre(dir1, dir2, out_fig_dir)