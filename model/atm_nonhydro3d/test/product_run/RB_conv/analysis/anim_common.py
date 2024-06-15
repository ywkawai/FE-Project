import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import ArtistAnimation
import os

def get_fpathlist(dir1,ftype,domainlabel,timelabel,
                  PRC_NUM_X,PRC_NUM_Y, PRC_NUM_XS,PRC_NUM_XE,PRC_NUM_YS,PRC_NUM_YE):
    return [ [ dir1 + "{0}{1}{2}.pe{3:06d}.nc".format( ftype, domainlabel, timelabel, i + j*PRC_NUM_X ) for i in range(PRC_NUM_XS-1,PRC_NUM_XE) ] for j in range(PRC_NUM_YS-1,PRC_NUM_YE) ]

def merge_xy(fpath, dim=["y","x"]):
    return xr.open_mfdataset(fpath, decode_times=False,combine="nested", concat_dim=dim)

def open_tmp_nc(dir, varname):
  print(f'{dir}/history.pe*.nc')
  return xr.open_mfdataset(f'{dir}/history.pe*.nc', decode_times=False, combine='by_coords')[varname]

def get_var_list_y0(varname, dir_list, pe_num_x, pe_nun_y, pe_xs, pe_xe, pe_ys, pe_ye):
  var_list = []
  for dir in dir_list:
    fpath = get_fpathlist(dir, "history", "","", pe_num_x, pe_nun_y, pe_xs, pe_xe, pe_ys, pe_ye)
    var = merge_xy(fpath)[varname].isel(y=0)
    var_list.append(var)
  return var_list    

def gen_anim( var_listlist, vmin, vmax, 
             anim_file):
    
  var_num = len(var_listlist)
  X_list = []; Z_list = []
  fig, axes = plt.subplots(1, var_num, figsize=(8*var_num,5))
  
  for v in range(0,var_num):
    axes[v].set_xlabel("x [m]")
    axes[v].set_ylabel("z [m]")
    var_list = var_listlist[v]
    x, z = np.meshgrid(var_list[0].x.values, var_list[0].z.values)
    X_list.append(x); Z_list.append(z)

  artists = []
  istart=0; time_offset = 0.0

  var_list_1 = var_listlist[0]
  for k in range(0,len(var_list_1)):
    time = var_list_1[k].time
    for time_ind in range(istart,len(time.values)):
        print("time="+f"{time_offset} + "+str(time.values[time_ind]))
        cax_list = []
        for v in range(0,var_num):
            var_list = var_listlist[v]
            cax = axes[v].pcolormesh(X_list[v], Z_list[v], var_list[k].isel(time=time_ind).values,
                vmin=vmin, vmax=vmax, 
                cmap = "jet"
            )
            cax_list.append(cax)
        
            if time_ind == 0:
                fig.colorbar(cax, 
                            orientation='horizontal', ax=axes[v], shrink=0.7)#, pad=0.015)                    
        
        title = axes[1].text(-500.0, 1700, "PT (time="+str(time_offset+time.values[time_ind]) + ") [s]", backgroundcolor="white", size="large")
        cax_list.append(title)
        artists.append(cax_list)
      
    istart=1; time_offset = time_offset + time.values[-1]
    
  ani = ArtistAnimation(fig, artists, interval=90)
  ani.save(anim_file, writer='imagemagick')


