import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import xarray as xr
import numpy as np
import os
import joblib


def open_tmp_nc(dir, varname):
  print(f'{dir}/history.pe*.nc')
  return xr.open_mfdataset(f'{dir}/history.pe*.nc', decode_times=False, combine='by_coords')[varname]

def hori_1Daxis_fmt(tick_val, pos):
  #val = int(tick_val/1000)
  val = int(tick_val)
  return f'{val}'
def vert_1Daxis_fmt(tick_val, pos):
  val = int(tick_val/1000)
  return f'{val}'

def set_fig_XY_axis(ax, xlabel):
  ax.tick_params(labelsize=18, length=10)  
  print("xlabel:", xlabel)
  if (xlabel):  
    ax.set_xlabel('longitude [deg]', fontsize=18)
  ax.xaxis.set_major_formatter(tick.FuncFormatter(hori_1Daxis_fmt))
  ax.xaxis.set_major_locator(tick.MultipleLocator(30.0))
  ax.xaxis.set_minor_locator(tick.MultipleLocator(10.0))
  ax.set_ylabel('latitude [deg]', fontsize=18)
  ax.yaxis.set_major_formatter(tick.FuncFormatter(hori_1Daxis_fmt))
  ax.yaxis.set_major_locator(tick.MultipleLocator(30.0))
  ax.yaxis.set_minor_locator(tick.MultipleLocator(5.0))

def set_fig_XZ_axis(ax, xlabel):
  ax.tick_params(labelsize=18, length=10)  
  print("xlabel:", xlabel)
  if (xlabel):  
    ax.set_xlabel('longitude [deg]', fontsize=18)
  ax.xaxis.set_major_formatter(tick.FuncFormatter(hori_1Daxis_fmt))
  ax.xaxis.set_major_locator(tick.MultipleLocator(30.0))
  ax.xaxis.set_minor_locator(tick.MultipleLocator(10.0))
  ax.set_ylabel('height [km]', fontsize=18)
  ax.yaxis.set_major_formatter(tick.FuncFormatter(vert_1Daxis_fmt))
  ax.yaxis.set_major_locator(tick.MultipleLocator(5000.0))
  ax.yaxis.set_minor_locator(tick.MultipleLocator(250.0))

def set_fig_YZ_axis(ax, ylabel):
  ax.tick_params(labelsize=18, length=10)  
  print("ylabel:", ylabel)
  if (ylabel):  
    ax.set_xlabel('latitude [deg]', fontsize=18)
  ax.xaxis.set_major_formatter(tick.FuncFormatter(hori_1Daxis_fmt))
  ax.xaxis.set_major_locator(tick.MultipleLocator(30.0))
  ax.xaxis.set_minor_locator(tick.MultipleLocator(10.0))
  ax.set_ylabel('height [km]', fontsize=18)
  ax.yaxis.set_major_formatter(tick.FuncFormatter(vert_1Daxis_fmt))
  ax.yaxis.set_major_locator(tick.MultipleLocator(5000.0))
  ax.yaxis.set_minor_locator(tick.MultipleLocator(250.0))

def set_fig_XZ_axis_1D(ax, xlabel):
  ax.tick_params(labelsize=18, length=10)  
#  print("ylabel:", xlabel)
  if (xlabel):  
    ax.set_xlabel('lon [deg]', fontsize=18)
  ax.xaxis.set_major_formatter(tick.FuncFormatter(hori_1Daxis_fmt))
  ax.xaxis.set_minor_formatter(plt.NullFormatter())  
  ax.xaxis.set_major_locator(tick.MultipleLocator(2))
  ax.xaxis.set_minor_locator(tick.MultipleLocator(.2))
  ax.set_ylabel('Error', fontsize=18)

def create_graph_numerror1D(data_list, vmin0, vmax0, scale_ord, lon_lim, width=13):
    print(f"mkgraph ..")
    row_num = len(data_list)
    col_num = row_num
    fig, axes = plt.subplots(1,col_num, figsize=(width*col_num,1*row_num), sharex=True)
    fig.subplots_adjust(bottom=0.2, wspace=0.02, hspace=0.05)
    for j in range(0,col_num):
      vmin=vmin0/float((j+1)**scale_ord); vmax=vmax0/float((j+1)**scale_ord)      
      for i in range(0,1):
        ax = axes[j]    
        set_fig_XZ_axis_1D(ax, i+1==row_num)

        data = data_list[j]        
        ax.plot(data.coords["lon"], data.values, marker="o")
        ax.set_xlim(lon_lim[0],lon_lim[1])
        ax.set_ylim(vmin,vmax)   

def create_graph_zoom(data_list, vmin, vmax, cnt_levels, cnt_zero_val, lon_lim, z_lim, width=13):
    print(f"mkgraph ..")
    row_num = len(data_list)

    fig, axes = plt.subplots(row_num,1, figsize=(width,5*row_num), sharex=True, sharey=True)
    fig.subplots_adjust(bottom=0.2)
    for i, ax in enumerate(axes):    
      ax.set_xlim(lon_lim[0],lon_lim[1])
      ax.set_ylim(z_lim[0],z_lim[1])
      set_fig_XZ_axis(ax, i+1==row_num)

      data = data_list[i]
      
      lon = data.coords["lon"]
      z = data.coords["z"]

      Lon, Z = np.meshgrid(lon,z)          
      pcm = ax.pcolormesh(Lon, Z, data.values, vmin=vmin, vmax=vmax, cmap='jet')
      if (cnt_levels is not None):
        levels = cnt_levels[cnt_levels != cnt_zero_val]    
        cont = ax.contour(Lon, Z, data.values, levels=levels, colors=['black'])  
    #    cont.clabel(fmt='%2.0f', fontsize=12)
        cont = ax.contour(Lon, Z, data.values, levels=[cnt_zero_val], colors=['darkgray'], linewidths=1.0)  
      
      fmt = tick.ScalarFormatter(useMathText=True)
      cbar = fig.colorbar(pcm, ax=ax, aspect=30.0, extend='both', orientation='vertical', shrink=0.8, format=fmt, pad=0.02)
      cbar.ax.tick_params(labelsize=14)
      cbar.ax.yaxis.get_offset_text().set_fontsize(14)
    
def create_graph_zoom2(data_list, vmin, vmax, cnt_levels, cnt_zero_val, lat_lim, z_lim, width=13):
    print(f"mkgraph ..")
    row_num = len(data_list)

    fig, axes = plt.subplots(row_num,1, figsize=(width,5*row_num), sharex=True, sharey=True)
    fig.subplots_adjust(bottom=0.2)
    for i, ax in enumerate(axes):    
      ax.set_xlim(lat_lim[0],lat_lim[1])
      ax.set_ylim(z_lim[0],z_lim[1])
      set_fig_XZ_axis(ax, i+1==row_num)

      data = data_list[i]
      
      lat = data.coords["lat"]
      z = data.coords["z"]

      Lat, Z = np.meshgrid(lat,z)          
      pcm = ax.pcolormesh(Lat, Z, data.values, vmin=vmin, vmax=vmax, cmap='jet')
      if (cnt_levels is not None):
        levels = cnt_levels[cnt_levels != cnt_zero_val]    
        cont = ax.contour(Lat, Z, data.values, levels=levels, colors=['black'])  
    #    cont.clabel(fmt='%2.0f', fontsize=12)
        cont = ax.contour(Lat, Z, data.values, levels=[cnt_zero_val], colors=['darkgray'], linewidths=1.0)  
      
      fmt = tick.ScalarFormatter(useMathText=True)
      cbar = fig.colorbar(pcm, ax=ax, aspect=30.0, extend='both', orientation='vertical', shrink=0.8, format=fmt, pad=0.02)
      cbar.ax.tick_params(labelsize=14)
      cbar.ax.yaxis.get_offset_text().set_fontsize(14)
      
def create_graph_zoom_hori(data_list, vmin, vmax, cnt_levels, cnt_zero_val, lon_lim, lat_lim, width=13):
    print(f"mkgraph ..")
    row_num = len(data_list)

    fig, axes = plt.subplots(row_num,1, figsize=(width,5*row_num), sharex=True, sharey=True)
    fig.subplots_adjust(bottom=0.2)
    for i, ax in enumerate(axes):    
      ax.set_xlim(lon_lim[0],lon_lim[1])
      ax.set_ylim(lat_lim[0],lat_lim[1])
      set_fig_XY_axis(ax, i+1==row_num)

      data = data_list[i]
      
      lon = data.coords["lon"]
      lat = data.coords["lat"]

      Lon, Lat = np.meshgrid(lon,lat)          
      pcm = ax.pcolormesh(Lon, Lat, data.values, vmin=vmin, vmax=vmax, cmap='jet')
      if (cnt_levels is not None):
        levels = cnt_levels[cnt_levels != cnt_zero_val]    
        cont = ax.contour(Lon, Lat, data.values, levels=levels, colors=['black'])  
    #    cont.clabel(fmt='%2.0f', fontsize=12)
        cont = ax.contour(Lon, Lat, data.values, levels=[cnt_zero_val], colors=['darkgray'], linewidths=1.0)  
      
      fmt = tick.ScalarFormatter(useMathText=True)
      cbar = fig.colorbar(pcm, ax=ax, aspect=30.0, extend='both', orientation='vertical', shrink=0.8, format=fmt, pad=0.02)
      cbar.ax.tick_params(labelsize=14)
      cbar.ax.yaxis.get_offset_text().set_fontsize(14)
    
#    plt.savefig(out_pngname)

def create_graph_numerror(data_list, vmin0, vmax0, scale_ord, lon_lim, z_lim, width=13):
    print(f"mkgraph ..")
    row_num = len(data_list)

    fig, axes = plt.subplots(row_num,3, figsize=(width*3,5*row_num), sharex=True, sharey=True)
    fig.subplots_adjust(bottom=0.2, wspace=0.02, hspace=0.05)
    for j in range(0,3):
      fac = 0.5**(j*scale_ord)
      vmin=vmin0*fac; vmax=vmax0*fac; contour_levs=np.linspace(vmin,vmax,21)      
      print(f"vmin: {vmin}, vmax: {vmax}, j={j}, fac={1.0/fac}")
            
      for i in range(0,row_num):
        ax = axes[i,j]    
        ax.set_xlim(lon_lim[0],lon_lim[1])
        ax.set_ylim(z_lim[0],z_lim[1])
        set_fig_XZ_axis(ax, i+1==row_num)

        data = data_list[i]
        
        lon = data.coords["lon"]
        z = data.coords["z"]

        Lon, Z = np.meshgrid(lon,z)          
        pcm = ax.pcolormesh(Lon, Z, data.values, vmin=vmin, vmax=vmax, cmap='jet')
        fmt = tick.ScalarFormatter(useMathText=True)
        cbar = fig.colorbar(pcm, ax=ax, aspect=30.0, extend='both', orientation='vertical', shrink=0.8, format=fmt, pad=0.02)
        cbar.ax.tick_params(labelsize=14)
        cbar.ax.yaxis.get_offset_text().set_fontsize(14)

def create_graph_numerror2(data_list, vmin0, vmax0, scale_ord, lon_lim, z_lim, width=13):
    print(f"mkgraph ..")
    row_num = len(data_list)

    fig, axes = plt.subplots(row_num,3, figsize=(width*3,5*row_num), sharex=True, sharey=True)
    fig.subplots_adjust(bottom=0.2, wspace=0.02, hspace=0.05)
    for j in range(0,3):
      fac = 0.5**(j*scale_ord)
      vmin=vmin0*fac; vmax=vmax0*fac; contour_levs=np.linspace(vmin,vmax,21)      
      print(f"vmin: {vmin}, vmax: {vmax}, j={j}, fac={1.0/fac}")
      
      for i in range(0,row_num):
        ax = axes[i,j]    
        ax.set_xlim(lon_lim[0],lon_lim[1])
        ax.set_ylim(z_lim[0],z_lim[1])
        set_fig_XZ_axis(ax, i+1==row_num)

        data = data_list[i]
        
        lat = data.coords["lat"]
        z = data.coords["z"]

        Lat, Z = np.meshgrid(lat,z)          
        pcm = ax.pcolormesh(Lat, Z, data.values, vmin=vmin, vmax=vmax, cmap='jet')
        fmt = tick.ScalarFormatter(useMathText=True)
        cbar = fig.colorbar(pcm, ax=ax, aspect=30.0, extend='both', orientation='vertical', shrink=0.8, format=fmt, pad=0.02)
        cbar.ax.tick_params(labelsize=14)
        cbar.ax.yaxis.get_offset_text().set_fontsize(14)

def create_graph_numerror_hori(data_list, vmin0, vmax0, scale_ord, lon_lim, lat_lim, width=13):
    print(f"mkgraph ..")
    row_num = len(data_list)

    fig, axes = plt.subplots(row_num,3, figsize=(width*3,5*row_num), sharex=True, sharey=True)
    fig.subplots_adjust(bottom=0.2, wspace=0.02, hspace=0.05)
    for j in range(0,3):
      for i in range(0,row_num):
        ax = axes[i,j]    
        ax.set_xlim(lon_lim[0],lon_lim[1])
        ax.set_ylim(lat_lim[0],lat_lim[1])
        set_fig_XY_axis(ax, i+1==row_num)

        vmin=vmin0/float((j+1)**scale_ord); vmax=vmax0/float((j+1)**scale_ord); contour_levs=np.linspace(vmin,vmax,21)

        data = data_list[i]
        
        lon = data.coords["lon"]
        lat = data.coords["lat"]

        Lon, Lat = np.meshgrid(lon,lat)          
        pcm = ax.pcolormesh(Lon, Lat, data.values, vmin=vmin, vmax=vmax, cmap='jet')
        fmt = tick.ScalarFormatter(useMathText=True)
        cbar = fig.colorbar(pcm, ax=ax, aspect=30.0, extend='both', orientation='vertical', shrink=0.8, format=fmt, pad=0.02)
        cbar.ax.tick_params(labelsize=14)
        cbar.ax.yaxis.get_offset_text().set_fontsize(14)


#------------------------------------

def output_merge_nc_y0_sub(outdata_dir, tmp_dir, varname, time_sec_list):
  os.makedirs(tmp_dir, exist_ok=True)
  var = open_tmp_nc(outdata_dir, varname).sel(lat=0).isel(lat=0)
  for time_sec in time_sec_list:
    print(f"Output: {tmp_dir}/{varname}.nc, time_sec={time_sec}")
    var.sel(time=time_sec).to_netcdf(f"{tmp_dir}/{varname}_y0_t{time_sec}.nc")
  
def output_merge_nc_y0(outdata_dir_list, tmp_dir_list, varname, time_sec_list):
  result = joblib.Parallel(n_jobs=4, verbose=2)(joblib.delayed(output_merge_nc_y0_sub)(outdata_dir, tmp_dir_list[i], varname, time_sec_list) for i, outdata_dir in enumerate(outdata_dir_list) )

def output_merge_nc_cutlon_sub(outdata_dir, tmp_dir, varname, time_sec_list, lon):
  os.makedirs(tmp_dir, exist_ok=True)
  var_ = open_tmp_nc(outdata_dir, varname)
  lonindx = get_index(var_, "lon", lon)
  var = var_.isel(lon=lonindx)
  
  for time_sec in time_sec_list:
    print(f"Output: {tmp_dir}/{varname}.nc, time_sec={time_sec}")
    var.sel(time=time_sec).to_netcdf(f"{tmp_dir}/{varname}_lon{lon}_t{time_sec}.nc")
  
def output_merge_nc_cutlon(outdata_dir_list, tmp_dir_list, varname, time_sec_list, lon):
  result = joblib.Parallel(n_jobs=4, verbose=2)(joblib.delayed(output_merge_nc_cutlon_sub)(outdata_dir, tmp_dir_list[i], varname, time_sec_list, lon) for i, outdata_dir in enumerate(outdata_dir_list) )


def get_index(var, coord_name, target_pos):
  # print(var)
  coord = var.coords[coord_name].values
  return np.argmin((coord-target_pos)**2)

def output_merge_nc_hori_sub(outdata_dir, tmp_dir, varname, time_sec_list, zlev):
  os.makedirs(tmp_dir, exist_ok=True)
  var_ = open_tmp_nc(outdata_dir, varname)
  zindx = get_index(var_, "z", zlev)
  var = var_.isel(z=zindx)
  
  for time_sec in time_sec_list:
    print(f"Output: {tmp_dir}/{varname}.nc, time_sec={time_sec}")
    var.sel(time=time_sec).to_netcdf(f"{tmp_dir}/{varname}_z{zlev}_t{time_sec}.nc")
  
def output_merge_nc_hori(outdata_dir_list, tmp_dir_list, varname, time_sec_list, zlev):
  result = joblib.Parallel(n_jobs=4, verbose=2)(joblib.delayed(output_merge_nc_hori_sub)(outdata_dir, tmp_dir_list[i], varname, time_sec_list, zlev) for i, outdata_dir in enumerate(outdata_dir_list) )


def open_merge_nc(tmp_dir, varname, postfix):
  print(f'{tmp_dir}/{varname}{postfix}.nc')
  return xr.open_mfdataset(f'{tmp_dir}/{varname}{postfix}.nc', decode_times=False, combine='by_coords')[varname]

