import pyshtools as pysh
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.ticker as tick
import numpy as np
import os
import linear_analytic_sol

EXP_list = ["Eh12Ez6P7"]

EQ_TINT_type = "rhot_heve"
OUTPUT_DATA_DIR="outdata_vis"
ANALYSIS_OUT_DIR="analysis_out/h25m"
TMP_DATA_DIR="./tmp_data/Eh12Ez6P7/"

h0=25.0
#-------------------------

def v_1Daxis_fmt(tick_val, pos):
  val = int(tick_val/1000)
  return f'{val}'

def set_fig_Xz_axis(ax):
  ax.tick_params(labelsize=18, length=8)
  ax.set_xlabel('longitude [deg]', fontsize=20)  
  ax.xaxis.set_major_locator(tick.MultipleLocator(30.0))
  ax.xaxis.set_minor_locator(tick.MultipleLocator(10.0))
  ax.set_ylabel('height [km]', fontsize=20)
  ax.yaxis.set_major_locator(tick.MultipleLocator(2000.0))
  ax.yaxis.set_minor_locator(tick.MultipleLocator(500.0))
  ax.yaxis.set_major_formatter(tick.FuncFormatter(v_1Daxis_fmt))

def set_fig_Xz_flat_axis(ax):
  ax.tick_params(labelsize=18, length=8)
  ax.set_xlabel('x [km]', fontsize=20)  
  ax.xaxis.set_major_locator(tick.MultipleLocator(10e3))
  ax.xaxis.set_minor_locator(tick.MultipleLocator(2e3))
  ax.xaxis.set_major_formatter(tick.FuncFormatter(v_1Daxis_fmt))  
  ax.set_ylabel('height [km]', fontsize=20)
  ax.yaxis.set_major_locator(tick.MultipleLocator(2000.0))
  ax.yaxis.set_minor_locator(tick.MultipleLocator(500.0))
  ax.yaxis.set_major_formatter(tick.FuncFormatter(v_1Daxis_fmt))

def set_fig_XY_axis(ax):
  ax.tick_params(labelsize=18, length=8)
  ax.set_xlabel('longitude [deg]', fontsize=20)  
  ax.xaxis.set_major_locator(tick.MultipleLocator(30.0))
  ax.xaxis.set_minor_locator(tick.MultipleLocator(10.0))
  ax.set_ylabel('latitude [deg]', fontsize=20)  
  ax.yaxis.set_major_locator(tick.MultipleLocator(30.0))
  ax.yaxis.set_minor_locator(tick.MultipleLocator(10.0))

def plot_var_xz(v, fig_title, vmin, vmax, vint, cnt_levels, cbar_int, lon_lim, zlim, 
                png_name):
  x = v.coords["lon"]
  z = v.coords["z"]

  X, Z = np.meshgrid(x,z)
  fig = plt.figure(figsize=(14,7)) 
  ax = fig.add_subplot(1,1,1)
  set_fig_Xz_axis(ax)
  ax.set_title(fig_title, fontsize=22)
  ax.set_xlim(lon_lim)  
  ax.set_ylim(zlim)
  
  print(f"output: {png_name} ..")
#  pcm = ax.pcolormesh(X, Z, v.values, vmin=vmin, vmax=vmax, cmap='jet')
  lv = np.linspace(vmin,vmax,int((vmax-vmin)/vint)+1)  
  pcm = ax.contourf(X, Z, v.values, levels=lv, cmap='jet', extend='both')  
  levels = cnt_levels[cnt_levels != 0.0]
  cont = ax.contour(X, Z, v.values, levels=levels, colors=['black'])  
#  cont.clabel(fmt='%2.0f', fontsize=12)
  cont = ax.contour(X, Z, v.values, levels=[0.0], colors=['lightgreen'], alpha=1.0) 

#  fmt = tick.ScalarFormatter(useMathText=True)
  cbar = plt.colorbar(pcm, aspect=50.0, ticks=cnt_levels, extend='both', shrink=1.0, orientation='vertical', pad=0.03)
  cbar.ax.ticklabel_format(style='sci', scilimits=(-3,3)) 
  cbar.ax.tick_params(labelsize=18)
  cbar.set_ticks(np.arange(vmin, vmax+1e-10, cbar_int))

  plt.savefig(png_name)
  
def plot_var_hori(v, fig_title, vmin, vmax, cnt_levels, lon_lim, lat_lim, 
                png_name):
  print(v)
  x = v.coords["lon"]
  y = v.coords["lat"]

  X, Y = np.meshgrid(x,y)
  fig = plt.figure(figsize=(14,7)) 
  ax = fig.add_subplot(1,1,1)
  set_fig_XY_axis(ax)
  ax.set_title(fig_title, fontsize=22)
  ax.set_xlim(lon_lim)  
  ax.set_ylim(lat_lim)
  
  print(f"output: {png_name} ..")
#  pcm = ax.pcolormesh(X, Z, v.values, vmin=vmin, vmax=vmax, cmap='jet')
  pcm = ax.contourf(X, Y, v.values, levels=cnt_levels, cmap='jet', extend='both')  
  levels = cnt_levels[cnt_levels != 0.0]
  cont = ax.contour(X, Y, v.values, levels=levels, colors=['black'])  
#  cont.clabel(fmt='%2.0f', fontsize=12)
  cont = ax.contour(X, Y, v.values, levels=[0.0], colors=['lightgreen'], alpha=1.0) 

#  fmt = tick.ScalarFormatter(useMathText=True)
  cbar = plt.colorbar(pcm, aspect=50.0, ticks=cnt_levels, extend='both', shrink=1.0, orientation='vertical', pad=0.03)
  cbar.ax.ticklabel_format(style='sci', scilimits=(-3,3)) 
  cbar.ax.tick_params(labelsize=18)

  plt.savefig(png_name)
  

def plot_var_xz_flat(v, fig_title, vmin, vmax, vint, cnt_levels, cbar_int, png_name):
  x = v.coords["x"]
  z = v.coords["z"]

  X, Z = np.meshgrid(x,z)
  fig = plt.figure(figsize=(14,7)) 
  ax = fig.add_subplot(1,1,1)
  set_fig_Xz_flat_axis(ax)
  ax.set_title(fig_title, fontsize=22)

  print(f"output: {png_name} ..")
  lv = np.linspace(vmin,vmax,int((vmax-vmin)/vint)+1)  
  print(lv)
  pcm = ax.contourf(X, Z, v.values, levels=lv, cmap='jet', extend="both")
  levels = cnt_levels[cnt_levels != 0.0]
  cont = ax.contour(X, Z, v.values, levels=levels, colors=['black'])  
#  cont.clabel(fmt='%2.0f', fontsize=12)
  cont = ax.contour(X, Z, v.values, levels=[0.0], colors=['lightgreen'], alpha=1.0) 

#  fmt = tick.ScalarFormatter(useMathText=True)
  cbar = plt.colorbar(pcm, aspect=50.0, extend='both', shrink=1.0, orientation='vertical', pad=0.03)
  cbar.ax.ticklabel_format(style='sci', scilimits=(-3,3)) 
  cbar.ax.tick_params(labelsize=18)
  cbar.set_ticks(np.arange(vmin, vmax+1e-10, cbar_int))
  
  plt.savefig(png_name)

def mkgraph(exp_name, time_list):
    os.makedirs(f"{ANALYSIS_OUT_DIR}/{exp_name}/", exist_ok=True)

    dir = f"{EQ_TINT_type}/{exp_name}/{OUTPUT_DATA_DIR}"
    
    print(f"{dir}/history.pe0*nc")
    w = xr.open_mfdataset(f"{dir}/history.pe0*nc", decode_times=False, combine='by_coords')["W"]

    #---
    w_lat0 = w.sel(lat=0).isel(lat=0)
    for time in time_list:
      print(f"plot: W time={time} sec ..")
      cnt_levels = np.arange(-1e-1, 1e-1+1e-2, 1e-2)  
      plot_var_xz(w_lat0.sel(time=time), "W", -.1, .1, .0025, cnt_levels, 2.5e-2, 
                  [140,220], [0, 10e3], 
                  f"{ANALYSIS_OUT_DIR}/{exp_name}/W_lat0_time{time}.pdf")    
            
    #---
    os.makedirs(TMP_DATA_DIR, exist_ok=True)
    w_hori = w.sel(z=8e3).isel(z=0).sel(time=time_list)    
    w_hori.to_netcdf(f"{TMP_DATA_DIR}/W_hori_z8km.nc")
    w_hori = xr.open_mfdataset(f"{TMP_DATA_DIR}/W_hori_z8km.nc", decode_times=False, combine='by_coords')["W"]
    for time in time_list:      
      print(f"plot: W time={time} sec ..")      
      cnt_levels = np.arange(-1e-1, 1e-1+1e-2, 1e-2)        
      plot_var_hori(w_hori.sel(time=time), "W", -1, 1, cnt_levels, [140,280], [-70,70], 
                  f"{ANALYSIS_OUT_DIR}/{exp_name}/W_z8km_time{time}.pdf")    

def mkgraph_linsol(u_lin, w_lin):
    os.makedirs(f"{ANALYSIS_OUT_DIR}/", exist_ok=True)
    w_ = w_lin.sel(x=slice(-26e3,26e3),z=slice(0,10e3))
    print(f"plot: W linear sol")
    cnt_levels = np.arange(-1e-1, 1e-1+1e-2, 1e-2)  
    plot_var_xz_flat(w_, "W", -.1, .1, .0025, cnt_levels, 2.5e-2, 
                     f"{ANALYSIS_OUT_DIR}/W_linsol.pdf")    

#---------------------------------

U0=20.0; TEMP0=300.0;
Lx=200e3; Lx_lin=Lx*8; Lz=30e3
Nx=4096*4; Nz=300
topo_params = {"name": "Schaer", "h0": h0, "a": 5e3, "lam":4e3, "xc": Lx/2}
u_lin, w_lin = linear_analytic_sol.gen_linsol(topo_params, Lx_lin, Nx, Lz, Nz, 
                                              U0, TEMP0)
mkgraph_linsol(u_lin, w_lin)

time_list=[600, 1200, 1800, 3600, 5400, 7200]
for exp_name in EXP_list:
    mkgraph(exp_name, time_list)