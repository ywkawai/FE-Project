import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import os

#--
run_dir = "rhot_hevi/Eh16Ez12P7"
init_dir="rhot_hevi/Eh16Ez12P7"
#--
# run_dir = "rhot_hevi/Eh32Ez24P7_dtx0.5"
# init_dir="rhot_hevi/Eh32Ez24P7"
#--
# run_dir = "rhot_hevi/Eh32Ez24P7_dtx0.25"
# init_dir="rhot_hevi/Eh32Ez24P7"
#--
# run_dir = "rhot_hevi/Eh64Ez48P7_dtx0.5"
# init_dir="rhot_hevi/Eh64Ez48P7"
#--
# run_dir = "rhot_hevi/Eh64Ez48P7_dtx0.25"
# init_dir="rhot_hevi/Eh64Ez48P7"
#--
# run_dir = "rhot_hevi/Eh6Ez4P11"
# init_dir="rhot_hevi/Eh6Ez4P11"
#--
# run_dir = "rhot_hevi/Eh12Ez8P11"
# init_dir="rhot_hevi/Eh12Ez8P11"

dist_dir = f"analysis_out/{run_dir}"

#------------------
def open_nc(dir, varname):
    return  xr.open_mfdataset(f"{dir}/history.pe000*.nc", decode_times=False, combine='by_coords')[varname]

def open_nc_bs(dir, varname):
    return  xr.open_mfdataset(f"{dir}/bs.pe000*.nc", decode_times=False, combine='by_coords')[varname]

def v_1Daxis_fmt(tick_val, pos):
  val = int(tick_val/1000)
  return f'{val}'

def set_fig_Xz_axis(ax):
  ax.tick_params(labelsize=16, length=8)
  ax.set_xlabel('longitude [deg]', fontsize=16)  
  ax.xaxis.set_major_locator(tick.MultipleLocator(30.0))
  ax.xaxis.set_minor_locator(tick.MultipleLocator(10.0))
  ax.set_ylabel('height [km]', fontsize=16)
  ax.yaxis.set_major_locator(tick.MultipleLocator(2000.0))
  ax.yaxis.set_minor_locator(tick.MultipleLocator(500.0))
  ax.yaxis.set_major_formatter(tick.FuncFormatter(v_1Daxis_fmt))

def plot_var_xz(v, fig_title, vmin, vmax, vint, cnt_levels, png_name):
  x = v.coords["lon"]
  z = v.coords["z"]

  X, Z = np.meshgrid(x,z)
  fig = plt.figure(figsize=(9,6)) 
  ax = fig.add_subplot(1,1,1)
  set_fig_Xz_axis(ax)
#  ax.set_title(fig_title, fontsize=24)

  print(f"output: {png_name} ..")
  lv = np.linspace(vmin,vmax,int((vmax-vmin)/vint)+1)
  print(lv)
  pcm = ax.contourf(X, Z, v.values, levels=lv, cmap='jet', extend="both")
  levels = cnt_levels[cnt_levels**2 > 5e-17**2]
#  print(levels)
  cont = ax.contour(X, Z, v.values, levels=levels, colors=['black'])  
#  cont.clabel(fmt='%2.0f', fontsize=12)
  cont = ax.contour(X, Z, v.values, levels=[0.0], colors=['lightgreen'], linewidths=2.0)  

#  fmt = tick.ScalarFormatter(useMathText=True)
  cbar = plt.colorbar(pcm, aspect=60.0, 
                      extend='both', shrink=0.9, orientation='horizontal', pad=0.14)
#  cbar.set_ticks(np.linspace(np.floor(vmin), np.ceil(vmax), num=5, endpoint=True))
  cbar.ax.ticklabel_format(style='sci', scilimits=(-3,3), useMathText=True)
  cbar.ax.xaxis.get_offset_text().set(size=14) 
  cbar.ax.tick_params(labelsize=14)

  plt.savefig(png_name)

  #----------------------------------------------------

os.makedirs(dist_dir, exist_ok=True)    

vars = {}
for varname in ["Umet", "W", "THERM", "DDENS"]:
    vars[varname] = open_nc(f"{run_dir}/outdata", varname)

for varname in ["PRES_hyd", "DENS_hyd"]:
    vars[varname] = open_nc_bs(f"{init_dir}/outdata", varname)

Rd=287.04
Cp=1004.64
Gam=Cp/(Cp-Rd)
PRE0 = 1000e2

rhot_hyd = PRE0 / Rd * (vars["PRES_hyd"] / PRE0)**(1.0/Gam) 
vars["PT_dash"] = rhot_hyd / vars["DENS_hyd"] * ( ( 1.0 + vars["THERM"] / rhot_hyd ) / ( 1.0 + vars["DDENS"] / vars["DENS_hyd"] ) - 1.0 )

for tsec in [43200]:

    levels = np.arange(-5e-3, 5e-3+5e-4, 5e-4)    
    plot_var_xz(vars["PT_dash"].sel(lat=0, time=tsec)[:,0,:], 'Potential temperature perturbation [K]', 
                -3.5e-3, 3.5e-3, 1e-4, levels,  f"{dist_dir}/PTdash_t{tsec}sec.pdf")
    
    levels = np.arange(-8e-3, 8e-3+1e-3, 1e-3)    
    plot_var_xz(vars["Umet"].sel(lat=0, time=tsec)[:,0,:], 'Zonal wind [m/s]', 
                -8e-3, 8e-3, 2e-4, levels,  f"{dist_dir}/Umet_t{tsec}sec.pdf")
    
    levels = np.arange(-3e-5, 3e-5+2e-6, 2e-6)    
    plot_var_xz(vars["W"].sel(lat=0, time=tsec)[:,0,:], 'Vertical wind [m/s]',
                -2.5e-5, 2.5e-5, 5e-7, levels,  f"{dist_dir}/W_t{tsec}sec.pdf")

for tsec in [86400]:
    
    levels = np.arange(-2.5e-3, 2.5e-3+2.5e-4, 2.5e-4)    
    plot_var_xz(vars["PT_dash"].sel(lat=0, time=tsec)[:,0,:], 'Potential temperature perturbation [K]', 
                -3e-3, 3e-3, 1e-4, levels,  f"{dist_dir}/PTdash_t{tsec}sec.pdf")
    
    levels = np.arange(-6e-3, 6e-3+1e-3, 1e-3)    
    plot_var_xz(vars["Umet"].sel(lat=0, time=tsec)[:,0,:], 'Zonal wind [m/s]', 
                -4e-3, 4e-3, 2e-4, levels,  f"{dist_dir}/Umet_t{tsec}sec.pdf")
    
    levels = np.arange(-2e-5, 2e-5+1e-6, 1e-6)    
    plot_var_xz(vars["W"].sel(lat=0, time=tsec)[:,0,:], 'Vertical wind [m/s]',
                -2e-5, 2e-5, 2e-7, levels,  f"{dist_dir}/W_t{tsec}sec.pdf")

for tsec in [129600, 172800]:    
    levels = np.arange(-2.5e-3, 2.5e-3+2.5e-4, 2.5e-4)    
    plot_var_xz(vars["PT_dash"].sel(lat=0, time=tsec)[:,0,:], 'Potential temperature perturbation [K]', 
                -3e-3, 3e-3, 1e-4, levels,  f"{dist_dir}/PTdash_t{tsec}sec.pdf")
  
    levels = np.arange(-6e-3, 6e-3+1e-3, 1e-3)    
    plot_var_xz(vars["Umet"].sel(lat=0, time=tsec)[:,0,:], 'Zonal wind [m/s]',
                -4e-3, 4e-3, 2e-4, levels,  f"{dist_dir}/Umet_t{tsec}sec.pdf")
    
    levels = np.arange(-1.8e-5, 1.8e-5+1e-6, 1e-6)        
    plot_var_xz(vars["W"].sel(lat=0, time=tsec)[:,0,:], 'Vertical wind [m/s]', 
                -1.8e-5, 1.8e-5, 5e-7, levels,  f"{dist_dir}/W_t{tsec}sec.pdf")
