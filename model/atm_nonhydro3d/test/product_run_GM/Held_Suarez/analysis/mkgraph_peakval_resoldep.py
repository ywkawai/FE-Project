import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.ticker as tick
import numpy as np
import os


EXP_list_p3 = ["Eh12Ez8P3", "Eh24Ez16P3", "Eh48Ez32P3"]
EXP_list_p7 = [ "Eh6Ez4P7", "Eh12Ez8P7", "Eh24Ez16P7", "Eh48Ez32P7"]
EXP_list_p11 = [ "Eh4Ez3P11", "Eh8Ez6P11", "Eh16Ez12P11" ]

EXP_listlist = {
    "P3": EXP_list_p3, "P7": EXP_list_p7, "P11": EXP_list_p11
}
EXP_list_DOF = {
    "Eh12Ez8P3": 48, "Eh24Ez16P3": 96, "Eh48Ez32P3": 192,
    "Eh6Ez4P7": 48, "Eh12Ez8P7": 96, "Eh24Ez16P7": 192, "Eh48Ez32P7": 384, 
    "Eh4Ez3P11": 48, "Eh8Ez6P11": 96, "Eh16Ez12P11": 192, 
}
EXP_list_color = {
    "P3": "blue", "P7": "red", "P11": "green",  
}
ANALYSIS_OUT_DIR="./analysis_out/peak_resoldep"

#--------------------------------------
PREVSTUDY_LCOLOR_LIST = {
  "TS2004": "cyan", 
  "Wan2013": "purple", 
  "Wan2008": "red",   
  "UJ2012": "orange",     
}
PREVSTUDY_LTYPE_LIST = {
  "TS2004": "--"
}

# PREVSTUDY_PEAKVAL_Umet = {
#   "TS2004": [[240e3, 120e3, 60e3], 
#              [32, 30, 30], [35, 33, 32.5]], 
#   "Wan2013": [[277e3, 138e3, 69e3], 
# }
PREVSTUDY_PEAKVAL_EDDYTEMP = {
   "Wan2013": [[277e3, 138e3, 69e3], 
               [21,38,46],[24,39,49]], 
   "Wan2008": [[180e3], 
               [44],[45]], 
   "UJ2012": [[208e3], 
               [40],[44]],    
   
}
PREVSTUDY_PEAKVAL_EDDYKinE = {
   "Wan2013": [[277e3, 138e3, 69e3], 
               [185,340,440],[195,380,450]], 
   "Wan2008": [[180e3], 
               [420],[430]], 
   "UJ2012": [[208e3], 
               [340],[370]],     
}
PREVSTUDY_PEAKVAL_EDDYMOMFlux = {
  "TS2004": [[240e3, 120e3, 60e3], 
             [65, 72.5, 75], [72.5, 77.5, 82.5]], 
  "Wan2013": [[277e3, 138e3, 69e3], 
               [52,75,80],[60,80,85]], 
  "Wan2008": [[180e3], 
              [72],[75]], 
  "UJ2012": [[208e3], 
              [62],[67]],        
     
}
PREVSTUDY_PEAKVAL_EDDYHEATFlux = {
  "TS2004": [[240e3, 120e3, 60e3],
             [19,22,22], [20, 23, 23] ], 
  "Wan2013": [[277e3, 138e3, 69e3], 
               [18,22,23],[20,23,25]], 
  "Wan2008": [[180e3], 
              [22],[23]], 
  "UJ2012": [[208e3], 
              [18],[19]],        
}

#---------------------------------------

def v_1Daxis_fmt(tick_val, pos):
  val = int(tick_val/100)
  return f'{val}'

def set_fig_Yp_axis(ax):
  ax.tick_params(labelsize=12, length=8)
  ax.set_xlabel('latitude [degrees]', fontsize=18) 
  ax.set_xlim([-89,89]) 
  ax.xaxis.set_major_locator(tick.MultipleLocator(30.0))
  ax.xaxis.set_minor_locator(tick.MultipleLocator(10.0))
  ax.set_ylabel('pressure [hPa]', fontsize=18)
  ax.invert_yaxis()
  ax.yaxis.set_major_locator(tick.MultipleLocator(100e2))
  ax.yaxis.set_major_formatter(tick.FuncFormatter(v_1Daxis_fmt))
  
def get_max_value(var):
  var_max = var.where(var==var.max(), drop=True).squeeze()
  var_max_ = np.array(var_max.values).flatten()
  max_loc_lat = np.array(var_max.coords["lat"].values).flatten()
  max_loc_plev = np.array(var_max.coords["p"].values).flatten()
  return var_max_[0], max_loc_lat[0], max_loc_plev[0]

def get_max(var_list, EXP_list, EXP_DOF_info):
  sh_var_data = []; sh_max_loc_lat_data = []; sh_max_loc_plev_data = []; 
  dof_list = []
  for exp_name in EXP_list:
    var = np.abs(var_list[exp_name].mean(["time"]))
    var_max, max_loc_lat, max_loc_plev = get_max_value(var)
    print(f"{var.name} SH: {var_max} (lat={max_loc_lat}, plev={max_loc_plev})")
    sh_var_data.append(var_max)
    sh_max_loc_lat_data.append(max_loc_lat)
    sh_max_loc_plev_data.append(max_loc_plev)
    dof_list.append(EXP_DOF_info[exp_name])
      
  sh_var = xr.DataArray(sh_var_data, dims=["DOF"], coords={'DOF': dof_list})  
  sh_max_loc_lat = xr.DataArray(sh_max_loc_lat_data, dims=["DOF"], coords={'DOF': dof_list})  
  sh_max_loc_plev = xr.DataArray(sh_max_loc_plev_data, dims=["DOF"], coords={'DOF': dof_list})  
  return sh_var, sh_max_loc_lat, sh_max_loc_plev

def mkgraph(var_list, vmin, vmax, EXP_list):
  fig, axs = plt.subplots(len(var_list), figsize=(10,12))
  fig.suptitle(var_list[EXP_list[0]].name)
  i = 0
  for exp_name in EXP_list:
    set_fig_Yp_axis(axs[i])
    
    var = var_list[exp_name]
    lat = var.coords["lat"]  
    plev = var.coords["p"] 
    
    Lat, Plev = np.meshgrid(lat,plev)   
    pcm = axs[i].pcolormesh(Lat, Plev, var.mean(["time"]), vmin=vmin, vmax=vmax, cmap="jet")
    fmt = tick.ScalarFormatter(useMathText=True)
    
    cbar = plt.colorbar(pcm, aspect=30.0, extend='both', orientation='vertical', shrink=0.8, format=fmt, pad=0.02, ax=axs[i])
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.yaxis.get_offset_text().set_fontsize(14)  
    
    i = i +1
    
def get_sh_nh_var_list(var_list, EXP_list):
  sh_var_list = {}; nh_var_list = {}
  for exp_name in EXP_list:
    sh_var_list[exp_name] = var_list[exp_name].sel(lat=slice(-90,0))
    nh_var_list[exp_name] = var_list[exp_name].sel(lat=slice(0,90))
  return sh_var_list, nh_var_list


def func_minor_tick(x, pos):
  ticks = [20, 40, 60, 80, 200]
  if np.any(np.isclose(np.ones_like(ticks)*x, ticks)):
      return f"{x:g}"
  else:
      return ""
      
def set_axis(ax, vmin, vmax):
  ax.tick_params(labelsize=14, length=8)  
  ax.set_xlim(20,300)
  ax.set_ylim(vmin,vmax)
  ax.invert_xaxis()
  ax.set_xscale("log")
  ax.xaxis.set_minor_locator(tick.LogLocator(numticks=999, subs="all"))
  ax.xaxis.set_minor_formatter(tick.FuncFormatter(func_minor_tick))
  ax.set_xlabel("Horizontal grid spacing [km]", fontsize=13)
  
def mkgraph_sh_nh_max_conv( 
                           var_listlist, EXP_listlist, EXP_DOF_info, EXP_list_color, vmin, vmax, lat_min, lat_max, plev_min, plev_max, 
                           prev_study_data, 
                           pngname ):
  fig = plt.figure(figsize=(18,6))
    
  ax1 = fig.add_subplot(1,3,1)
  set_axis(ax1, vmin,vmax)  

  ax2 = fig.add_subplot(1,3,2)
  set_axis(ax2, lat_min,lat_max)    

  ax3 = fig.add_subplot(1,3,3)
  set_axis(ax3, plev_min,plev_max) 
  ax3.invert_yaxis()
   
  
  for key, data in prev_study_data.items():
    doff = np.array(data[0])/1e3
    prev_data_sh = np.array(data[1])
    prev_data_nh = np.array(data[2])
    print(f"{key}: {doff} {prev_data_sh} {prev_data_nh}")
    ax1.errorbar(doff, 0.5*(prev_data_sh+prev_data_nh), yerr=0.5*np.abs(prev_data_sh-prev_data_nh), elinewidth=6, capthick=0, capsize=3, markersize=5, alpha=0.2, 
                 color=PREVSTUDY_LCOLOR_LIST[key], linestyle="None")    
  
  i = 0
  for explist_key, EXP_list in EXP_listlist.items():
    sh_var_list, nh_var_list = get_sh_nh_var_list(var_listlist[i], EXP_list)
    sh_var, sh_max_loc_lat, sh_max_loc_plev = get_max(sh_var_list, EXP_list, EXP_DOF_info)
    nh_var, nh_max_loc_lat, nh_max_loc_plev = get_max(nh_var_list, EXP_list, EXP_DOF_info)
    dof = 1e4/sh_var.DOF
      
    # ax1.plot(dof, sh_var, label="SH", color="blue", linestyle="--", marker="o", alpha=0.25)
    # ax1.plot(dof, nh_var, label="NH", color="cyan", linestyle="--", marker="o", alpha=0.25)
    ax1.plot(dof, (0.5*(sh_var+nh_var)), label=explist_key, color=EXP_list_color[explist_key], linestyle="-", marker="o")
    ax1.errorbar(dof, 0.5*(sh_var+nh_var), yerr=0.5*np.abs(sh_var-nh_var), capsize=5,  markersize=10, alpha=0.5, color=EXP_list_color[explist_key])    
    ax1.legend()
        
    # ax2.plot(dof, -sh_max_loc_lat, label="SH", linestyle="--", marker="o", alpha=0.25, color=EXP_list_color[explist_key])
    # ax2.plot(dof, nh_max_loc_lat, label="NH", linestyle="-.", marker="o", alpha=0.25, color=EXP_list_color[explist_key])
    ax2.plot(dof, (0.5*(-sh_max_loc_lat+nh_max_loc_lat)), label=explist_key, linestyle="-", marker="o", color=EXP_list_color[explist_key])
    ax2.errorbar(dof, 0.5*(-sh_max_loc_lat+nh_max_loc_lat), yerr=0.5*np.abs(sh_max_loc_lat+nh_max_loc_lat), capsize=5,  markersize=10, alpha=0.5, color=EXP_list_color[explist_key])        
    ax2.legend()  
    
    # ax3.plot(dof, sh_max_loc_plev, label="SH", linestyle="--", marker="o", alpha=0.25, color=EXP_list_color[explist_key])
    # ax3.plot(dof, nh_max_loc_plev, label="NH",  linestyle="-.", marker="o", alpha=0.25, color=EXP_list_color[explist_key])
    ax3.plot(dof, (0.5*(sh_max_loc_plev+nh_max_loc_plev)), label=explist_key, linestyle="-", marker="o", color=EXP_list_color[explist_key])
    ax3.errorbar(dof, 0.5*(sh_max_loc_plev+nh_max_loc_plev), yerr=0.5*np.abs(sh_max_loc_plev-nh_max_loc_plev), capsize=5,  markersize=10, alpha=0.5, color=EXP_list_color[explist_key])            
    ax3.legend()    
    i = i+1
                    
  plt.savefig(pngname)


umet_listlist = []
merid_eddy_hflx_listlist = []
merid_eddy_momflx_listlist = []
eddy_kinetic_energy_listlist = []
eddy_temp_variance_listlist = []

for EXP_list_key, EXP_list in EXP_listlist.items():
  umet_list = {}; 
  merid_eddy_hflx_list = {}
  merid_eddy_momflx_list = {}
  eddy_kinetic_energy_list = {}
  eddy_temp_variance_list = {}
  for exp_name in EXP_list:
    dir = f"analysis_out/{exp_name}/tmp_data"
    print(dir)
    umet_list[exp_name] = xr.open_mfdataset(f'{dir}/Umet.nc', decode_times=False, combine='by_coords')["Umet"]
    merid_eddy_hflx_list[exp_name] = xr.open_mfdataset(f'{dir}/merid_eddy_hflx.nc', decode_times=False, combine='by_coords')["merid_eddy_hflx"]
    merid_eddy_momflx_list[exp_name] = xr.open_mfdataset(f'{dir}/merid_eddy_momflx.nc', decode_times=False, combine='by_coords')["merid_eddy_momflx"]  
    eddy_kinetic_energy_list[exp_name] = xr.open_mfdataset(f'{dir}/eddy_kinetic_energy.nc', decode_times=False, combine='by_coords')["eddy_kinetic_energy"]  
    eddy_temp_variance_list[exp_name] = xr.open_mfdataset(f'{dir}/eddy_temp_variance.nc', decode_times=False, combine='by_coords')["eddy_temp_variance"]

  umet_listlist.append(umet_list)
  merid_eddy_hflx_listlist.append(merid_eddy_hflx_list)
  merid_eddy_momflx_listlist.append(merid_eddy_momflx_list)
  eddy_kinetic_energy_listlist.append(eddy_kinetic_energy_list)
  eddy_temp_variance_listlist.append(eddy_temp_variance_list)


os.makedirs(ANALYSIS_OUT_DIR, exist_ok=True)
# mkgraph_sh_nh_max_conv(umet_listlist, EXP_listlist, EXP_list_DOF, EXP_list_color, 29.5, 33.2, 40.0,46.5, 200e2, 260e2, 
#                        PREVSTUDY_PEAKVAL_Umet, 
#                        f"{ANALYSIS_OUT_DIR}/Umet_resoldep.pdf" )

mkgraph_sh_nh_max_conv(eddy_temp_variance_listlist, EXP_listlist, EXP_list_DOF, EXP_list_color, 20, 50, 35, 45, 780e2, 870e2, 
                       PREVSTUDY_PEAKVAL_EDDYTEMP, 
                       f"{ANALYSIS_OUT_DIR}/Tvariance_resoldep.pdf" )

mkgraph_sh_nh_max_conv(eddy_kinetic_energy_listlist, EXP_listlist, EXP_list_DOF, EXP_list_color, 180, 560, 40.5,46, 200e2, 350e2, 
                       PREVSTUDY_PEAKVAL_EDDYKinE, 
                       f"{ANALYSIS_OUT_DIR}/EddyKinEnergy_resoldep.pdf" )

mkgraph_sh_nh_max_conv(merid_eddy_momflx_listlist, EXP_listlist, EXP_list_DOF, EXP_list_color, 52, 85, 28, 38, 220e2, 280e2, 
                       PREVSTUDY_PEAKVAL_EDDYMOMFlux, 
                       f"{ANALYSIS_OUT_DIR}/EddyMomFlux_resoldep.pdf" )

mkgraph_sh_nh_max_conv(merid_eddy_hflx_listlist, EXP_listlist, EXP_list_DOF, EXP_list_color, 18.0,26, 33.5,43.5, 830e2, 890e2, 
                       PREVSTUDY_PEAKVAL_EDDYHEATFlux, 
                       f"{ANALYSIS_OUT_DIR}/EddyHeatFlux_resoldep.pdf" )
