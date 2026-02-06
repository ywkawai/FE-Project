import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname('__file__'), '../../analysis_lib/'))
import vstructure_common as common

#--
exp_dir = "../"
analysis_out_dir="./fig_vstructure/"
out_figext = "svg"


exp_list = {
   "Dx25m_P3": {"run_s": 3, "start_time":2*7200, "run_sm": 4,  "run_e": 24, "integ_time_per_run": 7200},     
   "Dx25m_P7": {"run_s": 1, "start_time":0, "run_sm": 4,  "run_e": 24, "integ_time_per_run": 7200},            
   "Dx27m_P11": {"run_s": 3, "start_time":2*7200, "run_sm": 4,  "run_e": 24, "integ_time_per_run": 7200},   
#--            
   "Dx12.5m_P3": {"run_s": 3, "start_time":2*7200, "run_sm": 8,  "run_e": 28, "integ_time_per_run": 7200},       
   "Dx12.5m_P7": {"run_s": 1, "start_time":0, "run_sm": 4,  "run_e": 24, "integ_time_per_run": 7200},    
   "Dx13m_P11": {"run_s": 3, "start_time":2*7200, "run_sm": 8,  "run_e": 28, "integ_time_per_run": 7200},                             
#--               
   "Dx6.3m_P3": {"run_s": 13, "start_time":19*3600, "run_sm": 13,  "run_e": 57, "integ_time_per_run": 3600},            
   "Dx6.3m_P7": {"run_s": 8, "start_time":7200*7, "run_sm": 8,  "run_e": 52, "integ_time_per_run": 3600},                  
   "Dx6.7m_P11": {"run_s": 13, "start_time":3600*19, "run_sm": 13,  "run_e": 34, "integ_time_per_run": 7200}, 
#--              
   "Dx3.1m_P7": {"run_s": 13, "start_time": 19*3600, "run_sm": 13,  "run_e": 100, "integ_time_per_run": 1800}, 
}

Exp_lsize_list = {
    "Dx25m_P3": 1,  "Dx25m_P3_MFweak": 1, 
    "Dx25m_P7": 2, "Dx25m_P7_MF": 1, 
    "Dx27m_P11": 1, "Dx27m_P11_MF": 1, 
    "Dx12.5m_P3": 1, "Dx12.5m_P3_MFweak": 2,
    "Dx12.5m_P7": 2, "Dx12.5m_P7_MFweak": 1,  
    "Dx13m_P11": 1, "Dx13m_P11_MFweak": 1,  
    "Dx6.3m_P3": 1, "Dx6.3m_P7": 2, "Dx6.7m_P11": 2, 
    "Dx3.1m_P7": 3
}
Exp_color_list = {
  "Dx25m_P3": "black", "Dx25m_P3_MF": "cyan",
  "Dx25m_P7": "black", "Dx25m_P7_MF": "gray", 
  "Dx27m_P11":  "black", "Dx27m_P11_MF":  "green",   
  "Dx12.5m_P3": "red", 'Dx12.5m_P3_MFweak': "orange", 
  "Dx12.5m_P7": "red", 'Dx12.5m_P7_MF': "yellow",    
  'Dx13m_P11': "red",
  'Dx6.3m_P3': "green", 'Dx6.3m_P7': "green",  'Dx6.7m_P11': "green", 
  'Dx3.1m_P7': "blue", 
}
Exp_ltype_list = {
  "Dx25m_P3": ":", "Dx25m_P3_MF": "--", 
  "Dx25m_P7": "--",  "Dx25m_P7_MF": "-",  
  "Dx27m_P11": "-",  "Dx27m_P11_MF":  "-",
  "Dx12.5m_P3": ":", 'Dx12.5m_P3_MFweak': ":", 
  "Dx12.5m_P7": "--", "Dx12.5m_P7_MF": "-", 
  'Dx13m_P11': "-", 
  'Dx6.3m_P3': ":", 'Dx6.3m_P7': "--", 'Dx6.7m_P11': "-", 
  'Dx3.1m_P7': "--"
}
Exp_label_list = {
  "Dx25m_P3":  "Δ=25m,P3",  
  "Dx25m_P7":  "Δ=25m,P7",   
  "Dx27m_P11":  "Δ=27m,P11",  
  "Dx12.5m_P3":  "Δ=13m,P3",
  "Dx12.5m_P7":  "Δ=13m,P7",     
  'Dx13m_P11': 'Δ=13m,P11', 
  'Dx6.3m_P3': 'Δ=6.3m,P3',
  'Dx6.3m_P7': 'Δ=6.3m,P7',  
  'Dx6.7m_P11': 'Δ=6.7m,P11', 
  'Dx3.1m_P7': 'Δ=3.1m,P7', 
}

#-------------------------------
# Functions
def set_tick(ax, xtick_major, xtick_minor, ytick_major, ytick_minor):
  ax.set_xticks(xtick_major)
  ax.set_xticks(xtick_minor, minor=True)
  ax.xaxis.set_ticks_position('both')
  ax.set_yticks(ytick_major)
  ax.set_yticks(ytick_minor, minor=True)

def set_ax(ax, vmin, vmax, tick_major, tick_minor, fig_title):
  ax.legend(fontsize=18)
  ax.set_xlim(vmin, vmax)
  ax.set_ylim(-10, 1610)
  ax.set_ylabel("height [m]", fontsize=24) 
  ax.tick_params(which="both", labelsize=18)    
  set_tick(ax, tick_major, tick_minor, 
               np.linspace(0,1500,4), np.linspace(0,1500,16) )
  ax.set_title(fig_title, fontsize=18)
  ax.legend()
  
def draw_lines(var_list, ax, vmin, vmax, tick_major, tick_minor, fig_title, exp_color_list, exp_ltype_style, exp_ltype_width, exp_label_list, do_uniform_vgrid_sampling):
  var_list_tmp0 = {}
  if do_uniform_vgrid_sampling:
    for exp_name, var in var_list.items():
      print(exp_name)
      zcoord_name = var.dims[0]
      z = var.coords[zcoord_name]
      dz = 1.6e3/len(z.values)
      z_ = np.arange(0.0+0.5*dz,1600+0.5*dz,dz)
      var_list_tmp0[exp_name] = var.groupby(zcoord_name).mean(zcoord_name).interp({zcoord_name: z_}) 
  
  var_list_tmp = {}  
  dz_ref = 3.125
  z_ref = np.arange(0.0,1600.0+dz_ref,dz_ref)
  for exp_name, var in var_list.items():
    print(exp_name)
    zcoord_name = var.dims[0]
    if do_uniform_vgrid_sampling:
      var_list_tmp[exp_name] = var_list_tmp0[exp_name].interp({zcoord_name: z_ref}, kwargs={"fill_value": "extrapolate"}).values
    else:
      var_list_tmp[exp_name] = var.groupby(zcoord_name).mean(zcoord_name).interp({zcoord_name: z_ref}) 

  var_s = np.vstack(list(var_list_tmp.values()))
  v_min = np.minimum.reduce(var_s)
  v_max = np.maximum.reduce(var_s)
  ax.fill_betweenx( z_ref, v_min, v_max, color="lightgray", alpha=0.6 )

  for exp_name, var in var_list.items():
    if "P7" in exp_name:
      print(exp_name)
      zcoord_name = var.dims[0]
      z = var.coords[zcoord_name]
      # ax.plot( var, z, 
      #       linestyle=exp_ltype_style[exp_name], color=exp_color_list[exp_name], label=exp_label_list[exp_name], linewidth=exp_ltype_width[exp_name], alpha=0.5 )
      ax.plot( var_list_tmp[exp_name], z_ref, 
            linestyle=exp_ltype_style[exp_name], color=exp_color_list[exp_name], label=exp_label_list[exp_name], linewidth=exp_ltype_width[exp_name], alpha=0.5 )

  set_ax(ax, vmin, vmax, tick_major, tick_minor, fig_title)
  
def create_fig_vertical_structure_fillmaxmin( var_list, vmin, vmax, tick_major, tick_minor, fig_title, exp_color_list, exp_ltype_style, exp_ltype_width, exp_label_list, out_pngname, do_uniform_vgrid_sampling=True ):
  fig, ax = plt.subplots(figsize=(10.0, 8.0))
  draw_lines(var_list, ax, vmin, vmax, tick_major, tick_minor, fig_title, exp_color_list, exp_ltype_style, exp_ltype_width, exp_label_list, do_uniform_vgrid_sampling)      
  plt.savefig(out_pngname)
  
def create_fig_vertical_structure_fillmaxmin_2( var_listlist, vmin, vmax, tick_major, tick_minor, fig_title, exp_color_list, exp_ltype_style, exp_ltype_width, exp_label_list, out_pngname, do_uniform_vgrid_sampling=True ):
  fig, ax = plt.subplots(figsize=(10.0, 8.0))
  for var_list in var_listlist:
    draw_lines(var_list, ax, vmin, vmax, tick_major, tick_minor, fig_title, exp_color_list, exp_ltype_style, exp_ltype_width, exp_label_list, do_uniform_vgrid_sampling)      
  plt.savefig(out_pngname)    


#--

os.makedirs(analysis_out_dir, exist_ok=True)

pt, eddy_momz_flux, sgs_momz_flux, mean_heat_flux, eddy_heat_flux, sgs_heat_flux, sgs_heat_bnd_stabflux = common.get_pbl_analysis(exp_dir, exp_list)
hwind_vari, vwind_vari = common.get_pbl_analysis_vel_vari(exp_dir, exp_list)

#- plot PT

create_fig_vertical_structure_fillmaxmin( pt, 299, 302, np.arange(299.0,302.0,0.5), np.arange(299.0,302.0,0.1), 
                                         "PT", Exp_color_list, Exp_ltype_list, Exp_lsize_list, Exp_label_list, 
                                         f"{analysis_out_dir}/PT.{out_figext}" )

#- plot heat flux
eddy_plus_sgs_heat_flux = {}
for key, ehf in eddy_heat_flux.items():
    eddy_plus_sgs_heat_flux[key] = ehf + sgs_heat_flux[key] + sgs_heat_bnd_stabflux[key]

create_fig_vertical_structure_fillmaxmin( eddy_plus_sgs_heat_flux, 0, 200.0, np.arange(0,220,20), np.arange(0,225,5), 
                                         "RESOLVED EDDY + SGS + SGS_BND_STAB", Exp_color_list, Exp_ltype_list, Exp_lsize_list, Exp_label_list, 
                                         f"{analysis_out_dir}/eddy_plus_SGS_plus_BNDstab_heat_flux.{out_figext}", False )

hwind_vari_tmp = {}
vwind_vari_tmp = {}
for key, hwind_vari_ in hwind_vari.items():
    hwind_vari_tmp[key] = hwind_vari_**0.5
    vwind_vari_tmp[key] = vwind_vari[key]**0.5

create_fig_vertical_structure_fillmaxmin( hwind_vari_tmp, 0, 6, np.arange(0,7,1), np.arange(0,7,0.1), 
                                         "HWind_vari", Exp_color_list, Exp_ltype_list, Exp_lsize_list, Exp_label_list, 
                                         f"{analysis_out_dir}/HWind_vari.{out_figext}" )
create_fig_vertical_structure_fillmaxmin( vwind_vari_tmp, 0, 6, np.arange(0,7,1), np.arange(0,7,0.1), 
                                         "VWind_vari", Exp_color_list, Exp_ltype_list, Exp_lsize_list, Exp_label_list, 
                                         f"{analysis_out_dir}/VWind_vari.{out_figext}" )
create_fig_vertical_structure_fillmaxmin_2( [hwind_vari_tmp, vwind_vari_tmp], 0, 5, np.arange(0,5,1), np.arange(0,5,0.1), 
                                         "VWind_vari", Exp_color_list, Exp_ltype_list, Exp_lsize_list, Exp_label_list, 
                                         f"{analysis_out_dir}/VWind_vari.{out_figext}" )