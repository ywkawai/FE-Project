import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from joblib import Parallel, delayed
import glob
import os
from matplotlib.lines import Line2D  


def create_tmp_nc(dir, tmp_dir, penum, zlev_list, time_list, n_bins, bin_min, bin_max, Nproc=4):  
  Parallel(n_jobs=Nproc)( [delayed(create_tmp_nc_sub1)(dir, pe, tmp_dir, zlev_list, time_list, n_bins, bin_min, bin_max) for pe in range(0,penum) ] )

def create_tmp_nc_sub1(dir, pe, tmp_dir, zlev_list, time_list, n_bins, bin_min, bin_max):
#   print(f'{dir}history.pe{pe:06}.nc')
  print(f"{dir}: pe={pe} vel_abs ..")

  for z in zlev_list:
    var_list = {}
    for varname in ["U", "V", "W"]:
        ds = xr.open_mfdataset(f"{tmp_dir}/tmp_{varname}_z{z}_interp.pe{pe:06}.nc", decode_times=False, combine='by_coords')
        var_list[varname] = ds[varname]
    #   print(f"Output {varname} z={z} pe={pe}..")
    #   var.to_netcdf(f"{tmp_dir}/tmp_{varname}_z{z}_interp.pe{pe:06}.nc")

    vel_abs = ( (var_list["U"]**2 + var_list["V"]**2 + var_list["W"]**2)**0.5 ).rename('vel_abs')
    vel_abs.to_netcdf(f"{tmp_dir}/tmp_vel_abs_z{z}_interp.pe{pe:06}.nc")
    hist, bin_edges = np.histogram(vel_abs.values.flatten(), bins=n_bins, range=(bin_min, bin_max))
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    # print(f"{dir}, {pe}: {hist}")
    hist_da = xr.DataArray(
        data=np.array(hist, dtype='int64'),
        coords={"bin_center": bin_centers},
        dims=["bin_center"],
        name="vel_abs_frequency")
    # print(hist_da)    
    hist_da.to_netcdf(f"{tmp_dir}/tmp_vel_abs_frequency_z{z}_interp.pe{pe:06}.nc")
    
def create_pdf_uvw_single(u, v, w, 
                          nbin, vel_range, 
                          figname):
    fig, ax = plt.subplots()
    xr.plot.hist(u, density=True, alpha=0.5, color="green", bins=nbin, range=vel_range, histtype='stepfilled', label="U")
    xr.plot.hist(v, density=True, alpha=0.5, color="blue", bins=nbin, range=vel_range, histtype='stepfilled', label="V")
    xr.plot.hist(w, density=True, alpha=0.5, color="red", bins=nbin, range=vel_range, histtype='stepfilled', label="W")
    plt.legend(fontsize=14)

    plt.ylim(0, 0.65)
    plt.xlabel("wind speed [m/s]", fontsize=14)
    plt.ylabel("probability density", fontsize=14)
    plt.tick_params(labelsize=14)
    plt.savefig(figname)
    
def create_pdf_uvw_overplot(varname, v_list, nbin, v_range, ylim_range, 
                            exp_color_list, exp_linetype_list, figname):
    exp_name_list = v_list.keys()

    fig, ax = plt.subplots()
    for exp_name in exp_name_list:
        xr.plot.hist(v_list[exp_name], density=True, color=exp_color_list[exp_name], bins=nbin, range=v_range, histtype='step', label=exp_name, linestyle=exp_linetype_list[exp_name])

    handles, labels = ax.get_legend_handles_labels()
    print(handles[0])
    new_handles = [Line2D([], [], c=h.get_edgecolor(), ls=h.get_linestyle()) for h in handles]        
    plt.legend(handles=new_handles, labels=labels, fontsize=10)

    plt.ylim(ylim_range[0], ylim_range[1])
    plt.xlabel(f"{varname} [m/s]", fontsize=14)
    plt.ylabel("probability density", fontsize=14)
    plt.tick_params(labelsize=14)
    plt.savefig(figname)

def create_pdf_uvw_overplot_log( 
                                varname, v_list, nbin, v_range, ylim_range, exp_color_list, exp_linetype_list, figname):
    exp_name_list = v_list.keys()

    fig, ax = plt.subplots(figsize=(8,10))
    for exp_name in exp_name_list:
        xr.plot.hist(v_list[exp_name], density=True, color=exp_color_list[exp_name], bins=nbin, range=v_range, histtype='step', label=exp_name, linestyle=exp_linetype_list[exp_name],
                     linewidth=3)

    # handles, labels = ax.get_legend_handles_labels()
    # print(handles[0])
    # new_handles = [Line2D([], [], c=h.get_edgecolor(), ls=h.get_linestyle()) for h in handles]        
    # plt.legend(handles=new_handles, labels=labels, fontsize=10)

    plt.yscale("log")
    plt.ylim(ylim_range[0], ylim_range[1])
    plt.xlabel(f"{varname} [m/s]", fontsize=14)
    plt.ylabel("probability density", fontsize=14)
    plt.tick_params(labelsize=14)
    plt.savefig(figname)
    

def get_var(tmp_data_dir, suffix, varname, time_list):
    return xr.open_mfdataset(f'{tmp_data_dir}/tmp_{varname}{suffix}.pe*.nc', decode_times=False, combine='by_coords')[varname].sel(time=time_list)

def merge_tmp_nc_sub(varname, tmpdir, exp_name, dir_ind, pe, suffix=""):    
    print(f"merge.. {varname} pe={pe}")    
    tmp_data_dir = f"{tmpdir}/tmp_{exp_name}_{dir_ind[0]}/"
    var_tmp = xr.open_mfdataset(f'{tmp_data_dir}/tmp_{varname}{suffix}.pe{pe:06}.nc')[varname]
    for i in range(1,len(dir_ind)):
        di = dir_ind[i]
        tmp_data_dir = f"{tmpdir}/tmp_{exp_name}_{di}/"
        var = xr.open_mfdataset(f'{tmp_data_dir}/tmp_{varname}{suffix}.pe{pe:06}.nc')[varname]
        var_tmp = var + var_tmp
    var_tmp.rename(varname).to_netcdf(f"{tmpdir}/tmp_{exp_name}/tmp_{varname}{suffix}.pe{pe:06}.nc")

def merge_tmp_nc(varname, tmpdir, exp_name, dir_ind, penum, suffix="", Nproc=4):  
    Parallel(n_jobs=Nproc)( [delayed(merge_tmp_nc_sub)(varname, tmpdir, exp_name, dir_ind, pe, suffix) for pe in range(0,penum) ] )

    tmp_data_dir = f"{tmpdir}/tmp_{exp_name}"
    var_tmp = xr.open_mfdataset(f'{tmp_data_dir}/tmp_{varname}{suffix}.pe{0:06}.nc')[varname]
    for pe in range(1,penum):
        var = xr.open_mfdataset(f'{tmp_data_dir}/tmp_{varname}{suffix}.pe{pe:06}.nc')[varname]
        var_tmp = var + var_tmp    
    var_tmp.rename(varname).to_netcdf(f"{tmp_data_dir}/tmp_{varname}{suffix}.nc")

def create_tmp_data( top_exp_dir, exp_name_list, tmpdir, dir_ind_list, runno_inidata,
                    zlev_listlist, time_listlist, pe_num_list, 
                    n_bins, bin_min, bin_max, nproc=4 ):
    for exp_name in exp_name_list:
        tmp_data_dir = f"{tmpdir}/tmp_{exp_name}"    
        os.makedirs(tmp_data_dir, exist_ok=True)   
        
        time_list = time_listlist[exp_name]
        zlev_list = zlev_listlist[exp_name]
        
        for di in dir_ind_list[exp_name]:
            target_dir = f"{top_exp_dir}/{exp_name}/run{di}/outdata/"
            target_dir_ini = f"{top_exp_dir}/{exp_name}/run{runno_inidata}/outdata/"

            tmpdir_di = tmp_data_dir + f"_{di}/"
            os.makedirs(tmpdir_di, exist_ok=True)        
            create_tmp_nc(target_dir, tmpdir_di, pe_num_list[exp_name], zlev_list, time_list[di], n_bins, bin_min, bin_max, nproc )

def merge_tmp_data( exp_name_list, tmpdir, dir_ind_list,
                    zlev_listlist, pe_num_list, 
                    nproc=4 ):
    for exp_name in exp_name_list:
        tmp_data_dir = f"{tmpdir}/tmp_{exp_name}"    
        os.makedirs(tmp_data_dir, exist_ok=True)   
        zlev_list = zlev_listlist[exp_name]
        for varname in ["vel_abs_frequency"]:
            for zlev in zlev_list:
                merge_tmp_nc(varname, tmpdir, exp_name, 
                            dir_ind_list[exp_name], pe_num_list[exp_name], 
                            f"_z{zlev}_interp", nproc)
    