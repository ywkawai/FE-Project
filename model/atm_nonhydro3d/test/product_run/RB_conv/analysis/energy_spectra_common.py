import numpy as np
from numpy.fft import fftn
import xarray as xr
import glob
import os
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


def create_tmp_zslice_nc_sub(dir, dir_ini, tmp_dir, pe, zlevel_list, time_list):
  print(f"{dir}history.pe{pe:06}.nc pe={pe} zslice..")
  ds = xr.open_mfdataset(f'{dir}/history.pe{pe:06}.nc', decode_times=False, combine='by_coords').sel(time=time_list, method="nearest")
  ds_basic = xr.open_mfdataset(f'{dir_ini}/basic_state.pe{pe:06}.nc', decode_times=False, combine='by_coords')

  for zlev in zlevel_list:
    for var_name in ["U", "V", "W"]:
      print(f"Output {var_name} z={zlev}m: pe={pe}..")
      ds[var_name].interp(z=zlev).to_netcdf(f"{tmp_dir}/tmp_{var_name}_z{zlev}_interp.pe{pe:06}.nc")

    print(f"Output DENS z={zlev}m : pe={pe}..")
    dens_hyd = ds_basic["DENS_hyd"].interp(z=zlev)
    ddens = ds["DDENS"].interp(z=zlev)
    var = dens_hyd + ddens
    var.name = "DENS"
    var.to_netcdf(f"{tmp_dir}/tmp_{var.name}_z{zlev}_interp.pe{pe:06}.nc")

def create_zslice_nc(dir, dir_ini, tmp_dir, penum, zlevel_list, time_list, Nproc=4):
  r = Parallel(n_jobs=Nproc)( [delayed(create_tmp_zslice_nc_sub)(dir, dir_ini, tmp_dir, pe, zlevel_list, time_list) for pe in range(0,penum) ] )

def calc_KE_spectra(dens, u, v, w, nx, ny, L):

  dx = L/float(nx)
  dy = dx
  freqx = np.fft.fftfreq(nx, dx)
  freqy = np.fft.fftfreq(nx, dy)
  kx = 2.0*np.pi*freqx
  ky = 2.0*np.pi*freqy
  print(f"horivel_power_spectra Nx={nx}, Ny={ny}, dx=dy={dx}")
  #print(f"kx={kx[0:-1]}")
  dkx = kx[2] - kx[1]

  dens_sqrt = dens**0.5
  N = nx
  s_ru = np.fft.fft2(dens_sqrt*u)/np.sqrt(N**2)
  s_rv = np.fft.fft2(dens_sqrt*v)/np.sqrt(N**2)
  s_rw = np.fft.fft2(dens_sqrt*w)/np.sqrt(N**2)
  s_f_abs_hvel = 0.5*(np.abs(s_ru)**2 + np.abs(s_rv)**2)/N**2
  s_f_abs_wvel = 0.5*np.abs(s_rw)**2/N**2
  s_f_abs_3dvel = s_f_abs_hvel + s_f_abs_wvel


  #s_f_abs = np.fft.fftshift(s_f_abs)
#  iend = int(N/2)+1
#  plt.plot(n[:iend],s_f_abs[:iend])
#  print(s_f_abs)

  Ncut = int(N/2)
  kxx, kyy = np.meshgrid(kx, ky, indexing="ij")
  nrr = np.sqrt(kxx**2 + kyy**2)
  s_f_abs_1d_hvel = np.zeros(Ncut)
  s_f_abs_1d_wvel = np.zeros(Ncut)

  print(kxx.shape)
  print(kyy.shape)
  print(s_f_abs_hvel.shape)

  for i in range(0,Ncut):
    a = np.where(( (nrr >= kx[i]-0.5*dkx ) & (nrr < kx[i]+0.5*dkx) ))
    s_f_abs_1d_hvel[i] = s_f_abs_hvel[a[1],a[0]].sum(axis=0)
    s_f_abs_1d_wvel[i] = s_f_abs_wvel[a[1],a[0]].sum(axis=0)
  
  s_f_abs_1d_3dvel = s_f_abs_1d_hvel + s_f_abs_1d_wvel
  #print(s_f_abs_1d_3dvel)
  
  check_total_KE = np.sum(0.5*dens*(u**2 + v**2 + w**2))/float(N**2)
  print(f"check total KE={check_total_KE}, {np.sum(s_f_abs_3dvel)}, {np.sum(s_f_abs_1d_3dvel)}")
  return kx[0:Ncut], s_f_abs_1d_hvel, s_f_abs_1d_wvel

def get_varnp(tmpdir, varname, time, zlev, ds_cache):

  ncpath = f"{tmpdir}/tmp_{varname}_z{zlev}_interp.pe*.nc"
  if  not (ncpath in ds_cache):
    print(f"Add {ncpath} to ds cache..")
    ds = xr.open_mfdataset( ncpath, decode_times=False, combine='by_coords')
    ds_cache[ncpath] = ds
  
  ds = ds_cache[ncpath]
  return ds[varname].sel(time=time, method="nearest").values

def gen_tmp_nc(
    exp_top_dir, exp_name, dir_ind, 
    zlevel_list, time_list, penum, 
    tmp_dir, runno_inidata, 
    nproc, skip_flag):
  time_list_ex = []
  for di in dir_ind:
    target_dir = f"{exp_top_dir}/{exp_name}/run{di}/outdata/"
    target_dir_ini = f"{exp_top_dir}/{exp_name}/run{runno_inidata}/outdata/"

    tmpdir = tmp_dir + f"_{di}/"  
    if skip_flag==False:
      os.makedirs(tmpdir, exist_ok=True)  
      create_zslice_nc(target_dir, target_dir_ini, tmpdir, penum, zlevel_list, time_list[di], nproc)

    time_list_ex.extend(time_list[di])

  return time_list_ex

def energy_spectra_analysis(exp_top_dir, exp_name, nx, penum, L, 
                            dir_ind, runno_inidata, 
                            TIME_list, ZLEVEL_list, 
                            OUTNC_suffix, 
                            nproc, gem_tmp_data_skip_flag):
  #----
  tmp_dir=f"{exp_top_dir}/tmp_data_energy_spectra/tmp_{exp_name}"  
  os.makedirs(tmp_dir, exist_ok=True)  

  time_list_ex = gen_tmp_nc( exp_top_dir, exp_name, dir_ind, ZLEVEL_list, TIME_list, penum, tmp_dir, runno_inidata, nproc, gem_tmp_data_skip_flag)
  r = Parallel(n_jobs=nproc)( [delayed(energy_spectra_analysis_core_sub)( tmp_dir, nx, L, di, TIME_list[di], ZLEVEL_list, OUTNC_suffix) for di in dir_ind ] )

  
def energy_spectra_analysis_core_sub( tmp_dir, nx, L, 
                            di, time_list, ZLEVEL_list, 
                            OUTNC_suffix ):
  #----
  Ncut = int(nx/2)
  freq = np.fft.fftfreq(nx, L/float(nx)) #Nnp.arange(-1.0/(2.0*dx), 1.0/(2.0*dx), 1.0/float(L))
  kx = 2.0*np.pi*freq

  ke_spectra_hvel = xr.DataArray( np.zeros((len(time_list), len(ZLEVEL_list), Ncut)), 
                                  coords=[time_list, np.array(ZLEVEL_list, dtype=float), kx[:Ncut]], dims=["time", "z", "k"])
  ke_spectra_hvel.name = "KE_spectra_momh"
  ke_spectra_hvel.attrs["units"] = "kg.m-3.m3.s-2" 
  ke_spectra_wvel = xr.DataArray( np.zeros((len(time_list), len(ZLEVEL_list), Ncut)), 
                                  coords=[time_list, np.array(ZLEVEL_list, dtype=float), kx[:Ncut]], dims=["time", "z", "k"])
  ke_spectra_wvel.name = "KE_spectra_momz"
  ke_spectra_wvel.attrs["units"] = "kg.m-3.m3.s-2" 

  tmpdir = tmp_dir + f"_{di}/"
  print(tmpdir)
  
  ds_cache = {} 
  for zind, zlev in enumerate(ZLEVEL_list):  
    dens_ = get_varnp(tmpdir, "DENS", time_list, zlev, ds_cache)
    u_ = get_varnp(tmpdir, "U", time_list, zlev, ds_cache)
    v_ = get_varnp(tmpdir, "V", time_list, zlev, ds_cache)
    w_ = get_varnp(tmpdir, "W", time_list, zlev, ds_cache)
      
    for tind, time in enumerate(time_list):
      print(f"calc_KE_spectra.. t={time} z={zlev}")
      dens = np.copy(dens_[:,:,tind])
      u = u_[tind,:,:]
      v = v_[tind,:,:]
      w = w_[tind,:,:]
      
      kx_, ke_spectra_hvel[tind,zind,:], ke_spectra_wvel[tind,zind,:] = \
        calc_KE_spectra( dens, u, v, w, nx, nx, L )

  for var in [ke_spectra_hvel, ke_spectra_wvel]: var.to_netcdf(f"{tmpdir}/{var.name}{OUTNC_suffix}.nc")

  
#---- mkgraph
def xaxis_txt(inv_lam, pos=None):
  l = str(int(1.0/inv_lam))
  return "  "+l+'$^{-1}$'

def xaxis_txt_minor(inv_lam, pos=None):
  l = str(int(1.0/inv_lam))
  if l=="20" or l=="50" or l=="200" or l=="500" or l=="2000" or l=="5000":
    return "  "+l+'$^{-1}$'
  else:
    return ""

def create_fig_energy_spectra_zoom(ke_spectra_list, zlev, figname, 
                                   exp_ltype_list, exp_ltype_width, exp_color_list, exp_label_list, 
                                   ylim_range, 
                                   K_EXP_NAME_LABEL, slope_m35_ampl = 1.7e-5):
  print(f"create_fig_energy_spectra: {figname}")

  fig, ax = plt.subplots(figsize=(5,10))
#pcm = plt.pcolormesh(kx[1:Ncut], ky[1:Ncut], np.log10(s_f_abs[1:Ncut,1:Ncut]), vmin=-6, vmax=6)
#plt.colorbar(pcm)
  ax.set_xlim(1.0/300.0,1.0/10.0)
  ax.set_ylim(ylim_range[0], 1e-2)
  #ax.set_ylim(1e-3, 1.5e-2)
  ax.set_yscale('log')
  ax.set_xscale('log')

  for exp_name in ke_spectra_list.keys(): 
    ke_spectra = ke_spectra_list[exp_name].sel(z=zlev)
    ax.plot(ke_spectra.k/(2.0*np.pi), ke_spectra, 
            linestyle=exp_ltype_list[exp_name], 
            linewidth=exp_ltype_width[exp_name], 
            color=exp_color_list[exp_name], 
            label=exp_label_list[exp_name] )

  kx = ke_spectra_list[K_EXP_NAME_LABEL].k
  slope_m35 = slope_m35_ampl * kx**(-5.0/3.0) 
  ax.plot(kx/(2.0*np.pi), slope_m35, color="grey", linestyle="-.", label="-5/3", linewidth=2)
  ax.set_xlabel("inverse of wavelength [m$^{-1}$]", fontsize=15)
  ax.set_ylabel("E(k)", fontsize=15)
  ax.xaxis.set_major_formatter(FuncFormatter(xaxis_txt))
  ax.xaxis.set_minor_formatter(FuncFormatter(xaxis_txt_minor))
  ax.tick_params(which="both", labelsize=15)

  ax.legend(loc='lower left')
  plt.savefig(figname)

def create_fig_energy_spectra(
    ke_spectra_list, zlev, figname, 
    exp_ltype_list, exp_ltype_width, exp_color_list, exp_label_list, 
    ylim_range, 
    K_EXP_NAME_LABEL, slope_m35_ampl = 1.7e-5):
  print(f"create_fig_energy_spectra: {figname}")

  fig, ax = plt.subplots(figsize=(8,6))
#pcm = plt.pcolormesh(kx[1:Ncut], ky[1:Ncut], np.log10(s_f_abs[1:Ncut,1:Ncut]), vmin=-6, vmax=6)
#plt.colorbar(pcm)
  ax.set_xlim(1.0/3000.0,1.0/10.0)
  ax.set_ylim(ylim_range[0], ylim_range[1])
  ax.set_yscale('log')
  ax.set_xscale('log')

  for exp_name in ke_spectra_list.keys(): 
    ke_spectra = ke_spectra_list[exp_name].sel(z=zlev)
    ax.plot(ke_spectra.k/(2.0*np.pi), ke_spectra, 
            linestyle=exp_ltype_list[exp_name], 
            linewidth=exp_ltype_width[exp_name], 
            color=exp_color_list[exp_name], 
            label=exp_label_list[exp_name])

#  kx = ke_spectra_list["reference"].k
  kx = ke_spectra_list[K_EXP_NAME_LABEL].k
  slope_m35 = slope_m35_ampl * kx**(-5.0/3.0) 
  ax.plot(kx/(2.0*np.pi), slope_m35, color="grey", linestyle="-.", label="-5/3", linewidth=2)
  ax.set_xlabel("inverse of wavelength [m$^{-1}$]", fontsize=15)
#  ax.set_ylabel("E(k)", fontsize=14)
  ax.xaxis.set_major_formatter(FuncFormatter(xaxis_txt))
  ax.xaxis.set_minor_formatter(FuncFormatter(xaxis_txt_minor))
  ax.tick_params(which="both", labelsize=13)

  ax.legend(loc='lower left')
  plt.savefig(figname)


def create_fig_energy_spectra_diffm53(ke_spectra_list, zlev, figname, 
                                      exp_ltype_list, exp_ltype_width, exp_color_list, exp_label_list, 
                                      ylim_range,                                       
                                      slope_m35_ampl = 1.7e-5):
  print(f"create_fig_energy_spectra: {figname}")

  fig, ax = plt.subplots(figsize=(10,7))
  ax.set_xlim(1.0/3000.0,1.0/18.0)
  ax.set_ylim(ylim_range[0], ylim_range[1])
  ax.set_yscale('log')
  ax.set_xscale('log')

  for exp_name in ke_spectra_list.keys(): 
    ke_spectra = ke_spectra_list[exp_name].sel(z=zlev)
    slope_m35 = slope_m35_ampl*ke_spectra.k**(-5.0/3.0) 
    ax.plot(ke_spectra.k/(2.0*np.pi), ke_spectra/slope_m35, 
      linestyle=exp_ltype_list[exp_name], 
      linewidth=exp_ltype_width[exp_name], 
      color=exp_color_list[exp_name], 
      label=exp_label_list[exp_name] )

  kx = ke_spectra_list[exp_name].k
  ax.plot(kx/(2.0*np.pi), slope_m35/slope_m35, color="grey", linestyle="-.", label="-5/3")

  ax.set_xlabel("inverse of wavelength [m$^{-1}$]", fontsize=20)
  ax.xaxis.set_major_formatter(FuncFormatter(xaxis_txt))
  ax.xaxis.set_minor_formatter(FuncFormatter(xaxis_txt_minor))
  ax.tick_params(which="both", labelsize=16, length=3)

  plt.savefig(figname)
  
def read_spectra_data(exp_name, tmp_dir, dir_ind, ke_spectra_3dvel_list, 
                      ncut, OUTNC_suffix=""):
  di = dir_ind[0]
  fac = 1.0 / float(len(dir_ind))
  ke_spectra_hvel = xr.open_mfdataset(f'{tmp_dir}_{di}/KE_spectra_momh{OUTNC_suffix}.nc', decode_times=False, combine='by_coords').KE_spectra_momh
  ke_spectra_wvel= xr.open_mfdataset(f'{tmp_dir}_{di}/KE_spectra_momz{OUTNC_suffix}.nc', decode_times=False, combine='by_coords').KE_spectra_momz
  ke_spectra_3dvel_tmp = fac * np.mean(ke_spectra_hvel[1:,:,1:ncut] + ke_spectra_wvel[1:,:,1:ncut], axis=0)
  for i in range(1,len(dir_ind)):
    di = dir_ind[i]
    ke_spectra_hvel = xr.open_mfdataset(f'{tmp_dir}_{di}/KE_spectra_momh{OUTNC_suffix}.nc', decode_times=False, combine='by_coords').KE_spectra_momh
    ke_spectra_wvel = xr.open_mfdataset(f'{tmp_dir}_{di}/KE_spectra_momz{OUTNC_suffix}.nc', decode_times=False, combine='by_coords').KE_spectra_momz
    ke_spectra_3dvel_tmp = ke_spectra_3dvel_tmp + fac * np.mean(ke_spectra_hvel[1:,:,1:ncut] + ke_spectra_wvel[1:,:,1:ncut], axis=0)
    
  ke_spectra_3dvel_list[exp_name] = ke_spectra_3dvel_tmp
  return

def read_spectra_data_hv(exp_name, tmp_dir, dir_ind, ke_spectra_hvel_list, ke_spectra_wvel_list, 
                      ncut, ncut_s=1,OUTNC_suffix=""):
  di = dir_ind[0]
  fac = 1.0 / float(len(dir_ind))
  ke_spectra_hvel = xr.open_mfdataset(f'{tmp_dir}_{di}/KE_spectra_momh{OUTNC_suffix}.nc', decode_times=False, combine='by_coords').KE_spectra_momh
  ke_spectra_wvel= xr.open_mfdataset(f'{tmp_dir}_{di}/KE_spectra_momz{OUTNC_suffix}.nc', decode_times=False, combine='by_coords').KE_spectra_momz
  ke_spectra_hvel_tmp = fac * np.mean(ke_spectra_hvel[1:,:,ncut_s:ncut], axis=0)
  ke_spectra_wvel_tmp = fac * np.mean(ke_spectra_wvel[1:,:,ncut_s:ncut], axis=0)  
  for i in range(1,len(dir_ind)):
    di = dir_ind[i]
    ke_spectra_hvel = xr.open_mfdataset(f'{tmp_dir}_{di}/KE_spectra_momh{OUTNC_suffix}.nc', decode_times=False, combine='by_coords').KE_spectra_momh
    ke_spectra_wvel = xr.open_mfdataset(f'{tmp_dir}_{di}/KE_spectra_momz{OUTNC_suffix}.nc', decode_times=False, combine='by_coords').KE_spectra_momz
    ke_spectra_hvel_tmp = ke_spectra_hvel_tmp + fac * np.mean(ke_spectra_hvel[1:,:,ncut_s:ncut], axis=0)
    ke_spectra_wvel_tmp = ke_spectra_wvel_tmp + fac * np.mean(ke_spectra_wvel[1:,:,ncut_s:ncut], axis=0)
    
  ke_spectra_hvel_list[exp_name] = ke_spectra_hvel_tmp
  ke_spectra_wvel_list[exp_name] = ke_spectra_wvel_tmp  
  return
