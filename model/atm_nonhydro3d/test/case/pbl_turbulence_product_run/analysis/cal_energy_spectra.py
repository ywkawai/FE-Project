import numpy as np
from numpy.fft import fftn
import xarray as xr
import glob
import os
from joblib import Parallel, delayed

ZLEVEL_list=[500]
OUTNC_suffix=""

#exp_name='E80P11'
#exp_name='E120P7'
#exp_name='E120P7_RK3'
#exp_name='E120P7_MF8OrdAlphx10'
#exp_name='E120P7_MF16Ord'
#exp_name='E120P7_MF16OrdAlphx10'
#exp_name='E120P7_MF32OrdAlphx10'
#exp_name='E120P7_hypupw'
#exp_name='E160P5'
#exp_name='E192P4'
#exp_name='E192P4_RK3'
# exp_name='E192P4_MF8OrdAlphx10'
#exp_name='E192P4_MF16OrdAlphx10'
#exp_name='E192P4_MF32OrdAlphx10'
#exp_name='E192P4_hypupw'
#exp_name='E240P3_RK3'
#exp_name='E240P3'
#exp_name='E240P3_MF8OrdAlphx10'
#exp_name='E240P3_MF16Ord'
exp_name='E240P3_MF16OrdAlphx10'
#exp_name='E240P3_MF32OrdAlphx10'
#exp_name='E240P3_hypupw'
#exp_name='E480P1'
#exp_name='E480P1_MF16Ord'
#exp_name='E480P1_MF32Alphx10'
#exp_name='E480P1_MFoff'
#exp_name='E480P1_hypupw'

dir=f"../run8_{exp_name}/outdata/"
dir_ind = [1]
TIME_list={ 
  1: np.arange(3.5*3600, 14440, 60),
  2: np.arange(13560, 4.0*3600+120, 120) }

ds_cache = {}

#-------------------------

tmp_dir=f"./tmp_data/tmp_{exp_name}"
out_dir=f"./{exp_name}"
Nproc = 6
PeNum = 144
if exp_name == 'E80P11':
  PeNum = 288
CpDry=1004.64
L = 9.6e3
Nx = 960
Ny = 960

def create_tmp_zslice_nc_sub(dir, tmp_dir, pe, zlevel_list, time_list):
  print(f"{dir}history.pe{pe:06}.nc pe={pe} zslice..")
  ds = xr.open_mfdataset(f'{dir}history.pe{pe:06}.nc', decode_times=False, combine='by_coords').sel(time=time_list, method="nearest")
  ds_basic = xr.open_mfdataset(f'{dir}basic_state.pe{pe:06}.nc', decode_times=False, combine='by_coords')

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

def calc_KE_spectra(dens, u, v, w):

  dx = L/float(Nx)
  dy = dx
  freqx = np.fft.fftfreq(Nx, dx)
  freqy = np.fft.fftfreq(Ny, dy)
  kx = 2.0*np.pi*freqx
  ky = 2.0*np.pi*freqy
  print(f"horivel_power_spectra Nx={Nx}, Ny={Ny}, dx=dy={dx}")
  #print(f"kx={kx[0:-1]}")
  dkx = kx[2] - kx[1]

  dens_sqrt = dens**0.5
  N = Nx
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

def get_varnp(tmpdir, varname, time, zlev):

  ncpath = f"{tmpdir}/tmp_{varname}_z{zlev}_interp.pe*.nc"
  if  not (ncpath in ds_cache):
    print(f"Add {ncpath} to ds cache..")
    ds = xr.open_mfdataset( ncpath, decode_times=False, combine='by_coords')
    ds_cache[ncpath] = ds
  
  ds = ds_cache[ncpath]
  return ds[varname].sel(time=time, method="nearest").values

def create_zslice_nc(dir, tmpdir, time_list):  
  r = Parallel(n_jobs=Nproc)( [delayed(create_tmp_zslice_nc_sub)(dir, tmpdir, pe, ZLEVEL_list, time_list) for pe in range(0,PeNum) ] )


#----

# make directory
os.makedirs(tmp_dir, exist_ok=True)
os.makedirs(out_dir, exist_ok=True)

#^^
time_list_ex = []
for di in dir_ind:
  create_zslice_nc(dir, tmp_dir, TIME_list[di])
  time_list_ex.extend(TIME_list[di])

print(time_list_ex)
Ncut = int(Nx/2)
freq = np.fft.fftfreq(Nx, L/float(Nx)) #Nnp.arange(-1.0/(2.0*dx), 1.0/(2.0*dx), 1.0/float(L))
kx = 2.0*np.pi*freq

ke_spectra_hvel = xr.DataArray( np.zeros((len(time_list_ex), len(ZLEVEL_list), Ncut)), 
                                coords=[time_list_ex, np.array(ZLEVEL_list, dtype=float), kx[:Ncut]], dims=["time", "z", "k"])
ke_spectra_hvel.name = "KE_spectra_momh"
ke_spectra_hvel.attrs["units"] = "kg.m-3.m3.s-2" 
ke_spectra_wvel = xr.DataArray( np.zeros((len(time_list_ex), len(ZLEVEL_list), Ncut)), 
                                coords=[time_list_ex, np.array(ZLEVEL_list, dtype=float), kx[:Ncut]], dims=["time", "z", "k"])
ke_spectra_wvel.name = "KE_spectra_momz"
ke_spectra_wvel.attrs["units"] = "kg.m-3.m3.s-2" 

tind_off = 0
for di in dir_ind:
  time_list = TIME_list[di]
  tmpdir = tmp_dir #+ f"_{di}/"
  print(tmpdir)
  for tind, time in enumerate(time_list):
    for zind, zlev in enumerate(ZLEVEL_list):
      print(f"calc_KE_spectra.. t={time} z={zlev}")
      kx_, ke_spectra_hvel[tind_off+tind,zind,:], ke_spectra_wvel[tind_off+tind,zind,:] = \
        calc_KE_spectra( 
          get_varnp(tmpdir, "DENS", time, zlev), 
          get_varnp(tmpdir, "U", time, zlev), 
          get_varnp(tmpdir, "V", time, zlev), 
          get_varnp(tmpdir, "W", time, zlev) )

  tind_off = len(time_list)

for var in [ke_spectra_hvel, ke_spectra_wvel]: var.to_netcdf(f"{tmp_dir}/{var.name}{OUTNC_suffix}.nc")

