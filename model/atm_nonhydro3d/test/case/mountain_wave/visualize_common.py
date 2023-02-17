import numpy as np
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
#from joblib import Parallel, delayed


Nx=4096*4
Nz=300
Cp=1004.0
Rd=287.0
Cv=Cp-Rd
Grav=9.81

#---------------------

def open_tmp_nc(dir, varname):
  print(f'{dir}/history.pe*.nc')
  return xr.open_mfdataset(f'{dir}/history.pe*.nc', decode_times=False, combine='by_coords')[varname]

def hori_1Daxis_fmt(tick_val, pos):
  val = int(tick_val/1000)
  return f'{val}'

def set_fig_YZ_axis(ax):
  ax.set_xlabel('x (km)', fontsize=12)  
  ax.xaxis.set_major_formatter(tick.FuncFormatter(hori_1Daxis_fmt))
  ax.xaxis.set_major_locator(tick.MultipleLocator(10e3))
  ax.xaxis.set_minor_locator(tick.MultipleLocator(2.5e3))
  ax.set_ylabel('z (km)', fontsize=12)
  ax.yaxis.set_major_formatter(tick.FuncFormatter(hori_1Daxis_fmt))
  ax.yaxis.set_major_locator(tick.MultipleLocator(1e3))
  ax.yaxis.set_minor_locator(tick.MultipleLocator(250.0))
  ax.tick_params(axis='both', labelsize=12)  

def plot_var_xz(title, png_name, v, vmin, vmax, levels, xlim, zlim, v_lin):

  x = v.coords["x"]
  z = v.coords["z"]
  X, Z = np.meshgrid(x,z)
  
  x2 = v_lin.coords["x"]
  z2 = v_lin.coords["z"]
  X2, Z2 = np.meshgrid(x2,z2)
  
  fig = plt.figure(figsize=(8,9)) 
  ax = fig.add_subplot(1,1,1)
  set_fig_YZ_axis(ax)
  ax.set_xlim(xlim)
  ax.set_ylim(zlim)  
  ax.set_title(title, fontsize=15)

  pcm = ax.pcolormesh(X, Z, v, vmin=vmin, vmax=vmax, cmap='jet', alpha=0.1)
  ax.contour(X, Z, v, levels, cmap='jet', linewidths=2, linestyles="-")
  ax.contour(X2, Z2, v_lin, levels, colors='gray', linewidths=1, linestyles="--")
  
  fmt = tick.ScalarFormatter(useMathText=True)
  fmt.set_powerlimits((0,0))
  cbar = plt.colorbar(pcm, aspect=40.0, extend='both', orientation='horizontal', shrink=0.5, format=fmt)
  cbar.ax.xaxis.get_offset_text().set_fontsize(7)

  plt.savefig(png_name)


def gen_linsol(topo_params, U0, PTEMP0, BruntFreq, Lx, Nx, Lz, Nz):
    
  S=BruntFreq**2/Grav
  
  x_ = np.linspace(0e3,Lx,Nx)
  z_ = np.linspace(0e3,Lz,Nz)
  k_ = np.fft.fftfreq(Nx) * 2.0*np.pi * Nx / Lx
  m2_ = (BruntFreq / U0)**2 - k_**2

  hs = gen_topo_data(x_, topo_params)
  w_hs = np.fft.fft(hs)
  v = np.zeros((len(w_hs)), dtype=complex)
  phi_np = np.zeros((Nz,len(w_hs)))
  u_np = np.zeros((Nz,len(w_hs)))
  w_np = np.zeros((Nz,len(w_hs)))

  trap_ind = np.where(m2_ < 0.0)
  prop_ind1 = np.where((m2_ >=0.0) & (k_ >=0.0))
  prop_ind2 = np.where((m2_ >=0.0) & (k_ < 0.0))
  v[trap_ind] = - np.sqrt(np.abs(m2_[trap_ind]))
  v[prop_ind1] = 1.0j * np.sqrt(np.abs(m2_[prop_ind1]))
  v[prop_ind2] = - 1.0j * np.sqrt(np.abs(m2_[prop_ind2]))

  for k in range(0,len(z_)):
      fac = (np.exp(-S*z_[k]) * (1.0 - Grav/(Cp*PTEMP0*S)*(1.0 - np.exp(-S*z_[k])) )**(Cv/Rd) )**(-0.5)
      w_hs_X_exp_imz = w_hs[:] * np.exp(v[:] * z_[k])
      
      w_w_z = 1.0j * U0 * k_[:] * w_hs_X_exp_imz[:]
      w_u_z = - U0 * v[:] * w_hs_X_exp_imz[:]  
      w_psi_z = U0**2 * v[:] * w_hs_X_exp_imz[:]
      
      phi_np[k,:] = fac*np.real(np.fft.ifft(w_psi_z))
      u_np[k,:] = fac*np.real(np.fft.ifft(w_u_z))  
      w_np[k,:] = fac*np.real(np.fft.ifft(w_w_z))
      
      mwt_mask = np.where(hs > z_[k])
      phi_np[k,mwt_mask] = np.nan        
      u_np[k,mwt_mask] = np.nan        
      w_np[k,mwt_mask] = np.nan        

  phi_lin = xr.DataArray(phi_np, coords={'x': x_, 'z': z_}, dims=['z', 'x'])
  u_lin = xr.DataArray(u_np, coords={'x': x_, 'z': z_}, dims=['z', 'x'])
  w_lin = xr.DataArray(w_np, coords={'x': x_, 'z': z_}, dims=['z', 'x'])
  return u_lin, w_lin

def gen_topo_data(x, params):
    h0 = params["h0"]
    if params["name"] =='bell_shape':
        a = params['a']
        xc = params['xc']        
        hs = h0 / ( ((x-xc)/a)**2 + 1.0 )
    if params["name"] =='Schaer':
        a = params['a']
        lam = params['lam']        
        xc = params['xc']        
        hs = h0 * np.exp( - ((x-xc)/a)**2 ) * np.cos(np.pi * (x-xc) / lam)**2
    
    return hs

#-----------------------------------