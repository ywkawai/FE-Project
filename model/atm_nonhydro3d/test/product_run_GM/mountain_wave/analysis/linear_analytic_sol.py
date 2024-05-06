import xarray as xr
import numpy as np

Cp=1004.64
Rd=287.04
Cv=Cp-Rd
Grav=9.80665

def gen_linsol(topo_params, Lx, Nx, Lz, Nz, U0, TEMP0):
    BruntFreq=(Grav**2 / (Cp*TEMP0))**0.5
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
        fac0 = Grav/(Cp*TEMP0*S)
        fac = (np.exp(-S*z_[k]) * (1.0 - fac0 * (1.0 - np.exp(-S*z_[k])) )**(Cv/Rd) )**(-0.5)
        fac_z = - 0.5 * fac**3 * ( - S * np.exp(-S*z_[k]) * (1.0 - fac0 * (1.0 - np.exp(-S*z_[k])) )**(Cv/Rd) 
                                  - S * np.exp(-2.0*S*z_[k]) * fac0 * Cv/Rd * ( 1.0 - fac0 * (1.0 - np.exp(-S*z_[k])) )**(Cv/Rd-1.0) )
        w_hs_X_exp_imz = w_hs[:] * np.exp(v[:] * z_[k])
        
        print(f"{k}: {fac}, {fac_z}, {fac * Grav / 350.0**2}")
        w_w_z = 1.0j * U0 * k_[:] * w_hs_X_exp_imz[:]
        w_u_z = - U0 * ( fac * v[:] + fac_z - fac * Grav / 350.0**2 ) * w_hs_X_exp_imz[:]
        w_psi_z = U0**2 * v[:] * w_hs_X_exp_imz[:]
        
        phi_np[k,:] = fac*np.real(np.fft.ifft(w_psi_z))
        u_np[k,:] = np.real(np.fft.ifft(w_u_z))  
        w_np[k,:] = fac*np.real(np.fft.ifft(w_w_z))

        mwt_mask = np.where(hs > z_[k])
        phi_np[k,mwt_mask] = np.nan        
        u_np[k,mwt_mask] = np.nan        
        w_np[k,mwt_mask] = np.nan        

    xc = topo_params["xc"]
    #lon_ = 180.0 * ( ( x_ - xc ) / ( np.pi * RPlanet ) + 1.0 )
    x = x_ - xc
    
    phi_lin = xr.DataArray(phi_np, coords={'x': x, 'z': z_}, dims=['z', 'x'])
    u_lin = xr.DataArray(u_np, coords={'x': x, 'z': z_}, dims=['z', 'x'])
    w_lin = xr.DataArray(w_np, coords={'x': x, 'z': z_}, dims=['z', 'x'])
    return u_lin, w_lin

def gen_topo_data(x, params):
    h0 = params["h0"]
    if params["name"] =='bell_shape':
        a = params['a']
        xc = params['xc']        
        hs = h0 / ( ((x-72e3)/a)**2 + 1.0 )
    if params["name"] =='Schaer':
        a = params['a']
        lam = params['lam']        
        xc = params['xc']        
        hs = h0 * np.exp( - ((x-xc)/a)**2 ) * np.cos(np.pi * (x-xc) / lam)**2
    
    return hs
