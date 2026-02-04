import numpy as np
import xarray as xr

import re
from collections import defaultdict
from scipy.optimize import root_scalar
from functools import partial


def eval_effective_resolution( ke_spectra_ratio, cr: float, dx: float): 
    rms =  (ke_spectra_ratio-cr)**2
    k = ke_spectra_ratio.k
    argmin_ind = np.argmin(rms.values)
    k_eff = k[argmin_ind].values

    eff_resol_grid_len = 2.0*np.pi / ( k_eff * dx )
    return k_eff, eff_resol_grid_len

def eval_energy_pile_error(Enum, Eref, k1: float, k2: float):
    kind_max = np.min([Enum.k.size, Eref.k.size])

    Eref_ = Eref.isel(k=slice(0,kind_max))
    Enum_ = Enum.isel(k=slice(0,kind_max))
    band0 = (Enum_.k >= k1) & (Enum_.k <= k2)    
    band = (Enum_.k >= k1) & (Enum_.k <= k2) & (Enum_.values > Eref_.values)

    diff = np.abs(Enum_ - Eref_)
    abs_error = diff.where(band, drop=True).integrate(coord='k')
    ref_energy = Eref_.where(band0, drop=True).integrate(coord='k')
    rel_error = abs_error / ref_energy
    return abs_error.values, rel_error.values

def diff(x, f_fit, y):
  return f_fit(x) - y

def eval_energy_pile_error_equiv_dx(energy_accum_list):
    grouped = defaultdict(dict)
    pattern = r'Dx([0-9.]+)m_P(\d+)'

    for key, value in energy_accum_list.items():
        match = re.match(pattern, key)
        if match:
            dx_str, p_str = match.groups()
            dx = float(dx_str)
            p = int(p_str)
            grouped[p][dx] = value
    dataarrays = {}
    for p, dx_values in grouped.items():
        dx_list = sorted(dx_values.keys())
        val_list = [dx_values[dx] for dx in dx_list]
        da = xr.DataArray(val_list, coords={"Dx": dx_list}, dims=["Dx"])
        dataarrays[p] = da    

    #-
    energy_accum_dx_list = {}
    for p, da_ in dataarrays.items():
    
        da = da_.where(da_.Dx > 6.0, drop=True)
        x = da.Dx
        coeffs = np.polyfit(x, da.values, deg=2)
        x_fit = np.linspace(x.min().values, x.max().values, 100)
        y_fit = np.polyval(coeffs, x_fit)
    
        da_P3 = dataarrays[3]
        f_fit = lambda x: np.polyval(coeffs, x)
        for dx in x.values:
            diff_partial = partial(diff, f_fit=f_fit, y=da_P3.sel(Dx=dx, method="nearest"))
            result = root_scalar(diff_partial, method="newton", x0=da.Dx[-1])
            print(f"Dx={dx}, p={p} Root: Dx={result.root:12.2f} (y={da_P3.sel(Dx=dx, method="nearest").values:12.4f})")
            energy_accum_dx_list[f"Dx{str(dx).strip('.0')}m_P{p}"] = result.root

    return energy_accum_dx_list

def output_energy_spectra_eval( eff_resol_list, energy_accum_list, num_energy_accum_equiv_dx_list, out_fname ):
    key_list = list(eff_resol_list.keys())

    line_list = []
    sect_str = "info"
    header_items = []
    header_items.append(f"{sect_str:<22}")
    for exp_name in key_list:
        header_items.append(f"{exp_name:>12}")
    line = ", ".join(header_items)
    print(line); line_list.append(line+"\n")

    eff_resol_str = []
    num_en_accum_str = []
    num_en_accum_dx_str = []

    for exp_name in key_list:
        eff_resol_str.append(f"{eff_resol_list[exp_name]:12.2f}")
        num_en_accum_str.append(f"{energy_accum_list[exp_name]:12.2f}")
        num_en_accum_dx_str.append(f"{num_energy_accum_equiv_dx_list[exp_name]:12.2f}")
        
    line = f"{'Effective Resolution':<22}, {', '.join(eff_resol_str)}"
    print(line); line_list.append(line+"\n")
    line = f"{'Num. Energy Accum.':<22}, {', '.join(num_en_accum_str)}"
    print(line); line_list.append(line+"\n")
    line = f"{'Num. Energy Accum. Dx':<22}, {', '.join(num_en_accum_dx_str)}"
    print(line); line_list.append(line+"\n")

    with open(out_fname, mode='w') as f:
        f.writelines(line_list)    