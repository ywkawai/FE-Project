import pandas as pd
import xarray as xr
import re
import numpy as np
import matplotlib.pyplot as plt

#---------------------------

def get_comp_resource_usage(fname):
  df = pd.read_csv(fname, skipinitialspace=True)
  df.columns = df.columns.str.strip()
  df["section"] = df["section"].str.strip()

  ser = df.set_index("section").stack()
  ser.index = ser.index.set_names(["section", "case"])
  da = ser.to_xarray()
  da.name = "NH"
  return da.to_dataset() 

def get_EffResol_NumEnAccum(fname):
  df = pd.read_csv(fname)
  df.columns = df.columns.str.strip()
  df["info"] = df["info"].str.strip()  
  df = df.set_index("info")
  ds = xr.Dataset(
      {
          "Effective_Resolution": (["case"], df.loc["Effective Resolution"].values),
          "Num_Energy_Accum": (["case"], df.loc["Num. Energy Accum."].values),
          "Num_Energy_Accum_Dx": (["case"], df.loc["Num. Energy Accum. Dx"].values),          
      },
      coords={"case": df.columns}
  )
  return ds.Effective_Resolution, ds.Num_Energy_Accum, ds.Num_Energy_Accum_Dx

def print_sect_nh(ds_comp_resource):
  print(ds_comp_resource.case.values)
  for sec in ds_comp_resource.section.values:
    nh_str = []
    for case in ds_comp_resource.case.values:
      nh = ds_comp_resource.sel(section=sec, case=case).NH
      nh_str.append(f"{nh:12.2f}")
    nh_str = " ,".join(nh_str)
    print(f"{sec:<45}: {nh_str}")

def get_rhs_eval_numinfo(dt_info, rkstage, integ_time_per_run):
  rhs_eval_num_info = {}
  for key, dt in dt_info.items():
    rhs_eval_num_info[key] = integ_time_per_run/dt * rkstage
  return rhs_eval_num_info


def get_comp_cost_per_eval_rhs(ds, rhs_eval_num_info, rhs_eval_num_info_dg, resol_table, comp_cost_per_eval_rhs_dg=None):
  case_list = ds.case.values
  
  comp_cost_per_eval_rhs = {}
  comp_cost_per_eval_rhs_normalized = {}
  
  normalized_fac = {}
  for key, eval_num in rhs_eval_num_info.items():
    normalized_fac[key] = rhs_eval_num_info_dg["P3"]/eval_num
  print(normalized_fac)
  
  for case in case_list:
    porder = case.split("_")[-1]
    nh = ds.sel(case=case).NH * normalized_fac[porder]
    comp_cost_per_eval_rhs[case] = float(nh)
    print(f"{case} {float(nh):12.2f} [NH]")
          
  for case in case_list:
    resol = case.split("_")[0]
    case_ref = resol_table[resol] + "_P3"
    if comp_cost_per_eval_rhs_dg:
      comp_cost_per_eval_rhs_ref = comp_cost_per_eval_rhs_dg[case_ref] 
    else:
      comp_cost_per_eval_rhs_ref = comp_cost_per_eval_rhs[case_ref]
    
    comp_cost_per_eval_rhs_normalized[case] = comp_cost_per_eval_rhs[case] / comp_cost_per_eval_rhs_ref
    print(f"{case} {float(comp_cost_per_eval_rhs_normalized[case]):12.2f}")
    if "P11" in case or "UD7" in case:
      print("")
    
  return comp_cost_per_eval_rhs, comp_cost_per_eval_rhs_normalized

def get_comp_cost_cost_with_eff_resol( ds, ds_dg, eff_resol_list, eff_resol_list_dg, cr_max_ratio,  cr_max_ratio_dg, resol_table):
  case_list = eff_resol_list.case.values
  print(case_list)
  
  comp_cost = {}
  comp_cost_normalized = {}
    
  for case in case_list:
    resol = case.split("_")[0]    
    porder = case.split("_")[-1]
    case_ref = resol_table[resol] + "_P3"

    nh = ds.sel(case=case).NH / cr_max_ratio[porder]
    comp_cost[case] = float(nh) * ( float(eff_resol_list.sel(case=case)) / float(eff_resol_list_dg.sel(case=case_ref)) )**4
    
    print(f"{case} {comp_cost[case]:12.2f}")
      
  for case in case_list:
    resol = case.split("_")[0]
    porder = case.split("_")[-1]
    case_ref = resol_table[resol] + "_P3"
    
    nh_ref = ds_dg.sel(case=case_ref).NH / cr_max_ratio_dg["P3"]
    comp_cost_normalized[case] = comp_cost[case] / float(nh_ref)
    print(f"{case} {float(comp_cost_normalized[case]):12.2f}")

    if "P11" in case or "UD7" in case:
      print("")
    
  return comp_cost, comp_cost_normalized

def get_comp_cost_cost_with_num_en_accum( ds, ds_dg, num_en_accum_dx, cr_max_ratio,  cr_max_ratio_dg, resol_table):
  case_list = num_en_accum_dx.case.values
  print(case_list)
  
  comp_cost = {}
  comp_cost_normalized = {}
    
  for case in case_list:
    resol = case.split("_")[0]
    porder = case.split("_")[-1]
    dx = float(re.findall(r'\d+\.?\d*', resol)[0])
    case_ref = resol_table[resol] + "_P3"
    
    nh = ds.sel(case=case).NH / cr_max_ratio[porder]
    comp_cost[case] = float(nh) * ( dx / float(num_en_accum_dx.sel(case=case)) )**4
    
    print(f"{case} {comp_cost[case]:12.2f}")
      
  for case in case_list:
    resol = case.split("_")[0]
    case_ref = resol_table[resol] + "_P3"

    nh_ref = ds_dg.sel(case=case_ref).NH / cr_max_ratio_dg["P3"]
    comp_cost_normalized[case] = comp_cost[case] / nh_ref
    print(f"{case} {float(comp_cost_normalized[case]):12.2f}")

    if "P11" in case or "UD7" in case:
      print("")
    
  return comp_cost, comp_cost_normalized


def mkgraph_compcost(comp_cost_dat, breakdown_sects, resol_table, 
                     EXP_groups, Exp_label_list, 
                     EXP_group_colors, Sect_color_list, 
                     cr_max_ratio, cr_max_ratio_dg, 
                     out_fname,
                     is_per_rhs=False, rhs_eval_num_info=None, rhs_eval_num_info_dg=None, 
                     comp_cost_dat_dg=None, ylim=[0.0,1.5]):

    offset = 0
    gap = 0.5
    labels = []; values = []; colors = []; positions = []
    scheme_num = len(EXP_groups["Group1"])

    comp_cost_dat_focused_sects = comp_cost_dat.sel(section="MAIN ATM_DYN_FOCUSED_SECTS")
    if comp_cost_dat_dg:
        comp_cost_dat_focused_sects_ref = comp_cost_dat_dg.sel(section="MAIN ATM_DYN_FOCUSED_SECTS").NH 
    else:
        comp_cost_dat_focused_sects_ref = comp_cost_dat.sel(section="MAIN ATM_DYN_FOCUSED_SECTS").NH 
    
    for group, keys in EXP_groups.items():
        for i, k in enumerate(keys):
            labels.append(Exp_label_list[k])
            values.append(comp_cost_dat_focused_sects.sel(case=k).NH)
            colors.append(EXP_group_colors[group])
            positions.append(offset + i)
        offset += len(keys) + gap

    values_sect = {}; colors_sect = {}; break_down_ratio_sect = {}
    for sect in breakdown_sects:
        comp_cost = comp_cost_dat.sel(section=sect)

        value_tmp = []; color_sect_tmp = []; break_down_ratio_tmp = []
        for group, keys in EXP_groups.items():
            for i, k in enumerate(keys):
                resol = k.split("_")[0]
                porder = k.split("_")[-1]
                case_ref = resol_table[resol] + "_P3"
                
                if is_per_rhs:
                    normalized_tot = (  comp_cost_dat_focused_sects.sel(case=k).NH / rhs_eval_num_info[porder] ) / ( comp_cost_dat_focused_sects_ref.sel(case=case_ref) / rhs_eval_num_info_dg["P3"] )
                else:
                    normalized_tot = (  comp_cost_dat_focused_sects.sel(case=k).NH / cr_max_ratio[porder] ) / ( comp_cost_dat_focused_sects_ref.sel(case=case_ref) / cr_max_ratio_dg["P3"] )
                
                break_down_ratio = comp_cost.sel(case=k).NH / comp_cost_dat_focused_sects.sel(case=k).NH

                break_down_ratio_tmp.append( break_down_ratio )
                value_tmp.append( break_down_ratio * normalized_tot )
                color_sect_tmp.append(Sect_color_list[sect])
                # print(f"{sect} {Exp_label_list[k]}: {comp_cost[k]} NH, {comp_cost[k]/comp_cost_dat_focused_sects[k]*100.0} %" )
        values_sect[sect] = value_tmp
        colors_sect[sect] = color_sect_tmp
        break_down_ratio_sect[sect] = break_down_ratio_tmp

    fig, ax = plt.subplots(1,1,figsize=(2*scheme_num+0.5,5))
    
    height = np.array(values)*0.0
    bottom = np.array(values)*0.0
    print(f"pos :{positions}")
    for section, sect_label in breakdown_sects.items():
        height = height + np.array(values_sect[section])
        print("")
        print(f"{section}: {bottom}")
        print(f"{section}: {height}")
        bars = ax.bar(positions, values_sect[section], bottom=bottom, color=colors_sect[section], 
                      edgecolor='white', linewidth=0.6, label=sect_label )
        
        if not is_per_rhs:
            for i, bar in enumerate(bars):
                if (break_down_ratio_sect[section][i] > 0.1):
                    plt.text(
                        bar.get_x() + bar.get_width() / 2,  
                        bottom[i] + 0.35 * values_sect[section][i],                             
                        f"{break_down_ratio_sect[section][i]*100:.0f}%",                    
                        ha="center", va="bottom", fontsize=10 )
    
        bottom = height

    for i, bar in enumerate(bars):
        plt.text(
            bar.get_x() + bar.get_width() / 2,  
            bottom[i] + 0.02,                             
            f"{bottom[i]:.2f}",                    
            ha="center", va="bottom", fontsize=12 )


    # ax.set_title(sect)
    ax.set_ylim(ylim)
    ax.legend(bbox_to_anchor=(1.02, 1.0), loc="upper left", borderaxespad=0.0)

    ax.set_xticks(positions, labels, rotation=30, ha="right")
    plt.savefig(out_fname, bbox_inches="tight")

def mkgraph_compcost_rhs_normalize(sect, comp_cost_dat, EXP_groups, Exp_label_list, EXP_group_colors, out_fname):

    offset = 0
    gap = 0.5
    labels = []; values = []; colors = []; positions = []
    scheme_num = len(EXP_groups["Group1"])

    for group, keys in EXP_groups.items():
        for i, k in enumerate(keys):
            labels.append(Exp_label_list[k])
            values.append(comp_cost_dat[k])
            colors.append(EXP_group_colors[group])
            positions.append(offset + i)
        offset += len(keys) + gap

    fig, ax = plt.subplots(1,1,figsize=(2*scheme_num+0.5+5,4))
    bars = ax.bar(positions, values, color=colors)
    ax.set_title(sect)
    ax.set_ylim([0.0,2.3])
    ax.set_xticks(positions, labels, rotation=30, ha="right")
    
    for bar in bars:
        height = bar.get_height()
        plt.text(
            bar.get_x() + bar.get_width() / 2,  
            height,                             
            f"{height:.2f}",                    
            ha="center", va="bottom", fontsize=12
        )

    plt.savefig(out_fname, bbox_inches="tight")


def mkgraph_compcost_eff_resol_en_accum( comp_cost_eff_resol_dg, comp_cost_num_en_accum_dg, comp_cost_eff_resol_fv, comp_cost_num_en_accum_fv, 
                                         Exp_label_list_dg, Exp_label_list_fv, out_fname, visible_text=False, ylim=[1e0,2e3] ):
  fig, ax = plt.subplots(1,1, figsize=(10,8))
  ax.set_xscale("log")
  ax.set_yscale("log")
  for key, comp_cost_eff_resol in comp_cost_eff_resol_dg.items():
      ax.plot( comp_cost_eff_resol, comp_cost_num_en_accum_dg[key], "o", color="darkblue")
      if visible_text:
        ax.text(comp_cost_eff_resol, comp_cost_num_en_accum_dg[key], 
                Exp_label_list_dg[key], fontsize=8, ha="left", va="bottom") 

  for key, comp_cost_eff_resol in comp_cost_eff_resol_fv.items():
      ax.plot( comp_cost_eff_resol, comp_cost_num_en_accum_fv[key], "o", color="orange")
      if visible_text:
        ax.text(comp_cost_eff_resol, comp_cost_num_en_accum_fv[key], 
                Exp_label_list_fv[key], fontsize=8, ha="left", va="bottom") 

  ax.set_xlim([1e0, 1e4])
  ax.set_ylim(ylim)
  ax.tick_params(which="both", labelsize=15)
  ax.grid()
  ax.grid(which="minor", linestyle=":", linewidth=0.7)
  plt.savefig(out_fname, bbox_inches="tight")

