import os
import sys
sys.path.append(os.path.join(os.path.dirname('__file__'), '../../analysis_lib/'))
import comp_cost_analysis_common as common
#---------------------------

DG_dir = "./"
FV_dir = "../../../../../../../../scale_rbconv/scale-rm/test/case/rbconv_productrun/DNS/"
COMP_COST_OUTDIR = "comp_cost_info/"
COMP_COST_FIG_OUTDIR = "fig_comp_cost_info/"
FIG_EXT = "svg"

dt_info_dg = { "P3": 5.0, "P7": 3.125, "P11": 2.5}
ratio_cr_max_dg = {"P3": 1.25, "P7": 1.28, "P11": 1.02}
rkstage_dg = 10

dt_info_fv = { "UD3": 5.0, "UD7": 5.0}
ratio_cr_max_fv = {"UD3": 1.20, "UD7": 1.20}
rkstage_fv = 4

correct_fac_P11 = 1.3

integ_time_per_run = 7200.0 # [s] # In collect_comp_time.py, the required computational resources for all cases are normalized by integ_time_per_run = 7200 [s]

sects_fv2dg = {
    'MAIN ATM_DYN_upadte_caltend_ex':  'MAIN ATM_DYN_update_caltend_ex + ebnd_flux',
    'MAIN ATM_DYN_applyBC_prgv': 'MAIN ATM_DYN_applyBC_prgv',
    'MAIN ATM_DYN_update_pre': 'MAIN ATM_DYN_update_pre', 
    'MAIN ATM_DYN_update_post': 'MAIN ATM_DYN_update_post',             
    'MAIN ATM_DYN_cal_pres': 'MAIN ATM_DYN_cal_pres',
    'MAIN ATM_DYN_exchange_prgv': 'MAIN ATM_DYN_exchange_prgv + wait',
    'MAIN ATM_DYN_FOCUSED_SECTS': 'MAIN ATM_DYN_FOCUSED_SECTS' }

Exp_label_list_dg = {
  "Dx25m_P3":  "Δ=25m, P3",  
  "Dx25m_P7":  "Δ=25m, P7",   
  "Dx27m_P11":  "Δ=27m, P11",  
  "Dx12.5m_P3":  "Δ=13m, P3",
  "Dx12.5m_P7":  "Δ=13m, P7",     
  'Dx13m_P11': 'Δ=13m, P11', 
  'Dx6.3m_P3': 'Δ=6.3m, P3',
  'Dx6.3m_P7': 'Δ=6.3m, P7',  
  'Dx6.7m_P11': 'Δ=6.7m, P11', 
  'Dx3.1m_P7': 'Δ=3.1m, P7', 
}
Exp_label_list_fv = {
  "Dx25m_UD3":  "Δ=25m, UD3",  
  "Dx25m_UD7":  "Δ=25m, UD7",   
  "Dx12.5m_UD3":  "Δ=13m, UD3",
  "Dx12.5m_UD7":  "Δ=13m, UD7",     
  'Dx6.25m_UD3': 'Δ=6.3m, UD3',
  'Dx6.25m_UD7': 'Δ=6.3m, UD7',  
}

EXP_groups_dg = {
    "Group1": ["Dx25m_P3", "Dx25m_P7", "Dx27m_P11"],
    "Group2": ["Dx12.5m_P3", "Dx12.5m_P7", "Dx13m_P11"],
    "Group3": ["Dx6.3m_P3", "Dx6.3m_P7", "Dx6.7m_P11"],
}
EXP_groups_fv = {
    "Group1": ["Dx25m_UD3", "Dx25m_UD7"],
    "Group2": ["Dx12.5m_UD3", "Dx12.5m_UD7"],
    "Group3": ["Dx6.25m_UD3", "Dx6.25m_UD7"],
}
EXP_group_colors = {"Group1": "gray", "Group2": "red", "Group3": "green"}


# DG
sects_label_info = { 
    'MAIN ATM_DYN_update_caltend_ex + ebnd_flux': "RHS eval.", 
    'MAIN ATM_DYN_exchange_prgv + wait': "Halo handling, Data comm.", 
    # 'MAIN ATM_DYN_update_advance': "RK scheme",
    'MAIN ATM_DYN_applyBC_prgv': "Boundary condition",
    'MAIN ATM_DYN_cal_pres': "Pressure calculation", 
    'MAIN ATM_DYN_update_modalfilter': "Modal filtering",
    'MAIN ATM_DYN_update_pre': "Pre proc.", 
    'MAIN ATM_DYN_update_post': "Post proc." }

sects_color_info = { 
    'MAIN ATM_DYN_update_caltend_ex + ebnd_flux': "orange", 
    'MAIN ATM_DYN_exchange_prgv + wait': "deepskyblue", 
    'MAIN ATM_DYN_update_advance': "darkviolet",
    'MAIN ATM_DYN_applyBC_prgv': "orangered",
    'MAIN ATM_DYN_cal_pres': "pink", 
    'MAIN ATM_DYN_update_modalfilter': "gray",
    'MAIN ATM_DYN_update_pre': "green", 
    'MAIN ATM_DYN_update_post': "lightgreen" }
# FV
sects_label_info_fv = { 
    'MAIN ATM_DYN_update_caltend_ex + ebnd_flux': "RHS eval.", 
    'MAIN ATM_DYN_exchange_prgv + wait': "Halo handling, Data comm.", 
    # 'MAIN ATM_DYN_update_advance': "RK scheme",
    'MAIN ATM_DYN_applyBC_prgv': "Boundary condition",
    'MAIN ATM_DYN_cal_pres': "Pressure calculation", 
    'MAIN ATM_DYN_update_modalfilter': "Modal filtering",
    'MAIN ATM_DYN_update_pre': "Pre proc.", 
    'MAIN ATM_DYN_update_post': "Post proc." }

sects_label_info_fv = { 
    'MAIN ATM_DYN_upadte_caltend_ex': "RHS eval.",
    'MAIN ATM_DYN_exchange_prgv': "Halo handling, Data comm.",  
    'MAIN ATM_DYN_applyBC_prgv': 'Boundary condition',
    'MAIN ATM_DYN_cal_pres': "Pressure calculation",    
    'MAIN ATM_DYN_update_pre': "Pre proc.", 
    'MAIN ATM_DYN_update_post': "Post proc.",}
sects_color_info_fv = { 
    'MAIN ATM_DYN_upadte_caltend_ex': "orange",
    'MAIN ATM_DYN_applyBC_prgv': 'orangered',
    'MAIN ATM_DYN_update_pre': "green", 
    'MAIN ATM_DYN_update_post': "lightgreen",             
    'MAIN ATM_DYN_cal_pres': "pink", 
    'MAIN ATM_DYN_exchange_prgv': "deepskyblue",}

#----------------------------------------
# Funcations
resol_table = {"Dx25m": "Dx25m", "Dx27m": "Dx25m",
                "Dx12.5m": "Dx12.5m", "Dx13m": "Dx12.5m",
                "Dx6.25m": "Dx6.3m", "Dx6.3m": "Dx6.3m", "Dx6.7m": "Dx6.3m",}

#---
ds_comp_resource_dg = common.get_comp_resource_usage(f"{DG_dir}/{COMP_COST_OUTDIR}/comp_resourse_usage.dat")
eff_resol_dg, num_en_accum_dg, num_en_accum_dx_dg = common.get_EffResol_NumEnAccum(f"{DG_dir}/{COMP_COST_OUTDIR}/EffResol_NumEnAccum.dat")

for exp_case in ds_comp_resource_dg.case.values:
    if "P11" in exp_case:
        for sect in ds_comp_resource_dg.section.values:
            ds_comp_resource_dg.NH.loc[dict(case=exp_case,section=sect)] = correct_fac_P11 * ds_comp_resource_dg.sel(case=exp_case,section=sect).NH.values


ds_comp_resource_fv = common.get_comp_resource_usage(f"{FV_dir}/{COMP_COST_OUTDIR}/comp_resourse_usage.dat")
eff_resol_fv, num_en_accum_fv, num_en_accum_dx_fv = common.get_EffResol_NumEnAccum(f"{FV_dir}/{COMP_COST_OUTDIR}/EffResol_NumEnAccum.dat")

common.print_sect_nh(ds_comp_resource_dg)
common.print_sect_nh(ds_comp_resource_fv)

#-
# DG
comp_cost_list_dg = {}
comp_cost_normalized_list_dg = {}

rhs_eval_num_info_dg = common.get_rhs_eval_numinfo(dt_info_dg, rkstage_dg, integ_time_per_run)
print(rhs_eval_num_info_dg)

for sect in ds_comp_resource_dg.section.values:
  print(sect)
  comp_cost_list_dg[sect], comp_cost_normalized_list_dg[sect] = common.get_comp_cost_per_eval_rhs(
      ds_comp_resource_dg.sel(section=sect), 
      rhs_eval_num_info_dg, rhs_eval_num_info_dg, resol_table )
  
# FV
comp_cost_list_fv = {}
comp_cost_normalized_list_fv = {}

rhs_eval_num_info_fv = common.get_rhs_eval_numinfo(dt_info_fv, rkstage_fv, integ_time_per_run)
print(rhs_eval_num_info_fv)

for sect in ds_comp_resource_fv.section.values:
  print(sect)
  sect_dg = sects_fv2dg[sect]
  comp_cost_list_fv[sect_dg], comp_cost_normalized_list_fv[sect_dg] = common.get_comp_cost_per_eval_rhs(
     ds_comp_resource_fv.sel(section=sect), 
     rhs_eval_num_info_fv, rhs_eval_num_info_dg, resol_table, comp_cost_list_dg[sect_dg])

#--
ds_comp_resource_dg_ = ds_comp_resource_dg.copy(deep=True)

for exp_case in Exp_label_list_dg.keys():
    if not "Dx3.1m_P7" in exp_case:
        selA = dict(case=exp_case, section="MAIN ATM_DYN_update_caltend_ex + ebnd_flux")
        selB = dict(case=exp_case, section="MAIN ATM_DYN_update_advance")
        ds_comp_resource_dg_["NH"].loc[selA] = ds_comp_resource_dg["NH"].sel(**selA) + ds_comp_resource_dg["NH"].sel(**selB)

if not os.path.exists(COMP_COST_FIG_OUTDIR):
    os.makedirs(COMP_COST_FIG_OUTDIR)

common.mkgraph_compcost( ds_comp_resource_dg_, sects_label_info, resol_table, 
                        EXP_groups_dg, Exp_label_list_dg, EXP_group_colors, sects_color_info, 
                        ratio_cr_max_dg, ratio_cr_max_dg, 
                        f"{COMP_COST_FIG_OUTDIR}/comp_cost_dg_breakdown_integ_period.{FIG_EXT}", 
                        ylim=[0.0, 3.4] )
common.mkgraph_compcost( ds_comp_resource_dg_, sects_label_info, resol_table, 
                        EXP_groups_dg, Exp_label_list_dg, EXP_group_colors, sects_color_info, 
                        ratio_cr_max_dg, ratio_cr_max_dg, 
                        f"{COMP_COST_FIG_OUTDIR}/comp_cost_dg_breakdown_perRKstage.{FIG_EXT}", 
                        True, rhs_eval_num_info_dg, rhs_eval_num_info_dg, 
                        ylim=[0.0,2.35] )

common.mkgraph_compcost( ds_comp_resource_fv, sects_label_info_fv, resol_table, 
                        EXP_groups_fv, Exp_label_list_fv, EXP_group_colors, sects_color_info_fv, 
                        ratio_cr_max_fv, ratio_cr_max_dg, 
                        f"{COMP_COST_FIG_OUTDIR}/comp_cost_fv_breakdown_integ_period.{FIG_EXT}", 
                        False, None, None, ds_comp_resource_dg, 
                        ylim=[0.0, 3.4] )
common.mkgraph_compcost( ds_comp_resource_fv, sects_label_info_fv, resol_table, 
                        EXP_groups_fv, Exp_label_list_fv, EXP_group_colors, sects_color_info_fv, 
                        ratio_cr_max_fv, ratio_cr_max_dg, 
                        f"{COMP_COST_FIG_OUTDIR}/comp_cost_fv_breakdown_perRKstage.{FIG_EXT}", 
                        True, rhs_eval_num_info_fv, rhs_eval_num_info_dg, ds_comp_resource_dg, 
                        ylim=[0.0, 2.35])

#--
comp_cost_eff_resol_dg, comp_cost_eff_resol_normalized_dg = common.get_comp_cost_cost_with_eff_resol( 
    ds_comp_resource_dg.sel(section="MAIN ATM_DYN_FOCUSED_SECTS"), 
    ds_comp_resource_dg.sel(section="MAIN ATM_DYN_FOCUSED_SECTS"), 
    eff_resol_dg, eff_resol_dg, ratio_cr_max_dg, ratio_cr_max_dg, resol_table )

comp_cost_num_en_accum_dg, comp_cost_num_en_accum_normalized_dg = common.get_comp_cost_cost_with_num_en_accum( 
    ds_comp_resource_dg.sel(section="MAIN ATM_DYN_FOCUSED_SECTS"), 
    ds_comp_resource_dg.sel(section="MAIN ATM_DYN_FOCUSED_SECTS"), 
    num_en_accum_dx_dg, ratio_cr_max_dg, ratio_cr_max_dg, resol_table )

comp_cost_eff_resol_fv, comp_cost_eff_resol_normalized_fv = common.get_comp_cost_cost_with_eff_resol( 
    ds_comp_resource_fv.sel(section="MAIN ATM_DYN_FOCUSED_SECTS"), 
    ds_comp_resource_dg.sel(section="MAIN ATM_DYN_FOCUSED_SECTS"), 
    eff_resol_fv, eff_resol_dg, ratio_cr_max_fv, ratio_cr_max_dg, resol_table )
comp_cost_num_en_accum_fv, comp_cost_num_en_accum_normalized_fv = common.get_comp_cost_cost_with_num_en_accum( 
    ds_comp_resource_fv.sel(section="MAIN ATM_DYN_FOCUSED_SECTS"), 
    ds_comp_resource_dg.sel(section="MAIN ATM_DYN_FOCUSED_SECTS"), 
    num_en_accum_dx_fv, ratio_cr_max_fv, ratio_cr_max_dg, resol_table )


common.mkgraph_compcost_eff_resol_en_accum( comp_cost_eff_resol_dg, comp_cost_num_en_accum_dg, comp_cost_eff_resol_fv, comp_cost_num_en_accum_fv, 
                                            Exp_label_list_dg, Exp_label_list_fv, 
                                            f"{COMP_COST_FIG_OUTDIR}/comp_cost_diagram.{FIG_EXT}" )
common.mkgraph_compcost_eff_resol_en_accum( comp_cost_eff_resol_dg, comp_cost_num_en_accum_dg, comp_cost_eff_resol_fv, comp_cost_num_en_accum_fv, 
                                            Exp_label_list_dg, Exp_label_list_fv, 
                                            f"{COMP_COST_FIG_OUTDIR}/comp_cost_diagram_with_lbl.{FIG_EXT}", True )  