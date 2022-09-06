import textwrap
import os

def gen_conf( out_conf_dir,  
            out_dir, 
            beta, 
            basis_type, form_type, 
            MF_order, MF_alph, MF_alph_lbl, 
            Cr_DG_dxeff, 
            tlev_list,                        
            tscheme_name ):

    s = f'''
      &PARAM_ANALYSIS
        OUT_DIR="{out_dir}", 
        ADV_VEL=1.0,     
        Helem=1.0,       
        BETA={beta},        
        basis_type="{basis_type}",  
        form_type="{form_type}",   
        MF_order={MF_order},    
        MF_alph={MF_alph},     
        MF_alph_lbl="{MF_alph_lbl}", 
        Cr_DG_dxeff={Cr_DG_dxeff},   
        tlev_num={len(tlev_list)},    
        tlev_slot={", ".join(map(str,tlev_list))},  
        tscheme_name="{tscheme_name}", 
      /   
    '''
#    print(textwrap.dedent(s)[1:-1])
    with open(f"{out_conf_dir}/analysis.conf", mode="w") as f:
      f.write(textwrap.dedent(s)[1:-1])
    
Cr_DG_dxeff = 5.0 * 0.0125 / 10.0
time_list = [ 1, 10, 50, 100, 200, 500, 1000, 3000 ]
beta      = 1.0

#out_conf_dir = 'stdupwind_semidiscrete/'; tscheme_name =""; MForder=1; MFalph=0.0; MFlbl=""
#out_conf_dir = 'stdupwind_RK3/'; tscheme_name ="RK3"; MForder=1; MFalph=0.0; MFlbl=""
#out_conf_dir = 'stdupwind_RK4/'; tscheme_name ="RK4"; MForder=1; MFalph=0.0; MFlbl=""
#out_conf_dir = 'stdupwind_RKo4s10/'; tscheme_name ="RKo4s10"; MForder=1; MFalph=0.0; MFlbl=""
out_conf_dir = 'stdupwind_MF32Alph1E-3_RKo4s10/'; tscheme_name ="RKo4s10"; MForder=32; MFalph=1e-3; MFlbl="MF32Alph1E-3"

# out_conf_dir = 'stdupwind_MF32Alph1E-3_RK4/'; tscheme_name ="RK4"; MForder=32; MFalph=1e-3; MFlbl="MF32Alph1E-3"
# out_conf_dir = 'stdupwind_MF32Alph1E-2_RK4/'; tscheme_name ="RK4"; MForder=32; MFalph=1e-2; MFlbl="MF32Alph1E-2"
# out_conf_dir = 'stdupwind_MF16Alph1E-2_RK4/'; tscheme_name ="RK4"; MForder=16; MFalph=1e-2; MFlbl="MF16Alph1E-2"
# out_conf_dir = 'stdupwind_MF08Alph1E-2_RK4/'; tscheme_name ="RK4"; MForder= 8; MFalph=1e-2; MFlbl="MF08Alph1E-2"

#out_conf_dir = 'stdupwind_RK3/'; tscheme_name ="RK3"; MForder=1; MFalph=0.0; MFlbl=""
#out_conf_dir = 'stdupwind_RK4/'; tscheme_name ="RK4"; MForder=1; MFalph=0.0; MFlbl=""
#out_conf_dir = 'stdupwind_RKo4s10/'; tscheme_name ="RKo4s10"; MForder=1; MFalph=0.0; MFlbl=""

# Cr_DG_dxeff = 50.0 * 5.0 * 0.0125 / 10.0
# out_conf_dir = 'stdupwind_Crx50_semidiscrete/'; tscheme_name =""; MForder=1; MFalph=0.0; MFlbl=""
# out_conf_dir = 'stdupwind_Crx50_RK3/'; tscheme_name ="RK3"; MForder=1; MFalph=0.0; MFlbl=""
# out_conf_dir = 'stdupwind_Crx50_RK4/'; tscheme_name ="RK4"; MForder=1; MFalph=0.0; MFlbl=""
# out_conf_dir = 'stdupwind_Crx50_RKo4s10/'; tscheme_name ="RKo4s10"; MForder=1; MFalph=0.0; MFlbl=""

# Cr_DG_dxeff = 190.0 * 5.0 * 0.0125 / 10.0
# out_conf_dir = 'stdupwind_Crx190_RKo4s10/'; tscheme_name ="RKo4s10"; MForder=1; MFalph=0.0; MFlbl=""

# beta     = 70.0
# out_conf_dir = 'ovupwind_RKo4s10/'; tscheme_name ="RKo4s10"; MForder=1; MFalph=0.0; MFlbl=""
# out_conf_dir = 'ovupwind_MF32Alph1E-3_RKo4s10/'; tscheme_name ="RKo4s10"; MForder=32; MFalph=1e-3; MFlbl="MF32Alph1E-3"


os.makedirs(out_conf_dir, exist_ok=True)
gen_conf( out_conf_dir, 
        "data/", beta, 'nodal', 'strong', 
        MForder, MFalph, MFlbl, Cr_DG_dxeff, time_list, tscheme_name )