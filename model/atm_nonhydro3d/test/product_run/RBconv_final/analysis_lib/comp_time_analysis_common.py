import re

def get_prof_rapreport(log_file):
    print(log_file)
    prof_lines = []    
    with open(log_file, "r") as f:
        lines = f.readlines()    
        s_lines = [line.strip() for line in lines]
        flag = False
        for line in s_lines:
            if flag:
                prof_lines.append(line)
            if "INFO  [PROF_rapreport]" in line:
                flag = True
    
    prof_info_dict = {}
    for prof_line in prof_lines:
        m = re.match(r'(\S+\s\S+)\s*lev=\s*(\d): T\(avg\)=(.*), T\(max\)=(.*)\[(.*)\], T\(min\)=(.*)\[(.*)\], N=(.*)', prof_line)        
#        print(m.groups())
        prof_info = { "sect": m.group(1), "lev": int(m.group(2)), "Tavg": float(m.group(3)), 
                     "Tmax": float(m.group(4)), "Tmax_prc": int(m.group(5)),
                     "Tmin": float(m.group(6)), "Tmin_prc": int(m.group(7)), 
                     "N": int(m.group(8)) }
        prof_info_dict[m.group(1)] = prof_info
        
    return prof_info_dict

def output_comp_time(prof_info_rn, fname = None):
    sorted_prof_list = sorted(prof_info_rn[0].values(), key=lambda x:x["Tavg"], reverse=True)
    sect_list = []
    
    lines = []
    for prof_list in sorted_prof_list:
        if ( "ATM_" in prof_list["sect"] and prof_list["lev"]==2):
            sect_list.append(prof_list["sect"])
    for sect in sect_list:
        tavg_list = []
        tavg_mean = 0.0; tmax = 0.0; tmin = 1e20
        for prof_info in prof_info_rn:
            tavg = prof_info[sect]['Tavg']
            tmin_ = prof_info[sect]['Tmin']
            tmax_ = prof_info[sect]['Tmax']
            if tmax <= tmax_:
                tmax = tmax_
            if tmin >= tmin_:
                tmin = tmin_
            tavg_mean = tavg_mean + tavg
            tavg_list.append(f"{tavg:>10}")
        tavg_mean = tavg_mean / float(len(prof_info_rn))
        line = f"{sect:<40}: {(', ').join(tavg_list)} | {tavg_mean:8.1f} [sec] {(tmin-tavg_mean)/tavg_mean*100.0:5.1f} [%]  {(tmax-tavg_mean)/tavg_mean*100.0:4.1f} [%]"
        lines.append(line+"\n")
        print(line)
        
    if fname:
        with open(fname, "w", encoding="utf-8") as f:
            f.writelines(lines)

def get_sect_tmin(prof_info_rn, sect):
    tmin = 1e20
    ind_min = -1
    for ind, prof_info in enumerate(prof_info_rn):
        tavg = prof_info[sect]['Tavg']
        if (tavg < tmin):
            tmin = tavg
            ind_min = ind
#    print(f"{sect:<40}: {(', ').join(tavg_list)} | {tavg_mean:8.1f} [sec]")        
    return tmin, ind_min

def get_sect_tavg(prof_info_rn, sect):
    tavg_list = []
    tavg_mean = 0.0
    for prof_info in prof_info_rn:
        tavg = prof_info[sect]['Tavg']
        tavg_mean = tavg_mean + tavg
        tavg_list.append(f"{tavg:>10}")
    tavg_mean = tavg_mean / float(len(prof_info_rn))
#    print(f"{sect:<40}: {(', ').join(tavg_list)} | {tavg_mean:8.1f} [sec]")        
    return tavg_mean

def output_comp_resource_mean_usage(prof_info_list, sect_list, node_info_list, is_overlap_comp_comm=False, fname=None):

    nhinfo_list_focused_sect = {}
    if is_overlap_comp_comm:
        dycore_focused_sects = [ 'MAIN ATM_DYN_update_caltend_ex + ebnd_flux', 
                                'MAIN ATM_DYN_exchange_prgv + wait', 
                                'MAIN ATM_DYN_update_advance',
                                'MAIN ATM_DYN_applyBC_prgv',
                                'MAIN ATM_DYN_cal_pres', 
                                'MAIN ATM_DYN_update_modalfilter',
                                'MAIN ATM_DYN_update_pre', 
                                'MAIN ATM_DYN_update_post',
                                ]
    else:
        dycore_focused_sects = [ 'MAIN ATM_DYN_update_caltend_ex', 
                                'MAIN ATM_DYN_exchange_prgv', 
                                'MAIN ATM_DYN_exchange_prgv_wait',                                                                 
                                'MAIN ATM_DYN_update_advance',
                                'MAIN ATM_DYN_applyBC_prgv',
                                'MAIN ATM_DYN_cal_pres', 
                                'MAIN ATM_DYN_update_modalfilter',
                                'MAIN ATM_DYN_update_pre', 
                                'MAIN ATM_DYN_update_post',
                                ]

    line_list = []
    sect_str = "section"
    header_items = []
    header_items.append(f"{sect_str:<47}")
    for exp_name in prof_info_list.keys():
        header_items.append(f"{exp_name:>12}")
        nhinfo_list_focused_sect[exp_name] = 0.0
    line = ", ".join(header_items)
    print(line); line_list.append(line+"\n")
    
    
    for sect in sect_list:
        tinfo_list = []
        nhinfo_list = []
        if is_overlap_comp_comm and 'MAIN ATM_DYN_exchange_prgv_wait'==sect:
            sect_out = "MAIN ATM_DYN_exchange_prgv + wait"
        elif is_overlap_comp_comm and 'MAIN ATM_DYN_update_caltend_ex'==sect:
            sect_out = "MAIN ATM_DYN_update_caltend_ex + ebnd_flux"            
        else:
            sect_out = sect
        
        for key, prof_info in prof_info_list.items():
            if 'MAIN ATM_DYN_exchange_prgv' == sect_out:
                tinfo, tmin_runno_id = get_sect_tmin(prof_info, sect)
            elif is_overlap_comp_comm and 'MAIN ATM_DYN_exchange_prgv_wait' == sect:
                tinfo, tmin_runno_id = get_sect_tmin(prof_info, sect)
                prof_info_ = prof_info[tmin_runno_id]
                tinfo2 = prof_info_['MAIN ATM_DYN_exchange_prgv']['Tavg']
                tinfo = tinfo + tinfo2
            elif is_overlap_comp_comm and 'MAIN ATM_DYN_update_caltend_ex'==sect:
                tinfo = get_sect_tavg(prof_info, sect)
                tinfo2 = get_sect_tavg(prof_info, 'MAIN ATM_DYN_ebnd_flux')
                tinfo = tinfo + tinfo2                
            else:
                tinfo = get_sect_tavg(prof_info, sect)
            
            tinfo_list.append(f"{tinfo:12.1f}")
            
            nh = tinfo/3600.0*node_info_list[key]
            nhinfo_list.append(f"{nh:12.2f}")
            
            if sect_out in dycore_focused_sects:
                nhinfo_list_focused_sect[key] = nhinfo_list_focused_sect[key] + nh
            
        tinfo_list_s = ", ".join(tinfo_list)
        nhinfo_list_s = ", ".join(nhinfo_list)
        line = f"  {sect_out:<45}, {nhinfo_list_s}"
        print(line); line_list.append(line+"\n")


    sect_out = "MAIN ATM_DYN_FOCUSED_SECTS"
    nhinfo_list_s = ", ".join(list([ f"{v:12.2f}" for k, v in nhinfo_list_focused_sect.items()]))
    line = f"  {sect_out:<45}, {nhinfo_list_s}"
    print(line); line_list.append(line+"\n")
    

    if fname:
        with open(fname, mode='w') as f:
            f.writelines(line_list)    
