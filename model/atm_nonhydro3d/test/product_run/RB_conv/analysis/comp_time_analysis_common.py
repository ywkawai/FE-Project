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

def output_comp_time(prof_info_rn):
    sorted_prof_list = sorted(prof_info_rn[0].values(), key=lambda x:x["Tavg"], reverse=True)
    sect_list = []
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
        print(f"{sect:<40}: {(', ').join(tavg_list)} | {tavg_mean:8.1f} [sec] {(tmin-tavg_mean)/tavg_mean*100.0:5.1f} [%]  {(tmax-tavg_mean)/tavg_mean*100.0:4.1f} [%]")    