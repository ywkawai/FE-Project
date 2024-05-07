import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../../common'))
import batch_job_common

dir_list = [
# * p=7
    "Eh64Ez34P7", 
# * p=4
    "Eh100Ez52P4",  
# * p=3    
    "Eh128Ez64P3",
# * p=7 (without no shallow atmosphere approximation)
    "Eh64Ez34P7_deepatm",     
# * p=4 (without no shallow atmosphere approximation)    
    "Eh100Ez52P4_deepatm",                 
# * p=3 (without no shallow atmosphere approximation)        
    "Eh128Ez64P3_deepatm",         
]

num_run_list = {
  "Eh64Ez34P7": 8,  "Eh100Ez52P4": 8, "Eh128Ez64P3": 8,
  "Eh64Ez34P7_deepatm": 8,  "Eh100Ez52P4_deepatm": 8, "Eh128Ez64P3_deepatm": 8,
}

start_subdir_list = {
  "Eh64Ez34P7": "run1", "Eh128Ez64P3": "run7", "Eh100Ez52P4": "run7",
  "Eh64Ez34P7_deepatm": "run1", "Eh128Ez64P3_deepatm": "run7", "Eh100Ez52P4_deepatm": "run7",
}

    
def run_job(exp_dir, num_run, start_subdir):
    old_path = os.getcwd()
    
    subdir_list_cand = []
    for runno in range(1,num_run+1):
        subdir_list_cand.append(f"run{runno}")

    add_flag = False; subdir_list = []
    for subdir in subdir_list_cand:
        if add_flag == False and subdir == start_subdir:
            add_flag = True
        if add_flag:        
            subdir_list.append(subdir)

    job_id = -1; subjob_ind = -1            
    for subdir in subdir_list:
        os.chdir(f"{exp_dir}/{subdir}")
        print(os.getcwd())
        
        job_id, subjob_ind = batch_job_common.pjsub(job_id, subjob_ind, "job_run.sh")
#        job_id, subjob_ind = batch_job_common.pjsub(job_id, subjob_ind, "job_regrid.sh")    
                
        os.chdir(old_path)

for dir in dir_list:
    run_job(dir, num_run_list[dir], start_subdir_list[dir])
