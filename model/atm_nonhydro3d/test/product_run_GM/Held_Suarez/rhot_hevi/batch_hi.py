import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../../common'))
import batch_job_common

dir_list = [
  "Eh48Ez32P3",
  "Eh48Ez32P7", 
  "Eh16Ez12P11",   
]

num_run_list = {
  "Eh48Ez32P3": 20, 
  "Eh48Ez32P7": 30,  
  "Eh16Ez12P11": 20, 
}

start_subdir_list = {
  "Eh48Ez32P3": "run1",
  "Eh48Ez32P7": "run1",  
  "Eh16Ez12P11": "run1", 
}

#------------------------------------------------------------------
    
def run_job(exp_dir, num_spinup, num_run, start_subdir):
    old_path = os.getcwd()
        
    subdir_list_cand = []
    for runno in range(1,num_spinup+1):
        subdir_list_cand.append(f"spinup{runno}")        
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
        job_id, subjob_ind = batch_job_common.pjsub(job_id, subjob_ind, "job_regrid.sh")            
        if "run" in subdir:
            job_id, subjob_ind = batch_job_common.pjsub(job_id, subjob_ind, "job_regrid_p.sh")            
        
        os.chdir(old_path)

for dir in dir_list:
    run_job(dir, 0, num_run_list[dir], start_subdir_list[dir])
