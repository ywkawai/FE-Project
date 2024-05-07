import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../../common'))
import batch_job_common

dir_list = [
#* p=3
   "./Eh12Ez8P3",  "./Eh24Ez16P3",
#* p=7
   "./Eh6Ez4P7", "./Eh12Ez8P7", "./Eh24Ez16P7", 
#* p=11
   "./Eh4Ez3P11", "./Eh8Ez6P11", 
]

num_spinup_run_list = {
  "./Eh12Ez8P3":2, "./Eh24Ez16P3":2,       
  "./Eh6Ez4P7":2, "./Eh12Ez8P7":2, "./Eh24Ez16P7":4,   
  "./Eh4Ez3P11":2, "./Eh8Ez6P11":2,     
}

num_run_list = {
  "./Eh12Ez8P3": 8, "./Eh24Ez16P3": 16,   
  "./Eh6Ez4P7": 8, "./Eh12Ez8P7": 16, "./Eh24Ez16P7": 20, 
  "./Eh4Ez3P11": 8, "./Eh8Ez6P11": 16,   
}

start_subdir_list = {
  "./Eh12Ez8P3": "spinup1", "./Eh24Ez16P3": "spinup1", 
  "./Eh6Ez4P7": "spinup1", "./Eh12Ez8P7": "spinup1", "./Eh24Ez16P7": "spinup1",
  "./Eh4Ez3P11": "spinup1", "./Eh8Ez6P11": "spinup1",
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
    run_job(dir, num_spinup_run_list[dir], num_run_list[dir], start_subdir_list[dir])
