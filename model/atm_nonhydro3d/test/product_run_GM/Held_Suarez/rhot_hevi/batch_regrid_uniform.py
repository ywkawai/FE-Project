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

num_run_list = {
  "./Eh12Ez8P3":3, "./Eh24Ez16P3":2,   
  "./Eh6Ez4P7":4, "./Eh12Ez8P7":8, "./Eh24Ez16P7":2, "./Eh48Ez32P7":2,  
  "./Eh4Ez3P11":4, "./Eh8Ez6P11":6,   
}


#------------------------------------------------------------------
def run_job(exp_dir, num_run):
    old_path = os.getcwd()
    
    subdir_list_cand = []
    for runno in range(1,num_run+1):
        subdir_list_cand.append(f"run{runno}")

    add_flag = False; subdir_list = []
    for subdir in subdir_list_cand:
        if add_flag == False:
            add_flag = True
        if add_flag:        
            subdir_list.append(subdir)
    
    job_id = -1; subjob_ind = -1    
    for subdir in subdir_list:
        os.chdir(f"{exp_dir}/{subdir}")
        print(os.getcwd())
        
        job_id, subjob_ind = batch_job_common.pjsub(job_id, subjob_ind, "job_regrid_p_uniform.sh")        
        os.chdir(old_path)

for dir in dir_list:
    run_job(dir, num_run_list[dir])
