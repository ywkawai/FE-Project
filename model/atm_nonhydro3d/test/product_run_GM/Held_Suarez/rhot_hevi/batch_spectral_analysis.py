import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../../common'))
import batch_job_common

dir_list = [
#* p=3
   "./Eh12Ez8P3",  "./Eh24Ez16P3", "./Eh48Ez32P3", 
#* p=7
   "./Eh6Ez4P7", "./Eh12Ez8P7",  "./Eh24Ez16P7", "./Eh48Ez32P7", 
#* p=11    
   "./Eh4Ez3P11", "./Eh8Ez6P11", "./Eh16Ez12P11"    
]

num_run_list = {
  "./Eh12Ez8P3":4, "./Eh24Ez16P3": 8, "./Eh48Ez32P3": 20,   
  "./Eh6Ez4P7":4, "./Eh12Ez8P7":8, "./Eh24Ez16P7":20, "./Eh48Ez32P7":30,  
  "./Eh4Ez3P11":4, "./Eh8Ez6P11":8,  "./Eh16Ez12P11":20,    
}

start_subdir_list = {
  "./Eh12Ez8P3": "run1", "./Eh24Ez16P3": "run1", "./Eh48Ez32P3": "run1",   
  "./Eh6Ez4P7": "run1", "./Eh12Ez8P7": "run1", "./Eh24Ez16P7": "run1", "./Eh48Ez32P7": "run1", 
  "./Eh4Ez3P11": "run1", "./Eh8Ez6P11": "run1", "./Eh16Ez12P11":"run1",     
}


#------------------------------------------------------------------
    
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
        
        job_id, subjob_ind = batch_job_common.pjsub(job_id, subjob_ind, "job_regrid_p_spectra.sh")
        job_id, subjob_ind = batch_job_common.pjsub(job_id, subjob_ind, "job_spectral_analysis.sh")        
        os.chdir(old_path)

for dir in dir_list:
    run_job(dir, num_run_list[dir], start_subdir_list[dir])
