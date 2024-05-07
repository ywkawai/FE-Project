import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../../common'))
import batch_job_common

dir_list = [
# * p=3
    "./Eh10Ez8P3", "./Eh20Ez8P3", "./Eh40Ez8P3", "./Eh80Ez8P3", 
# * p=3 (control)
    "./Eh160Ez8P3", 
# * p=7    
   "./Eh5Ez4P7", "./Eh10Ez4P7", "./Eh20Ez4P7", "./Eh40Ez4P7", 
# * p=7 (control)   
    "./Eh80Ez4P7", 
# * p=11
    "./Eh3Ez3P11", "./Eh6Ez3P11", "./Eh12Ez3P11", "./Eh24Ez3P11", 
# * p=11 (control)    
    "./Eh48Ez3P11", 
]
    
def run_job(exp_dir):
    old_path = os.getcwd()
    
    job_id = -1; subjob_ind = -1            
    for runno in [1, 2]:
        os.chdir(f"{exp_dir}_{runno}")
        print(os.getcwd())
        
        job_id, subjob_ind = batch_job_common.pjsub(job_id, subjob_ind, "job_run.sh")
        job_id, subjob_ind = batch_job_common.pjsub(job_id, subjob_ind, "job_regrid.sh")    
        job_id, subjob_ind = batch_job_common.pjsub(job_id, subjob_ind, "job_regrid_p.sh")    
        
        os.chdir(old_path)

for dir in dir_list:
    run_job(dir)
