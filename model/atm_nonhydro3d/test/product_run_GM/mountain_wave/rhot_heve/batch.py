import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../../common'))
import batch_job_common

dir_list = [
#* p=3
    "Eh24Ez12P3", "Eh48Ez20P3",  "Eh96Ez36P3", 
#* p=7 
    "Eh12Ez6P7", "Eh24Ez12P7", "Eh48Ez20P7", 
#* p=7 (Reference experiment)
    "Eh96Ez36P7", 
#* p=11
    "Eh8Ez5P11", "Eh16Ez8P11", "Eh32Ez14P11",      
]

def run_job(exp_dir):
    old_path = os.getcwd()

    os.chdir(f"{exp_dir}")
    print(os.getcwd())

    job_id = -1; subjob_ind = -1        
    job_id, subjob_ind = batch_job_common.pjsub(job_id, subjob_ind, "job_run.sh")
    job_id, subjob_ind = batch_job_common.pjsub(job_id, subjob_ind, "job_regrid.sh")    
    job_id, subjob_ind = batch_job_common.pjsub(job_id, subjob_ind, "job_regrid_compari.sh")
        
    os.chdir(old_path)

for dir in dir_list:
    run_job(dir)
