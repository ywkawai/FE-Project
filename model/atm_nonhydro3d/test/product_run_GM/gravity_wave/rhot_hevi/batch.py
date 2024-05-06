import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../../common'))
import batch_job_common

dir_list = [
# * p=1
    "Eh16Ez12P1", "Eh32Ez24P1", "Eh64Ez48P1", "Eh128Ez96P1", 
# * p=3    
    "Eh8Ez6P3", "Eh16Ez12P3", "Eh32Ez24P3", "Eh64Ez48P3", 
# * p=7
    "Eh4Ez3P7", "Eh8Ez6P7",  "Eh16Ez12P7", "Eh32Ez24P7", "Eh64Ez48P7", 
    "Eh4Ez3P7_dtx0.5", "Eh8Ez6P7_dtx0.5",  "Eh16Ez12P7_dtx0.5", "Eh32Ez24P7_dtx0.5", "Eh64Ez48P7_dtx0.5", 
    "Eh4Ez3P7_dtx0.25", "Eh8Ez6P7_dtx0.25",  "Eh16Ez12P7_dtx0.25", "Eh32Ez24P7_dtx0.25", "Eh64Ez48P7_dtx0.25",     
# * p=11
    "Eh3Ez2P11", "Eh6Ez4P11", "Eh12Ez8P11", 
    "Eh3Ez2P11_dtx0.5", "Eh6Ez4P11_dtx0.5", "Eh12Ez8P11_dtx0.5",     
    "Eh3Ez2P11_dtx0.25", "Eh6Ez4P11_dtx0.25", "Eh12Ez8P11_dtx0.25",        
]
    
def run_job(exp_dir):
    old_path = os.getcwd()
    
    os.chdir(f"{exp_dir}")
    print(os.getcwd())
        
    job_id = -1; subjob_ind = -1        
    job_id, subjob_ind = batch_job_common.pjsub(job_id, subjob_ind, "job_run.sh")
    job_id, subjob_ind = batch_job_common.pjsub(job_id, subjob_ind, "job_regrid.sh")    
        
    os.chdir(old_path)

for dir in dir_list:
    run_job(dir)
