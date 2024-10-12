import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../../common'))
import batch_job_common

dir_list = [
   "Eh16Ez12P1", "Eh32Ez12P1", "Eh64Ez12P1", "Eh128Ez12P1",         
   "Eh8Ez6P3", "Eh16Ez6P3", "Eh32Ez6P3",  "Eh64Ez6P3", 
#    "Eh4Ez3P7", "Eh8Ez3P7", "Eh16Ez3P7", 
   "Eh32Ez3P7", 
#    "Eh4Ez3P7_check", "Eh8Ez3P7_check", "Eh16Ez3P7_check", "Eh32Ez3P7_check",     
   "Eh2Ez3P11", "Eh4Ez3P11", "Eh8Ez3P11", "Eh16Ez3P11", 
]

def run_job(exp_dir):
    old_path = os.getcwd()
    
    os.chdir(f"{exp_dir}")
    print(os.getcwd())

    job_id = -1; subjob_ind = -1        
    job_id, subjob_ind = batch_job_common.pjsub(job_id, subjob_ind, "job_regrid_analysis.sh")
    job_id, subjob_ind = batch_job_common.pjsub(job_id, subjob_ind, "job_analysis.sh")    
        
    os.chdir(old_path)

for dir in dir_list:
    run_job(dir)