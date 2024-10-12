import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../common'))
import batch_job_common

dir_list = [
    "Eh8Ez6P3_MFAlpha1D-3", "Eh16Ez6P3_MFAlpha1D-3", "Eh32Ez6P3_MFAlpha1D-3", "Eh64Ez6P3_MFAlpha1D-3", 
    "Eh4Ez3P7_MFAlpha1D-3", "Eh8Ez3P7_MFAlpha1D-3", "Eh16Ez3P7_MFAlpha1D-3", "Eh32Ez3P7_MFAlpha1D-3", 
    "Eh2Ez3P11_MFAlpha1D-3", "Eh4Ez3P11_MFAlpha1D-3", "Eh8Ez3P11_MFAlpha1D-3", "Eh16Ez3P11_MFAlpha1D-3", 
#--
    "Eh8Ez6P3_MFAlpha1D-2", "Eh16Ez6P3_MFAlpha1D-2", "Eh32Ez6P3_MFAlpha1D-2", "Eh64Ez6P3_MFAlpha1D-2", 
    "Eh4Ez3P7_MFAlpha1D-2", "Eh8Ez3P7_MFAlpha1D-2", "Eh16Ez3P7_MFAlpha1D-2", "Eh32Ez3P7_MFAlpha1D-2", 
    "Eh2Ez3P11_MFAlpha1D-2", "Eh4Ez3P11_MFAlpha1D-2", "Eh8Ez3P11_MFAlpha1D-2", "Eh16Ez3P11_MFAlpha1D-2", 
#---    
    "Eh8Ez6P3_MFAlpha1D-1", "Eh16Ez6P3_MFAlpha1D-1", "Eh32Ez6P3_MFAlpha1D-1", "Eh64Ez6P3_MFAlpha1D-1", 
    "Eh4Ez3P7_MFAlpha1D-1", "Eh8Ez3P7_MFAlpha1D-1", "Eh16Ez3P7_MFAlpha1D-1", "Eh32Ez3P7_MFAlpha1D-1", 
    "Eh2Ez3P11_MFAlpha1D-1", "Eh4Ez3P11_MFAlpha1D-1", "Eh8Ez3P11_MFAlpha1D-1", "Eh16Ez3P11_MFAlpha1D-1", 
#---            
    "Eh8Ez6P3_MFAlpha1D0", "Eh16Ez6P3_MFAlpha1D0", "Eh32Ez6P3_MFAlpha1D0", "Eh64Ez6P3_MFAlpha1D0", 
    "Eh4Ez3P7_MFAlpha1D0", "Eh8Ez3P7_MFAlpha1D0", "Eh16Ez3P7_MFAlpha1D0", "Eh32Ez3P7_MFAlpha1D0", 
    "Eh2Ez3P11_MFAlpha1D0", "Eh4Ez3P11_MFAlpha1D0", "Eh8Ez3P11_MFAlpha1D0", "Eh16Ez3P11_MFAlpha1D0", 
#---        
    "Eh8Ez6P3_MFAlpha1D1", "Eh16Ez6P3_MFAlpha1D1", "Eh32Ez6P3_MFAlpha1D1", "Eh64Ez6P3_MFAlpha1D1", 
    "Eh4Ez3P7_MFAlpha1D1", "Eh8Ez3P7_MFAlpha1D1", "Eh16Ez3P7_MFAlpha1D1", "Eh32Ez3P7_MFAlpha1D1", 
    "Eh2Ez3P11_MFAlpha1D1", "Eh4Ez3P11_MFAlpha1D1", "Eh8Ez3P11_MFAlpha1D1", "Eh16Ez3P11_MFAlpha1D1", 
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