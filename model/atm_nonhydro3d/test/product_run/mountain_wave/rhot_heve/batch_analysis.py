import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../../common'))
import batch_job_common

dir_list = [
    "Eh96Ez15P3_MFoff", "Eh192Ez30P3_MFoff", "Eh384Ez60P3_MFoff", "Eh768Ez120P3_MFoff", "Eh1536Ez240P3_MFoff", 
    "Eh96Ez15P3", "Eh192Ez30P3", "Eh384Ez60P3", "Eh768Ez120P3", "Eh1536Ez240P3",     
    "Eh48Ez6P7", "Eh96Ez12P7", "Eh192Ez20P7", 
#    "Eh384Ez60P7",  
    "Eh32Ez5P11", "Eh64Ez10P11", "Eh128Ez20P11", #"Eh256Ez40P11", 
]

def run_job(exp_dir):
    old_path = os.getcwd()

    os.chdir(f"{exp_dir}")
    print(os.getcwd())

    job_id = -1; subjob_ind = -1        
    job_id, subjob_ind = batch_job_common.pjsub(job_id, subjob_ind, "job_analysis.sh")
        
    os.chdir(old_path)

for dir in dir_list:
    run_job(dir)
