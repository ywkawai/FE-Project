import os
import re
import subprocess

dir_list = [
   "./Eh10Ez6P3",  "./Eh20Ez6P3", "./Eh40Ez6P3",
#   "./Eh5Ez3P7",  "./Eh10Ez3P7", "./Eh20Ez3P7", 
#    "./Eh40Ez3P7", 
#    "./Eh4Ez3P11", "./Eh8Ez3P11", "./Eh16Ez3P11", 
]
def extract_job_id(ret_stdout):
    ret = re.match(r".+Job ([0-9]+)_[0-9]+ submitted.", ret_stdout)
    return ret.group(1)
    
def run_job(exp_dir, runno):
    old_path = os.getcwd()
    
    job_id = -1
    subjob_ind = -1
    
    os.chdir(f"{exp_dir}_{runno}")
    print(os.getcwd())
        
    res_regrid = subprocess.run(["pjsub", "--step job_regrid_analysis.sh"], capture_output=True, text=True )    
    print(f"regrid_analysis.sh: {res_regrid.stdout}")
    job_id = extract_job_id(res_regrid.stdout)
    subjob_ind = subjob_ind+1

    res_analysis= subprocess.run(["pjsub", "--step", f"--sparam jid={job_id},sd=ec!=0:after:{subjob_ind}", "job_analysis.sh"], capture_output=True, text=True )        
    print(f"analysis.sh: {res_analysis.stdout}")
    subjob_ind = subjob_ind+1
        
    os.chdir(old_path)

for dir in dir_list:
    run_job(dir, 1)
    run_job(dir, 2)
