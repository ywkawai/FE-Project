import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../common'))
import batch_job_common

dir_list = [
  "Dx3.1m_P7", 
]
num_run_list = {
  "Dx3.1m_P7": 100,
}
start_subdir_list = {
  "Dx3.1m_P7": "run13", 
}

#-------------   
def analysis_job(exp_dir, num_run, start_subdir):
  old_path = os.getcwd()

  job_id = -1; subjob_ind = -1
  for subdir in batch_job_common.get_subdir_list(num_run, start_subdir):
    os.chdir(f"{exp_dir}/{subdir}")
    print(os.getcwd())
    job_id, subjob_ind = batch_job_common.pjsub(job_id, subjob_ind, "job_analysis.sh")        
    os.chdir(old_path)

def regrid_job(exp_dir, num_run, start_subdir):
  old_path = os.getcwd()

  job_id = -1; subjob_ind = -1
  for subdir in batch_job_common.get_subdir_list(num_run, start_subdir):
    os.chdir(f"{exp_dir}/{subdir}")
    print(os.getcwd())
    job_id, subjob_ind = batch_job_common.pjsub(job_id, subjob_ind, "job_regrid.sh")        
    os.chdir(old_path)

for dir in dir_list:
#  regrid_job(dir, num_run_list[dir], start_subdir_list[dir])
 analysis_job(dir, num_run_list[dir], start_subdir_list[dir])
