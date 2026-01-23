import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../common'))
import batch_job_common

dir_list = [
  "Dx12.5m_P3",   
  "Dx12.5m_P7", 
  "Dx13m_P11",   
]
num_run_list = {
  "Dx12.5m_P3": 24, 
  "Dx12.5m_P7": 24,
  "Dx13m_P11": 24, 
}
start_subdir_list = {
  "Dx12.5m_P3": "run8", 
  "Dx12.5m_P7": "run4", 
  "Dx13m_P11": "run8",  
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
  # regrid_job(dir, num_run_list[dir], start_subdir_list[dir])
  analysis_job(dir, num_run_list[dir], start_subdir_list[dir])
