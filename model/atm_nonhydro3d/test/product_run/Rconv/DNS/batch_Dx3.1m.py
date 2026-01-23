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

def run_job(exp_dir, num_run, start_subdir):
    old_path = os.getcwd()

    job_id = -1; subjob_ind = -1      
    subdir_list_cand = []
    for runno in range(1,num_run+1):
      subdir_list_cand.append(f"run{runno}")

    add_flag = False; subdir_list = []
    for subdir in subdir_list_cand:
      if add_flag == False and subdir == start_subdir:
          add_flag = True
      if add_flag:        
          subdir_list.append(subdir)
    
    for subdir in subdir_list:
      os.chdir(f"{exp_dir}/{subdir}")
      print(os.getcwd())
        
      job_id, subjob_ind = batch_job_common.pjsub(job_id, subjob_ind, "job_run.sh")
      os.chdir(old_path)

for dir in dir_list:
    run_job(dir, num_run_list[dir], start_subdir_list[dir])
