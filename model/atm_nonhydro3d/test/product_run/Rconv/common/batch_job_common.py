import re
import subprocess
import math

FUGAKU_LANG_ENV="tcsds-1.2.41"
SPACK_LIB_PATH="/vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.10.0/snappy-1.1.10-j7x2zk63lmbm36rc724j3ik66antzv5k/lib:/vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.10.0/zlib-ng-2.1.4-m4nsqgrityww4eqljits3svzi4wynus3/lib:/vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.10.0/zstd-1.5.5-o7dcxyae4kjgruht44yuloa27ddisxfe/lib:/vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.10.0/c-blosc-1.21.5-66ynsnfbl2evrc6id523jcamtzvu2ugx/lib64:/vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.10.0/hdf5-1.14.3-yhazdvld6vknkhmbcqrbl34ifsac2hao/lib:/vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.10.0/libaec-1.0.6-5mds5a4x6u7m2ysxqxxgm2l6tpax77ds/lib64:/vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.10.0/netcdf-c-4.9.2-cczsjh6lmalmjqhr72yatvxzmc3iwdl2/lib:/vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.10.0/netcdf-fortran-4.6.1-kjm7jl5naxulm7neknsgxnsbkf7iet5j/lib:/vol0004/apps/oss/spack-v0.21/opt/spack/linux-rhel8-a64fx/fj-4.10.0/parallel-netcdf-1.12.3-zs75olgkg547zcyxaoubsyevmvtkby5o/lib"
LLIO_GFSCACHE="/vol0005:/vol0004"
GROUPID="ra000005"

#------------------------------------------------------------------

def get_subdir_list(num_run, start_subdir):
  subdir_list_cand = []
  for runno in range(1,num_run+1):
    subdir_list_cand.append(f"run{runno}")

  add_flag = False; subdir_list = []
  for subdir in subdir_list_cand:
    if add_flag == False and subdir == start_subdir:
        add_flag = True
    if add_flag:        
        subdir_list.append(subdir)
  return subdir_list

def extract_job_id(ret_stdout):
    ret = re.match(r".+Job ([0-9]+)_[0-9]+ submitted.", ret_stdout)
    return ret.group(1)
    
def pjsub(job_id, subjob_ind, job_script):
    if subjob_ind < 0:
        res = subprocess.run(["pjsub", f"--step {job_script}"], capture_output=True, text=True )    
        print(f"{job_script}: {res.stdout}")
        job_id = extract_job_id(res.stdout)
    else:
        res = subprocess.run(["pjsub", "--step", f"--sparam jid={job_id},sd=ec!=0:after:{subjob_ind}", job_script], capture_output=True, text=True )    
        print(f"{job_script}: {res.stdout}")
        
    return job_id, subjob_ind + 1

def get_job_header( job_name, nprc, elapse_time, node="" 
    ):
  node_num = math.ceil(nprc/4)
  if node_num > 384:
    rscgrp = "large"
  elif node_num == 384:
    rscgrp = "large"
    node_num = 385    
  else:
    rscgrp = "small"
  
  if len(node) > 0:
    node_num = node
    
  jobshell_s = f"""################################################################################
#
# for Fugaku
#
################################################################################
#PJM --rsc-list "rscunit=rscunit_ft01"
#PJM --name "{job_name}"
#PJM -x PJM_LLIO_GFSCACHE={LLIO_GFSCACHE}
#PJM --rsc-list "rscgrp={rscgrp}"
#PJM --rsc-list "node={node_num}"
#PJM --rsc-list "elapse={elapse_time}"
#PJM -g {GROUPID}
#PJM --mpi "max-proc-per-node=4"
#PJM -S

module purge
module load lang/{FUGAKU_LANG_ENV}

llio_transfer /home/apps/oss/scale/llio.list
export LD_LIBRARY_PATH=/lib64:/opt/FJSVxtclanga/tcsds-mpi-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:`cat /home/apps/oss/scale/llio.list | sed 's:\\(.*/lib\\(\\|64\\)\\)/.*:\\1:' | uniq | sed -z 's/\\n/:/g'`

#export XOS_MMM_L_ARENA_FREE=1
export FORT90L="-Wl,-T"
#export OMPI_MCA_plm_ple_memory_allocation_policy=bind_local
export PLE_MPI_STD_EMPTYFILE="off"
export PARALLEL=12
export OMP_NUM_THREADS=12
#export fu11bf=1
  """
  return jobshell_s  


#---
def mkconf_regrid( conf_path,
                nprcx, nprcy, nex, ney, nez, porder, 
                regrid_nprcx, regrid_nprcy, 
                regrid_nex, regrid_ney, regrid_nez, 
                regrid_porder, 
                in_basename, vars_list, out_basename ): 
    conf_run_s = f"""#--- Configuration file for a test case of RB convection  -------
&PARAM_IO
 IO_LOG_BASENAME = "regrid_LOG"
! IO_LOG_ALLNODE  = .true., 
/
&PARAM_CONST
  CONST_OHM = 0.0D0, 
/
&PARAM_REGRID_MESH
 in_MeshType  = "STRUCTURED3D", 
 out_MeshType = "STRUCTURED3D", 
/
&PARAM_REGRID_INTERP_FIELD
  !- input --------------------
  in_basename="{in_basename}",      
  vars = {vars_list}, ! "U", "V", "W", "DDENS", "PT", 
  !out_tinterval = 5,
/
&PARAM_REGRID_FILE
  !-- output ----------------
  out_basename="{out_basename}", 
  out_UniformGrid=.true., 
/
&PARAM_REGRID_INMESH3D_STRUCTURED
  dom_xmin         = -1.6D3, 
  dom_xmax         =  1.6D3, 
  isPeriodicX      = .true.,
  dom_ymin         = -1.6D3,  
  dom_ymax         =  1.6D3,  
  isPeriodicY      = .true.,  
  dom_zmin         = 0.0D0, 
  dom_zmax         = 1.6D3,   
  NprcX            = {nprcx}, 
  NeX              = {nex},
  NprcY            = {nprcy},   
  NeY              = {ney},
  NeGZ             = {nez},
  PolyOrder_h      = {porder},
  PolyOrder_v      = {porder},    
/
&PARAM_REGRID_OUTMESH3D_STRUCTURED
  dom_xmin         = -1.6D3, 
  dom_xmax         =  1.6D3, 
  isPeriodicX      = .true.,
  dom_ymin         = -1.6D3,  
  dom_ymax         =  1.6D3,  
  isPeriodicY      = .true.,  
  dom_zmin         = 0.0D0, 
  dom_zmax         = 1.6D3, 
  NprcX       = {regrid_nprcx},       
  NeX         = {regrid_nex},           
  NprcY       = {regrid_nprcy}, 
  NeY         = {regrid_ney},    
  NeGZ        = {regrid_nez}, 
  PolyOrder_h = {regrid_porder}, 
  PolyOrder_v = {regrid_porder}, 
/
    """
    
    with open(conf_path, 'w') as f:
        f.write(conf_run_s)

#--