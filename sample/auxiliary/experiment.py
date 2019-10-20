#!/usr/bin/env python
from abc import ABCMeta, abstractmethod
import os
import textwrap
import shutil

class Experiment:
  def __init__(self, exp_topdir, name, org_pe_path):
    self._expcase_list = {}
    self._name = name
    self._topdir = exp_topdir
    self._org_pe_path = os.path.abspath(org_pe_path)

    pe_basename = os.path.basename(self._org_pe_path)
    self._pe_path = os.path.abspath(self.topdir + "/" + pe_basename)

  def add_expcase(self, expcase):
    self._expcase_list[expcase.name] = expcase
  
  @property
  def name(self):
    return self._name
  
  @property
  def topdir(self):
    return self._topdir

  @property
  def expcase_list(self): 
    return self._expcase_list

  @property
  def pe_path(self):
    return self._pe_path

  def print_info(self):
    print('===== Experiment Information ================')
    print('name: '+self.name)
    print('top dir: '+self.topdir)  
    print('original pe: '+self._org_pe_path)      
    print(str(len(self._expcase_list)) +' cases described below is contained..')
    for k, v in self.expcase_list.items():
      v.print_info()

  def prepair(self):
    shutil.copyfile(self._org_pe_path, self._pe_path)
    os.chmod(self.pe_path, 0o777)

    print("prepair & generate configure files and job scripts ..")
    for case in self.expcase_list.values():
      print("case= " + case.name + " path=" + case.path + "...")
      case.prepair()
      case.generate_configfile()
      case.generate_jobscript(self._topdir + case.name+".sh")      

class ExpCase(metaclass=ABCMeta):
  def __init__(self, exp, path_from_exptopdir, name, confname, 
               duration_sec, dt_sec, 
               pename, nproc=1, nthread=1):
     self.exp = exp    
     self._name = name
     self._path = exp.topdir + "/" + path_from_exptopdir
     self._conf_name = confname
     self._pename = pename
     self._nproc = nproc
     self._nthread = nthread
     self._dt_sec = dt_sec
     self._duration_sec = duration_sec

  @property
  def path(self):
    return self._path

  @property
  def name(self):
    return self._name
  
  @property
  def conf_name(self):
    return self._conf_name

  def print_info(self):
    print('++ ExpCase Information ++++++')
    print('name: '+self.name)
    print('path: '+self.path)

  def prepair(self):
    os.makedirs(self.path, exist_ok='True')
    link_dist=self.path + "/" + self._pename
    if os.path.exists(link_dist): os.unlink(link_dist) 
    os.symlink(self.exp.pe_path, link_dist)    

  @abstractmethod
  def generate_configfile(self):
    pass

  def generate_jobscript(self, jobshell_path):

    conf = textwrap.dedent('''\
    #!/bin/bash
    #********************************************************************************************
    #PBS -q s
    #PBS -l nodes=1:ppn=1
    #PBS -N {expname}
    #PBS -m ae
    #PBS -l walltime=24:00:00
    ##################################################

    MPIRUN=mpiexec.hydra

    #---------
    export I_MPI_PIN_DOMAIN=omp
    export KMP_AFFINITY=compact
    source /etc/profile.d/modules.sh
    module unload mpt/2.12
    module load intelmpi/5.1.2.150

    export OMP_NUM_THREADS={nthread}
    PBS_O_WORKDIR={expdir}
    ## End of setting *******************************************************************************

    cd $PBS_O_WORKDIR
    $MPIRUN -n {nproc} {pename} {confname}    
    ''').format( 
      expname=self.name, pename=self._pename, 
      expdir=os.path.abspath(self.path), confname=self.conf_name, 
      nproc=self._nproc, nthread=self._nthread )

    with open(jobshell_path, mode='w') as f:
      f.write( conf ) 




