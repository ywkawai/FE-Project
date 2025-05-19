################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../../..)
TESTDIR     = ../../..


# user-defined source files
CODE_DIR    = ..
ORG_SRCS    = mod_user.F90

# parameters for run
INITCONF    = init.conf
RUNCONF     = run.conf
TPROC       = 4
POSTCONF    = regrid_topo.conf,regrid.conf
TPROC_POST  = 4,4
OMP_NUM_THREADS = 1

# required data (parameters,distributed files)
DATPARAM    =
DATDISTS    =

# user-defined targets
USER_target = make_expdir

# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common

# Make dirs
make_expdir:
	mkdir -p reg_outdata
