################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../..


# user-defined source files
CODE_DIR    = .
ORG_SRCS    = mod_user.F90

# parameters for run
INITCONF    = init.conf
RUNCONF     = run.conf
TPROC       = 24
OMP_NUM_THREADS = 1
# parameters for post-processing
POSTCONF    = regrid.conf
TPROC_POST  = 4

# user-defined targets
USER_target = make_outdir


# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common

# Make dirs
make_outdir:
	mkdir -p outdata
