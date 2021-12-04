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
TPROC       = 6

# required data (parameters,distributed files)
DATPARAM    =
DATDISTS    =

# user-defined targets
USER_target = make_expdir

# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common

# Make dirs
make_expdir:
	mkdir -p run_{spinup1,spinup2,1,2}/outdata_p

