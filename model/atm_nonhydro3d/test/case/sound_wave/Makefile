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
INITCONF        = init.conf
RUNCONF         = run.conf
TPROC           = 1
OMP_NUM_THREADS = 4


# required data (parameters,distributed files)
DATPARAM    =
DATDISTS    =


POSTFIX="_test"
# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
