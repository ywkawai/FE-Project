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
#INITCONF    = init.conf
RUNCONF     = run.conf
TPROC       = 1

# required data (parameters,distributed files)
DATPARAM    =
DATDISTS    =

#--
BINNAME = scale-dg_analytic_metric
FFLAG_TESTCASE = "-DSCALE_PRODUCT_RUN_RM_MOUNTAIN_WAVE -DSCALE_PRODUCT_RUN_RM_MOUNTAIN_WAVE_ANALYTIC_METRIC"

# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
