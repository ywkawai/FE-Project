################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../..)
SAMPLE_DIR  = ../

BINNAME = test_euler3d_hevi
OBJS = 

RUNCONF         = test.conf
TPROC           = 1
OMP_NUM_THREADS = 2

#-------------------------------------------------------
include $(SAMPLE_DIR)/Makefile.common

vis:
	bash ./visualize/visualize.sh
