################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../..)
SAMPLE_DIR  = ../

BINNAME = test_advect3d
OBJS = 

#-------------------------------------------------------
include $(SAMPLE_DIR)/Makefile.common

run:
	OMP_NUM_THREADS=2 mpirun -n 4 ./$(BINNAME) test.conf

vis:
	bash ./visualize/visualize.sh
