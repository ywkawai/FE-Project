################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../..)
SAMPLE_DIR  = ../

BINNAME = test_advect1d
OBJS = 

#-------------------------------------------------------
include $(SAMPLE_DIR)/Makefile.common

run:
	OMP_NUM_THREADS=1 mpirun -n 1 ./$(BINNAME) test.conf

vis:
	bash ./visualize/visualize.sh
