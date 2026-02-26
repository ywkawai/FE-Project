################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../..)
SAMPLE_DIR  = ../

BINNAME = test_advect3d
OBJS = mod_advect3d_kernel.o    \
			 mod_advect3d_numerror.o

RUNCONF         = test.conf
TPROC           = 4
OMP_NUM_THREADS = 2

#-------------------------------------------------------
include $(SAMPLE_DIR)/Makefile.common

vis:
	bash ./visualize/visualize.sh

$(BUILD_DIR)/test_advect3d.o : $(BUILD_DIR)/mod_advect3d_kernel.o \
															 $(BUILD_DIR)/mod_advect3d_numerror.o