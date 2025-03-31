################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../..)
SAMPLE_DIR  = ../

BINNAME = test_advect2d
OBJS = mod_advect2d_numerror.o

RUNCONF         = test.conf
TPROC           = 1
OMP_NUM_THREADS = 4

#-------------------------------------------------------
include $(SAMPLE_DIR)/Makefile.common

vis:
	bash ./visualize/visualize.sh

clean_analysis:
	rm -f analysis

$(BUILD_DIR)/test_advect2d.o : $(BUILD_DIR)/mod_advect2d_numerror.o