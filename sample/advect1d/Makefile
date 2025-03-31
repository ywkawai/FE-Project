################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../..)
SAMPLE_DIR  = ../

BINNAME = test_advect1d
OBJS    = mod_advect1d_numerror.o

RUNCONF         = test.conf
TPROC           = 1
OMP_NUM_THREADS = 1

#------------------------------------------------------
include $(SAMPLE_DIR)/Makefile.common

vis:
	bash ./visualize/visualize.sh

clean_analysis:
	rm -f analysis

$(BUILD_DIR)/test_advect1d.o : $(BUILD_DIR)/mod_advect1d_numerror.o