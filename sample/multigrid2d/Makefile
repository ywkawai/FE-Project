################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../..)
SAMPLE_DIR  = ../

BINNAME = test_multigrid2d
OBJS    = mod_poisson2d_mg.o \
		  mod_poisson2d_smoother.o

RUNCONF         = test.conf
TPROC           = 4
OMP_NUM_THREADS = 1

#------------------------------------------------------
include $(SAMPLE_DIR)/Makefile.common

vis:
	bash ./visualize/visualize.sh

clean_analysis:
	rm -f analysis

$(BUILD_DIR)/test_multigrid2d.o : $(BUILD_DIR)/mod_poisson2d_mg.o
$(BUILD_DIR)/mod_poisson2d_mg.o : $(BUILD_DIR)/mod_poisson2d_smoother.o
