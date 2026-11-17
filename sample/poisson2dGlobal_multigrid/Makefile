################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../..)
SAMPLE_DIR  = ../

BINNAME = test_poisson2dGlobal_multigrid
OBJS    = mod_poisson2dGlobal_smoother.o \
		  mod_poisson2dGlobal_numerror.o

RUNCONF         = test.conf
TPROC           = 6
OMP_NUM_THREADS = 1

#------------------------------------------------------
include $(SAMPLE_DIR)/Makefile.common

vis:
	bash ./visualize/visualize.sh

clean_analysis:
	rm -f analysis

$(BUILD_DIR)/test_poisson2dGlobal_multigrid.o : \
	$(BUILD_DIR)/mod_poisson2dGlobal_smoother.o \
	$(BUILD_DIR)/mod_poisson2dGlobal_numerror.o
