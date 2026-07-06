################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../..)
SAMPLE_DIR  = ../

BINNAME = test_poisson3d_multigrid
OBJS    = mod_poisson3d_smoother.o

RUNCONF         = test.conf
TPROC           = 4
OMP_NUM_THREADS = 1

#------------------------------------------------------
include $(SAMPLE_DIR)/Makefile.common

vis:
	bash ./visualize/visualize.sh

clean_analysis:
	rm -f analysis

$(BUILD_DIR)/test_multigrid3d.o : \
	$(BUILD_DIR)/mod_poisson3d_smoother.o

