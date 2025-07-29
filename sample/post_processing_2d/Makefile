################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../..)
SAMPLE_DIR  = ../

BINNAME = test_post_processing_2d
OBJS    = 

RUNCONF         = test.conf
TPROC           = 1
OMP_NUM_THREADS = 1

#------------------------------------------------------
include $(SAMPLE_DIR)/Makefile.common

vis:
	bash ./visualize/visualize.sh

clean_analysis:
	rm -f analysis
