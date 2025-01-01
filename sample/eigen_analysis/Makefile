################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../..)
SAMPLE_DIR  = ../

BINNAME = test_eigen_analysis
OBJS    = 

RUNCONF         = test.conf
TPROC           = 1
OMP_NUM_THREADS = 1

#-------------------------------------------------------
run: prep_run
prep_run:
	mkdir -p data

include $(SAMPLE_DIR)/Makefile.common

vis:
	python ./visualize/dispersion_dissipation_graph.py

clean_analysis: 
	rm -rf ./data
	rm -rf ./visualize/*.png
