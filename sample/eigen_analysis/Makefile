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

#-------------------------------------------------------
include $(SAMPLE_DIR)/Makefile.common

run:
	export OMP_NUM_THREADS=1
	mkdir -p data
	mpirun -n 1 ./$(BINNAME) test.conf

vis:
	python ./visualize/dispersion_dissipation_graph.py

clean_analysis: 
	rm -rf ./data
	rm -rf ./visualize/*.png
