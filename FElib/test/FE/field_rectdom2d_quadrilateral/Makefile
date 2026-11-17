################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../..)
FELIB_TEST_DIR = ../../

BINNAME = test_field_rectdom2d_quadrilateral
OBJS    = mod_test_meshfield_util.o

RUNCONF         = test.conf
TPROC           = 1
OMP_NUM_THREADS = 1

#------------------------------------------------------
include $(FELIB_TEST_DIR)/Makefile.common