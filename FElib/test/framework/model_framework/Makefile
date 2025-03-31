################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../..)
FELIB_TEST_DIR = ../../

BINNAME = test_model_framework
OBJS = mod_atmos_component.o \
       mod_atmos_vars.o         \
			 mod_atmos_mesh.o   \
			 mod_atmos_dyn.o

RUNCONF         = test.conf
TPROC           = 1
OMP_NUM_THREADS = 1

#------------------------------------------------------
include $(FELIB_TEST_DIR)/Makefile.common

$(BUILD_DIR)/test_model_framework.o: $(BUILD_DIR)/mod_atmos_component.o
$(BUILD_DIR)/mod_atmos_component.o:  \
	$(BUILD_DIR)/mod_atmos_vars.o      \
	$(BUILD_DIR)/mod_atmos_dyn.o
$(BUILD_DIR)/mod_atmos_vars.o: $(BUILD_DIR)/mod_atmos_mesh.o