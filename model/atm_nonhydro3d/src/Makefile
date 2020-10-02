################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../..)
BUILD_DIR   = ./.libs
SYSDEP_DIR  = $(TOPDIR)/sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_FE_SYS)
include $(TOPDIR)/Mkinclude

BINNAME = scale-dg

BINS = $(BINDIR)/$(BINNAME)

VERSION = $(shell git rev-parse --short HEAD 2> /dev/null)
ifeq ($(VERSION),)
  VERSION  = $(shell cat VERSION)
else
  VERSION := $(VERSION)
endif

VPATH = \
	$(BUILD_DIR): \
	admin:        \
	atmos:        \
	user:         

OBJS = \
	mod_dg_driver.o              \
	mod_atmos_component.o        \
	mod_atmos_dyn.o              \
	mod_atmos_dyn_vars.o         \
	mod_atmos_mesh.o             \
	mod_atmos_vars.o             \
	mod_atmos_bnd.o              \
	mod_user.o                   \
	mod_exp.o

LIBS = $(LIBDIR)/libScaleFECore.a

SAMPLE_AUX_DIR = $(TOPDIR)/sample/auxiliary
SAMPLE_AUX_INCLUDE = $(SAMPLE_AUX_DIR)/.libs
SAMPLE_AUX_LIB = $(SAMPLE_AUX_DIR)/libScaleFESampleAux.a


all:
	$(MAKE) envlog
	$(MAKE) makedir
	$(MAKE) make_sample_aux
	@echo;echo "Entering atm_nonhydro3d ...";echo "Current version is " $(VERSION)
	$(MAKE) makebin
	@echo "Complete making."

makedir:
	mkdir -p $(BINDIR)
	mkdir -p $(BUILD_DIR)

makebin: $(BINS)

run:
	mpirun -n 1 ./$(BINNAME) $(CONF_NAME)

vis:
	bash ./visualize/visualize.sh

make_sample_aux:
	make -C $(SAMPLE_AUX_DIR)

$(BINDIR)/$(BINNAME): $(BUILD_DIR)/$(BINNAME).o $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(SAMPLE_AUX_LIB) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(LIBS) $(CONTRIB_LIBS) $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)

$(BUILD_DIR)/$(BINNAME).o: $(BINNAME).F90 $(patsubst %,$(BUILD_DIR)/%,$(OBJS))

distclean: clean
	rm -f $(BINNAME)

clean:
	rm -rf $(BUILD_DIR)

dependency:
	./makedepend .

.SUFFIXES:
.SUFFIXES: .o .F90 .mod

%.mod: %.F90
	$(MAKE) $(patsubst %.F90,%.o,$<)

$(BUILD_DIR)/%.o: %.F90
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(MODDIR) -I../include $(CONTRIB_INCLUDE) -I$(SCALEFELIBDIR)/include -I$(SAMPLE_AUX_INCLUDE) $(SCALE_NETCDF_INCLUDE) $(MODDIROPT) $(BUILD_DIR) -o $@ -c $<

.PHONY : clean distclean allclean

include $(TOPDIR)/utils/make/Make_environments
-include depend

$(BUILD_DIR)/mod_user.o: $(BUILD_DIR)/mod_exp.o
