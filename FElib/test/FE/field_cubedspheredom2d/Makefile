################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../..)
BUILD_DIR   = ./.libs
SYSDEP_DIR  = $(TOPDIR)/sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_FE_SYS)
include $(TOPDIR)/Mkinclude

BINNAME = test_field_cubedspheredom2d

LIBS = $(LIBDIR)/libScaleFECore.a

OBJS =

all:
	$(MAKE) envlog
	$(MAKE) makedir
	$(MAKE) makebin
	$(MAKE) run

makedir:
	mkdir -p $(BUILD_DIR)

makebin: $(BINNAME)
	@echo "Complete making."

run:
	mpirun -n 1 ./$(BINNAME) test.conf

vis:
	bash ./visualize/visualize.sh

$(BINNAME): $(BUILD_DIR)/$(BINNAME).o $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(LIBS) $(CONTRIB_LIBS)  $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)

$(BUILD_DIR)/$(BINNAME).o: $(BINNAME).f90 $(patsubst %,$(BUILD_DIR)/%,$(OBJS))

distclean: clean
	rm -f $(BINNAME)

clean:
	rm -rf $(BUILD_DIR)
	rm -f *.nc

.SUFFIXES:
.SUFFIXES: .o .f90 .mod

%.mod: %.f90
	$(MAKE) $(patsubst %.f90,%.o,$<)

$(BUILD_DIR)/%.o: %.f90
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(MODDIR) $(CONTRIB_INCLUDE) -I$(SCALEFELIBDIR)/include  $(SCALE_NETCDF_INCLUDE) $(MODDIROPT) $(BUILD_DIR) -o $@ -c $<

.PHONY : clean distclean allclean

include $(TOPDIR)/utils/make/Make_environments
