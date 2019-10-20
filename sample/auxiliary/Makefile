################################################################################
#
# Makefile for auxiliary
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../..)
BUILD_DIR   = ./.libs
SYSDEP_DIR  = $(TOPDIR)/sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_FE_SYS)
include $(TOPDIR)/Mkinclude

SAMPLE_AUX_LIBS = libScaleFESampleAux.a
SAMPLE_AUX_OBJS = \
	mod_operator_fvm.o \
	mod_output_fvm.o   \
	mod_fieldutil.o

all:
	$(MAKE) envlog
	$(MAKE) makedir
	$(MAKE) makelib

makedir:
	mkdir -p $(BUILD_DIR)

makelib: $(SAMPLE_AUX_LIBS)
	@echo "Complete making."

$(SAMPLE_AUX_LIBS): $(patsubst %,$(BUILD_DIR)/%,$(SAMPLE_AUX_OBJS))
	ar -r $@ $?
	ranlib $@

distclean: clean
	rm -f $(SAMPLE_AUX_LIBS)

clean:
	rm -rf $(BUILD_DIR)

.SUFFIXES:
.SUFFIXES: .o .f90 .mod


%.f90 : %.f90.erb
	erb $< > $@

%.mod: %.f90
	$(MAKE) $(patsubst %.f90,%.o,$<)

$(BUILD_DIR)/%.o: %.f90
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(MODDIR) $(CONTRIB_INCLUDE) -I$(SCALEFELIBDIR)/include  $(SCALE_NETCDF_INCLUDE) $(MODDIROPT) $(BUILD_DIR) -o $@ -c $<

.PHONY : clean distclean allclean

include $(TOPDIR)/utils/make/Make_environments
