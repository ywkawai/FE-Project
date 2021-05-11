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

BINNAME  = global_shallow_water
INITNAME = global_shallow_water_init

BINS = $(BINDIR)/$(BINNAME)$(POSTFIX)  \
       $(BINDIR)/$(INITNAME)$(POSTFIX) 

LIBS = $(LIBDIR)/libScaleFECore.a

VERSION = $(shell git rev-parse --short HEAD 2> /dev/null)
ifeq ($(VERSION),)
  VERSION  = $(shell cat VERSION)
else
  VERSION := $(VERSION)
endif

VPATH = \
	$(BUILD_DIR): \
	util: 

OBJS = \
	mod_user.o                  \
	mod_exp.o                   \
	mod_globalsw_component.o    \
	mod_sw_dyn.o                \
	mod_sw_dyn_vars.o           \
	mod_sw_vars.o               \
	mod_sw_bnd.o                \
	mod_sw_mesh.o


all:
	$(MAKE) build
	$(MAKE) install

build:
	$(MAKE) info
	$(MAKE) makelib
	@echo;echo "Entering scale-dg global shallow water ...";echo "Current version is " $(VERSION)
	mkdir -p $(BUILD_DIR)
	$(MAKE) $(subst $(BINDIR),$(BUILD_DIR),$(BINS))
	@echo "Complete making scale-dg-globalsw ."

install:
	mkdir -p $(BINDIR)
	$(MAKE) $(BINS)

info:
	@$(MAKE) -C $(SCALEFELIBDIR)/src --no-print-directory info
	@echo
	@echo "SCALE-DG"
	@$(MAKE) --no-print-directory conflog

makelib:
	$(MAKE) -C $(SCALEFELIBDIR)/src

$(BUILD_DIR)/$(BINNAME)$(POSTFIX): $(BUILD_DIR)/$(BINNAME).o $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(LIBS) $(CONTRIB_LIBS) $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)

$(BUILD_DIR)/$(INITNAME)$(POSTFIX) : $(BUILD_DIR)/$(INITNAME).o $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(LIBS) $(CONTRIB_LIBS) $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)

$(BUILD_DIR)/$(BINNAME).o  : $(LIBS) $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o
$(BUILD_DIR)/$(INITNAME).o : $(LIBS) $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o
$(BUILD_DIR)/mod_user.o    : $(LIBS) $(patsubst %,$(BUILD_DIR)/%,$(SUB_OBJS)) $(BUILD_DIR)/mod_exp.o

$(BINDIR)/$(BINNAME)$(POSTFIX)  : $(BUILD_DIR)/$(BINNAME)$(POSTFIX)
	install $< $@

$(BINDIR)/$(INITNAME)$(POSTFIX) : $(BUILD_DIR)/$(INITNAME)$(POSTFIX)
	install $< $@

allclean: distclean
	$(MAKE) -C $(SCALEFELIBDIR)/src allclean
	rm -f $(BINDIR)/$(BINNAME)*
	rm -f $(BINDIR)/$(PPNAME)*
	rm -rf $(TOPDIR)/bin


distclean: clean
	rm -f $(BINDIR)/$(BINNAME)

clean:
	rm -rf $(BUILD_DIR)
	rm -f *.o *.mod *.lst *~

.SUFFIXES:
.SUFFIXES: .o .F90 .mod

%.mod: %.F90
	$(MAKE) $(patsubst %.F90,%.o,$<)

$(BUILD_DIR)/%.o: %.F90
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(MODDIR) $(CONTRIB_INCLUDE) -I$(SCALEFELIBDIR)/include $(SCALE_NETCDF_INCLUDE) $(MODDIROPT) $(BUILD_DIR) -o $@ -c $<

.PHONY : clean distclean allclean

include $(TOPDIR)/utils/make/Make_environments


$(BUILD_DIR)/$(BINNAME).o: \
	$(BUILD_DIR)/mod_globalsw_component.o \
	$(BUILD_DIR)/mod_user.o 

$(BUILD_DIR)/$(INITNAME).o: \
	$(BUILD_DIR)/mod_globalsw_component.o \
	$(BUILD_DIR)/mod_user.o 

$(BUILD_DIR)/mod_globalsw_component.o: \
	$(BUILD_DIR)/mod_sw_dyn.o

$(BUILD_DIR)/mod_sw_dyn.o:       \
	$(BUILD_DIR)/mod_sw_dyn_vars.o \
	$(BUILD_DIR)/mod_sw_vars.o     \
	$(BUILD_DIR)/mod_sw_bnd.o

$(BUILD_DIR)/mod_user.o: \
  mod_user.F90 \
	$(BUILD_DIR)/mod_globalsw_component.o \
	$(BUILD_DIR)/mod_sw_vars.o \
	$(BUILD_DIR)/mod_exp.o

$(BUILD_DIR)/mod_exp.o: \
	$(BUILD_DIR)/mod_sw_vars.o

$(BUILD_DIR)/mod_sw_vars.o: \
	$(BUILD_DIR)/mod_sw_mesh.o

$(BUILD_DIR)/mod_sw_dyn_vars.o: \
	$(BUILD_DIR)/mod_sw_mesh.o

$(BUILD_DIR)/mod_sw_bnd.o: \
	$(BUILD_DIR)/mod_sw_mesh.o
