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

BINNAME  = scale-dg
INITNAME = scale-dg_init
PPNAME   = scale-dg_pp

BINS = $(BINDIR)/$(BINNAME)$(POSTFIX)  \
       $(BINDIR)/$(INITNAME)$(POSTFIX) \
       $(BINDIR)/$(PPNAME)$(POSTFIX) 
LIBS = $(LIBDIR)/libScaleFECore$(POSTFIX).a

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
	user:         \
	preprocess:

OBJS = \
	mod_dg_launcher.o            \
	mod_dg_driver.o              \
	mod_atmos_component.o        \
	mod_atmos_dyn.o              \
	mod_atmos_dyn_vars.o         \
	mod_atmos_phy_sfc.o          \
	mod_atmos_phy_sfc_vars.o     \
	mod_atmos_phy_mp.o           \
	mod_atmos_phy_mp_vars.o      \
	mod_atmos_phy_tb.o           \
	mod_atmos_phy_tb_vars.o      \
	mod_atmos_mesh.o             \
	mod_atmos_mesh_rm.o          \
	mod_atmos_mesh_gm.o          \
	mod_atmos_vars.o             \
	mod_user_base.o              \
	mod_user.o                   \
	mod_experiment.o             \
	\
	mod_dg_prep.o                \
	mod_mktopo.o                 \
	mod_mkinit_util.o            \
	mod_mkinit.o


all:
	$(MAKE) build
	$(MAKE) install

build:
	$(MAKE) info
	$(MAKE) makelib
	@echo;echo "Entering scale-dg...";echo "Current version is " $(VERSION)
	mkdir -p $(BUILD_DIR)
	$(MAKE) makebin
	@echo "Complete making scale-dg."

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

makebin:  $(BUILD_DIR)/$(BINNAME)$(POSTFIX) $(BUILD_DIR)/$(PPNAME)$(POSTFIX)

$(BUILD_DIR)/$(BINNAME)$(POSTFIX)  : $(BUILD_DIR)/$(BINNAME).o  $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(LIBS) $(CONTRIB_LIBS) $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)

$(BUILD_DIR)/$(PPNAME)$(POSTFIX) : $(BUILD_DIR)/$(PPNAME).o $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(LIBS) $(CONTRIB_LIBS) $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)

$(BUILD_DIR)/$(BINNAME).o  : $(LIBS) $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o
$(BUILD_DIR)/$(PPNAME).o : $(LIBS) $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(BUILD_DIR)/mod_user.o
$(BUILD_DIR)/mod_user.o    : $(LIBS) $(patsubst %,$(BUILD_DIR)/%,$(SUB_OBJS)) $(BUILD_DIR)/mod_experiment.o $(BUILD_DIR)/mod_mkinit_util.o

$(BINDIR)/$(BINNAME)$(POSTFIX)  : $(BUILD_DIR)/$(BINNAME)$(POSTFIX)
	install $< $@

$(BINDIR)/$(PPNAME)$(POSTFIX)   : $(BUILD_DIR)/$(PPNAME)$(POSTFIX)
	install $< $@

$(BINDIR)/$(INITNAME)$(POSTFIX) : $(BINDIR)/$(PPNAME)$(POSTFIX)
	(cd $(BINDIR); ln -sf `basename $<` $@)


distclean: clean
	rm -f $(BINDIR)/$(BINNAME)$(POSTFIX)
	rm -f $(BINDIR)/$(INITNAME)$(POSTFIX)
	rm -f $(BINDIR)/$(PPNAME)$(POSTFIX)
	rm -rf ./.libs*

clean:
	rm -rf $(BUILD_DIR)
	rm -f *.o *.mod *.lst *~

dependency:
	./makedepend .

.SUFFIXES:
.SUFFIXES: .o .F90 .mod

%.mod: %.F90
	$(MAKE) $(patsubst %.F90,%.o,$<)

$(BUILD_DIR)/%.o: %.F90
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(MODDIR) -I../include $(CONTRIB_INCLUDE) -I$(SCALEFELIBDIR)/include $(SCALE_NETCDF_INCLUDE) $(MODDIROPT) $(BUILD_DIR) -o $@ -c $<

$(BUILD_DIR)/%.o : %.c
	$(CC) $(CFLAGS) -o $@ -c $<

allclean: distclean
	$(MAKE) -C $(SCALEFELIBDIR)/src allclean
	rm -f $(BINDIR)/$(BINNAME)*
	rm -f $(BINDIR)/$(INITNAME)*
#        rm -f $(BINDIR)/$(PPNAME)*
	rm -rf $(TOPDIR)/bin


distclean: clean
	rm -f $(BINDIR)/$(BINNAME)$(POSTFIX)
	rm -f $(BINDIR)/$(INITNAME)$(POSTFIX)
#        rm -f $(BINDIR)/$(PPNAME)$(POSTFIX)
	rm -rf ./.libs*

.PHONY : clean distclean allclean depend dependency

include $(TOPDIR)/utils/make/Make_environments
-include depend
