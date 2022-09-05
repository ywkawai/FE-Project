################################################################################
#
# Makefile
#
################################################################################

TOPDIR     = $(abspath ../../../../../../../../)
BUILD_DIR  = ./.libs
SYSDEP_DIR = $(TOPDIR)/sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_FE_SYS)
include $(TOPDIR)/Mkinclude

BINNAME = pbltb_analysis_fr

LIBS = $(LIBDIR)/libScaleFECore.a

VPATH = plugin:

OBJS = \
	mod_diag_tb.o \
	mod_common.o


all:
	$(MAKE) build
	$(MAKE) install

build:
	mkdir -p $(BUILD_DIR)
	$(MAKE) $(BUILD_DIR)/$(BINNAME)

install:
	mkdir -p $(BINDIR)
	$(INSTALL) $(BUILD_DIR)/$(BINNAME) $(BINDIR)/$(BINNAME)

$(OBJECT): $(MODS)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)
	@mv *.mod inc/; mv *.o inc/;

$(BUILD_DIR)/$(BINNAME) : $(BUILD_DIR)/prg_$(BINNAME).o  $(patsubst %,$(BUILD_DIR)/%,$(OBJS))
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(LIBS) $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS) $(CONTRIB_LIBS)

$(BUILD_DIR)/prg_$(BINNAME).o : prg_$(BINNAME).F90 $(patsubst %,$(BUILD_DIR)/%,$(OBJS))

distclean: clean
	rm -f $(BINDIR)/$(BINNAME)

clean:
	rm -rf $(BUILD_DIR)



.SUFFIXES:
.SUFFIXES: .o .F90 .mod

%.mod: %.F90
	make $(patsubst %.F90,%.o,$<)

$(BUILD_DIR)/%.o : %.F90
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(MODDIR) -I../include $(CONTRIB_INCLUDE) -I$(SCALEFELIBDIR)/include -I$(SAMPLE_AUX_INCLUDE) $(SCALE_NETCDF_INCLUDE) $(MODDIROPT) $(BUILD_DIR) -o $@ -c $<

.PHONY : clean distclean allclean

include $(TOPDIR)/utils/make/Make_environments

$(BUILD_DIR)/prg_vertical_structure.o: \
	$(BUILD_DIR)/mod_diag_tb.o $(BUILD_DIR)/mod_common.o
$(BUILD_DIR)/mod_diag_tb.o: \
	$(BUILD_DIR)/mod_common.o
