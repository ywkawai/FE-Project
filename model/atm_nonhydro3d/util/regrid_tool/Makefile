################################################################################
#
# Makefile
#
################################################################################

TOPDIR     = $(abspath ../../../..)
BUILD_DIR  = ./.libs
SYSDEP_DIR = ../../../../sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_FE_SYS)
include $(TOPDIR)/Mkinclude

BINNAME = regrid_tool

LIBS = $(LIBDIR)/libScaleFECore.a

VPATH = plugin:

OBJS = \
	mod_regrid_mesh.o          \
	mod_regrid_mesh_base.o     \
	mod_regrid_nodemap.o       \
	mod_regrid_interp_vcoord.o \
	mod_regrid_interp_field.o  \
	mod_regrid_file.o


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

$(BUILD_DIR)/mod_regrid_mesh.o : $(BUILD_DIR)/mod_regrid_nodemap.o $(BUILD_DIR)/mod_regrid_mesh_base.o mod_regrid_mesh.F90
$(BUILD_DIR)/mod_regrid_nodemap.o : $(BUILD_DIR)/mod_regrid_mesh_base.o mod_regrid_nodemap.F90
$(BUILD_DIR)/mod_regrid_mesh_base.o : mod_regrid_mesh_base.F90
$(BUILD_DIR)/mod_regrid_file.o     : $(BUILD_DIR)/mod_regrid_interp_field.o $(BUILD_DIR)/mod_regrid_mesh_base.o mod_regrid_file.F90
$(BUILD_DIR)/mod_regrid_interp_field.o    : $(BUILD_DIR)/mod_regrid_mesh_base.o $(BUILD_DIR)/mod_regrid_nodemap.o mod_regrid_interp_field.F90
$(BUILD_DIR)/mod_regrid_interp_vcoord.o : $(BUILD_DIR)/mod_regrid_interp_field.o $(BUILD_DIR)/mod_regrid_mesh_base.o $(BUILD_DIR)/mod_regrid_nodemap.o mod_regrid_interp_vcoord.F90

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
