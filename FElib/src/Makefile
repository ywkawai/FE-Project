################################################################################
#
# Makefile for scale library
#
################################################################################

TOPDIR      =../../
SYSDEP_DIR :=../../sysdep

BUILD_DIR = "undefined"
include $(SYSDEP_DIR)/Makedef.$(SCALE_FE_SYS)
include $(TOPDIR)/Mkinclude

LIBNAME=libScaleFECore$(POSTFIX).a

VPATH = \
	common:           \
	element:          \
	mesh:             \
	data:             \
	file:             \
	fluid_dyn_solver: \
	turbulence:       \
	microphysics:     \
	model_framework:

LIBVERSION = $(shell git rev-parse --short HEAD 2> /dev/null)
ifeq ($(LIBVERSION),)
  LIBVERSION  = $(shell cat VERSION)
else
  LIBVERSION := $(LIBVERSION)
endif


OBJS_NAME_COMMON = \
	scale_polynominal.o            \
	scale_polygon.o                \
	scale_sparsemat.o              \
	scale_linalgebra.o             \
	scale_time_manager.o           \
	scale_quicksort.o              \
	scale_timeint_rk_butcher_tab.o \
	scale_timeint_rk.o             \
	scale_linkedlist.o             \
	scale_gmres.o                  \
	scale_geographic_coord_cnv.o   \
	scale_cubedsphere_coord_cnv.o  \
	scale_coriolis_param.o

OBJS_NAME_ELEMENT = \
	scale_element_base.o          \
	scale_element_line.o          \
	scale_element_quadrilateral.o \
	scale_element_hexahedral.o    \
	scale_element_modalfilter.o

OBJS_NAME_MESH = \
	scale_meshutil_1d.o            \
	scale_meshutil_2d.o            \
	scale_meshutil_3d.o            \
	scale_meshutil_vcoord.o        \
	scale_meshutil_cubedsphere2d.o \
	scale_meshutil_cubedsphere3d.o \
	scale_mesh_bndinfo.o           \
	scale_localmesh_base.o         \
	scale_localmesh_1d.o           \
	scale_localmesh_2d.o           \
	scale_localmesh_3d.o           \
	scale_mesh_base.o	           \
	scale_mesh_base1d.o	           \
	scale_mesh_base2d.o	           \
	scale_mesh_base3d.o	           \
	scale_mesh_linedom1d.o         \
	scale_mesh_rectdom2d.o         \
	scale_mesh_cubedom3d.o         \
	scale_mesh_cubedspheredom2d.o  \
	scale_mesh_cubedspheredom3d.o  \
	scale_mesh_topography.o

OBJS_NAME_DATA = \
	scale_variableinfo.o                     \
	scale_localmeshfield_base.o              \
	scale_meshfield_base.o                   \
	scale_meshfieldcomm_base.o               \
	scale_meshfieldcomm_1d.o                 \
	scale_meshfieldcomm_rectdom2d.o          \
	scale_meshfieldcomm_cubedspheredom2d.o   \
	scale_meshfieldcomm_cubedom3d.o          \
	scale_meshfieldcomm_cubedspheredom3d.o   \
	scale_meshfield_statistics.o             \
	scale_meshfield_analysis_numerror_base.o \
	scale_meshfield_analysis_numerror.o

OBJS_NAME_FILE = \
	scale_file_common_meshfield.o  \
	scale_file_base_meshfield.o    \
	scale_file_history_meshfield.o \
	scale_file_restart_meshfield.o \
	scale_file_monitor_meshfield.o

OBJS_NAME_FLUID_DYN_SOLVER = \
	scale_atm_dyn_dgm_modalfilter.o                  \
	scale_atm_dyn_dgm_spongelayer.o                  \
	scale_atm_dyn_dgm_hydrostatic.o                  \
	scale_atm_dyn_dgm_nonhydro3d_numdiff.o           \
	scale_atm_dyn_dgm_nonhydro2d.o                   \
	scale_atm_dyn_dgm_nonhydro3d_common.o            \
	scale_atm_dyn_dgm_nonhydro3d_etot_heve_numflux.o \
	scale_atm_dyn_dgm_nonhydro3d_rhot_heve_numflux.o \
	scale_atm_dyn_dgm_nonhydro3d_etot_heve.o         \
	scale_atm_dyn_dgm_nonhydro3d_rhot_heve.o         \
	scale_atm_dyn_dgm_nonhydro3d_etot_hevi_common.o  \
	scale_atm_dyn_dgm_nonhydro3d_rhot_hevi_common.o  \
	scale_atm_dyn_dgm_nonhydro3d_etot_hevi_numflux.o \
	scale_atm_dyn_dgm_nonhydro3d_rhot_hevi_numflux.o \
	scale_atm_dyn_dgm_nonhydro3d_etot_hevi.o         \
	scale_atm_dyn_dgm_nonhydro3d_rhot_hevi.o         \
	scale_atm_dyn_dgm_nonhydro3d_rhot_heve_splitform.o \
	scale_atm_dyn_dgm_nonhydro3d_rhot_hevi_splitform.o \
	scale_atm_dyn_dgm_nonhydro3d_hevi_gmres.o        \
	scale_atm_dyn_dgm_globalnonhydro3d_rhot_heve.o   \
	scale_atm_dyn_dgm_globalnonhydro3d_etot_heve.o   \
	scale_atm_dyn_dgm_globalnonhydro3d_rhot_hevi.o   \
	scale_atm_dyn_dgm_globalnonhydro3d_etot_hevi.o   \
	scale_atm_dyn_dgm_bnd.o                          \
	scale_atm_dyn_dgm_driver_base.o                  \
	scale_atm_dyn_dgm_driver_nonhydro3d.o            \
	scale_atm_dyn_dgm_globalsw.o                     \
	scale_atm_dyn_dgm_driver_trcadv3d.o	             \
	scale_atm_dyn_dgm_trcadvect3d_heve.o	

OBJS_NAME_TURBULENCE = \
	scale_atm_phy_tb_dgm_driver.o    \
	scale_atm_phy_tb_dgm_dns.o       \
	scale_atm_phy_tb_dgm_smg.o       \
	scale_atm_phy_tb_dgm_globalsmg.o \
	scale_atm_phy_tb_dgm_common.o 

OBJS_NAME_MICROPHYS = \
	scale_atm_phy_mp_dgm_common.o

OBJS_NAME_MODEL_FRAMEWORK = \
	scale_model_meshbase_manager.o \
	scale_model_mesh_manager.o     \
	scale_model_component.o        \
	scale_model_component_proc.o   \
	scale_model_var_manager.o

OBJS_NAME = \
	$(OBJS_NAME_COMMON)           \
	$(OBJS_NAME_ELEMENT)          \
	$(OBJS_NAME_MESH)             \
	$(OBJS_NAME_DATA)		      \
	$(OBJS_NAME_FILE)             \
	$(OBJS_NAME_FLUID_DYN_SOLVER) \
	$(OBJS_NAME_TURBULENCE)       \
	$(OBJS_NAME_MICROPHYS)        \
	$(OBJS_NAME_MODEL_FRAMEWORK)


OBJS_NAME_SCALELIB = \
	scale_precision.o  \
	scale_const.o      \
	scale_io.o         \
	scale_time.o       \
	scale_prc_cartesC.o                \
	scale_comm_cartesC.o               \
	scale_atmos_thermodyn.o  \
	scale_atmos_saturation.o \
	scale_atmos_grid_cartesC_index.o   \
	scale_atmos_grid_cartesC.o         \
	\
	scale.o	\
	scale_sys_sigvars.o \
	scale_file_netcdf.o \
	scale_file_grads.o	\
	scale_calendar.o	\
	scale_debug.o \
	scale_fpm.o \
	scale_index.o \
	scale_mapprojection.o \
	scale_matrix.o	\
	scale_sort.o	\
	scale_hash.o	\
	scale_prc.o	       \
	scale_prof.o	\
	scale_random.o	\
	scale_sigvars.o \
	scale_statistics.o	\
	scale_tracer.o	\
	scale_vector.o	\
	scale_monitor.o	\
	scale_file.o	\
	scale_file_cartesC.o	\
	scale_file_h.o	\
	scale_file_history.o	\
	scale_file_netcdf.o \
	\
	scale_atmos_adiabat.o \
	scale_atmos_aerosol.o \
	scale_atmos_hydrometeor.o \
	\
	scale_atmos_phy_mp_common.o \
	scale_atmos_phy_mp_kessler.o	\
	scale_atmos_phy_mp_sn14.o	\
	scale_atmos_phy_mp_suzuki10.o \
	scale_atmos_phy_mp_tomita08.o \
	\
	scale_atmos_phy_sf_bulk.o \
	scale_atmos_phy_sf_const.o \
	\
	scale_ocean_grid_cartesC.o	\
	scale_ocean_grid_cartesC_index.o	\
	scale_ocean_grid_cartesC_real.o \
	scale_ocean_grid_icoA.o \
	scale_ocean_grid_icoA_index.o \
	scale_ocean_dyn_slab.o	\
	scale_ocean_dyn_offline.o \
	scale_ocean_phy_ice_simple.o	\
	scale_ocean_phy_albedo.o	\
	scale_ocean_phy_albedo_nakajima00.o \
	scale_ocean_phy_roughness.o \
	scale_ocean_phy_roughness_miller92.o	\
	scale_ocean_phy_roughness_moon07.o	\
	scale_ocean_phy_tc.o	\
	scale_land_grid_cartesC.o \
	scale_land_grid_cartesC_index.o \
	scale_land_grid_cartesC_real.o	\
	scale_land_grid_icoA.o	\
	scale_land_grid_icoA_index.o	\
	scale_land_dyn_bucket.o \
	scale_land_phy_matsiro.o	\
	scale_land_phy_snow_diagnos.o \
	scale_land_phy_snow_ky90.o	\
	scale_urban_grid_cartesC.o	\
	scale_urban_grid_cartesC_index.o	\
	scale_urban_grid_cartesC_real.o \
	scale_urban_grid_icoA.o \
	scale_urban_grid_icoA_index.o \
	scale_urban_dyn_kusaka01.o	\
	scale_bulkflux.o	\
	scale_landuse.o \
	scale_cpl_sfc_index.o \
	scale_cpl_phy_sfc_fixed_temp.o	\
	scale_cpl_phy_sfc_skin.o \
	scale_spnudge.o \
	scale_dft.o	

MODS_SCALELIB = $(addprefix $(SCALE)/include/, \
	$(shell echo $(patsubst %.o,%.mod, \
	$(filter-out scale_sys_sigvars.o, \
	$(filter-out scale_file_netcdf.o, \
	$(filter-out read_toshiba.o, \
	$(filter-out read_toshiba_mpr.o, \
	$(OBJS_NAME_SCALELIB) ))))) | tr A-Z a-z) )

#------------------------------------------
OBJS=$(patsubst %, $(BUILD_DIR)/%, $(OBJS_NAME))

all:
	$(MAKE) prep_scalelib
	$(MAKE) build
	$(MAKE) install

build:
	@echo;echo "Entering scaleFElib..."; echo "Current version is " $(LIBVERSION)
	mkdir -p ${BUILD_DIR}
	$(MAKE) $(LIBNAME)
	@echo "Complete making scaleFElib."

install:
	$(MAKE) modules
	$(MAKE) $(LIBDIR)/$(LIBNAME)


info:
	@$(MAKE) --no-print-directory envlog
	@echo
	@echo "SCALEFElib"
	@$(MAKE) --no-print-directory conflog


prep_scalelib: 
	echo "Prepair libscale_sub.a and include files for FE-Project .."
	$(MAKE) $(LIBDIR)/libscale_sub.a

$(LIBDIR)/libscale_sub.a: $(SCALE)/lib/libscale.a
	mkdir -p $(LIBDIR)
	$(AR) -x $< $(OBJS_NAME_SCALELIB)
	$(AR) $(ARFLAGS) $(@F) $(OBJS_NAME_SCALELIB)
	install $(@F) $(LIBDIR)
	rm -f $(OBJS_NAME_SCALELIB)

	mkdir -p $(MODDIR)/scalelib
	install $(MODS_SCALELIB) $(MODDIR)/scalelib

$(LIBDIR)/$(LIBNAME): $(LIBNAME)
	mkdir -p $(LIBDIR)
	install $< $@

$(LIBNAME): $(OBJS)
	echo "build $(LIBNAME).."
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

modules: $(OBJS)
	mkdir -p $(MODDIR)
	install $(MODS) $(MODDIR)

clean:
	rm -rf $(BUILD_DIR)
	rm -f *.o *.mod *.lst

distclean: clean
	rm -rf $(LIBNAME)
	rm -rf libscale_sub.a

allclean: distclean

dependency:
	./makedepend .


.SUFFIXES:
.SUFFIXES: .o .F90 .c .erb .mod

%.F90 : %.F90.erb
	erb $< > $@

common/scale_quicksort.F90: scale_quicksort.F90.erb
common/scale_timeint_rk.F90: scale_timeint_rk.F90.erb

$(BUILD_DIR)/%.o : %.F90
	$(FC) $(FFLAGS) $(SCALE_INCLUDE) $(INCLUDE) -I$(SCALEFELIBDIR)/include -I$(BUILD_DIR) $(MODDIROPT) $(BUILD_DIR) -o $@ -c $<

.PHONY : clean distclean allclean depend modules dependency

include $(TOPDIR)/utils/make/Make_environments

-include depend

# DO NOT DELETE