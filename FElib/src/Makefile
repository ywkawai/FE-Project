TOPDIR      =../../
SYSDEP_DIR :=../../sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_FE_SYS)
include $(TOPDIR)/Mkinclude

VPATH = \
	common:           \
	element:          \
	mesh:             \
	data:             \
	file:             \
	fluid_dyn_solver: \
	model_framework:  

OBJS_NAME_COMMON = \
	scale_polynominal.o  \
	scale_sparsemat.o    \
	scale_linalgebra.o   \
	scale_time_manager.o \
	scale_quicksort.o    \
	scale_timeint_rk.o   \
	scale_linkedlist.o   \
	scale_gmres.o

OBJS_NAME_ELEMENT = \
	scale_element_base.o          \
	scale_element_line.o          \
	scale_element_quadrilateral.o \
	scale_element_hexahedral.o 

OBJS_NAME_MESH = \
	scale_meshutil_1d.o    \
	scale_meshutil_2d.o    \
	scale_meshutil_3d.o    \
	scale_mesh_bndinfo.o   \
	scale_localmesh_base.o \
	scale_localmesh_1d.o   \
	scale_localmesh_2d.o   \
	scale_localmesh_3d.o   \
	scale_mesh_base.o	     \
	scale_mesh_base1d.o	   \
	scale_mesh_base2d.o	   \
	scale_mesh_base3d.o	   \
	scale_mesh_linedom1d.o \
	scale_mesh_rectdom2d.o \
	scale_mesh_cubedom3d.o

OBJS_NAME_DATA = \
	scale_variableinfo.o             \
	scale_localmeshfield_base.o      \
	scale_meshfield_base.o           \
	scale_meshfieldcomm_base.o       \
	scale_meshfieldcomm_1d.o         \
	scale_meshfieldcomm_rectdom2d.o  \
	scale_meshfieldcomm_cubedom3d.o

OBJS_NAME_FILE = \
	scale_file_common_meshfield.o  \
	scale_file_history_meshfield.o \
	scale_file_restart_meshfield.o

OBJS_NAME_FLUID_DYN_SOLVER = \
	scale_atm_dyn_modalfilter.o          \
	scale_atm_dyn_nonhydro3d_numdiff.o   \
	scale_atm_dyn_nonhydro2d.o           \
	scale_atm_dyn_nonhydro3d_heve.o      \
	scale_atm_dyn_nonhydro3d_hevi.o      \
	scale_atm_dyn_nonhydro3d_hevi_gmres.o 


OBJS_NAME_MODEL_FRAMEWORK = \
	scale_model_mesh_manager.o   \
	scale_model_component.o      \
	scale_model_component_proc.o \
	scale_model_var_manager.o

OBJS_NAME = \
	$(OBJS_NAME_COMMON)           \
	$(OBJS_NAME_ELEMENT)          \
	$(OBJS_NAME_MESH)             \
	$(OBJS_NAME_DATA)		          \
	$(OBJS_NAME_FILE)             \
	$(OBJS_NAME_FLUID_DYN_SOLVER) \
	$(OBJS_NAME_MODEL_FRAMEWORK)
	


OBJ_DIR= .
MOD_DIR=$(TOPDIR)/include/
LIBNAME=libScaleFECore.a

#------------------------------------------
OBJS=$(patsubst %, $(OBJ_DIR)/%, $(OBJS_NAME))

all:
	$(MAKE) makedir
	@echo;echo "Entering scaleFElib..."; echo "Current version is " $(LIBVERSION)
	$(MAKE) $(LIBDIR)/$(LIBNAME)
	$(MAKE) modules
	@echo "Complete making scaleFElib."
	make modules

makedir:
	mkdir -p $(OBJ_DIR)
	mkdir -p $(MODDIR)
	mkdir -p $(LIBDIR)

$(LIBDIR)/$(LIBNAME): $(LIBNAME)
	install $< $@

$(LIBNAME): $(OBJS)
	echo "build $(LIBNAME).."
	ar -r $@ $?
	ranlib $@

modules: $(OBJS)
	install $(MODS) $(MOD_DIR)

clean:
	rm -rf $(OBJS) *.mod

distclean: clean
	rm -rf $(LIBNAME)

allclean: distclean

dependency:
	./makedepend .


.SUFFIXES:
.SUFFIXES: .o .F90 .c .erb .mod

%.F90 : %.F90.erb
	erb $< > $@

common/scale_quicksort.F90: scale_quicksort.F90.erb
common/scale_timeint_rk.F90: scale_timeint_rk.F90.erb

$(OBJ_DIR)/%.o : %.F90
	$(FC) $(FFLAGS) $(SCALE_INCLUDE) $(INCLUDE) -I$(SCALEFELIBDIR)/include -o $@ -c $<

-include depend

	# DO NOT DELETE
# DO NOT DELETE
