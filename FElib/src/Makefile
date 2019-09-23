TOPDIR      =../../
SYSDEP_DIR :=../../sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_FE_SYS)
include $(TOPDIR)/Mkinclude

VPATH = \
	common:  \
	element: \
	mesh:    \
	data:    \
	file:

OBJS_NAME_COMMON = \
	scale_polynominal.o  \
	scale_sparsemat.o    \
	scale_linalgebra.o   \
	scale_time_manager.o \
	scale_quicksort.o

OBJS_NAME_ELEMENT = \
	scale_element_base.o          \
	scale_element_line.o          \
	scale_element_quadrilateral.o \
	scale_element_hexahedral.o 

OBJS_NAME_MESH = \
	scale_meshutil_1d.o    \
	scale_meshutil_2d.o    \
	scale_meshutil_3d.o    \
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
	scale_localmeshfield_base.o      \
	scale_meshfield_base.o           \
	scale_meshfieldcomm_base.o       \
	scale_meshfieldcomm_1d.o         \
	scale_meshfieldcomm_rectdom2d.o  \
	scale_meshfieldcomm_cubedom3d.o

OBJS_NAME_FILE = \
	scale_file_history_meshfield.o

#OBJS_NAME_MESH = \
	MeshCubedSphere3D_mod.o \
	MeshCubedSphere2D_mod.o \
	Mesh3D_mod.o            \
	Mesh2D_mod.o            \
	LocalMesh3D_mod.o       \
	LocalMesh2D_mod.o       \
	MeshBase_mod.o          \
	CubedSphereUtil_mod.o   \
	MeshUtil3D_mod.o        \
	MeshUtil2D_mod.o

#OBJS_NAME_DATA = \
	LocalMeshField_mod.o     \
	LocalMeshFieldList_mod.o \
	MeshField_mod.o          \
	CSMeshField_mod.o        \
	CSMeshFieldList_mod.o    \
	MeshFieldComm_mod.o      \
	CSMeshFieldComm_mod.o 

#OBJS_NAME_IO = \
	LocalMeshFieldWriter_mod.o    \
	CSLonLatMeshFieldWriter_mod.o


OBJS_NAME = \
	$(OBJS_NAME_COMMON)  \
	$(OBJS_NAME_ELEMENT) \
	$(OBJS_NAME_MESH)    \
	$(OBJS_NAME_DATA)		 \
	$(OBJS_NAME_FILE)


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

$(OBJ_DIR)/%.o : %.F90
	$(FC) $(FFLAGS) $(SCALE_INCLUDE) $(INCLUDE) -I$(SCALEFELIBDIR)/include -o $@ -c $<

-include depend

	# DO NOT DELETE
# DO NOT DELETE
