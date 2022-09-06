################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../..)
BUILD_DIR   = ./.libs
SYSDEP_DIR  = $(TOPDIR)/sysdep

include $(SYSDEP_DIR)/Makedef.$(SCALE_FE_SYS)
include $(TOPDIR)/Mkinclude

BINNAME = test_eigen_analysis

OBJS = mod_dg_eigen_analysis.o

LIBS = $(LIBDIR)/libScaleFECore.a

SAMPLE_AUX_DIR = ../auxiliary
SAMPLE_AUX_INCLUDE = $(SAMPLE_AUX_DIR)/.libs
SAMPLE_AUX_LIB = $(SAMPLE_AUX_DIR)/libScaleFESampleAux.a

all:
	$(MAKE) envlog
	$(MAKE) makedir
	$(MAKE) make_sample_aux
	$(MAKE) makebin

makedir:
	mkdir -p $(BUILD_DIR)

makebin: $(BINNAME)
	@echo "Complete making."

run:
	export OMP_NUM_THREADS=1
	mkdir -p data
	mpirun -n 1 ./$(BINNAME) test.conf

vis:
#	python ./visualize/eigen_graph.py
#	python ./visualize/ampl_phase_error_graph.py
#	python ./visualize/Rdiff_Rdisp_graph.py
	python ./visualize2/Rdiff_Rdisp_graph.py
#	python ./visualize/dispersion_dissipation_graph.py

make_sample_aux:
	make -C $(SAMPLE_AUX_DIR)

$(BINNAME): $(BUILD_DIR)/$(BINNAME).o $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(SAMPLE_AUX_LIB) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(LIBS) $(CONTRIB_LIBS) $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS)

$(BUILD_DIR)/$(BINNAME).o: $(BINNAME).f90 $(patsubst %,$(BUILD_DIR)/%,$(OBJS))

distclean: clean
	rm -f $(BINNAME)
	rm -rf ./data
	rm -rf ./visualize/*.png

clean:
	rm -rf $(BUILD_DIR)

.SUFFIXES:
.SUFFIXES: .o .f90 .mod

%.mod: %.f90
	$(MAKE) $(patsubst %.f90,%.o,$<)

$(BUILD_DIR)/%.o: %.f90
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(MODDIR) $(CONTRIB_INCLUDE) -I$(SCALEFELIBDIR)/include -I$(SAMPLE_AUX_INCLUDE) $(SCALE_NETCDF_INCLUDE) $(MODDIROPT) $(BUILD_DIR) -o $@ -c $<

.PHONY : clean distclean allclean

include $(TOPDIR)/utils/make/Make_environments

