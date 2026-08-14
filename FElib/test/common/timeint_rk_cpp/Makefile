################################################################################
#
# Makefile for the C++ unit test of scale_timeint_rk / scale_timeint_rk_butcher_tab
# (FElib/src_cpp)
#
# Like sparsemat_cpp, this test has no MPI/scale_io dependency, so it does not
# use FElib/test/Makefile.common's jobshell-based `run` target. Like
# polynomial_cpp, scale_timeint_rk_butcher_tab calls LAPACK (dgetrf/dgetrs)
# directly, so the final link additionally needs $(SCALE_MATHLIB_LIBS)
# (-llapack -lblas).
#
################################################################################

TOPDIR           = $(abspath ../../../..)
SYSDEP_DIR       := $(TOPDIR)/sysdep
SCALEFELIBDIR    := $(TOPDIR)/FElib

BUILD_DIR = "undefined"
include $(SYSDEP_DIR)/Makedef.$(SCALE_FE_SYS)
include $(TOPDIR)/Mkinclude

BINNAME = test_timeint_rk_cpp
LIBCPP  = $(LIBDIR)/libScaleFECoreCpp$(POSTFIX).a

all: build

build: $(BINNAME)
	@echo "Complete making $(BINNAME)."

$(BINNAME): $(BUILD_DIR)/$(BINNAME).o $(LIBCPP)
	$(CXX) $(CXXFLAGS) -o $@ $(BUILD_DIR)/$(BINNAME).o $(LIBCPP) $(SCALE_MATHLIB_LIBS)

$(BUILD_DIR)/$(BINNAME).o: $(BINNAME).cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -I$(SCALEFELIBDIR)/src_cpp -I$(SCALEFELIBDIR)/src_cpp/common -MMD -MP -o $@ -c $<

$(BUILD_DIR):
	mkdir -p $@

$(LIBCPP):
	@echo "$(LIBCPP) not found -- build and install it first:"
	@echo "  \$$(MAKE) -C $(SCALEFELIBDIR)/src_cpp build install"
	@exit 1

run: $(BINNAME)
	./$(BINNAME)

check: run

clean:
	rm -rf $(BUILD_DIR)

distclean: clean
	rm -f $(BINNAME)

.SUFFIXES:
.SUFFIXES: .o .cpp .d

.PHONY : all build run check clean distclean

-include $(BUILD_DIR)/$(BINNAME).d
