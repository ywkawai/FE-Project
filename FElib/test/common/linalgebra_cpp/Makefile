################################################################################
#
# Makefile for the C++ unit test of scale_linalgebra (FElib/src_cpp).
# See FElib/test/Makefile.common_cpp for the shared build recipe. This test
# calls LAPACK (dgetrf/dgetrs) directly, hence NEEDS_LAPACK.
#
################################################################################

TOPDIR           = $(abspath ../../../..)
SCALEFELIBDIR    := $(TOPDIR)/FElib

BINNAME      = test_linalgebra_cpp
NEEDS_LAPACK = T

include $(SCALEFELIBDIR)/test/Makefile.common_cpp
