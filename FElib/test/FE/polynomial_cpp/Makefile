################################################################################
#
# Makefile for the C++ unit test of scale_polynomial (FElib/src_cpp).
# See FElib/test/Makefile.common_cpp for the shared build recipe. This test
# calls LAPACK (dstev) directly, hence NEEDS_LAPACK.
#
################################################################################

TOPDIR           = $(abspath ../../../..)
SCALEFELIBDIR    := $(TOPDIR)/FElib

BINNAME      = test_polynomial_cpp
NEEDS_LAPACK = T

include $(SCALEFELIBDIR)/test/Makefile.common_cpp
