################################################################################
#
# Makefile for the C++ unit test of scale_element_base / scale_element_line
# (FElib/src_cpp). See FElib/test/Makefile.common_cpp for the shared build
# recipe. This test (via scale_element_base -> scale_linalgebra) calls
# LAPACK directly, hence NEEDS_LAPACK.
#
################################################################################

TOPDIR           = $(abspath ../../../..)
SCALEFELIBDIR    := $(TOPDIR)/FElib

BINNAME      = test_element_line_cpp
NEEDS_LAPACK = T

include $(SCALEFELIBDIR)/test/Makefile.common_cpp
