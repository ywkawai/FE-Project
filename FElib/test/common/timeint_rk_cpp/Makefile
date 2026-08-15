################################################################################
#
# Makefile for the C++ unit test of scale_timeint_rk / scale_timeint_rk_butcher_tab
# (FElib/src_cpp). See FElib/test/Makefile.common_cpp for the shared build
# recipe. This test (via scale_timeint_rk_butcher_tab -> scale_linalgebra)
# calls LAPACK directly, hence NEEDS_LAPACK.
#
################################################################################

TOPDIR           = $(abspath ../../../..)
SCALEFELIBDIR    := $(TOPDIR)/FElib

BINNAME      = test_timeint_rk_cpp
NEEDS_LAPACK = T

include $(SCALEFELIBDIR)/test/Makefile.common_cpp
