################################################################################
#
# Makefile for FElib C++ implementation (libScaleFECpp)
#
################################################################################

TOPDIR      =../../
SYSDEP_DIR :=../../sysdep

BUILD_DIR = "undefined"
include $(SYSDEP_DIR)/Makedef.$(SCALE_FE_SYS)
include $(TOPDIR)/Mkinclude

LIBNAME = libScaleFECoreCpp$(POSTFIX).a

VPATH = \
	common:

OBJS_NAME_COMMON = \
	scale_sparsemat.o   \
	scale_polynomial.o

OBJS_NAME = \
	$(OBJS_NAME_COMMON)

#------------------------------------------
OBJS = $(patsubst %, $(BUILD_DIR)/%, $(OBJS_NAME))
DEPS = $(OBJS:.o=.d)

all:
	$(MAKE) build
	$(MAKE) install

build:
	@echo;echo "Entering FElib/src_cpp..."
	mkdir -p $(BUILD_DIR)
	$(MAKE) $(LIBNAME)
	@echo "Complete making $(LIBNAME)."

install:
	mkdir -p $(LIBDIR)
	install $(LIBNAME) $(LIBDIR)

$(LIBNAME): $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

$(BUILD_DIR):
	mkdir -p $@

clean:
	rm -rf $(BUILD_DIR)
	rm -f *.o *.d

distclean: clean
	rm -f $(LIBNAME)

allclean: distclean
	rm -f $(LIBDIR)/$(LIBNAME)

.SUFFIXES:
.SUFFIXES: .o .cpp .d

$(BUILD_DIR)/%.o : %.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -I$(SCALEFELIBDIR)/src_cpp -MMD -MP -o $@ -c $<

.PHONY : all build install clean distclean allclean

-include $(DEPS)
