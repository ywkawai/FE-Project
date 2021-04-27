################################################################################
#
# Makefile for the reference manual
#
################################################################################

TOPDIR = ../..

all: html tex

dox:
				$(MAKE) -C $(TOPDIR)/FElib/doc dox
				$(MAKE) -C $(TOPDIR)/model/atm_nonhydro3d/doc dox

html: dox
				doxygen FE-Project.cfg

tex: latex/refman.pdf

latex/refman.pdf: latex/refman.tex
				$(MAKE) -C latex

clean:
				rm -rf html latex

allclean: clean
				$(MAKE) -C $(TOPDIR)/FElib/doc clean
				$(MAKE) -C $(TOPDIR)/model/atm_nonhydro3d/doc clean

.PHONY: doxygen clean
