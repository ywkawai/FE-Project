###########################################################
#
# Makefile for compiling tex document
#
###########################################################

# target file name
TEXNAME  = scale_dg_users_guide
VERFNAME = 00_version.tex
VERTEMPL = 00_version_template.tex

SCALE_VERSION = 5.5.5

# commands to compile documents
LATEX    = pdflatex
BIBTEX   = bibtex
HTMLBIN  = pandoc --to=html
CURL     = curl -f -s -S --create-dirs
IMAGE_URL = https://scale.riken.jp/archives/figure
IMAGE_URL2 = https://r-ccs-climate.riken.jp/members/kawai/SCALE-DG_doc/archives
IMAGE_PATH = $(IMAGE_URL)/$(if $(shell git branch --contains 2>/dev/null|grep -v master|grep -v release|grep -v hotfix),devel,$(SCALE_VERSION))
IMAGE_PATH2 = $(IMAGE_URL2)

VERSION  = $(shell cat ../../model/atm_nonhydro3d/VERSION)
VERSION_SCALE  = 5.5.5

IMAGE_BASENAME = scale_logo_final_ULWB.pdf
IMAGE_BASENAME2 = $(shell grep includegraphics *tex | ruby -n -e '/^\s*\\includegraphics.*\{(.+)\}/ =~ $$_.split(":")[1] && print(/\A([^.]+)(\..+)?\Z/=~File.basename($$1) && $$1+($$2 || ".eps"),"\n") ') title_wallpaper.png

FIG_DIR = ../figure
IMAGES = $(patsubst %,$(FIG_DIR)/%,$(IMAGE_BASENAME))
IMAGES2 = $(patsubst %,$(FIG_DIR)/%,$(IMAGE_BASENAME2))

SRCS = $(wildcard *.tex) $(IMAGES)


all: update_images pdf

pdf: $(TEXNAME).pdf

html: $(TEXNAME).html

update_images:
	$(if $(IMAGE_URL),,$(error IMAGE_URL is empty))
	for image in $(IMAGES) ; do \
	  $(CURL) $(IMAGE_PATH)/`basename $$image` -z $$image -o $$image; \
	done
	for image in $(IMAGES2) ; do \
	  echo $$image; \
	  $(CURL) $(IMAGE_PATH2)/`basename $$image` -z $$image -o $$image; \
	done

clean:
	rm -f *.dvi *.aux *.bbl *.blg *.log *.out *.toc $(VERFNAME)

allclean: clean
	rm -rf *.pdf $(FIG_DIR)

$(TEXNAME).pdf: $(VERFNAME) $(SRCS) ../reference.bib
	$(LATEX)  $(TEXNAME)
	$(BIBTEX) $(TEXNAME)
	$(LATEX)  $(TEXNAME)

$(VERFNAME): $(VERTEMPL)
	$(shell cat $< | sed s/"#VERSION#"/$(VERSION)/g > $@)

$(TEXNAME).html: $(TEXNAME).tex
	$(HTMLBIN) $< -o $@


.SUFFIXES:
.SUFFIXES: .tex .bbl .dvi .pdf .bib .html .png .eps

# $(FIG_DIR)/%.png:
# 	$(if $(IMAGE_URL),,$(error IMAGE_URL is empty))
# 	$(CURL) $(IMAGE_PATH)/`basename $@` -o $@
# $(FIG_DIR)/%.pdf:
# 	$(if $(IMAGE_URL),,$(error IMAGE_URL is empty))
# 	$(CURL) $(IMAGE_PATH)/`basename $@` -o $@
# $(FIG_DIR)/%.eps:
# 	$(if $(IMAGE_URL),,$(error IMAGE_URL is empty))
# 	$(CURL) $(IMAGE_PATH)/`basename $@` -o $@

.PHONY : pdf clean html
