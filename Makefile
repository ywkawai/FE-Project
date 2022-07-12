FELIBDIR   = FElib/src
MODEL_GLOBALSW_DIR = model/global_shallow_water/src
MODEL_ATM2D_DIR = model/atm_nonhydro2d/src
MODEL_ATM3D_DIR = model/atm_nonhydro3d/src
MODEL_ATM3D_UTIL_REGRID = model/atm_nonhydro3d/util/regrid_tool

all: build

build: 
	$(MAKE) -C $(FELIBDIR) build
	$(MAKE) -C $(FELIBDIR) install
	$(MAKE) build-model-globalsw build-model-atm2d build-model-atm3d build-model-atm3d-util

install: install-model-globalsw install-model-atm2d install-model-atm3d install-model-atm3d-util

clean:
	$(MAKE) -C $(MODEL_GLOBALSW_DIR) clean
	$(MAKE) -C $(MODEL_ATM2D_DIR) clean
	$(MAKE) -C $(MODEL_ATM3D_DIR) clean
	$(MAKE) -C $(MODEL_ATM3D_UTIL_REGRID) clean
	$(MAKE) -C $(FELIBDIR) clean

distclean:
	$(MAKE) -C $(MODEL_GLOBALSW_DIR) distclean
	$(MAKE) -C $(MODEL_ATM2D_DIR) distclean
	$(MAKE) -C $(MODEL_ATM3D_DIR) distclean
	$(MAKE) -C $(MODEL_ATM3D_UTIL_REGRID) distclean
	$(MAKE) -C $(FELIBDIR) distclean

#--
build-model-globalsw:
	$(MAKE) -C $(MODEL_GLOBALSW_DIR) build

build-model-atm2d:
	$(MAKE) -C $(MODEL_ATM2D_DIR) build

build-model-atm3d:
	$(MAKE) -C $(MODEL_ATM3D_DIR) build

build-model-atm3d-util:
	$(MAKE) -C $(MODEL_ATM3D_UTIL_REGRID) build

#--
install-model-globalsw:
	$(MAKE) -C $(MODEL_GLOBALSW_DIR) install

install-model-atm2d:
	$(MAKE) -C $(MODEL_ATM2D_DIR) install

install-model-atm3d:
	$(MAKE) -C $(MODEL_ATM3D_DIR) install

install-model-atm3d-util:
	$(MAKE) -C $(MODEL_ATM3D_UTIL_REGRID) install

