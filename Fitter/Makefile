##############################################################################
#
#  Makefile for Hydra library libKineRefit.so
#
#  This makefile contains all definitions local to this module. All
#  general definitions are included from makefiles named "hydra.*.mk".
#
##############################################################################


LIB_NAME := KineRefit

USES_ORACLE= : no

SOURCE_FILES := hrefitcand.cc \
		hfitter.cc


include $(HADDIR)/hades.def.mk

# set this, while debugging
SO_CXX_FLAGS += -O0

include $(HADDIR)/hades.module.mk
