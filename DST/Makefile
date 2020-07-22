APP_NAME      := analysisDST

SOURCE_FILES  := analysisDST.cc
USES_RFIO     := no
USES_ORACLE   := yes
USES_GFORTRAN := yes

include $(HADDIR)/hades.def.mk

#INC_DIRS += $(PLUTO_DIR)/include
#LIB_DIRS += $(PLUTO_DIR)
#SO_LOADER_DIRS  += $(PLUTO_DIR)
#APP_LOADER_DIRS += $(PLUTO_DIR)
#LIBS          += -lPluto

HYDRA_LIBS    += -lDst
#HYDRA_LIBS    += -lParticle


.PHONY:  default
default: clean build install

include $(HADDIR)/hades.app.mk

