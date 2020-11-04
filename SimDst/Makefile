APP_NAME      := analysisDST

SOURCE_FILES  := analysisDST.cc
USES_RFIO     := no
USES_ORACLE   := yes
USES_GFORTRAN := yes

include $(HADDIR)/hades.def.mk

HYDRA_LIBS    += -lDst

PLUTODIR := ""

.PHONY:  default
default: clean build install

include $(HADDIR)/hades.app.mk

