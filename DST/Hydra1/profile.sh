#!/bin/bash

. /cvmfs/hades.gsi.de/install/root-5.34.34/bin/thisroot.sh
. /cvmfs/hades.gsi.de/install/5.34.34/old/hydra-11042016/defall.sh
. /cvmfs/hades.gsi.de/install/5.34.34/old/hgeant-11042016/defall.sh

export HADDIR=/cvmfs/hades.gsi.de/install/5.34.34/old/hydra-11042016/
#export CERN_ROOT=/cvmfs/hades.gsi.de/install/cernlib_gfortran/2005/
export HGEANT_DIR=/cvmfs/hades.gsi.de/install/5.34.34/old/hgeant-11042016/

export ORA_USER=hades_ana/hades@db-hades
export ORACLE_HOME=/cvmfs/it.gsi.de/oracle/product/12.1.2/client_x86_64_1

export LC_ALL=C
#export ROOTLOGON=$HADDIR/macros/rootlogon.C
export PATH=${INSTALL_DIR}/bin:${MYHADDIR}/bin:${PATH}
export LD_LIBRARY_PATH=${INSTALL_DIR}/lib:${MYHADDIR}/lib:${HADDIR}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=/cvmfs/it.gsi.de/oracle/product/12.1.2/client_x86_64_1/lib/:${LD_LIBRARY_PATH}
