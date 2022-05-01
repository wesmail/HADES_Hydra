#!/bin/bash

#jobscript=$1
jobarray=$1
outprefix=$2

jobscript=$(pwd)/submitJob.sh

singularity exec \
-B /cvmfs/hadessoft.gsi.de/install/debian8/install:/cvmfs/hades.gsi.de/install \
-B /cvmfs/hadessoft.gsi.de/install/debian8/oracle:/cvmfs/it.gsi.de/oracle \
-B /lustre \
/cvmfs/vae.gsi.de/debian8/containers/user_container-production.sif  ${jobscript} ${jobarray} ${outprefix}
