#!/bin/bash

if [ $# -lt 3 ]; then
  echo -e "\nJob script for submission of analysis jobs on KRONOS.\n"
  echo -e "USAGE: . sendScript.sh <array_stop> <in_list> <out_prefix>\n"
  sleep 10
  exit 1
fi

stop=$1
inlist=$2
outprefix=$3

jobscript=$(pwd)/submitJob.sh
jobarray=$(pwd)/$inlist

wrap=$(pwd)/wrap.sh
command="--array=1-${stop} -- ${wrap} ${jobarray} ${outprefix}"
sbatch $command
