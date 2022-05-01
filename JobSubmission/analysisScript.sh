#!/bin/bash

input=$1
name=$2

root -l -b -q "Analysis.C(\"$input\", \"$name\")"
