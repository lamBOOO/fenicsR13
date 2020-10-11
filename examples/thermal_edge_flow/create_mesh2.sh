#!/bin/bash

for exp1 in $(seq -3 -1 -3)
do
  for exp2 in $(seq -6 -1 -6)
  do
    for exp3 in $(seq -6 -1 -10)
    do
      outname=study2_"$exp1"_"$exp2"_"$exp3".h5
      echo $outname
      geoToH5 beamchamber_unstructured_nonuniform.geo "$outname" "-setnumber exp1 $exp1 -setnumber exp2 $exp2 -setnumber exp3 $exp3"
    done
  done
done
