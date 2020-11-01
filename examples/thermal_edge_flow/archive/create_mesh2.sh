#!/bin/bash

for p1 in $(seq -3 -1 -3)
do
  for p2 in $(seq -6 -1 -6)
  do
    for p3 in $(seq -6 -1 -11)
    do
      outname=study2_"$p1"_"$p2"_"$p3".h5
      echo $outname
      geoToH5 beamchamber_unstructured_nonuniform.geo "$outname" "-setnumber p1 $p1 -setnumber p2 $p2 -setnumber p3 $p3"
    done
  done
done
