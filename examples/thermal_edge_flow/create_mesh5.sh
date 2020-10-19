#!/bin/bash


for split in $(seq 0 1 3)
do
  outname=study5_"$split".h5
  echo $outname
  geoToH5 beamchamber_hybrid_curved.geo "$outname" "-setnumber split $split"
done
