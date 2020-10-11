#!/bin/bash

for split in $(seq 0 1 3)
do
  outname=study3_"$split".h5
  echo $outname
  geoToH5 beamchamber_structured_split.geo "$outname" "-setnumber split $split"
done
