#!/bin/bash


for split in $(seq 0 1 3)
do
  outname=study9_"$split".h5
  echo $outname
  geoToH5 study9.geo "$outname" "-setnumber split $split"
done
