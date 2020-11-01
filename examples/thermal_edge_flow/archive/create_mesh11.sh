#!/bin/bash


for split in $(seq 0 1 3)
do
  for exp5 in $(seq 12 1 12)
  do
    outname=study11_"$exp5"_"$split".h5
    echo $outname
    geoToH5 study11.geo "$outname" "-setnumber split $split -setnumber exp5 $exp5"
  done
done
