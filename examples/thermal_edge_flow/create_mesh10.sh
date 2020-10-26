#!/bin/bash


for split in $(seq 3 1 3)
do
  for exp5 in $(seq 8 1 11)
  do
    outname=study10_"$exp5"_"$split".h5
    echo $outname
    geoToH5 study10.geo "$outname" "-setnumber split $split -setnumber exp5 $exp5"
done
