#!/bin/bash


for exp3 in $(seq 10 1 15)
do
  outname=study6_"$exp3".h5
  echo $outname
  geoToH5 study6.geo "$outname" "-setnumber exp3 $exp3"
done
