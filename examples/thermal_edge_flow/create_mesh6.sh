#!/bin/bash


for exp1 in $(seq 4 1 5)
do
  for exp2 in $(seq 5 1 6)
  do
    for exp3 in $(seq 10 1 14)
    do
      outname=study6_"$exp1"_"$exp2"_"$exp3".h5
      echo $outname
      geoToH5 study6.geo "$outname" "-setnumber exp1 $exp1 -setnumber exp2 $exp2 -setnumber exp3 $exp3"
    done
  done
done
