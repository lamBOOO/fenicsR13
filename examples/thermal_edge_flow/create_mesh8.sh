#!/bin/bash


for exp1 in $(seq 6 1 6)
do
  for exp2 in $(seq 6 1 6)
  do
    for exp3 in $(seq 8 1 13)
    do
      for split in $(seq 0 1 0)
      do
        outname=study8_"$exp1"_"$exp2"_"$exp3"_"$split".h5
        echo $outname
        geoToH5 study8.geo "$outname" "-setnumber exp1 $exp1 -setnumber exp2 $exp2 -setnumber exp3 $exp3 -setnumber split $split"
      done
    done
  done
done
