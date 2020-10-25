#!/bin/bash


for exp1 in $(seq 4 1 4)
do
  for exp2 in $(seq 4 1 4)
  do
    for exp3 in $(seq 11 1 11)
    do
      for split in $(seq 0 1 2)
      do
        outname=study7_"$exp1"_"$exp2"_"$exp3"_"$split".h5
        echo $outname
        geoToH5 study7.geo "$outname" "-setnumber exp1 $exp1 -setnumber exp2 $exp2 -setnumber exp3 $exp3 -setnumber split $split"
      done
    done
  done
done
