#!/bin/bash


for p in $(seq 0 1 3)
do
  outname=unstructured_"$p".h5
  echo $outname
  geoToH5 unstructured.geo "$outname" "-setnumber p $p"

  outname=structured_left_"$p".h5
  echo $outname
  geoToH5 structured_left.geo "$outname" "-setnumber p $p"

  outname=structured_right_"$p".h5
  echo $outname
  geoToH5 structured_right.geo "$outname" "-setnumber p $p"

  outname=structured_sym_"$p".h5
  echo $outname
  geoToH5 structured_sym.geo "$outname" "-setnumber p $p"
done
