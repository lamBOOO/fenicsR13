#!/bin/bash


for p in $(seq 0 1 5)
do
  outname=simple_mesh_"$p".h5
  echo $outname
  geoToH5 simple_mesh.geo "$outname" "-setnumber p $p"
done
