#!/bin/bash


for p in $(seq 0 1 6)
do
  outname=simple_mesh_"$p".h5
  echo $outname
  geoToH5 simple_mesh.geo "$outname" "-setnumber p $p"

  outname=structured_simple_mesh_"$p".h5
  echo $outname
  geoToH5 structured_simple_mesh.geo "$outname" "-setnumber p $p"

  outname=structured_simple_mesh_flipped_"$p".h5
  echo $outname
  geoToH5 structured_simple_mesh_flipped.geo "$outname" "-setnumber p $p"

  outname=structured_simple_mesh_symmetric_"$p".h5
  echo $outname
  geoToH5 structured_simple_mesh_symmetric.geo "$outname" "-setnumber p $p"
done
