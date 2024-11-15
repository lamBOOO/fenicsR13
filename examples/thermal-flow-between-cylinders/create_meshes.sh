#!/bin/bash

for p in 0 1 2 3 4 5 6 7
do
  for name in mesh
  do
    geoToH5 "$name".geo "$name""$p".h5 "-setnumber p $p"
  done
done

for p in 0 1 2 3 4 5 6 7
do
  for name in mesh_adaptive
  do
    geoToH5 "$name".geo "$name""$p".h5 "-setnumber p $p"
  done
done
