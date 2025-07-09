#!/bin/bash

for pair in "0 2" "1 4" "2 8" "3 16" "4 32" "5 64" "6 128" "7 256" "8 512" "9 1024"; do
# for pair in "0 2"; do
  read -r i n <<< "$pair"   # splits on IFS (space by default)
  mpirun -n $n fenicsR13 ./input3d_$i.yml
done

# cleanup for sending
rm out3d_*/*.h5
rm out3d_*/*.xdmf
rm out3d_*/solve*
rm out3d_*/assemble*
rm out3d_*/avgvel*
