#!/bin/bash

for pair in "0 2" "1 4" "2 8" "3 16" "4 32" "5 64" "6 128" "7 256" "8 512" "9 1024"; do
  read -r i n <<< "$pair"   # splits on IFS (space by default)
  mpirun -n $n fenicsR13 ./input_$i.yml
done

# cleanup for sending
rm out_*/*.h5
rm out_*/*.xdmf
rm out_*/solve*
rm out_*/assemble*
rm out_*/avgvel*
