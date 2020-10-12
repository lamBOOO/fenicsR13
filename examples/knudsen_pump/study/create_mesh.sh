#!/bin/bash
for p in $(seq 4 9)
do
  geoToH5 ../knudsen_pump.geo knudsen_pump"$p".h5 "-setnumber p $p"
done
