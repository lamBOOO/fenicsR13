#!/bin/bash
for p in $(seq 4 8)
do
  geoToH5 ../knudsen_pump.geo knudsen_pump"$p".h5 "-setnumber p $p"
done
