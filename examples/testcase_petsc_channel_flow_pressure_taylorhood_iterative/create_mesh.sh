#!/bin/bash

for p in 8 9 10 11 12 13 14 15 16 17
do
  geoToH5 channel.geo channel_"$p".h5 "-setnumber p $p"
  geoToH5 channel3d.geo channel3d_"$p".h5 "-3 -setnumber p $p"
done
