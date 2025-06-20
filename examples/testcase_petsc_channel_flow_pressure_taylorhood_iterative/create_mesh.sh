#!/bin/bash

for p in 8 9 10 11 12 13 14 15
do
  geoToH5 channel.geo channel"$p".h5 "-setnumber p $p"
done
