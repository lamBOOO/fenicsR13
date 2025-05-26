#!/bin/bash

for p in 0 1 2 3 4 5 6
do
  geoToH5 channel.geo channel"$p".h5 "-setnumber p $p"
done
