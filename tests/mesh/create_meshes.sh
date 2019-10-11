#!/bin/bash

for p in 0 1 2 3 4 5 6 7
do
  python3 ../../src/geoToH5.py ring.geo ring"$p".h5 "-setnumber p $p"
done
