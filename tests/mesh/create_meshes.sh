#!/bin/bash

for p in 0 1 2 3 4 5 6 7
do
  for name in ring ring_antisym
  do
    python3 ../../src/geoToH5.py "$name".geo "$name""$p".h5 "-setnumber p $p"
  done
done
