#!/bin/bash

for p in 1 2 3 4 5
do
  for name in shell
  do
    geoToH5 "$name".geo "$name""$p".h5 "-3 -setnumber p $p"
  done
done
