#!/bin/bash

for p in 3 4 5 6 7 8
do
  geoToH5 lid_u.geo lid_u"$p".h5 "-setnumber p $p"
  geoToH5 lid.geo lid"$p".h5 "-setnumber p $p"
done
