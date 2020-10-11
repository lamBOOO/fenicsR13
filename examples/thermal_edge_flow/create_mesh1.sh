#!/bin/bash
#geoToH5 beamchamber.geo beamchamber4.h5 "-setnumber p 4"
# geoToH5 beamchamber_structured.geo beamchamber_structured4.h5 "-setnumber p 4"
# geoToH5 beamchamber_structured_2.geo beamchamber_structured_2.h5 "-setnumber p 4"
# geoToH5 beamchamber_structured.geo beamchamber_structured.h5 "-setnumber exp 7 nnodes 7"

# for exp in 2 3 4 5 6 7
# do
#   geoToH5 beamchamber_structured.geo "$exp".h5 "-setnumber exp $exp -setnumber nnodes 15"
# done

# for exp in $(seq 7 1 7)
# do
#   for nnodes in $(seq 22 2 24)
#   do
#     geoToH5 beamchamber_structured.geo "$exp"_"$nnodes".h5 "-setnumber exp $exp -setnumber nnodes $nnodes"
#   done
# done

for exp in $(seq 2 1 7)
do
  for nnodes in $((10 + 7 - $exp))
  do
    geoToH5 beamchamber_structured.geo "$exp"_"$nnodes".h5 "-setnumber exp $exp -setnumber nnodes $nnodes"
  done
done
