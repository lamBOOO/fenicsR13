#!/bin/bash
#geoToH5 beamchamber.geo beamchamber4.h5 "-setnumber p 4"
# geoToH5 beamchamber_structured.geo beamchamber_structured4.h5 "-setnumber p 4"
# geoToH5 beamchamber_structured_2.geo beamchamber_structured_2.h5 "-setnumber p 4"
# geoToH5 beamchamber_structured.geo beamchamber_structured.h5 "-setnumber exp 7 nnodes 7"

for exp in 2 3 4 5 6 7
do
  geoToH5 beamchamber_structured.geo "$exp".h5 "-setnumber exp $exp nnodes 5"
done
