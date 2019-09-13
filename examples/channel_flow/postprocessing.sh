#!/bin/bash

for KNUDSEN in 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0
do
  MASSFLOW=$(cat channel_flow$KNUDSEN/massflow_3002)
  echo "$KNUDSEN, $MASSFLOW" >> table.csv
done

cat table.csv