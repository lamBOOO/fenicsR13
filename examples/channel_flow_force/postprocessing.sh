#!/bin/bash

rm table.csv

for KNUDSEN in 0.03125 0.0625 0.125 0.25 0.5 1.0 2.0
do
  MASSFLOW=$(cat results_channel_flow_force$KNUDSEN/massflow_3002)
  echo "$KNUDSEN, $MASSFLOW" >> table.csv
done

cat table.csv
