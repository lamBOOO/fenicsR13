#!/bin/bash

# generate meshes
cd ../mesh
python3 create_meshes.py

# run simulation
cd ../heat
python3 ../../src/fenicsR13.py

# compare errors
diff -u errors.csv errors_ref.csv # returns exit code 1 if change. echo $?