#!/bin/bash

# Settings
gmshPath=/Applications/gmsh/Gmsh.app/Contents/MacOS/gmsh

${gmshPath} -2 ring.geo
dolfin-convert ring.msh ring.xml

# gzip ring.xml # not necessary needed