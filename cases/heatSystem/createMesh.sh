#!/bin/bash

# Settings
gmshPath=/Applications/gmsh/Gmsh.app/Contents/MacOS/gmsh

${gmshPath} -setnumber p 2 -2 ring.geo
# ${gmshPath} -setnumber p 2 -2 -order 2 ring.geo
dolfin-convert ring.msh ring.xml

# gzip ring.xml # not necessary needed