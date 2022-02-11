#!/bin/sh

cp -r 0.orig 0
blockMesh -dict system
setFields
decomposePar
mpirun -np 4 porousGasificationFoam -parallel > log&
