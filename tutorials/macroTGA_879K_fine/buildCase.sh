#!/bin/bash

foamCleanPolyMesh
rm -r 0
cp -r save879 0
blockMesh

#first refinement level
setSet -batch setSet.c0
refineHexMesh -overwrite c0

#second refinement level
setSet -batch setSet.c1
refineHexMesh -overwrite c1

setFields
renumberMesh -overwrite
checkMesh
