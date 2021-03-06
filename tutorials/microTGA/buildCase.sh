#!/bin/bash

foamCleanPolyMesh
rm -r 0
cp -r save  0
blockMesh
setFields
renumberMesh -overwrite
checkMesh
