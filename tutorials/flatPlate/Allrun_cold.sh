#!/bin/sh

cp -r 0.orig 0
cp initialData/pyrolysisProperties_cold constant/pyrolysisProperties
cp initialData/controlDict_cold system/controlDict
cp initialData/setFieldsDict_cold system/setFieldsDict

blockMesh -dict system/blockMeshDict
setSet -batch makeFaceSet.setSet
createPatch -overwrite
setFields

decomposePar

# Run
# porousGasificationFoam
mpirun -np 4 porousGasificationFoam -parallel > log&

# -----------------------------------------------------------------------------
