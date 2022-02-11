#!/bin/sh

cp -r 0.orig 0
cp initialData/pyrolysisProperties_hot constant/pyrolysisProperties
cp initialData/controlDict_hot system/controlDict
cp initialData/controlDict_hot_stat system/controlDict
cp initialData/setFieldsDict_hot system/setFieldsDict

blockMesh -dict system/blockMeshDict
setSet -batch makeFaceSet.setSet
createPatch -overwrite
setFields
cp -r initialData/0.pre/* 0/

decomposePar

# Run
# porousGasificationFoam
mpirun -np 4 porousGasificationFoam -parallel > log&

# -----------------------------------------------------------------------------
