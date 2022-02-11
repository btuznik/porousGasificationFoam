#!/bin/bash

cp Df_and_porosityProperties/Df2.56e+10 p04/0/Df
cp Df_and_porosityProperties/Df2.56e+10 p08/0/Df
cp Df_and_porosityProperties/Df2.56e+10 p10/0/Df
cp Df_and_porosityProperties/Df2.56e+10 p12/0/Df
cp Df_and_porosityProperties/Df2.56e+10 p16/0/Df

cp Df_and_porosityProperties/porosityProperties_f3.21e+6 p04/constant/porosityProperties
cp Df_and_porosityProperties/porosityProperties_f3.21e+6 p08/constant/porosityProperties
cp Df_and_porosityProperties/porosityProperties_f3.21e+6 p10/constant/porosityProperties
cp Df_and_porosityProperties/porosityProperties_f3.21e+6 p12/constant/porosityProperties
cp Df_and_porosityProperties/porosityProperties_f3.21e+6 p16/constant/porosityProperties

cd p04
porousGasificationFoam > log&
cd ../p08
porousGasificationFoam > log&
cd ../p10
porousGasificationFoam > log&
cd ../p12
porousGasificationFoam > log&
cd ../p16
porousGasificationFoam > log&

