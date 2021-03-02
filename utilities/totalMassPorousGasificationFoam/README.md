# totalMassPorousGasificationFoam -- User manual.

1. run totalMassPorousGasificationFoam inside the case directory:
    
    `totalMassPorousGasificationFoam`
    
2. the file 'totalMass.xy' will be created in the case directory. It contains two
   columns: time and the solid state mass over the whole computational domain.
   
3. the data can be visualised e.g. with gnuplot. Run gnuplot in the terminal
   and type:
   
    `"plot "totalMass.xy" using 1:2`
