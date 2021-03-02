# Total mass -- User manual.

1. run totalMass inside the case directory.
2. the file 'totalMass.xy' will be created in the directory. It contains two
   columns: time and total mass of the solid species.
3. the data can be visualised e.g. with gnuplot. Run gnuplot in the terminal
   and type:
   
    `"plot "totalMass.xy" using 1:2`
