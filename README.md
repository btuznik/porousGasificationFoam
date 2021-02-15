# porousGasificationFoam

Table of content:
1. [Installation guide](#installation)
2. [Source guide](#source)
3. [Documentation](#doc)

<a name="installation"></a>
## Installation guide

The installation guide is prepared under assumption 
that OpenFOAM is installed in standard location: `/opt/OpenFoam-8/`

To install the 'porousGasificationMedia':
1. Set the OpenFOAM environmental paths by typing (modify the commend if your
   OpenFOAM is installed elsewhere):

    `$ source /opt/OpenFoam-8/etc/bashrc`

2. check the environmental settings
    * run any OpenFOAM solver, e.g. icoFoam -help
    * type: `$ echo $WM_PROJECT_USER_DIR`
   

3. Optionally change the destination path. The default path is: `$WM_PROJECT_USER_DIR/`.
   To change destination path edit file porousGasificationMediaDirectories located
   in the library installation folder.

4. Set the package environment by typing:

    `$ source porousGasificationMediaDirectories`

5. Run the install script:

   `$ . ./install`

6. Test the installation by running the solver:

    `$ porousGasificationFoam`

In case of errors:
1. Make sure the OpenFoam is correctly installed,
   by testing one of the tutorial cases provided with OpenFOAM. 
   For more information on installing OpenFoam vist [OpenFOAM Foundation website](https://openfoam.org/version/8/).
   
2. Open new terminal and again set the nessesary paths:

    `$ source source /opt/OpenFoam-8/etc/bashrc`\
    `$ source <PATH TO porousGasificationMediaDirectories>/porousGasificationMediaDirectiories`

3. Check the setting:

    run: `icoFoam -help` (output should be 'usage info')\
    `echo $WM_PROJECT_USER_DIR (typically $HOME/OpenFoam/<user_name>-8/`\
    `echo $FOAM_HGS` (typically `$HOME/OpenFoam/<user_name>-7/porousGasificationMedia`)\
    
    If missing set these variables manually.
4. Make sure all required packages have been installed by typing:

    `$ sudo apt-get install build-essential flex bison cmake zlib1g-dev qt4-dev-tools libqt4-dev gnuplot libreadline-dev libncurses-dev libxt-dev`
    
5. Test the ld version: `ld -v`
   The 'GNU ld' is required. If the 'GNU gold' is installed remove it typing:
   
    `$ sudo apt-get remove --purge binutils-gold`
   
6. Run install script:

    `./install > log_install &`

7. Other possibility:
    * clean up the library settings: unset `$LD_LIBRARY_PATH`
      and repeat steps 1-3 and 6.

<a name="source"></a>
## Source guide

###  Installation part

Files for installation and sourcing paths:

* `./README.md` -- readme file

* `./porousGasificationMediaDirectories` -- file adding paths to installation folder and libraries

* `./install` -- installation script

###  setPorosity

Utility for creating porosity fields:
*  porosityF - porosity field
*  Df - Darcy porous resitance tensor


```
# Initialises the fields reads if needed.
./setPorosity/setPorosity.C 
./setPorosity/createFields.H

# User editable part for introducing porosity fields.
./setPorosity/medium.H

# Wmake files.
./setPorosity/Make/options
./setPorosity/Make/files
```

###  porousGasificationFoam solver

 porousGasificationFoam solver main code, that uses porousGasificationMedia library.
 All calculations are scheduled here.

```
# Main solver file.
./porousGasificationFoam/porousGasificationFoam.C

# The solver executable source code.
./porousGasificationFoam/createFields.H
./porousGasificationFoam/readPyrolysisTimeControls.H
./porousGasificationFoam/setMultiRegionDeltaT.H
./porousGasificationFoam/createPyrolysisModel.H
./porousGasificationFoam/createPorosity.H
./porousGasificationFoam/createHeterogeneousRadiationModel.H
./porousGasificationFoam/chemistry.H
./porousGasificationFoam/readChemistryProperties.H
./porousGasificationFoam/UEqn.H
./porousGasificationFoam/rhoEqn.H
./porousGasificationFoam/solidRegionDiffusionNo.H
./porousGasificationFoam/pEqn.H
./porousGasificationFoam/hsEqn.H
./porousGasificationFoam/YEqn.H
./porousGasificationFoam/radiation.H

# Wmake files.
./porousGasificationFoam/Make/options
./porousGasificationFoam/Make/files

```

### porousGasificationMedia library

porousGasificationMedia library inculding three major parts:
1. fieldPorosityModel -- implementation of mechanical properties of porous medium.
2. thermophysicalModels -- implementation of thermophysical and chemical properties od porous medium.
3. pyrolysisModels -- classes that evaluate porous medium state and properties.
4. radiationModels


####  fieldPorosityModel

Mechanical properties of porous medium, can be used as standalone library

```
# Base class for mechanical properties of porous zones.
./porousGasificationMedia/fieldPorosityModel/fieldPorosityModelTemplates.C

# Darcy law porous medium implementation.
./porousGasificationMedia/fieldPorosityModel/fieldPorosityModel.H
./porousGasificationMedia/fieldPorosityModel/fieldPorosityModel.C

# wmake files.
./porousGasificationMedia/fieldPorosityModel/Make/options
./porousGasificationMedia/fieldPorosityModel/Make/files
```


####  thermophyscialModels

Script building all files in thermophysicalModels library:

`./porousGasificationMedia/thermophysicalModels/Allwmake`

##### solid

```
# Base class for heat transfer model.
./porousGasificationMedia/thermophysicalModels/solid/heatTransfer/heatTransferModel/heatTransferModel.C
./porousGasificationMedia/thermophysicalModels/solid/heatTransfer/heatTransferModel/heatTransferModel.H

# Cylinder pores based heat transfer model
./porousGasificationMedia/thermophysicalModels/solid/heatTransfer/cylinder/cylinder.C
./porousGasificationMedia/thermophysicalModels/solid/heatTransfer/cylinder/cylinder.H

# Constant heat transfer model
./porousGasificationMedia/thermophysicalModels/solid/heatTransfer/const/const.H
./porousGasificationMedia/thermophysicalModels/solid/heatTransfer/const/const.C
```

####  pyrolysisModels
Libraries where properites of solid phase are calculated

```
# Base class for heterogenous pyrolysis models.
./porousGasificationMedia/pyrolysisModels/pyrolysisModel/heterogeneousPyrolysisModel/heterogeneousPyrolysisModel.C
./porousGasificationMedia/pyrolysisModels/pyrolysisModel/heterogeneousPyrolysisModel/heterogeneousPyrolysisModelI.H
./porousGasificationMedia/pyrolysisModels/pyrolysisModel/heterogeneousPyrolysisModel/heterogeneousPyrolysisModel.H
./porousGasificationMedia/pyrolysisModels/pyrolysisModel/heterogeneousPyrolysisModel/heterogeneousPyrolysisModelNew.C

# Implementation of volumetric prolysis model volPyrolysis.
./porousGasificationMedia/pyrolysisModels/pyrolysisModel/volPyrolysis/volPyrolysisI.H
./porousGasificationMedia/pyrolysisModels/pyrolysisModel/volPyrolysis/volPyrolysis.H
./porousGasificationMedia/pyrolysisModels/pyrolysisModel/volPyrolysis/volPyrolysis.C

# wmake files.
./porousGasificationMedia/pyrolysisModels/Make/options
./porousGasificationMedia/pyrolysisModels/Make/files
```

##### radiationModels

<a name="doc"></a>
# Documentation

The documentation of the model can be generated with Doxygen.
To generate documentation:\
`$ doxygen doc/Doxyfile.in`

In order to view to view the documentation open `./doc/html/index.html` in a desired internet browser.
