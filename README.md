# porousGasificationFoam

Table of content:
1. [Installation guide](#installation)
2. [Source guide](#source)
3. [Documentation](#doc)

<a name="installation"></a>
## Installation guide

The installation guide is prepared under the assumption 
that OpenFOAM is installed in standard location: `/opt/OpenFoam-8/`

To install the 'porousGasificationFoam':
1. Set the OpenFOAM environmental paths by typing (modify the command if your
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

1. Make sure the OpenFoam and all required packaches have been correctly installed,
   by testing one of the tutorial cases provided with OpenFOAM. 
   For more information on installing OpenFoam vist [OpenFOAM Foundation website](https://openfoam.org/version/8/).
   
2. Open new terminal and again set the nessesary paths:

    `$ source source /opt/OpenFoam-8/etc/bashrc`\
    `$ source <PATH TO porousGasificationMediaDirectories>/porousGasificationMediaDirectiories`

3. Check the setting:

    run: `icoFoam -help` (output should be 'usage info')\
    `echo $WM_PROJECT_USER_DIR (typically $HOME/OpenFoam/<user_name>-8/`\
    `echo $FOAM_HGS` (typically `$HOME/OpenFoam/<user_name>-8/porousGasificationMedia`)\
    
    If missing set these variables manually.
    
4. Run instalaltion script and check for errors:

    `./install > log_install &`

<a name="source"></a>

## Source guide

###  Installation part

Files for installation and sourcing paths:

* `./README.md` -- readme file

* `./porousGasificationMediaDirectories` -- file with enviromental variables
                                            needed for the installation

* `./install` -- installation script

###  porousGasificationFoam solver

 porousGasificationFoam's main code, that uses porousGasificationMedia library.
 All calculations are scheduled here.

### porousGasificationMedia library

porousGasificationMedia library inculding four major parts:
1. fieldPorosityModel -- implementation of mechanical properties of porous medium.
2. thermophysicalModels -- implementation of thermophysical and chemical properties od porous medium.
3. pyrolysisModels -- classes that evaluate porous medium state and properties.
4. radiationModels

###  Utilities

1. setPorosity -- utility for creating porosity fields:
    *  porosityF - porosity field
    *  Df - Darcy porous resitance tensor
2. totalMass -- calculates the total mass loss over time and writes to a file.

<a name="doc"></a>
# Documentation

The documentation of the model can be generated with Doxygen software. To build
the documentation the doxygen and graphviz packages are required.
For Ubuntu users the packages can be obtained with the following command:
     `sudo apt-get install doxygen graphviz`

* Running Doxygen
  In the `$POROUS_DOC_SRC/doc/Doxygen` directory type:
    `./Allwmake`
  
In order to view to view the documentation open
`$WM_PROJECT_DIR/doc/Doxygen/html/index.html` in a desired internet browser.
