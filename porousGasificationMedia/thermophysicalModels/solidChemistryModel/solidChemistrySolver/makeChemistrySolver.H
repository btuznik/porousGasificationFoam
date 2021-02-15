/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef makeChemistrySolver_H
#define makeChemistrySolver_H

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineChemistrySolver(Model, SolidThermo, SolidThermoPhysics, GasThermoPhysics)    \
                                                                               \
    typedef Model<SolidThermo, SolidThermoPhysics, GasThermoPhysics>                       \
        Model##SolidThermo##SolidThermoPhysics##GasThermoPhysics;                          \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Model##SolidThermo##SolidThermoPhysics##GasThermoPhysics,                          \
        (                                                                      \
            word(Model##SolidThermo##SolidThermoPhysics##GasThermoPhysics::typeName_()) + "<"  \
          + SolidThermo::typeName + "," + SolidThermoPhysics::typeName() + ","            \
          + GasThermoPhysics::typeName() + ">"                                    \
        ).c_str(),                                                             \
        0                                                                      \
    );


#define makeChemistrySolver(Solver, Model, SolidThermo, SolidThermoPhysics, GasThermoPhysics) \
                                                                               \
    typedef Solver<Model<SolidThermo, SolidThermoPhysics, GasThermoPhysics>>               \
        Solver##Model##SolidThermo##SolidThermoPhysics##GasThermoPhysics;                  \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Solver##Model##SolidThermo##SolidThermoPhysics##GasThermoPhysics,                  \
        (                                                                      \
            word(Solver##Model##SolidThermo##SolidThermoPhysics##GasThermoPhysics::typeName_()) \
          + "<" + word(Model<SolidThermo, SolidThermoPhysics, GasThermoPhysics>::typeName_())   \
          + "<" + word(SolidThermo::typeName_()) + "," + SolidThermoPhysics::typeName() + "," \
          + GasThermoPhysics::typeName() + ">>"                                   \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    BasicSolidChemistryModel<SolidThermo>::                                    \
        addthermoConstructorToTable                                            \
        <Solver##Model##SolidThermo##SolidThermoPhysics##GasThermoPhysics>                 \
        add##Solver##Model##SolidThermo##SolidThermoPhysics##GasThermoPhysics##\
thermoConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
