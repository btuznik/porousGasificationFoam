/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

Class
    Foam::BasicSolidChemistryModel

Description
    Basic chemistry model templated on thermodynamics

SourceFiles
    BasicChemistryModelI.H
    BasicSolidChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef BasicSolidChemistryModel_H
#define BasicSolidChemistryModel_H

#include "basicSolidChemistryModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     class BasicSolidChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class BasicSolidChemistryModel
:
    public basicSolidChemistryModel
{
protected:

    // Protected data

        //- Thermo
        const ReactionThermo& thermo_;


public:

    //- Runtime type information
    TypeName("BasicSolidChemistryModel");


    //- Thermo type
    typedef ReactionThermo reactionThermo;


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        BasicSolidChemistryModel,
        thermo,
        (const ReactionThermo& thermo),
        (thermo)
    );


    // Constructors

        //- Construct from thermo
        BasicSolidChemistryModel(const ReactionThermo& thermo);


    //- Selector
    static autoPtr<BasicSolidChemistryModel<ReactionThermo>> New
    (
        const ReactionThermo& thermo
    );


    //- Destructor
    virtual ~BasicSolidChemistryModel();


    // Member Functions

        //- Return const access to the thermo package
        inline const ReactionThermo& thermo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BasicSolidChemistryModel.C"
#endif

#include "BasicSolidChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
