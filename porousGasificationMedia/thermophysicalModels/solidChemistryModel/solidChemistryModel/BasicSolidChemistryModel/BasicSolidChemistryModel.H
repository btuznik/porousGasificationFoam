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

template<class SolidThermo>
class BasicSolidChemistryModel
:
    public basicSolidChemistryModel
{
protected:

    // Protected data

        //- Thermo
        const SolidThermo& solidThermo_;


public:

    //- Runtime type information
    TypeName("BasicSolidChemistryModel");


    //- Thermo type
    typedef SolidThermo reactionThermo;


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        BasicSolidChemistryModel,
        thermo,
        (
            const SolidThermo& thermo,
            PtrList<volScalarField>& gasPhaseGases
        ),
        (thermo, gasPhaseGases)
    );


    // Constructors

        //- Construct from thermo
        BasicSolidChemistryModel(const SolidThermo& thermo);

    //- Selector
    static autoPtr<BasicSolidChemistryModel<SolidThermo>> New
    (
        const SolidThermo& thermo,
        PtrList<volScalarField>& gasPhaseGases,
        const word thermoName
    );


    //- Destructor
    virtual ~BasicSolidChemistryModel();


    // Member Functions

        //- Return access to the solid thermo package
        //inline SolidThermo& solidThermo();

        //- Return access to the solid's thermo package
        inline const SolidThermo& solidThermo() const;

        //- Return total solids mass source term [kg/m3/s]
        virtual tmp<volScalarField> RRs() const = 0;

        //- Return total gases mass source term [kg/m3/s]
        virtual tmp<volScalarField> RRg() const = 0;

        //- Return chemical source terms for solids [kg/m3/s]
        virtual tmp<volScalarField> RRs(const label i) const = 0;

        //- Return chemical source terms for gases [kg/m3/s]
        virtual tmp<volScalarField> RRg(const label i) const = 0;

        //- Return sensible enthalpy for gas i [J/Kg]
        virtual tmp<volScalarField> gasHs
        (
            const volScalarField& p,
            const volScalarField& T,
            const label i
        ) const = 0;

        //- Return source for enthalpy equation [kg/m/s3]
        virtual tmp<volScalarField> Sh() const = 0;

        //- Return specie Table for gases
        virtual const speciesTable& gasTable() const = 0;

        //- Set reacting status of cell, cellI
        virtual void setCellReacting(const label cellI, const bool active) = 0;

        //- Calculates the reaction rates
        virtual void calculate() = 0;

        //- Return source for porosity equation [1/s]
        virtual tmp<volScalarField> RRpor(const volScalarField T) const = 0;

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
