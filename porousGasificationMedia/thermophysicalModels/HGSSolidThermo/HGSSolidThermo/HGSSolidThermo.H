/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::HGSSolidThermo

Description
    Abstract base-class for fluid and solid thermodynamic properties

SourceFiles
    HGSSolidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef HGSSolidThermo_H
#define HGSSolidThermo_H

#include "volFields.H"
#include "runTimeSelectionTables.H"
#include "fvMesh.H"
#include "typeInfo.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "wordIOList.H"
#include "basicSolidMixture.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class HGSSolidThermo Declaration
\*---------------------------------------------------------------------------*/

class HGSSolidThermo
:
    public IOdictionary
{

protected:

    // Protected data

                 
        //- Mesh reference
        const fvMesh& mesh_;
        
        // Fields

             //- Temperature [K]
            volScalarField T_;

            //- Density [kg/m3]
            volScalarField rho_;

            //- Absorption coefficient [1/m]
            volScalarField kappa_;

            //- Scatter coeffcient [1/m]
            volScalarField sigmaS_;

            //- Emissivity coeffcient []
            volScalarField emissivity_;
            
            
    // Protected Member Functions

        //- Construct as copy (not implemented)
        HGSSolidThermo(const HGSSolidThermo&);

       
public:

    //- Runtime type information
    TypeName("HGSSolidThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        HGSSolidThermo,
        mesh,
        (const fvMesh& mesh),
        (mesh)
    );
    
    declareRunTimeSelectionTable
    (
        autoPtr,
        HGSSolidThermo,
        dictionary,
        (const fvMesh& mesh, const dictionary& dict),
        (mesh, dict)
    );
    
    declareRunTimeSelectionTable
    (
        autoPtr,
        HGSSolidThermo,
        gasPhase,
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const PtrList<volScalarField>& gasPhaseGases
        ),
        (mesh, dict, gasPhaseGases)
    );

    // Constructors

        //- Construct from mesh 
        HGSSolidThermo(const fvMesh&);

        //- Construct from mesh and dictionary and phase name
        HGSSolidThermo(const fvMesh&, const dictionary&);
        
        //- Construct from mesh, dictionary and phase name
        HGSSolidThermo
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const PtrList<volScalarField>& gasPhaseGases
        );

    // Selectors
    
        //- Return a pointer to a new basicHGSSolidThermo created from
        //  the solidThermophysicalProperties dictionary
        static autoPtr<HGSSolidThermo> New(const fvMesh&);

        //- Return a pointer to a new basicHGSSolidThermo created from
        //  local dictionary
        static autoPtr<HGSSolidThermo> New(const fvMesh&, const dictionary&);
        
        static autoPtr<HGSSolidThermo> New
        (            const fvMesh& mesh,
            const dictionary& dict,
            const PtrList<volScalarField>& gasPhaseGases
        );

    //- Destructor
    virtual ~HGSSolidThermo();


    // Member Functions

        //- Update properties
        virtual void correct() = 0;
              
        //- Return the composition of the solid mixture
        virtual basicSolidMixture& composition();

        //- Return the composition of the solid mixture
        virtual const basicSolidMixture& composition() const;


        // Access to thermodynamic state variables

           //- mesh
            const fvMesh& mesh() const;

            //- Density [kg/m^3]
            const volScalarField& rho() const;

            //- Thermal conductivity [W/m/K]
            virtual const volScalarField& K() const;

            //- Thermal conductivity [W/m/K]
            virtual const volSymmTensorField& directionalK() const;

            //- Absorption coefficient [1/m]
            const volScalarField& kappa() const;

            //- Emissivity coefficient
            const volScalarField& sigmaS() const;

            //- Emissivity coefficient [1/m]
            const volScalarField& emissivity() const;

            //- Temperature [K]
            const volScalarField& T() const;

            //- non-const access for T
            volScalarField& T();

            //- non-const access for rho
            volScalarField& rho();


        // Derived thermal properties

            //- Specific heat capacity [J/kg/K]
            virtual tmp<volScalarField> Cp() const = 0;

            //- Sensible enthalpy [J/kg]
            virtual tmp<volScalarField> hs() const;

           //- Heat of formation [J/kg]
            virtual tmp<volScalarField> Hf() const = 0;


        // Per patch calculation

            //- Density [kg/m3]
            virtual tmp<scalarField> rho(const label patchI) const = 0;

            //- Specific heat capacity [J/kg/K)]
            virtual tmp<scalarField> Cp(const label patchI) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual tmp<scalarField> hs(const label patchI) const;

            //- Thermal conductivity [W//m/K]
            virtual tmp<scalarField> K(const label patchI) const;

            //- Thermal conductivity [W//m/K]
            virtual tmp<symmTensorField> directionalK(const label) const;

            //- Heat of formation [J/kg]
            virtual tmp<scalarField> Hf(const label patchI) const = 0;

            //- Scatter coefficient [1/m]
            virtual tmp<scalarField> sigmaS(const label) const = 0;

            //- Absorption coefficient [1/m]
            virtual tmp<scalarField> kappa(const label) const = 0;

            //- Emissivity []
            virtual tmp<scalarField> emissivity(const label) const = 0;

    // I-O

        //- Write the HGSSolidThermo properties
        virtual bool writeData(Ostream& os) const = 0;

        //- Read thermophysicalProperties dictionary
        virtual bool read() = 0;

        //- Ostream Operator
        friend Ostream& operator<<(Ostream& os, const HGSSolidThermo& s);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
