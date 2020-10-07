/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Foam::constHGSSolidThermo

Description
    The thermophysical properties of a constHGSSolidThermo

SourceFiles
    constHGSSolidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef constHGSSolidThermo_H
#define constHGSSolidThermo_H

#include "HGSSolidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class constHGSSolidThermo;

Ostream& operator<<(Ostream&, const constHGSSolidThermo&);

/*---------------------------------------------------------------------------*\
                      Class constHGSSolidThermo Declaration
\*---------------------------------------------------------------------------*/

class constHGSSolidThermo
:
    public HGSSolidThermo
{
private:

    // Private data

        //- Dictionary
        dictionary dict_;

        //- Constant thermal conductivity [W/(m.K)]
        dimensionedScalar constK_;

        //- Thermal conductivity field[W/(m.K)]
        volScalarField K_;

        //- Density [kg/m3]
        dimensionedScalar constRho_;

        //- Specific heat capacity [J/(kg.K)]
        dimensionedScalar constCp_;

        //- Heat of formation [J/kg]
        dimensionedScalar constHf_;

        //- Emissivity
        dimensionedScalar constEmissivity_;

        //- Absorptivity [1/m]
        dimensionedScalar constKappa_;

        //- Scatter [1/m]
        dimensionedScalar constSigmaS_;


public:

    //- Runtime type information
    TypeName("constHGSSolidThermo");


    // Constructors

        //- Construct from mesh
        constHGSSolidThermo(const fvMesh& mesh);

        //- Construct from mesh and dict
        constHGSSolidThermo(const fvMesh& mesh, const dictionary& dict);

	//- Construct from mesh, dict and gasPhaseGases
        constHGSSolidThermo(const fvMesh& mesh, const dictionary& dict, const PtrList<volScalarField>& gasPhaseGases);

    //- Destructor
    virtual ~constHGSSolidThermo();


    // Member Functions

        //- Update properties
        virtual void correct();


    // Acces functions

        //- Constant access to K
        virtual const volScalarField& K() const;

	//- zmiana non-const acces to K
	volScalarField& K();

        //- Thermal conductivity [W/(m.K)]
        virtual const volSymmTensorField& directionalK() const;


        // Derived properties

            //- Specific heat capacity [J/(kg.K)]
            virtual tmp<volScalarField> Cp() const;

            //- Heat of formation [J/kg]
            virtual tmp<volScalarField> Hf() const;


        // Per patch calculation

            //- Density [kg/m3]
            virtual tmp<scalarField> rho(const label patchI) const;
            virtual tmp<scalarField> rhoT(const label patchI);

            //- Specific heat capacity [J/kg/K)]
            virtual tmp<scalarField> Cp(const label patchI) const;

            //- Thermal conductivity [W//m/K]
            virtual tmp<scalarField> K(const label patchI) const;

            //- Thermal conductivity [W//m/K]
            virtual tmp<symmTensorField>directionalK(const label) const;

            //- Heat of formation [J/kg]
            virtual tmp<scalarField> Hf(const label patchI) const;

            //- Scatter coefficient [1/m]
            virtual tmp<scalarField> sigmaS(const label) const;

            //- Absorption coefficient [1/m]
            virtual tmp<scalarField> kappa(const label) const;

             //- Emissivity []
            virtual tmp<scalarField> emissivity(const label) const;


        // I-O

            //- Write the constHGSSolidThermo properties
            virtual bool writeData(Ostream& os) const;

            //- Read solidThermophysicalProperties dictionary
            virtual bool read();

            //- Read solidThermophysicalProperties dictionary
            bool read(const dictionary&);

            //- Ostream Operator
            friend Ostream& operator<<(Ostream&, const constHGSSolidThermo&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
