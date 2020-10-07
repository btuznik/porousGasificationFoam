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
    Foam::interpolatedHGSSolidThermo

Description
    Table interpolated solid thermo

SourceFiles
    interpolatedHGSSolidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef interpolatedHGSSolidThermo_H
#define interpolatedHGSSolidThermo_H

#include "HGSSolidThermo.H"
#include "interpolateSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class interpolatedHGSSolidThermo;

Ostream& operator<<(Ostream&, const interpolatedHGSSolidThermo&);

/*---------------------------------------------------------------------------*\
                  Class interpolatedHGSSolidThermo Declaration
\*---------------------------------------------------------------------------*/

class interpolatedHGSSolidThermo
:
    public HGSSolidThermo,
    public interpolateSolid
{
private:

    // Private data

        //- Dictionary
        const dictionary dict_;


public:

    // Constructors

        //- Construct from mesh
        interpolatedHGSSolidThermo(const fvMesh& mesh, const word);

        //- Construct from mesh and dictionary
        interpolatedHGSSolidThermo
        (
            const fvMesh& mesh,
            const word,
            const dictionary& dict
        );



    //- Destructor
    virtual ~interpolatedHGSSolidThermo();


    // Member Functions

        //- Calculate properties
        void calculate();


        // Derived properties

            //- Specific heat capacity [J/(kg.K)]
            virtual tmp<volScalarField> Cp() const;

            //- Heat of formation [J/kg]
            virtual tmp<volScalarField> Hf() const;


        // Per patch calculation

            //- Density [kg/m3]
            virtual tmp<scalarField> rho(const label patchI) const;

            //- Specific heat capacity [J/kg/K)]
            virtual tmp<scalarField> Cp(const label patchI) const;

            //- Heat of formation [J/kg]
            virtual tmp<scalarField> Hf(const label patchI) const;

            //- Scatter coefficient [1/m]
            virtual tmp<scalarField> sigmaS(const label) const;

            //- Absorption coefficient [1/m]
            virtual tmp<scalarField> kappa(const label) const;

             //- Emissivity []
            virtual tmp<scalarField> emissivity(const label) const;


        // I-O

            //- Write the interpolatedHGSSolidThermo properties
            virtual bool writeData(Ostream& os) const;

            //- Read the interpolatedHGSSolidThermo properties
            virtual bool read();

            //- Read the interpolatedHGSSolidThermo properties
            bool read(const dictionary& dict);

            //- Ostream Operator
            friend Ostream& operator<<
            (
                Ostream&,
                const interpolatedHGSSolidThermo&
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
