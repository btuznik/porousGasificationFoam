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
    Foam::radiation::heterogeneousP1

Description
    Works well for combustion applications where optical thickness, tau is
    large, i.e. tau = a*L > 3 (L = distance between objects)

    Assumes
     - all surfaces are diffuse
     - tends to over predict radiative fluxes from sources/sinks
       *** SOURCES NOT CURRENTLY INCLUDED ***

SourceFiles
    heterogeneousP1.C

\*---------------------------------------------------------------------------*/

#ifndef heterogeneousRadiationModelheterogeneousP1_H
#define heterogeneousRadiationModelheterogeneousP1_H

#include "heterogeneousRadiationModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{

/*---------------------------------------------------------------------------*\
                           Class heterogeneousP1 Declaration
\*---------------------------------------------------------------------------*/

class heterogeneousP1
:
    public heterogeneousRadiationModel
{
    // Private data

        //- Incident radiation / [W/m2]
        volScalarField G_;

        //- Total radiative heat flux [W/m2]
        volScalarField qr_;

        //- Absorption coefficient
        volScalarField a_;

        //- Absorption coefficient
        volScalarField as_;

        //- Absorption coefficient
        volScalarField borderAs_;

        //- Emission contribution
        volScalarField E_;

        //- Emission contribution
        dimensionedScalar borderL_;

        //- porosity field
        const volScalarField& porosityF_;

        //- porosity field
        const volScalarField& surfFI_;

        //- porosity field
        volScalarField surfF_;

	    //- field where porosity is
	    volScalarField whereIs_;

	    //- field where is no porosity
	    volScalarField whereIsNot_;

        //- enthalpy source term for solid
	    volScalarField solidSh_;


public:

    //- Runtime type information
    TypeName("heterogeneousP1");


    // Constructors

        //- Construct from components
        heterogeneousP1
        (
            const volScalarField& T,
            const volScalarField& porosityF,
            const volScalarField& surfF,
            const volScalarField& Ts
        );

        //- Disallow default bitwise copy construction
        heterogeneousP1(const heterogeneousP1&) = delete;


    //- Destructor
    virtual ~heterogeneousP1();


    // Member functions

        // Edit

            //- Solve radiation equation(s)
            void calculate();

            //- Read radiation properties dictionary
            bool read();


        // Access

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const;

            //- Source term component (constant)
            virtual tmp<DimensionedField<scalar, volMesh> > Ru() const;

            //- Source term component (constant)
            virtual tmp<volScalarField> solidSh() const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const heterogeneousP1&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
