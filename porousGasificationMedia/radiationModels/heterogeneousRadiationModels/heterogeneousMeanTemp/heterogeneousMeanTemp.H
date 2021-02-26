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
    Foam::radiation::heterogeneousMeanTemp

Description
    Works well for combustion applications where optical thickness, tau is
    large, i.e. tau = a*L > 3 (L = distance between objects)

    Assumes
     - all surfaces are diffuse
     - tends to over predict radiative fluxes from sources/sinks
       *** SOURCES NOT CURRENTLY INCLUDED ***

SourceFiles
    heterogeneousMeanTemp.C

\*---------------------------------------------------------------------------*/

#ifndef heterogeneousRadiationModelheterogeneousMeanTemp_H
#define heterogeneousRadiationModelheterogeneousMeanTemp_H

#include "heterogeneousRadiationModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{

/*---------------------------------------------------------------------------*\
                           Class heterogeneousMeanTemp Declaration
\*---------------------------------------------------------------------------*/

class heterogeneousMeanTemp
:
    public heterogeneousRadiationModel
{
        // Private data

            //- Incident radiation / [W/m2]
            volScalarField G_;

            //- Total radiative heat flux [W/m2]
            volScalarField Qr_;

            //- Absorption coefficient
            volScalarField borderAs_;

            //- Emission contribution
            dimensionedScalar borderL_;

            //- porosity field
            const volScalarField& porosityF_;

            //- porosity field
            const volScalarField& surfFI_;

            //- porosity field
            volScalarField surfF_;

            //- magnitude of G_ field gradient
            volScalarField solidSh_;


public:

    //- Runtime type information
    TypeName("heterogeneousMeanTemp");


    // Constructors

        //- Construct from components
        heterogeneousMeanTemp
        (
            const volScalarField& T,
            const volScalarField& porosityF,
            const volScalarField& surfF,
            const volScalarField& Ts
        );

        //- Disallow default bitwise copy construction
        heterogeneousMeanTemp(const heterogeneousMeanTemp&) = delete;


    //- Destructor
    virtual ~heterogeneousMeanTemp();


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
        void operator=(const heterogeneousMeanTemp&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
