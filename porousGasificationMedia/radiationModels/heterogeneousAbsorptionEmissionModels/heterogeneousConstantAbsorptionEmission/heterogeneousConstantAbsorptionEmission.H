/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
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
    Foam::radiationModels::heterogeneousConstantAbsorptionEmission

Description
    Constant radiation absorption and emission coefficients for continuous
    phase

SourceFiles
    heterogeneousConstantAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef heterogeneousConstantAbsorptionEmission_H
#define heterogeneousConstantAbsorptionEmission_H

#include "heterogeneousAbsorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace heterogeneousAbsorptionEmissionModels
{

/*---------------------------------------------------------------------------*\
           Class heterogeneousConstantAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class heterogeneousConstantAbsorptionEmission
:
    public heterogeneousAbsorptionEmissionModel
{
    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Absorption coefficient / [1/m]
        dimensionedScalar a_;

        //- Absorption coefficient / [1/m]
        dimensionedScalar as_;

        //- Absorption coefficient / [1/m]
        dimensionedScalar borderAs_;

        //- Emission contribution / [kg/(m s^3)]
        dimensionedScalar E_;

        //- Emission contribution / [1/m]
        dimensionedScalar borderL_;

public:

    //- Runtime type information
    TypeName("heterogeneousConstantAbsorptionEmission");


    // Constructors

        //- Construct from components
        heterogeneousConstantAbsorptionEmission
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Destructor
    virtual ~heterogeneousConstantAbsorptionEmission();


    // Member Operators

        // Access

            // Absorption coefficient

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> aCont(const label bandI = 0) const;

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> asCont(const label bandI = 0) const;

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> borderAsCont(const label bandI = 0) const;

            // Emission contribution

                //- Emission contribution for continuous phase
                tmp<volScalarField> ECont(const label bandI = 0) const;

                // Border thickness layer
                virtual dimensionedScalar borderL(const label bandI = 0) const;

    // Member Functions

        inline bool isGrey() const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heterogeneousAbsorptionEmissionModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
