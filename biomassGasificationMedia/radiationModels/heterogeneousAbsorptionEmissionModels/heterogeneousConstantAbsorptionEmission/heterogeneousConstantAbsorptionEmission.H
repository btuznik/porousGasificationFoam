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
    Foam::radiationModels::heterogeneousAbsorptionEmissionModels::heterogeneousConstant

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
                          Class heterogeneousConstant Declaration
\*---------------------------------------------------------------------------*/

class heterogeneousConstant
:
    public heterogeneousAbsorptionEmissionModel
{
    // Private Data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Absorption coefficient / [1/m]
        dimensionedScalar a_;

        //- Absorption coefficient / [1/m]
        dimensionedScalar as_;

        //- Absorption coefficient / [1/m]
        dimensionedScalar borderAs_;

        //- Emission coefficient / [1/m]
        dimensionedScalar e_;

        //- Emission coefficient / [1/m]
        dimensionedScalar es_;

        //- Emission coefficient / [1/m]
        dimensionedScalar borderEs_;

        //- Emission contribution / [kg/(m s^3)]
        dimensionedScalar E_;


public:

    //- Runtime type information
    TypeName("heterogeneousConstant");


    // Constructors

        //- Construct from components
        heterogeneousConstant(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~heterogeneousConstant();


    // Member Functions

        // Access

            // Absorption coefficient

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> aCont(const label bandI = 0) const;

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> asCont(const label bandI = 0) const;

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> borderAsCont(const label bandI = 0) const;


            // Emission coefficient

                //- Emission coefficient for continuous phase
                tmp<volScalarField> eCont(const label bandI = 0) const;

                //- Emission coefficient for continuous phase
                tmp<volScalarField> esCont(const label bandI = 0) const;

                //- Emission coefficient for continuous phase
                tmp<volScalarField> borderEsCont(const label bandI = 0) const;


            // Emission contribution

                //- Emission contribution for continuous phase
                tmp<volScalarField> ECont(const label bandI = 0) const;


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
