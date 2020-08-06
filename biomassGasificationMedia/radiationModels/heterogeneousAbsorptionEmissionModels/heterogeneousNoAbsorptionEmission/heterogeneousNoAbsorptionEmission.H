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
    Foam::radiationModels::heterogeneousAbsorptionEmissionsModels::heterogeneousNoAbsorptionEmission

Description
    Dummy absorption-emission model for 'none'

SourceFiles
    heterogeneousNoAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef heterogeneousNoAbsorptionEmission_H
#define heterogeneousNoAbsorptionEmission_H

#include "heterogeneousAbsorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace heterogeneousAbsorptionEmissionModels
{         

/*---------------------------------------------------------------------------*\
                  Class heterogeneousNoAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class heterogeneousNoAbsorptionEmission
:
    public heterogeneousAbsorptionEmissionModel
{
public:

    //- Runtime type information
    TypeName("heterogeneousNoAbsorptionEmission");


    // Constructors

        //- Construct from components
        heterogeneousNoAbsorptionEmission(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~heterogeneousNoAbsorptionEmission();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heterogeneousAbsorptionEmissionsModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
