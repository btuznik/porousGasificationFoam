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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "heterogeneousAbsorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::radiationModels::heterogeneousAbsorptionEmissionModel>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::New
(
    const dictionary& dict,
    const fvMesh& mesh
)
{
    word heterogeneousAbsorptionEmissionModelType(dict.lookup("heterogeneousAbsorptionEmissionModel"));

    Info<< "Selecting heterogeneousAbsorptionEmissionModel "
        << heterogeneousAbsorptionEmissionModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(heterogeneousAbsorptionEmissionModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "heterogeneousAbsorptionEmissionModel::New(const dictionary&, const fvMesh&)"
        )   << "Unknown heterogeneousAbsorptionEmissionModelType type "
            << heterogeneousAbsorptionEmissionModelType
            << ", constructor not in hash table" << nl << nl
            << "    Valid heterogeneousAbsorptionEmissionModel types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << exit(FatalError);
    }

    return autoPtr<heterogeneousAbsorptionEmissionModel>(cstrIter()(dict, mesh));
}


// ************************************************************************* //
