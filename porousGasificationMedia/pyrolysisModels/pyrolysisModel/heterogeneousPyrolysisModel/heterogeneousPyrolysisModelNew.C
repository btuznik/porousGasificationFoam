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

\*---------------------------------------------------------------------------*/

#include "heterogeneousPyrolysisModel.H"
#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heterogeneousPyrolysisModels
{

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<heterogeneousPyrolysisModel> heterogeneousPyrolysisModel::New
(
    const fvMesh& mesh,
    HGSSolidThermo& solidThermo,
    psiReactionThermo& gasThermo,
    volScalarField& whereIs
)
{
    // Get model name, but do not register the dictionary.
    const word modelType
    (
        IOdictionary
        (
            IOobject
            (
                "pyrolysisProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("heterogeneousPyrolysisModel")
    );

    Info<< "Selecting heterogeneousPyrolysisModel " << modelType << endl;

    noRadiationConstructorTable::iterator cstrIter =
        noRadiationConstructorTablePtr_->find(modelType);

    if (cstrIter == noRadiationConstructorTablePtr_->end())
    {
        FatalErrorIn("heterogeneousPyrolysisModel::New(const fvMesh&)")
            << "Unknown heterogeneousPyrolysisModel type " << modelType
            << nl << nl << "Valid heterogeneousPyrolisisModel types are:" << nl
            << noRadiationConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<heterogeneousPyrolysisModel>(cstrIter()(modelType, mesh, solidThermo, gasThermo, whereIs));
}

autoPtr<heterogeneousPyrolysisModel> heterogeneousPyrolysisModel::New
(
    const fvMesh& mesh,
    HGSSolidThermo& solidThermo,
    psiReactionThermo& gasThermo,
    volScalarField& whereIs,
    volScalarField& radiation
)
{
    // Get model name, but do not register the dictionary.
    const word modelType
    (
        IOdictionary
        (
            IOobject
            (
                "pyrolysisProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("heterogeneousPyrolysisModel")
    );

    Info<< "Selecting heterogeneousPyrolysisModel " << modelType << endl;

    radiationConstructorTable::iterator cstrIter =
        radiationConstructorTablePtr_->find(modelType);

    if (cstrIter == radiationConstructorTablePtr_->end())
    {
        FatalErrorIn("heterogeneousPyrolysisModel::New(const fvMesh&)")
            << "Unknown heterogeneousPyrolysisModel type " << modelType
            << nl << nl << "Valid heterogeneousPyrolisisModel types are:" << nl
            << radiationConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<heterogeneousPyrolysisModel>(cstrIter()(modelType, mesh, solidThermo, gasThermo, whereIs, radiation));
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heterogeneousPyrolysisModels
} // End namespace Foam

// ************************************************************************* //
