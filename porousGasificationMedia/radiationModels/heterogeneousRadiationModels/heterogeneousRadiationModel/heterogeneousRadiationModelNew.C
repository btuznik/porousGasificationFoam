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

#include "heterogeneousRadiationModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::heterogeneousRadiationModel>
Foam::heterogeneousRadiationModel::New
(
    const volScalarField& T
)
{
    // get model name, but do not register the dictionary
    const word modelType
    (
        IOdictionary
        (
            IOobject
            (
                "radiationProperties",
                T.time().constant(),
                T.mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("heterogeneousRadiationModel")
    );

    Info<< "Selecting heterogeneousRadiationModel " << modelType << endl;

    TConstructorTable::iterator cstrIter =
        TConstructorTablePtr_->find(modelType);

    if (cstrIter == TConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "heterogeneousRadiationModel::New(const volScalarField&)"
        )   << "Unknown heterogeneousRadiationModel type "
            << modelType << nl << nl
            << "Valid heterogeneousRadiationModel types are:" << nl
            << TConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<heterogeneousRadiationModel>(cstrIter()(T));
}


Foam::autoPtr<Foam::heterogeneousRadiationModel>
Foam::heterogeneousRadiationModel::New
(
    const volScalarField& T,
    const volScalarField& porosityF,
    const volScalarField& surfF,
    const volScalarField& Ts
)
{
    // get model name, but do not register the dictionary
    const word modelType
    (
        IOdictionary
        (
            IOobject
            (
                "radiationProperties",
                T.time().constant(),
                T.mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("heterogeneousRadiationModel")
    );

    Info<< "Selecting heterogeneousRadiationModel " << modelType << endl;

    porosityConstructorTable::iterator cstrIter =
        porosityConstructorTablePtr_->find(modelType);

    if (cstrIter == porosityConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "heterogeneousRadiationModel::New(const volScalarField&, const volScalarField&)"
        )   << "Unknown heterogeneousRadiationModel type "
            << modelType << nl << nl
            << "Valid heterogeneousRadiationModel types are:" << nl
            << porosityConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<heterogeneousRadiationModel>(cstrIter()(T,porosityF,surfF,Ts));
}

// ************************************************************************* //
