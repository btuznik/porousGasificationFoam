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

#include "heatTransferModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(heatTransferModel, 0);
    defineRunTimeSelectionTable(heatTransferModel, porosity);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModel::heatTransferModel
(
    const volScalarField& porosity,
    const volScalarField& initialPorosity
)
:
    IOdictionary
    (
        IOobject
        (
            "heatTransferProperties",
            porosity.mesh().time().constant(),
            porosity.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    runTime_(porosity.time()),
    mesh_(porosity.mesh()),
    porosity_(porosity),
    initialPorosity_(initialPorosity)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::heatTransferModel> Foam::heatTransferModel::New
(
    const volScalarField& porosity,
    const volScalarField& initialPorosity
)
{
    // Get model name, but do not register the dictionary
    // otherwise it is registered in the database twice
    const word modelType
    (
        IOdictionary
        (
            IOobject
            (
                "heatTransferProperties",
                porosity.time().constant(),
                porosity.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("heatTransferModel")
    );

    Info<< "Selecting heatTransfer model type " << modelType << endl;

    porosityConstructorTable::iterator cstrIter =
        porosityConstructorTablePtr_->find(modelType);

    if (cstrIter == porosityConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "heatTransferModel::New(const volVectorField&, "
            "const surfaceScalarField&)"
        )   << "Unknown heatTransferModel type "
            << modelType << nl << nl
            << "Valid heatTransferModel types:" << endl
            << porosityConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

 
    return autoPtr<heatTransferModel>
    (
        cstrIter()(porosity,initialPorosity)
    );

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
