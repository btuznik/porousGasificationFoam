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

#include "specieTransferModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(specieTransferModel, 0);
defineRunTimeSelectionTable(specieTransferModel, porosity);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

specieTransferModel::specieTransferModel
(
    const volScalarField& por,
    const volScalarField& por0
)
:
    IOdictionary
    (
        IOobject
        (
            "specieTransferProperties",
            por.mesh().time().constant(),
            por.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    runTime_(por.time()),
    mesh_(por.mesh()),
    por_(por),
    por0_(por0)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //


autoPtr<specieTransferModel> specieTransferModel::New
(
    const volScalarField& por,
    const volScalarField& por0
)
{
    // get model name, but do not register the dictionary
    // otherwise it is registered in the database twice

    porosityConstructorTable::iterator cstrIter;

    if (por.db().lookupObject<dictionary>("chemistryProperties").lookupOrDefault("diffusionLimitedReactions",false))
    {
        const word modelType
        (
            IOdictionary
            (
                IOobject
                (
                    "specieTransferProperties",
                    por.time().constant(),
                    por.db(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            ).lookup("specieTransferModel")
        );

        Info<< "Selecting specieTransfer model type " << modelType << endl;

        cstrIter = porosityConstructorTablePtr_->find(modelType);

        if (cstrIter == porosityConstructorTablePtr_->end())
        {
            FatalErrorIn
            (
                "specieTransferModel::New(const volVectorField&, "
                "const surfaceScalarField&)"
            )   << "Unknown specieTransferModel type "
                << modelType << nl << nl
                << "Valid specieTransferModel types:" << endl
                << porosityConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        } 
    }
    else
    {
        cstrIter = porosityConstructorTablePtr_->find("emptyST");
    }

    return autoPtr<specieTransferModel>
    (
        cstrIter()(por,por0)
    );

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
