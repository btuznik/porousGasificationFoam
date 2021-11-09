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

#include "emptyST.H"
#include "Time.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
//namespace specieTransfer
//{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(emptyST, 0);
addToRunTimeSelectionTable(specieTransferModel, emptyST, porosity);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

emptyST::emptyST
(
    const volScalarField& por,
    const volScalarField& por0
)
:
    specieTransferModel(por,por0)
{
    read(); 
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<emptyST> emptyST::New
(
    const volScalarField& por,
    const volScalarField& por0
)
{
    return autoPtr<emptyST>
    (
        new emptyST( por,por0)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



tmp<volScalarField> emptyST::ST() const
{
// eqZx2uHGn006
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "STempty",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "ST", dimless/dimTime, 0
            )
        )
    );
}

bool emptyST::read()
{
    return true;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//} // End namespace specieTransfer
} // End namespace Foam

// ************************************************************************* //
