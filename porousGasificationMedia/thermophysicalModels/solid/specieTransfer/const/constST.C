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

#include "constST.H"
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

defineTypeNameAndDebug(constST, 0);
addToRunTimeSelectionTable(specieTransferModel, constST, porosity);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constST::constST
(
    const volScalarField& por,
    const volScalarField& por0
)
:
    specieTransferModel(por,por0),
    hCoeff_(0.0),
    SAV_(0.0)
{
    read(); 
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<constST> constST::New
(
    const volScalarField& por,
    const volScalarField& por0
)
{
    return autoPtr<constST>
    (
        new constST( por,por0)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



tmp<volScalarField> constST::ST() const
{
// eqZx2uHGn006
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "STconst",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "ST", dimless/dimTime, hCoeff_*SAV_
            )
        )
    );
}

bool constST::read()
{

	IOdictionary dict
    (
        IOobject
        (
            "specieTransferProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const dictionary& params = dict.subDict("Parameters");
    params.lookup("SAV") >> SAV_;    
    params.lookup("h") >> hCoeff_;

    return true;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//} // End namespace specieTransfer
} // End namespace Foam

// ************************************************************************* //
