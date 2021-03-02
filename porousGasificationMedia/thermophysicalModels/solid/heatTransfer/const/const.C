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

#include "const.H"
#include "Time.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constCONV, 0);
    addToRunTimeSelectionTable(heatTransferModel, constCONV, porosity);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constCONV::constCONV
(
    const volScalarField& por,
    const volScalarField& por0
)
:
  heatTransferModel(por,por0),
  hCoeff_(0.0),
  SAV_(0.0)
{
   read(); 
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::constCONV> Foam::constCONV::New
(
    const volScalarField& por,
    const volScalarField& por0
)
{
    return autoPtr<constCONV>
    (
        new constCONV( por,por0)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::constCONV::CONV() const
{
// eqZx2uHGn006
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "CONVconst",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "CONV", dimEnergy/dimTime/dimTemperature/dimVolume, hCoeff_ * SAV_
            )
        )
    );
}



bool Foam::constCONV::read()
{

    IOdictionary dict
    (
        IOobject
        (
            "heatTransferProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const dictionary& params = dict.subDict("Parameters");
    params.lookup("h") >> hCoeff_;
    params.lookup("SAV") >> SAV_;

    return true;
}

// ************************************************************************* //
