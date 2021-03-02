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

#include "cylinder.H"
#include "Time.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cylinderCONV, 0);
    addToRunTimeSelectionTable(heatTransferModel, cylinderCONV, porosity);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::cylinderCONV::Re(const label cellI) const
{
    return 2 * cylinderRadius_ * rhop_[cellI] * mag(Up_[cellI]) / mup_[cellI];
}

Foam::scalar Foam::cylinderCONV::Pr(const label cellI) const
{
    return mup_[cellI] / alphap_[cellI];
}

Foam::scalar Foam::cylinderCONV::Nu(const label cellI) const //eqZx2uHGn019
{
    return 2. + 1.1 * pow(Re(cellI), 0.6) * cbrt(Pr(cellI));
}

Foam::scalar Foam::cylinderCONV::kf(const label cellI) const
{
    return thermop_.Cp().ref()[cellI] * alphap_[cellI] * rhop_[cellI];
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cylinderCONV::cylinderCONV
(
    const volScalarField& porosity,
    const volScalarField& initialPorosity
)
:
    heatTransferModel(porosity,initialPorosity),
    hCoeff_(1.0),
    cylinderRadius_(1.0),
    Up_(db().lookupObject<volVectorField>("U")),
    rhop_(db().lookupObject<volScalarField>("rho")),
    alphap_(db().lookupObject<volScalarField>("thermo:alpha")),
    mup_(db().lookupObject<volScalarField>("thermo:mu")),
    thermop_(db().lookupObject<fluidThermo>("thermophysicalProperties"))
{
   read();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::cylinderCONV> Foam::cylinderCONV::New
(
    const volScalarField& porosity,
    const volScalarField& initialPorosity
)
{
    return autoPtr<cylinderCONV>
    (
        new cylinderCONV(porosity,initialPorosity)
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::cylinderCONV::CONV() const
{

    Foam::tmp<Foam::volScalarField> CONVloc_ = Foam::tmp<Foam::volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "CONVloc",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "zero", dimEnergy/dimTime/dimTemperature/dimVolume, 0.0
            )
        )
    );

    forAll (CONVloc_(), cellI)
    {
        // Surface area to volume ratio.
        scalar SAV = 2.0 * sqrt(1 - porosity()[cellI]) * sqrt(1 - initialPorosity()[cellI])
                     / cylinderRadius_; // eqZx2uHGn007

        scalar h_conv = Nu(cellI) * kf(cellI) / (2 * cylinderRadius_); //eqZx2uHGn020

        CONVloc_.ref()[cellI] = SAV * h_conv;
    }

    return CONVloc_;
}

bool Foam::cylinderCONV::read()
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
    params.lookup("cylinderRadius") >> cylinderRadius_;


    return true;
}

// ************************************************************************* //
