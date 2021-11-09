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

#include "pipe.H"
#include "Time.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(pipeCONV, 0);
    addToRunTimeSelectionTable(heatTransferModel, pipeCONV, porosity);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::pipeCONV::kf(const label cellI, const scalar Cp) const
{
    return Cp * alphap_[cellI] * rhop_[cellI];
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pipeCONV::pipeCONV
(
    const volScalarField& porosity,
    const volScalarField& initialPorosity
)
:
  heatTransferModel(porosity,initialPorosity),
  pipeRadius_(1.0),
  Up_(db().lookupObject<volVectorField>("U")),
  rhop_(db().lookupObject<volScalarField>("rho")),
  alphap_(db().lookupObject<volScalarField>("thermo:alpha")),
  mup_(db().lookupObject<volScalarField>("thermo:mu")),
  thermop_(db().lookupObject<fluidThermo>("thermophysicalProperties"))
{
   read(); 
}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::pipeCONV> Foam::pipeCONV::New
(
    const volScalarField& porosity,
    const volScalarField& initialPorosity
)
{
    return autoPtr<pipeCONV>
    (
        new pipeCONV(porosity,initialPorosity)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



Foam::tmp<Foam::volScalarField> Foam::pipeCONV::CONV() const
{
// eqZx2uHGn007
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
    static const scalar Nu = 3.66;
    const volScalarField& Cp = thermop_.Cp();

    forAll (CONVloc_(),cellI)
    {
        // Surface area to volume ratio.
        scalar SAV =  2.0 * sqrt(porosity()[cellI]) * sqrt(initialPorosity()[cellI])
                     / pipeRadius_; // eqZx2uHGn007

        scalar h_conv = Nu * kf(cellI, Cp[cellI]) / (2 * pipeRadius_); //eqZx2uHGn020

        CONVloc_.ref()[cellI] = SAV * h_conv;
    }

    return CONVloc_;
}



bool Foam::pipeCONV::read()
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

    params.lookup("pipeRadius") >> pipeRadius_;

    return true;
}

// ************************************************************************* //
