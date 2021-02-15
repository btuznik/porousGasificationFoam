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

namespace Foam
{
//namespace heatTransfer
//{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(constCONV, 0);
addToRunTimeSelectionTable(heatTransferModel, constCONV, porosity);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constCONV::constCONV
(
    const volScalarField& por,
    const volScalarField& por0
)
:
  heatTransferModel(por,por0),
  hCoeff_(0.0),
  SAV_(0.0),
  cylinderRadius_(1.0),
  constHTC_(true),
  Up_(db().lookupObject<volVectorField>("U")),
  rhop_(db().lookupObject<volScalarField>("rho")),
  alphap_(db().lookupObject<volScalarField>("thermo:alpha")),
  mup_(db().lookupObject<volScalarField>("thermo:mu")),
  thermop_(db().lookupObject<fluidThermo>("thermophysicalProperties"))
{
   read(); 
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<constCONV> constCONV::New
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



tmp<volScalarField> constCONV::CONV() const
{
// eqZx2uHGn006
    if (constHTC_)
    {
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
                    "CONV", dimEnergy/dimTime/dimTemperature/dimVolume, hCoeff_*SAV_
                )
            )
        );
    }
    else
    {
        Foam::tmp<Foam::volScalarField> CONVloc_ = Foam::tmp<Foam::volScalarField>
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
                    "zero", dimEnergy/dimTime/dimTemperature/dimVolume, 0.0
                )
            )
        );
        volScalarField& Cp = thermop_.Cp().ref();
        forAll (CONVloc_(),cellI)
        {
            CONVloc_.ref()[cellI] = SAV_*(1.
                    + 0.55
                    *Foam::pow(2*cylinderRadius_*rhop_[cellI]*mag(Up_[cellI])/mup_[cellI],0.6)
                    *Foam::pow(mup_[cellI]/alphap_[cellI],0.33333333333))
                    *Cp[cellI]*alphap_[cellI]*rhop_[cellI]/cylinderRadius_;  //eqZx2uHGn019 eqZx2uHGn020
        }
        return CONVloc_;
    }
}



bool constCONV::read()
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
    dict.lookup("constHTC") >> constHTC_;
    params.lookup("SAV") >> SAV_;
    if (constHTC_)
    {
        params.lookup("h") >> hCoeff_;
    }
    else
    {
        params.lookup("h") >> hCoeff_;
        params.lookup("cylinderRadius") >> cylinderRadius_;
    }

    return true;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//} // End namespace heatTransfer
} // End namespace Foam

// ************************************************************************* //
