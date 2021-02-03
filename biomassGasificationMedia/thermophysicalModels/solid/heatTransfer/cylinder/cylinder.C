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



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
//namespace heatTransfer
//{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(cylinderCONV, 0);
addToRunTimeSelectionTable(heatTransferModel, cylinderCONV, porosity);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

cylinderCONV::cylinderCONV
(
    const volScalarField& por,
    const volScalarField& por0
)
:
    heatTransferModel(por,por0),
    CONVCoeff_(1.0),
    borderCONVCoeff_(1.0),
    pipeRadius_(1.0)
{
   read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> cylinderCONV::CONV()
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

     CONVloc_ = pow(1 - porosity(),0.5)*pow(1-initialPorosity(),0.5)*2.0/pipeRadius_*CONVCoeff_;

     return CONVloc_;

}


tmp<volScalarField> cylinderCONV::borderCONV()
{
    Foam::tmp<Foam::volScalarField> borderCONVloc_ = Foam::tmp<Foam::volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "CONVBorder",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "CONV", dimEnergy/dimTime/dimTemperature/dimVolume, 0.0
            )
        )
    );

    borderCONVloc_ = pow(1 - porosity(),0.5)*pow(1-initialPorosity(),0.5)*2.0/pipeRadius_*borderCONVCoeff_;

return borderCONVloc_;
}


void cylinderCONV::correct()
{
    heatTransferModel::correct();
}


bool cylinderCONV::read()
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

    params.lookup("CONV") >> CONVCoeff_;
    params.lookup("borderCONV") >> borderCONVCoeff_;
    params.lookup("poreRadius") >> pipeRadius_;

    return true;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//} // End namespace heatTransfer
} // End namespace Foam

// ************************************************************************* //
