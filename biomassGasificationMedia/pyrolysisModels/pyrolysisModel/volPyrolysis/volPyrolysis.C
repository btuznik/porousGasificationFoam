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

#include "volFields.H"
#include "volPyrolysis.H"
#include "dimensionSet.H"
#include "addToRunTimeSelectionTable.H"
#include "mapDistribute.H"
#include "zeroGradientFvPatchFields.H"
#include "surfaceInterpolate.H"
#include "fvm.H"
#include "fvcDiv.H"
#include "fvcVolumeIntegrate.H"
#include "fvMatrices.H"
#include "fvCFD.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heterogeneousPyrolysisModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(volPyrolysis, 0);

addToRunTimeSelectionTable(heterogeneousPyrolysisModel, volPyrolysis, noRadiation);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void volPyrolysis::readReactingOneDimControls()
{
    const dictionary& solution = this->solution().subDict("PIMPLE");
    solution.lookup("nNonOrthogonalCorrectors") >> nNonOrthCorr_;
    time_.controlDict().lookup("maxDi") >> maxDiff_;
    time_.controlDict().lookup("maxdT") >> maxDT_;
    equilibrium_.readIfPresent("equilibrium",coeffs_);
}

bool volPyrolysis::read()
{
    if (heterogeneousPyrolysisModel::read())
    {
        readReactingOneDimControls();
        return true;
    }
    else
    {
        return false;
    }
}

bool volPyrolysis::read(const dictionary& dict)
{
    if (heterogeneousPyrolysisModel::read(dict))
    {
        readReactingOneDimControls();
        return true;
    }
    else
    {
        return false;
    }
}

Foam::tmp<Foam::volScalarField> volPyrolysis::HTC() const
{
    Foam::tmp<Foam::volScalarField> HTCloc_ = Foam::tmp<Foam::volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "HTC",
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimTemperature/dimTime/dimVolume, 0.0)
        )
    );

    volScalarField borderHTC=HTmodel_->borderHTC();

    if(equilibrium_)
    {}
    else
    {
        HTCloc_ = HTmodel_->HTC()*whereIs_;
       // forAll(surfF_,cellI)
       // {
       //     HTCloc_.internalField()[surfF_[cellI]] = borderHTC[surfF_[cellI]];
       // }
    }
    return HTCloc_;
}

Foam::tmp<Foam::volScalarField> volPyrolysis::heatTransfer()
{
    Foam::tmp<Foam::volScalarField> Sh_ = Foam::tmp<Foam::volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "pyrolysisSh",
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );

    if(equilibrium_)
    {}
    else
    {
        volScalarField Tgas = gasThermo_.T();
        volScalarField HT(HTC_ * (Ts_-Tgas));
        Sh_ = HT * whereIs_;
    }
    return Sh_;
}

// TODO
/*
Foam::tmp<Foam::volScalarField> volPyrolysis::heatUpGasCalc() const
{

    Foam::tmp<Foam::volScalarField> hSh_ = Foam::tmp<Foam::volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "hSh_",
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );

    if (active_)
    {
        volScalarField Tgas = gasThermo_.T();

        if (equilibrium_)
        {}
        else
        {
            volScalarField tempSh = hSh_();
            volScalarField tempDiff = ((whereIs_*Ts_)+(whereIsNot_*Tgas));
            const speciesTable& gasTable = gasThermo_.composition().gasTable();
            forAll(gasTable,gasI)
            {
                tempSh += (solidChemistry_().gasHs(tempDiff,gasI) - solidChemistry_().gasHs(Tgas,gasI))*Srho(gasI);
            }
            hSh_= tempSh * whereIs_;
        }
    }
    return hSh_;
}
*/
// TODO
void volPyrolysis::evolvePorosity()
{
    if (active_)
    {}
    else
    {}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

volPyrolysis::volPyrolysis
(
    const word& modelType,
    const fvMesh& mesh,
    psiReactionThermo& gasThermo,
    volScalarField& whereIs
)
:
    heterogeneousPyrolysisModel(modelType, mesh),
    gasThermo_(gasThermo),
    equilibrium_(false),
    rho_
    (
        IOobject
        (
         "rhos",
         mesh_.time().timeName(),
         mesh_,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimMass/dimVolume, 0.0)
     ),
    Ts_
    (
         IOobject
         (
             "Ts",
             time_.timeName(),
             mesh_,
             IOobject::MUST_READ,
             IOobject::AUTO_WRITE
         ),
         mesh
    ),
    nNonOrthCorr_(-1),
    maxDiff_(10),
    porosity_(whereIs),
    initialPorosity_(whereIs),
    porosityArch_
    (
        IOobject
        (
            "porosityF0",
            time_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    voidFraction_(whereIs),

    whereIs_
    (
        IOobject
        (
            "whereIs",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        scalar(1.0)
    ),
    whereIsNot_
    (
        IOobject
        (
            "whereIsNot",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        scalar(0.0)
    ),
    HTmodel_(heatTransferModel::New(porosity_,porosityArch_)),
    HTC_
    (
        IOobject
        (
            "HTCvol",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",dimEnergy/dimTime/dimTemperature/dimVolume , 0.0)

    ),
    surfF_(),
    maxDT_(1.0)
{

    HTC_ = HTC();

    forAll(rho_, cellI)
    {
    	if (porosity_[cellI] == 1)
        {
            whereIsNot_[cellI] = 1;
            whereIs_[cellI] = 0;
        }
    }
    if (active_)
    {
        read();
    }

    DynamicList<label> dSurfF;
    forAll(whereIs_,cellI)
    {
    if (whereIs_[cellI] == 1)
    {
            bool surfC = false;
            forAll(mesh_.cellCells()[cellI],cellJ)
            {
                    if (whereIs_[mesh_.cellCells()[cellI][cellJ]] == 0) surfC = true;
            }
            if (surfC) dSurfF.append(cellI);
    }
    }
    surfF_ = dSurfF.shrink();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

volPyrolysis::~volPyrolysis()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Switch volPyrolysis::equilibrium() const
{
    return equilibrium_;
}

void volPyrolysis::preEvolveRegion()
{}


void volPyrolysis::evolveRegion()
{
    voidFraction_ = porosity_;
    if (equilibrium_)
    {}
    else
    {
        HTC_ = HTC();
    }

    evolvePorosity();
    info();
}

/// Return const temperature of solid phase [K]
const volScalarField& volPyrolysis::T() const
{
    return Ts_;
}

void volPyrolysis::info() const
{
    Info<< "\nPyrolysis: " << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace heterogeneousPyrolysisModels

// ************************************************************************* //
