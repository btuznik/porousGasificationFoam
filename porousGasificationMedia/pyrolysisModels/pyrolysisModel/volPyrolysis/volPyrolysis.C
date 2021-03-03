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
addToRunTimeSelectionTable(heterogeneousPyrolysisModel, volPyrolysis, radiation);

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

void volPyrolysis::solveSpeciesMass()
{
// eqZx2uHGn045
// eqZx2uHGn046

    if (debug)
    {
        Info<< "volPyrolysis:+solveSpeciesMass()" << endl;
    }

    if (active_)
    {
        volScalarField Mt(0.0 * Ym_[0]);
        volScalarField rhoLoc
        (
            max(rho_ * (1. - porosity_), dimensionedScalar("minRho",dimMass/dimVolume, SMALL))
        );

        for (label i = 0; i < Ys_.size(); ++i)
        {

            volScalarField& Yi = Ys_[i];
            Yi.ref() *= whereIs_;
            volScalarField sRhoSi = solidChemistry_->RRs(i);

            fvScalarMatrix YsEqn
            (
                fvm::ddt(rhoLoc,Yi)
             ==
                sRhoSi

            );

            YsEqn.relax();
            YsEqn.solve("Ys");

            Yi.max(0.0);
            Mt += Yi * rho_;

            Info<< "solid "<< Ys_[i].name()
                << " equation solved. Sources min/max   = " << gMin(sRhoSi)
                << ", " << gMax(sRhoSi)
                << "; values min Y = " << gMin(Ys_[i])
                << " max Y = " << gMax(Ys_[i]) << endl;
        }

        for (label i = 0; i < Ys_.size(); ++i)
        {
            Ym_[i] = whereIs_ * Ys_[i] * rho_;
            Ys_[i] = Ym_[i] / max(Mt,dimensionedScalar("minMass", dimMass/dimVolume, SMALL));
            Ym_[i] *= (1. - porosity_);
        }
    }
}


void volPyrolysis::solveEnergy()
{
// eqZx2uHGn047

    if (debug)
    {
        Info<< "volPyrolysis::solveEnergy()" << endl;
    }

    if (active_)
    {

        volTensorField composedK(K_ * (1 - porosity_) * anisotropyK_);
        radiationSh_ = radiation_;

        if (equilibrium_)
        {}
        else
        {
            volScalarField rhoCp
            (
                max
                (
                    rho_ * solidThermo_.Cp() * (1 - porosity_),
                    dimensionedScalar("minRhoCp",dimEnergy/dimTemperature/dimVolume,SMALL)
                )
            );

            // Simplistic immersed boundary for heat transport in solid phase.
            fvScalarMatrix TLap
            (
                fvm::laplacian(composedK, T_)
            );

            // Setting face fluxes otutside porous media to 0.
            surfaceScalarField  whereIsPatch  = fvc::snGrad(whereIs_);

            forAll(whereIsPatch,faceI)
            {
                if (whereIsPatch[faceI] != 0.)
                {
                    TLap.upper()[faceI] = 0.;
                }
            }

            TLap.diag() = 0;
            TLap.negSumDiag();

            // Correct on orthogonal meshes
            // For non-orthogona meshes
            // TLap.source() = 0. cancels the non-orthogonal correction from the
            // divergence of composedK, which is spurious on the edges of porous
            // media and negiligble inside porous media in most cases.
            TLap.source() = 0.;

            fvScalarMatrix TEqn
            (
                fvm::ddt(rhoCp, T_)
              - TLap
            ==
                chemistrySh_ // eqZx2uHGn004, eqZx2uHGn017
              - heatTransfer()() // eqZx2uHGn005
              - heatUpGas_
              + radiationSh_
            );

            TEqn.relax();
            TEqn.solve();
        }

        scalar minTemp = GREAT;
        scalar maxTemp = -GREAT;
        scalar areThere = 0;

        forAll(T_,cellI)
        {
            if (whereIs_[cellI] == 1 )
            {
                areThere = 1;
                if (T_[cellI] < minTemp)
                {
                    minTemp = T_[cellI];
                }
                if (T_[cellI] > maxTemp)
                {
                    maxTemp = T_[cellI];
                }
            }
        }
        reduce(maxTemp, maxOp<scalar>());
        reduce(minTemp, minOp<scalar>());
        reduce(areThere, maxOp<scalar>());

        if (areThere == 1)
        {
            Info<< " pyrolysis min/max(T) = " << minTemp << ", " << maxTemp << endl;
        }
        else
        {
            Info<< " no solid phase " << endl;
        }
    }
}

void volPyrolysis::calculateMassTransfer()
{
    if (infoOutput_)
    {
        totalGasMassFlux_ = fvc::domainIntegrate(solidChemistry_->RRg());
        totalHeatRR_ = fvc::domainIntegrate(chemistrySh_);

        addedGasMass_ +=
            fvc::domainIntegrate(solidChemistry_->RRg()) * time_.deltaT();
        lostSolidMass_ +=
            fvc::domainIntegrate(solidChemistry_->RRs()) * time_.deltaT();
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

volPyrolysis::volPyrolysis
(
    const word& modelType,
    const fvMesh& mesh,
    HGSSolidThermo& solidThermo,
    psiReactionThermo& gasThermo,
    volScalarField& whereIs
)
:
    heterogeneousPyrolysisModel(modelType, mesh),
    gasThermo_(gasThermo),
    Ygas_(gasThermo.composition().Y()),
    solidThermo_(solidThermo),
    solidChemistry_
    (
        BasicSolidChemistryModel<HGSSolidThermo>::New
        (
            solidThermo_,
            Ygas_,
            gasThermo.thermoName()
        )
    ),
    kappa_(solidThermo_.kappa()),
    K_(solidThermo_.K()),
    rho_(solidThermo_.rho()),
    rho0_
    (
        IOobject
        (
            "rhos0",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass/dimVolume, 0.0)
    ),
    Ys_(solidThermo_.composition().Y()),
    Ym_(Ys_.size()),
    Msolid_(Ys_.size()),
    Msolidtotal_(
        IOobject
        (
            "Msolidtotal",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass, 0.0)
    ),
    T_(solidThermo_.T()),
    equilibrium_(false),
    subintegrateSwitch_(false),
    nNonOrthCorr_(-1),
    maxDiff_(10),
    porosity_(whereIs),
    porosity0_(whereIs),
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
    radiation_(whereIs),
    phiHsGas_
    (
        IOobject
        (
            "phiHsGass",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
    ),
    chemistrySh_
    (
        IOobject
        (
            "solidChemistrySh",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),
    porositySource_
    (
        IOobject
        (
            "porosityS",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless/dimTime, 0.0)
    ),
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
    heatUpGas_
    (
        IOobject
        (
            "heatUpGas",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),
    HTmodel_(heatTransferModel::New(porosity_,porosityArch_)),
    CONV_
    (
        IOobject
        (
            "CONVvol",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",dimEnergy/dimTime/dimTemperature/dimVolume , 0.0)

    ),
    radiationSh_
    (
        IOobject
        (
            "radiationSh",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
    ),
    anisotropyK_
    (
        IOobject
        (
            "anisotropyK",
            time_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor("one", dimless, tensor(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
    ),
    surfF_(whereIs_),
    lostSolidMass_(dimensionedScalar("zero", dimMass, 0.0)),
    addedGasMass_(dimensionedScalar("zero", dimMass, 0.0)),
    totalGasMassFlux_(dimensionedScalar("zero", dimMass/dimTime, 0.0)),
    totalHeatRR_(dimensionedScalar("zero", dimEnergy/dimTime, 0.0)),
    timeChem_(1.0),
    maxDT_(1.0)
{
    mesh.setFluxRequired(T_.name());

    CONV_ = CONV();
    rho0_.ref() = rho_.ref();

    subintegrateSwitch_ = coeffs().lookupOrDefault("subintegrateHeatTransfer",false);

    forAll(Ys_, fieldI)
    {
        Ym_.set
        (
            fieldI,
            new volScalarField
            (
                IOobject
                (
                    Ys_[fieldI].name() + "m",
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimMass/dimVolume, 0.0)
            )
        );
        Ym_[fieldI].ref() = Ys_[fieldI].ref() * rho_.ref();

        // Msolid is used for renormalization of mass fractions of solid species.
        Msolid_.set
        (
            fieldI,
            new volScalarField
            (
                IOobject
                (
                    Ys_[fieldI].name() + "m2",
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimMass, 0.0)
            )
        );
        forAll(rho_,cellI)
        {
            Msolid_[fieldI][cellI] =
                Ys_[fieldI][cellI]
               * rho_[cellI]
               * mesh_.V()[cellI]
               * (1 - porosity_[cellI]);
        }
    }

    forAll(rho_,cellI)
    {
        if (porosity_[cellI] == 1)
        {
             whereIsNot_[cellI] = 1;
             whereIs_[cellI] = 0;
             rho0_[cellI] = 3.14;
        }
    }
    if (active_)
    {
        read();
    }

    surfF_ = 0;
    forAll(whereIs_,cellI)
    {
        if (whereIs_[cellI] == 1)
        {
            bool surfC = false;
            forAll(mesh_.cellCells()[cellI],cellJ)
            {
                    if (whereIs_[mesh_.cellCells()[cellI][cellJ]] == 0)
                    {
                        surfC = true;
                    }
            }
            if (surfC)
            {
                surfF_[cellI] = 1;
            }
        }
    }

    forAll(rho_,cellI)
    {
        Msolidtotal_[cellI] = 0.0;

        for (label i = 0; i < Ys_.size(); ++i)
        {
             Msolid_[i][cellI]=
                Ys_[i][cellI]
               *rho_[cellI]
               *mesh_.V()[cellI]
               *(1 - porosity_[cellI]);

             Msolidtotal_[cellI] += Msolid_[i][cellI];
        }

        for (label i = 0; i < Ys_.size(); ++i)
        {
            if(Msolidtotal_[cellI] > 0.0)
            {
               Ys_[i][cellI] = (Msolid_[i][cellI] / (Msolidtotal_[cellI]));
            }
        }
    }
}

volPyrolysis::volPyrolysis
(
    const word& modelType,
    const fvMesh& mesh,
    HGSSolidThermo& solidThermo,
    psiReactionThermo& gasThermo,
    volScalarField& whereIs,
    volScalarField& radiation
)
:
    heterogeneousPyrolysisModel(modelType, mesh),
    gasThermo_(gasThermo),
    Ygas_(gasThermo.composition().Y()),
    solidThermo_(solidThermo),
    solidChemistry_
    (
        BasicSolidChemistryModel<HGSSolidThermo>::New
        (
            solidThermo_,
            Ygas_,
            gasThermo.thermoName()
        )
    ),
    kappa_(solidThermo_.kappa()),
    K_(solidThermo_.K()),
    rho_(solidThermo_.rho()),
    rho0_
    (
        IOobject
        (
            "rhos0",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass/dimVolume, 0.0)
    ),
    Ys_(solidThermo_.composition().Y()),
    Ym_(Ys_.size()),
    Msolid_(Ys_.size()),
    Msolidtotal_(
        IOobject
        (
            "Msolidtotal",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass, 0.0)
    ),
    T_(solidThermo_.T()),
    equilibrium_(false),
    subintegrateSwitch_(false),
    nNonOrthCorr_(-1),
    maxDiff_(10),
    porosity_(whereIs),
    porosity0_(whereIs),
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
    radiation_(radiation),
    phiHsGas_
    (
        IOobject
        (
            "phiHsGass",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
    ),
    chemistrySh_
    (
        IOobject
        (
            "solidChemistrySh",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),
    porositySource_
    (
        IOobject
        (
            "porosityS",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless/dimTime, 0.0)
    ),
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
        scalar(0.0)
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
        scalar(1.0)
    ),
    heatUpGas_
    (
        IOobject
        (
            "heatUpGas",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),
    HTmodel_(heatTransferModel::New(porosity_,porosityArch_)),
    CONV_
    (
        IOobject
        (
            "CONVvol",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",dimEnergy/dimTime/dimTemperature/dimVolume , 0.0)
    ),
    radiationSh_
    (
        IOobject
        (
            "radiationSh",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
    ),
    anisotropyK_
    (
        IOobject
        (
            "anisotropyK",
            time_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor("one", dimless, tensor(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
    ),
    surfF_
    (
        IOobject
        (
            "surfF",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        scalar(0.0)
    ),
    lostSolidMass_(dimensionedScalar("zero", dimMass, 0.0)),
    addedGasMass_(dimensionedScalar("zero", dimMass, 0.0)),
    totalGasMassFlux_(dimensionedScalar("zero", dimMass/dimTime, 0.0)),
    totalHeatRR_(dimensionedScalar("zero", dimEnergy/dimTime, 0.0)),
    timeChem_(1.0),
    maxDT_(1.0)
{
    mesh.setFluxRequired(T_.name());

    CONV_ = CONV();
    rho0_.ref() = rho_.ref();
    subintegrateSwitch_ = coeffs().lookupOrDefault("subintegrateHeatTransfer",false);

    forAll(Ys_, fieldI)
    {
          Ym_.set
          (
                fieldI,
                new volScalarField
                (
                    IOobject
                    (
                        Ys_[fieldI].name() + "m",
                        time_.timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("zero",dimMass/dimVolume,0.0)
                )
          );

          Ym_[fieldI].ref() = Ys_[fieldI].ref() * rho_;
    }

    forAll(rho_,cellI)
    {
        if (porosity_[cellI] < 1.)
        {
             whereIsNot_[cellI] = 0.;
             whereIs_[cellI] = 1.;
        }
        else
        {
             whereIsNot_[cellI] = 1.;
             whereIs_[cellI] = 0.;
         rho0_[cellI] = 3.14;
        }
    }

    if (active_)
    {
        read();
    }

    surfF_ = 0;
    surfaceScalarField  whereIsPatch  = fvc::interpolate(whereIs_);

    forAll(whereIs_, cellI)
    {
        if (whereIs_[cellI] == 1)
        {
            bool surfC = false;
            forAll(mesh_.cells()[cellI],faceI)
            {
                label faceIl = mesh_.cells()[cellI][faceI];
                if (mesh_.isInternalFace(faceIl))
                {
                    if (whereIs_[mesh_.faceOwner()[faceIl]] != whereIs_[mesh_.faceNeighbour()[faceIl]])
                    {
                        surfC = true;
                    }
                }
                else
                {}
            }

            if (surfC) surfF_[cellI] = 1;
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

volPyrolysis::~volPyrolysis()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
volScalarField& volPyrolysis::rhoConst() const
{
    return rho_;
}

volScalarField& volPyrolysis::rho()
{
    return rho_;
}

const volScalarField& volPyrolysis::T() const
{
    return T_;
}

const tmp<volScalarField> volPyrolysis::Cp() const
{
    return solidThermo_.Cp();
}


const volScalarField& volPyrolysis::kappa() const
{
    return kappa_;
}

const volScalarField& volPyrolysis::K() const
{
    return K_;
}

const volScalarField& volPyrolysis::surf() const
{
    return surfF_;
}

scalar volPyrolysis::maxDiff() const
{
    return maxDiff_;
}

scalar volPyrolysis::solidRegionDiffNo() const
{
    scalar NumCprho = 0.0;
    scalar DiNum = 0.0;

    if (mesh_.nInternalFaces() > 0)
    {
        surfaceScalarField KrhoCpbyDelta
        (
            mesh_.surfaceInterpolation::deltaCoeffs()
          * fvc::interpolate(K_)
        );

        surfaceScalarField Cprho
        (
            fvc::interpolate(Cp()*rho_)
        );

        NumCprho = max(Cprho.ref()).value();
        reduce(NumCprho, maxOp<scalar>());
        if (NumCprho != 0.)
        {
            DiNum = gMax(KrhoCpbyDelta.internalField())*time_.deltaTValue()/NumCprho;
        }
        else
        {
            DiNum = SMALL;
        }
    }

    return DiNum;
}

scalar volPyrolysis::maxTime() const
{
return timeChem_;
}

Switch volPyrolysis::equilibrium() const
{
    return equilibrium_;
}

void volPyrolysis::preEvolveRegion() {
    // Added for completeness
    heterogeneousPyrolysisModel::preEvolveRegion();

    // Iterates over every cell and sets cells containing solid phase
    // as reacting cell.
    forAll(T_, cellI)
    {
        if ( active_ && whereIs_[cellI] != 0)
        {
            solidChemistry_->setCellReacting(cellI, true);
        }
        else
        {
            solidChemistry_->setCellReacting(cellI, false);
        }
    }

}

void volPyrolysis::evolveRegion()
{
    voidFraction_ = porosity_;

    if (equilibrium_)
    {}
    else
    {
        CONV_ = CONV();
    }

    timeChem_ = solidChemistry_->solve
    (
        time_.value() - time_.deltaTValue(),
        time_.deltaTValue()
    );

    heatUpGas_ = heatUpGasCalc()();

    chemistrySh_ = solidChemistry_->Sh()(); // eqZx2uHGn004

    solveSpeciesMass();

    evolvePorosity();

    solidThermo_.correct(); // eqZx2uHGn046

    if (equilibrium_)
    {}
    else
    {
        for (int nonOrth = 0; nonOrth <= nNonOrthCorr_; ++nonOrth)
        {
            solveEnergy();
        }
    }

    calculateMassTransfer();

    info();
}

void volPyrolysis::evolvePorosity()
{
    if (active_)
    {
        porositySource_ = solidChemistry_->RRpor(T_)();

        volScalarField& por = porosity_;

        fvScalarMatrix porosityEqn
        (
            fvm::ddt(por)
         ==
            porositySource_
        );

        porosityEqn.solve("porosity");

        Info<< "prosity equation solved. Sources min/max   = " << gMin(porositySource_)
            << ", " << gMax(porositySource_);

        Info<< "; values min Y = " << gMin(por)
            <<" max Y = " << gMax(por) << endl;

        forAll(porosity_,cellI)
        {
            if (porosity_[cellI] > 0.9999)
            {
                porosity_[cellI] = 1.0;
                T_[cellI] = gMin(T_);
            }
            if (porosity_[cellI] < 0.0001)
            {
                porosity_[cellI] = 0.0;
            }
            if (porosity_[cellI] < 1.0)
            {
                whereIs_[cellI] = 1.0;
                whereIsNot_[cellI] = 0.0;
            }
            else
            {
                whereIs_[cellI] = 0.0;
                whereIsNot_[cellI] = 1.0;
            }
        }
    }
    else
    {}
}

Foam::tmp<Foam::volScalarField> volPyrolysis::Srho() const
{
    Foam::tmp<Foam::volScalarField> tSrho
    (
        new volScalarField
        (
            IOobject
            (
                "thermoSingleLayer::Srho",
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    if (active_)
    {
        const speciesTable& gasTable = solidChemistry_->gasTable();

        forAll(gasTable,gasI)
        {
            tSrho = tSrho + Srho(gasI);
        }
    }

    return tSrho;
}


Foam::tmp<Foam::volScalarField> volPyrolysis::Srho(const label i) const
{

    if (active_)
    {
        const speciesTable& gasTable = solidChemistry_->gasTable();
        label j=-1;

        forAll(gasTable,gasI)
        {
            if (gasTable[gasI] == Ygas_[i].name())
            {
                j = gasI;
            }
        }

        if (j > -1)
        {

            tmp<volScalarField> tRRiGas = solidChemistry_->RRg(j);
            return tRRiGas;
        }
        else
        {
            return Foam::tmp<Foam::volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "volPyrolysis::Srho(" + Foam::name(i) + ")",
                        time_.timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE,
                        false
                    ),
                    mesh_,
                    dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
                )
            );

        }
    }
    else
    {
        return Foam::tmp<Foam::volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "volPyrolysis::Srho(" + Foam::name(i) + ")",
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                mesh_,
                dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
            )
        );
    }
}

Foam::tmp<Foam::volScalarField> volPyrolysis::heatTransfer()
{
// eqZx2uHGn005
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
        if (subintegrateSwitch_)
        {
            volScalarField rhoCpG(gasThermo_.rho() * gasThermo_.Cp() * porosity_);
            volScalarField Tgas = gasThermo_.T();
            volScalarField rhoCpS
                (
                        max
                        (
                            rho_ * solidThermo_.Cp() * (1 - porosity_),
                            dimensionedScalar("minRhoCp",dimEnergy/dimTemperature/dimVolume,SMALL)
                        )
                );
            volScalarField deltaTemp(T_ * 0);
            scalar deltaTime = time_.deltaTValue();

            forAll(deltaTemp,cellI)
            {
                if (whereIs_[cellI] == 1.0 && CONV_[cellI] > 0.0)
                {
                    deltaTemp[cellI] =
                        (Tgas[cellI]
                         * (rhoCpG[cellI]
                           +rhoCpS[cellI]
                           *exp(-(CONV_[cellI] * deltaTime*(rhoCpS[cellI]+rhoCpG[cellI]))
                                /(rhoCpS[cellI]*rhoCpG[cellI])))
                           +rhoCpS[cellI]*T_[cellI]
                           *(1 - exp( -(CONV_[cellI]*deltaTime*(rhoCpS[cellI]+rhoCpG[cellI]))/(rhoCpS[cellI]*rhoCpG[cellI])))
                        )
                        / (rhoCpS[cellI]+rhoCpG[cellI]) - Tgas[cellI];
                }
                else
                {
                    deltaTemp[cellI] = 0;
                }
            }

            forAll(Sh_(),cellI)
            {
                Sh_.ref()[cellI] = deltaTemp[cellI] * rhoCpG[cellI] * whereIs_[cellI] / deltaTime;
            }

            volScalarField HT(CONV_*(T_-Tgas));
            Info << "The heat transfer subintegration info:" << nl
                 << " no subintegration min, max:" << gMin(HT) << " " << gMax(HT) << nl
                 << " subintegration    min, max:" << gMin(Sh_()) << " " << gMax(Sh_()) << endl;
        }
        else
        {
            // This works only for small CONV otherwise oscillations appear.
            volScalarField Tgas = gasThermo_.T();
            volScalarField HT(CONV_*(T_-Tgas));
            Sh_ = HT * whereIs_;
        }
    }
    return Sh_;
}

Foam::tmp<Foam::volScalarField> volPyrolysis::CONV() const
{
    Foam::tmp<Foam::volScalarField> CONVloc_ = Foam::tmp<Foam::volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "CONV",
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

    if(equilibrium_)
    {}
    else
    {
        CONVloc_ = HTmodel_->CONV()*whereIs_;
    }
    return CONVloc_;
}

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
        volScalarField gasCp = gasThermo_.Cp();

        if (equilibrium_)
        {}
        else // eqZx2uHGn018
        {
            volScalarField tempSh = hSh_();
            tempSh = gasThermo_.Cp() * (T_ - gasThermo_.T()) * Srho();
            forAll(tempSh,cellI)
            hSh_ = tempSh * whereIs_ * (1 - porosity_);
        }
    }

    return hSh_;
}

Foam::tmp<Foam::volScalarField> volPyrolysis::heatUpGas() const
{
    return heatUpGas_;
}

Foam::tmp<Foam::volScalarField> volPyrolysis::solidChemistrySh() const
{
    return chemistrySh_;
}

void volPyrolysis::info() const
{
    Info<< "\nPyrolysis: " << endl;

    Info<< indent << "Total gas mass produced  [kg] = " << addedGasMass_.value() << nl
        << indent << "Total solid mass lost    [kg] = " << lostSolidMass_.value() << nl
        << indent << "Realese rate of pyrolysis gases  [kg/s] = " << totalGasMassFlux_.value() << nl
        << indent << "Heat release rate [J/s] = " << totalHeatRR_.value() << nl;

        Info << "!!!!!TimeCHEM: " << timeChem_ << endl;
        if (timeChem_ < GREAT )
        {
            Info << indent << "Suggested chemical time step from heterogeneous reactions [s] = "
                 << timeChem_ << nl;
        }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace heterogeneousPyrolysisModels

// ************************************************************************* //
