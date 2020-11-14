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

    if (debug)
    {
        Info<< "volPyrolysis:+solveSpeciesMass()" << endl;
    }

    if (active_)
    {
        volScalarField Mt(0.0*Ym_[0]);
        volScalarField rhoLoc(max(rho_,dimensionedScalar("minRho",dimMass/dimVolume,SMALL)));

        for (label i=0; i<Ys_.size(); i++)
        {

            volScalarField& Yi = Ys_[i];
            Yi.ref() *= whereIs_;

            fvScalarMatrix YsEqn
            (
                fvm::ddt(rhoLoc,Yi)
             ==
                solidChemistry_->RRs(i)
            );

            YsEqn.relax();
            YsEqn.solve("Ys");

            Yi.max(0.0);
            Mt += Yi*rho_;
        }

        for (label i=0; i<Ys_.size(); i++)
        {
            Ym_[i] = whereIs_ * Ys_[i] * rho_;
            Ys_[i] = Ym_[i]/max(Mt,dimensionedScalar("minMass", dimMass/dimVolume, SMALL));
            Ym_[i] *= (1. - porosity_);
        }
    }
}


void volPyrolysis::solveEnergy()
{
    if (debug)
    {
        Info<< "volPyrolysis::solveEnergy()" << endl;
    }

    if (active_)
    {

        volTensorField composedK(K_*(1-porosity_)*anisotropyK_);
        //radiationSh_ = radiation_;

        if (equilibrium_)
        {}
        else
        {
            volScalarField rhoCp(max(rho_*solidThermo_.Cp()*(1-porosity_),dimensionedScalar("minRhoCp",dimEnergy/dimTemperature/dimVolume,SMALL)));
            forAll(surfF_,cellI)
            {
                forAll(mesh_.cellCells()[surfF_[cellI]],cellJ)
                {
                    if (whereIs_[mesh_.cellCells()[surfF_[cellI]][cellJ]] == 0) T_[mesh_.cellCells()[surfF_[cellI]][cellJ]]=T_[surfF_[cellI]];
                }
            }
            fvScalarMatrix TEqn
            (
                fvm::ddt(rhoCp, T_)
              - fvm::laplacian(composedK, T_)
             ==
                chemistrySh_
              - heatTransfer()()
             // + radiationSh_
            );
            TEqn.solve();
        }

        T_ *= whereIs_;

        scalar minTemp = GREAT;
        scalar maxTemp = -GREAT;
        forAll(T_,cellI)
        {
            if (whereIs_[cellI] == 1 )
            {
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
        Info<< " pyrolysis min/max(T) = " << minTemp << ", " << maxTemp << endl;
    }
}

void volPyrolysis::calculateMassTransfer()
{
    if (infoOutput_)
    {
        totalGasMassFlux_ = fvc::domainIntegrate(solidChemistry_->RRg());
        totalHeatRR_ = fvc::domainIntegrate(chemistrySh_);

        addedGasMass_ +=
            fvc::domainIntegrate(solidChemistry_->RRg())*time_.deltaT();
        lostSolidMass_ +=
            fvc::domainIntegrate(solidChemistry_->RRs())*time_.deltaT();
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
    solidChemistry_(BasicSolidChemistryModel<HGSSolidThermo>::New(solidThermo_, Ygas_, gasThermo.thermoName())),
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
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass, 0.0)
    ),
    T_(solidThermo_.T()),
    equilibrium_(false),
    Tsolid_
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
    viscosityDropFactor_(1.0),
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
    viscosityDrop_
    (
        IOobject
        (
            "viscosityDrop",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        scalar(1.0)
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
    surfF_(),
    lostSolidMass_(dimensionedScalar("zero", dimMass, 0.0)),
    addedGasMass_(dimensionedScalar("zero", dimMass, 0.0)),
    totalGasMassFlux_(dimensionedScalar("zero", dimMass/dimTime, 0.0)),
    totalHeatRR_(dimensionedScalar("zero", dimEnergy/dimTime, 0.0)),
    timeChem_(1.0),
    maxDT_(1.0)
{

    HTC_ = HTC();
    rho0_.ref() = rho_.ref();
    viscosityDropFactor_ = coeffs().lookupOrDefault("viscosityDropFactor",1.0);
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
        Ym_[fieldI].ref() = Ys_[fieldI].ref()*rho_.ref();

        // Msolid is used for renormalization of mass fractions of solid species
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
                dimensionedScalar("zero",dimMass,0.0)
            )
        );
        forAll(rho_,cellI)
        {
        Msolid_[fieldI][cellI] = Ys_[fieldI][cellI]*rho_[cellI]*mesh_.V()[cellI]*(1-porosity_[cellI]);
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

    forAll(rho_,cellI)
    {
        Msolidtotal_[cellI]=0.0;

        for (label i=0; i<Ys_.size(); i++)
        {
             Msolid_[i][cellI]= Ys_[i][cellI]*rho_[cellI]*mesh_.V()[cellI]*(1-porosity_[cellI]);
             Msolidtotal_[cellI] += Msolid_[i][cellI];
        }

        for (label i=0; i<Ys_.size(); i++)
        {
            if(Msolidtotal_[cellI] > 0.0)
            {
               Ys_[i][cellI]=(Msolid_[i][cellI]/(Msolidtotal_[cellI]));
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
    solidThermo_(HGSSolidThermo::New
        (
         mesh,
         IOdictionary
             (
                 IOobject
                 (
                     "solidThermophysicalProperties",
                     mesh.time().constant(),
                     mesh,
                     IOobject::MUST_READ,
                     IOobject::NO_WRITE,
                     false
                 )
             ),
            Ygas_
                )()),
    solidChemistry_(BasicSolidChemistryModel<HGSSolidThermo>::New(solidThermo_, Ygas_, gasThermo.thermoName())),
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
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass, 0.0)
    ),
    T_(solidThermo_.T()),
    equilibrium_(false),
    Tsolid_
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
    viscosityDropFactor_(1.0),
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
    viscosityDrop_
    (
        IOobject
        (
            "viscosityDrop",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        scalar(1.0)
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
    surfF_(),
    lostSolidMass_(dimensionedScalar("zero", dimMass, 0.0)),
    addedGasMass_(dimensionedScalar("zero", dimMass, 0.0)),
    totalGasMassFlux_(dimensionedScalar("zero", dimMass/dimTime, 0.0)),
    totalHeatRR_(dimensionedScalar("zero", dimEnergy/dimTime, 0.0)),
    timeChem_(1.0),
    maxDT_(1.0)
{

    HTC_ = HTC();
    rho0_.ref() = rho_.ref();
    viscosityDropFactor_ = coeffs().lookupOrDefault("viscosityDropFactor",1.0);

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
          Ym_[fieldI].ref() = Ys_[fieldI].ref()*rho_;
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

const List<label>& volPyrolysis::surf() const
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
        DiNum = max(KrhoCpbyDelta.ref()).value()*time_.deltaTValue()/NumCprho;

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

// Prepares computational domain for evolving pyrolysis.
void volPyrolysis::preEvolveRegion() {
    // Added for completeness
    heterogeneousPyrolysisModel::preEvolveRegion();
    // Iterates over every cell and sets cells containing solid phase
    // as reacting cell.
    forAll(Tsolid_, cellI)
    {
        if ( active_ && whereIs_[cellI] != 0)
        {
            // sets this cell as reacting.
            solidChemistry_->setCellReacting(cellI, true);
        }
        else
        {
            // sets this cell as non-reacting
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
        HTC_ = HTC();
    }

    timeChem_ = solidChemistry_->solve
    (
        time_.value() - time_.deltaTValue(),
        time_.deltaTValue()
    );

    heatUpGas_ = heatUpGasCalc()();
    chemistrySh_ = solidChemistry_->Sh()();
    solveSpeciesMass();
    evolvePorosity();
    solidThermo_.correct();

    if (equilibrium_)
    {}
    else
    {
        for (int nonOrth=0; nonOrth<=nNonOrthCorr_; nonOrth++)
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

        Info << "  Evolve porosity porosity Source max/min " << max(porositySource_).value() << " / "<< min(porositySource_).value() << endl;
        volScalarField& por = porosity_;

        fvScalarMatrix porosityEqn
        (
            fvm::ddt(por)
         ==
            porositySource_
        );

        porosityEqn.solve("porosity");

        forAll(porosity_,cellI)
        {
            if (porosity_[cellI] > 0.9999)
            {
                porosity_[cellI] = 1.0;
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

        if (j>-1)
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
        volScalarField HT(HTC_ * (Tsolid_-Tgas));
        Sh_ = HT * whereIs_;
    }
    return Sh_;
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
        volScalarField pGas = gasThermo_.p();

        if (equilibrium_)
        {}
        else
        {
            volScalarField tempSh = hSh_();
            volScalarField tempDiff = ((whereIs_*T_)+(whereIsNot_*Tgas));
            const speciesTable& gasTable = solidChemistry_->gasTable();
            const dimensionedScalar stala("one",dimless,1.0);
            forAll(gasTable,gasI)
            {
                tempSh += (solidChemistry_().gasHs(pGas, tempDiff, gasI) - solidChemistry_().gasHs(pGas, Tgas,gasI))*Srho(gasI);
            }
            hSh_ = tempSh * whereIs_;
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

Foam::tmp<Foam::volScalarField> volPyrolysis::viscosityDrop()
{
    viscosityDrop_=1.0;

    forAll(radiation_,cellI)
    {
        if (radiation_[cellI]*whereIs_[cellI] > 0)
        {
            viscosityDrop_[cellI] = viscosityDropFactor_;
            forAll(mesh_.cellCells()[cellI],cellJ)
            {
                if (whereIsNot_[mesh_.cellCells()[cellI][cellJ]] > 0)
                {
                    viscosityDrop_[mesh_.cellCells()[cellI][cellJ]] = viscosityDropFactor_;
                }
            }
        }
    }
    return viscosityDrop_;
}

void volPyrolysis::info() const
{
    Info<< "\nPyrolysis: " << endl;

    Info<< indent << "Total gas mass produced  [kg] = "
        << addedGasMass_.value() << nl
        << indent << "Total solid mass lost    [kg] = "
        << lostSolidMass_.value() << nl
        << indent << "Total realese of pyrolysis gases  [kg/s] = "
        << totalGasMassFlux_.value() << nl
        << indent << "Total heat release rate [J/s] = "
        << totalHeatRR_.value() << nl;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace heterogeneousPyrolysisModels

// ************************************************************************* //
