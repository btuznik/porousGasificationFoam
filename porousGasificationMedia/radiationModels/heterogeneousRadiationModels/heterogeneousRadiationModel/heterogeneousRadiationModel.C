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

#include "heterogeneousRadiationModel.H"
#include "heterogeneousAbsorptionEmissionModel.H"
#include "scatterModel.H"
#include "sootModel.H"
#include "fvmSup.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heterogeneousRadiationModel, 0);
    defineRunTimeSelectionTable(heterogeneousRadiationModel, T);
    defineRunTimeSelectionTable(heterogeneousRadiationModel, porosity);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject
Foam::heterogeneousRadiationModel::createIOobject(const fvMesh& mesh) const
{
    IOobject io
    (
        "radiationProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}


void Foam::heterogeneousRadiationModel::initialise()
{
    solverFreq_ = max(1, lookupOrDefault<label>("solverFreq", 1));

    heterogeneousAbsorptionEmission_.reset
    (
        radiationModels::heterogeneousAbsorptionEmissionModel::New(*this, mesh_).ptr()
    );

    scatter_.reset(radiationModels::scatterModel::New(*this, mesh_).ptr());

    soot_.reset(radiationModels::sootModel::New(*this, mesh_).ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heterogeneousRadiationModel::heterogeneousRadiationModel(const volScalarField& T)
:
    IOdictionary
    (
        IOobject
        (
            "radiationProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    Ts_(T),
    coeffs_(dictionary::null),
    solverFreq_(0),
    firstIter_(true),
    heterogeneousAbsorptionEmission_(nullptr),
    scatter_(nullptr),
    soot_(nullptr)
{
    initialise();
}


Foam::heterogeneousRadiationModel::heterogeneousRadiationModel
(
    const word& type,
    const volScalarField& T
)
:
    IOdictionary
    (
        IOobject
        (
            "radiationProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    Ts_(T),
    coeffs_(dictionary::null),
    solverFreq_(readLabel(lookup("solverFreq"))),
    firstIter_(true),
    heterogeneousAbsorptionEmission_(nullptr),
    scatter_(nullptr),
    soot_(nullptr)
{
    initialise();
}


Foam::heterogeneousRadiationModel::heterogeneousRadiationModel
(
    const volScalarField& T,
    const volScalarField& porosityF,
    const List<label>& surfF,
    const volScalarField& Ts
)
:
    IOdictionary
    (
        IOobject
        (
            "radiationProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    Ts_(Ts),
    coeffs_(subOrEmptyDict(type() + "Coeffs")),
    solverFreq_(1),
    firstIter_(true),
    heterogeneousAbsorptionEmission_(nullptr),
    scatter_(nullptr),
    soot_(nullptr)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::heterogeneousRadiationModel::~heterogeneousRadiationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::heterogeneousRadiationModel::correct()
{
    if (firstIter_ || (time_.timeIndex() % solverFreq_ == 0))
    {
        calculate();
        firstIter_ = false;
    }

    if (!soot_.empty())
    {
        soot_->correct();
    }
}


bool Foam::heterogeneousRadiationModel::read()
{
    if (regIOobject::read())
    {
        coeffs_ = subOrEmptyDict(type() + "Coeffs");

        solverFreq_ = lookupOrDefault<label>("solverFreq", 1);
        solverFreq_ = max(1, solverFreq_);

        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::fvScalarMatrix> Foam::heterogeneousRadiationModel::Sh
(
    const basicThermo& thermo,
    const volScalarField& he
) const
{
    const volScalarField Cpv(thermo.Cpv());
    const volScalarField T3(pow3(T_));

    return
    (
        Ru()
      - fvm::Sp(4.0*Rp()*T3/Cpv, he)
      - Rp()*T3*(T_ - 4.0*he/Cpv)
    );
}


Foam::tmp<Foam::fvScalarMatrix> Foam::heterogeneousRadiationModel::ST
(
    const dimensionedScalar& rhoCp,
    volScalarField& T
) const
{
    return
    (
        Ru()/rhoCp
      - fvm::Sp(Rp()*pow3(T)/rhoCp, T)
    );
}


const Foam::radiationModels::heterogeneousAbsorptionEmissionModel&
Foam::heterogeneousRadiationModel::heterogeneousAbsorptionEmission() const
{
    if (!heterogeneousAbsorptionEmission_.valid())
    {
        FatalErrorInFunction
            << "Requested radiation heterogeneousAbsorptionEmission model, but model is "
            << "not active" << abort(FatalError);
    }

    return heterogeneousAbsorptionEmission_;
}


const Foam::radiationModels::sootModel&
Foam::heterogeneousRadiationModel::soot() const
{
    if (!soot_.valid())
    {
        FatalErrorInFunction
            << "Requested radiation sootModel model, but model is "
            << "not active" << abort(FatalError);
    }

    return soot_();
}


// ************************************************************************* //
