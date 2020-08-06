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
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(heterogeneousRadiationModel, 0);
        defineRunTimeSelectionTable(heterogeneousRadiationModel, dictionary);
        defineRunTimeSelectionTable(heterogeneousRadiationModel, porosity);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::heterogeneousRadiationModel::heterogeneousRadiationModel(const volScalarField& T)
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
    radiation_(false),
    coeffs_(dictionary::null),
    solverFreq_(0),
    firstIter_(true),
    absorptionEmission_(NULL),
    scatter_(NULL)
{}


Foam::radiation::heterogeneousRadiationModel::heterogeneousRadiationModel
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
    radiation_(lookup("radiation")),
    coeffs_(dictionary::null),
    solverFreq_(readLabel(lookup("solverFreq"))),
    firstIter_(true),
    absorptionEmission_(heterogeneousAbsorptionEmissionModel::New(*this, mesh_)),
    scatter_(scatterModel::New(*this, mesh_))
{
    solverFreq_ = max(1, solverFreq_);
}

Foam::radiation::heterogeneousRadiationModel::heterogeneousRadiationModel(const volScalarField& T, const volScalarField& porosityF,  const List<label>& surfF,const volScalarField& Ts)
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
    radiation_(false),
    coeffs_(dictionary::null),
    solverFreq_(0),
    firstIter_(true),
    absorptionEmission_(NULL),
    scatter_(NULL)
{}



// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiation::heterogeneousRadiationModel::~heterogeneousRadiationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::heterogeneousRadiationModel::read()
{
    if (regIOobject::read())
    {
        lookup("radiation") >> radiation_;

        lookup("solverFreq") >> solverFreq_,
        solverFreq_ = max(1, solverFreq_);

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiation::heterogeneousRadiationModel::correct()
{
    if (!radiation_)
    {
        return;
    }

    if (firstIter_ || (time_.timeIndex() % solverFreq_ == 0))
    {
        calculate();
        firstIter_ = false;
    }
}


Foam::tmp<Foam::fvScalarMatrix> Foam::radiation::heterogeneousRadiationModel::Sh
(
    basicThermo& thermo
) const
{
    volScalarField& h = thermo.h();
    const volScalarField Cp(thermo.Cp());
    const volScalarField T3(pow3(T_));

    return
    (
        Ru()
      - fvm::Sp(4.0*Rp()*T3/Cp, h)
      - Rp()*T3*(T_ - 4.0*h/Cp)
    );
}


Foam::tmp<Foam::fvScalarMatrix> Foam::radiation::heterogeneousRadiationModel::Shs
(
    basicThermo& thermo
) const
{
    volScalarField& hs = thermo.hs();
    const volScalarField Cp(thermo.Cp());
    const volScalarField T3(pow3(T_));

    return
    (
        Ru()
      - fvm::Sp(4.0*Rp()*T3/Cp, hs)
      - Rp()*T3*(T_ - 4.0*hs/Cp)
    );
}


// ************************************************************************* //
