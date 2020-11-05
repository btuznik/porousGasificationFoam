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

/** @file
 * Base class for pyrolysis models.
 * */
#include "heterogeneousPyrolysisModel.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heterogeneousPyrolysisModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(heterogeneousPyrolysisModel, 0);
defineRunTimeSelectionTable(heterogeneousPyrolysisModel, mesh);
defineRunTimeSelectionTable(heterogeneousPyrolysisModel, noRadiation);
defineRunTimeSelectionTable(heterogeneousPyrolysisModel, radiation);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
// Something that reads from pyrolysis dictionary.
void heterogeneousPyrolysisModel::readPyrolysisControls()
{
    reactionDeltaMin_ =
        coeffs_.lookupOrDefault<scalar>("reactionDeltaMin", 0.0);
}

bool heterogeneousPyrolysisModel::read()
{
        readPyrolysisControls();
        return true;
}

bool heterogeneousPyrolysisModel::read(const dictionary& dict)
{
        readPyrolysisControls();
        return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

heterogeneousPyrolysisModel::heterogeneousPyrolysisModel
(
    const word& modelType,
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "pyrolysisProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    time_(mesh.time()),
    active_(lookup("active")),
    infoOutput_(true),
    coeffs_(subOrEmptyDict("pyrolysisCoeffs"))
{
    if (active_)
    {
        read();
    }
}

heterogeneousPyrolysisModel::heterogeneousPyrolysisModel
(
    const word& modelType,
    const fvMesh& mesh,
    HGSSolidThermo& solidThermo,
    psiReactionThermo& gasThermo,
    volScalarField& whereIs
)
:
    IOdictionary
    (
        IOobject
        (
            "pyrolysisProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    time_(mesh.time()),
    active_(lookup("active")),
    infoOutput_(true),
    coeffs_(subOrEmptyDict("pyrolysisCoeffs"))
{
    if (active_)
    {
        read();
    }
}

heterogeneousPyrolysisModel::heterogeneousPyrolysisModel
(
    const word& modelType,
    const fvMesh& mesh,
    HGSSolidThermo& solidThermo,
    psiReactionThermo& gasThermo,
    volScalarField& whereIs,
    volScalarField& radiation
)
:
    IOdictionary
    (
        IOobject
        (
            "pyrolysisProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    time_(mesh.time()),
    active_(lookup("active")),
    infoOutput_(true),
    coeffs_(subOrEmptyDict("pyrolysisCoeffs"))
{
    if (active_)
    {
        read();
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

heterogeneousPyrolysisModel::~heterogeneousPyrolysisModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void heterogeneousPyrolysisModel::evolve()
{
    if (active_)
    {

        // Pre-evolve
        preEvolveRegion();
        // Increment the region equations up to the new time level
        evolveRegion();
        // Provide some feedback
        if (infoOutput_)
        {
            Info<< incrIndent;
            info();
            Info<< endl << decrIndent;
        }
    }
}


void heterogeneousPyrolysisModel::preEvolveRegion()
{
}

void heterogeneousPyrolysisModel::evolveRegion()
{
}

scalar heterogeneousPyrolysisModel::solidRegionDiffNo() const
{
    return VSMALL;
}

scalar heterogeneousPyrolysisModel::maxTime() const
{
    return VSMALL;
}


scalar heterogeneousPyrolysisModel::maxDiff() const
{
    return GREAT;
}

Switch heterogeneousPyrolysisModel::equilibrium() const
{
    return false;
}

Foam::tmp<Foam::volScalarField> heterogeneousPyrolysisModel::Srho() const
{
    notImplemented
    (
        "tmp<DimensionedField<scalar, volMesh> > heterogeneousPyrolysisModel::Srho() const"
    )

    return Foam::tmp<Foam::volScalarField>(NULL);
}


Foam::tmp<Foam::volScalarField> heterogeneousPyrolysisModel::Srho(const label) const
{
    notImplemented
    (
        "tmp<DimensionedField<scalar, volMesh> > heterogeneousPyrolysisModel::Srho"
        "(const label) const"
    )

    return Foam::tmp<Foam::volScalarField>(NULL);
}


Foam::tmp<Foam::volScalarField> heterogeneousPyrolysisModel::heatTransfer()
{
    notImplemented
    (
        "tmp<DimensionedField<scalar, volMesh> > heterogeneousPyrolysisModel::heatTransfer()"
    )

    return Foam::tmp<Foam::volScalarField>(NULL);
}

Foam::tmp<Foam::volScalarField> heterogeneousPyrolysisModel::heatUpGas() const
{
    notImplemented
    (
        "tmp<DimensionedField<scalar, volMesh> > heterogeneousPyrolysisModel::heatUpGas() const"
    )

    return Foam::tmp<Foam::volScalarField>(NULL);
}

Foam::tmp<Foam::volScalarField> heterogeneousPyrolysisModel::solidChemistrySh() const
{
    notImplemented
    (
        "tmp<DimensionedField<scalar, volMesh> > heterogeneousPyrolysisModel::solidChemistrySh() const"
    )

    return Foam::tmp<Foam::volScalarField>(NULL);
}

Foam::tmp<Foam::volScalarField> heterogeneousPyrolysisModel::viscosityDrop()
{
    notImplemented
    (
        "tmp<DimensionedField<scalar, volMesh> > heterogeneousPyrolysisModel::viscosityDrop()"
    )

    return Foam::tmp<Foam::volScalarField>(NULL);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heterogeneousPyrolysisModels
} // End namespace Foam

// ************************************************************************* //
