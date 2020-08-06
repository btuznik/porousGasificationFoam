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

#include "heterogeneousNoRadiation.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "radiationConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(heterogeneousNoRadiation, 0);

        addToRunTimeSelectionTable
        (
            heterogeneousRadiationModel,
            heterogeneousNoRadiation,
            porosity
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::heterogeneousNoRadiation::heterogeneousNoRadiation(const volScalarField& T,const volScalarField& porosityF,const List<label>& surfF, const volScalarField& Ts)
:
    heterogeneousRadiationModel(T)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::heterogeneousNoRadiation::~heterogeneousNoRadiation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::heterogeneousNoRadiation::read()
{
    return heterogeneousRadiationModel::read();
}


void Foam::radiation::heterogeneousNoRadiation::calculate()
{
    // Do nothing
}


Foam::tmp<Foam::volScalarField> Foam::radiation::heterogeneousNoRadiation::Rp() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Rp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "Rp",
                radiation::sigmaSB.dimensions()/dimLength,
                0.0
            )
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::radiation::heterogeneousNoRadiation::Ru() const
{
    return tmp<DimensionedField<scalar, volMesh> >
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "Ru",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "Ru", dimMass/dimLength/pow3(dimTime), 0.0
            )
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::radiation::heterogeneousNoRadiation::solidSh() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "solidSh",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "solidSh",
                dimMass/dimLength/pow3(dimTime),
                0.0
            )
        )
    );
}

// ************************************************************************* //
