/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
#include "physicoChemicalConstants.H"
#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
    defineTypeNameAndDebug(heterogeneousNoRadiation, 0);
    addToRadiationRunTimeSelectionTables(heterogeneousNoRadiation);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::heterogeneousNoRadiation::heterogeneousNoRadiation(const volScalarField& T)
:
    heterogeneousRadiationModel(T)
{}


Foam::radiationModels::heterogeneousNoRadiation::heterogeneousNoRadiation
(
        	const volScalarField& T,
        	const volScalarField& porosityF,
        	const List<label>& surfF,
        	const volScalarField& Ts
)
:
    heterogeneousRadiationModel(T, porosityF, surfF, Ts)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::heterogeneousNoRadiation::~heterogeneousNoRadiation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiationModels::heterogeneousNoRadiation::correct()
{}


void Foam::radiationModels::heterogeneousNoRadiation::calculate()
{}


bool Foam::radiationModels::heterogeneousNoRadiation::read()
{
    return heterogeneousRadiationModel::read();
}


Foam::tmp<Foam::volScalarField> 
Foam::radiationModels::heterogeneousNoRadiation::Rp() const
{
    return volScalarField::New
    (
        "Rp",
        mesh_,
        dimensionedScalar
        (
            constant::physicoChemical::sigma.dimensions()/dimLength,
            0
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiationModels::heterogeneousNoRadiation::Ru() const
{
    return volScalarField::Internal::New
    (
        "Ru",
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
    );
}


Foam::tmp<Foam::volScalarField> 
Foam::radiationModels::heterogeneousNoRadiation::solidSh() const
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
                dimMass/dimLength/pow3(dimTime),
                0.0
            )
        )
    );
}


// ************************************************************************* //
