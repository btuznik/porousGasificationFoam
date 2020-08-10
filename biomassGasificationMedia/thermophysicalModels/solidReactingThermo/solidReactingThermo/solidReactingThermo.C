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

#include "solidReactingThermo.H"
#include "fvMesh.H"


/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(solidReactingThermo, 0);
    defineRunTimeSelectionTable(solidReactingThermo, fvMesh);
    defineRunTimeSelectionTable(solidReactingThermo, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidReactingThermo::solidReactingThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicSolidThermo(mesh, phaseName),
    rho_
    (
        IOobject
        (
            phasePropertyName("rho"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimDensity
    )
{}


Foam::solidReactingThermo::solidReactingThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    basicSolidThermo(mesh, dict, phaseName),
    rho_
    (
        IOobject
        (
            phasePropertyName("rho"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimDensity
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidReactingThermo> Foam::solidReactingThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicSolidThermo::New<solidReactingThermo>(mesh, phaseName);
}


Foam::autoPtr<Foam::solidReactingThermo> Foam::solidReactingThermo::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
{
    return basicSolidThermo::New<solidReactingThermo>(mesh, dict, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidReactingThermo::~solidReactingThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::solidReactingThermo::rho() const
{
    return rho_;
}


Foam::tmp<Foam::scalarField> Foam::solidReactingThermo::rho(const label patchi) const
{
    return rho_.boundaryField()[patchi];
}


Foam::volScalarField& Foam::solidReactingThermo::rho()
{
    return rho_;
}


bool Foam::solidReactingThermo::read()
{
    return regIOobject::read();
}


// ************************************************************************* //
