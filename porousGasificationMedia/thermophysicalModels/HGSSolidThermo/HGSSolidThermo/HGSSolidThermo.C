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

#include "HGSSolidThermo.H"


/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(HGSSolidThermo, 0);
    defineRunTimeSelectionTable(HGSSolidThermo, mesh);
    defineRunTimeSelectionTable(HGSSolidThermo, dictionary);
    defineRunTimeSelectionTable(HGSSolidThermo, gasPhase);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HGSSolidThermo::HGSSolidThermo
(
    const fvMesh& mesh
)
:

    IOdictionary
    (
        IOobject
        (
            "solidThermophysicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    T_
    (
        IOobject
        (
            "Ts",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    rho_
    (
        IOobject
        (
            "rhos",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimMass/dimVolume
    ),
    kappa_
    (
        IOobject
        (
            "kappas",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless/dimLength
    ),
    sigmaS_
    (
        IOobject
        (
            "sigmaS",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless/dimLength
    ),
    emissivity_
    (
        IOobject
        (
            "emissivity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless
    )
{}


Foam::HGSSolidThermo::HGSSolidThermo
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
        
    IOdictionary
    (
        IOobject
        (
            "solidThermophysicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(mesh),
    T_
    (
        IOobject
        (
            "Ts",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    rho_
    (
        IOobject
        (
            "rhos",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimMass/dimVolume
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless/dimLength
    ),
    sigmaS_
    (
        IOobject
        (
            "sigmaS",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless/dimLength
    ),
    emissivity_
    (
        IOobject
        (
            "emissivity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless
    )
{}

Foam::HGSSolidThermo::HGSSolidThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const PtrList<volScalarField>& gasPhaseGases
)
:

      
    IOdictionary
    (
        IOobject
        (
            "solidThermophysicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(mesh),
    T_
    (
        IOobject
        (
            "Ts",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    rho_
    (
        IOobject
        (
            "rhos",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimMass/dimVolume
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless/dimLength
    ),
    sigmaS_
    (
        IOobject
        (
            "sigmaS",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless/dimLength
    ),
    emissivity_
    (
        IOobject
        (
            "emissivity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless
    )
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::HGSSolidThermo::~HGSSolidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fvMesh& Foam::HGSSolidThermo::mesh() const
{
    return mesh_;
}

Foam::volScalarField& Foam::HGSSolidThermo::T()
{
    return T_;
}


const Foam::volScalarField& Foam::HGSSolidThermo::T() const
{
    return T_;
}


const Foam::volScalarField& Foam::HGSSolidThermo::rho() const
{
    return rho_;
}


Foam::volScalarField& Foam::HGSSolidThermo::rho()
{
    return rho_;
}


const Foam::volScalarField& Foam::HGSSolidThermo::kappa() const
{
    return kappa_;
}


const Foam::volScalarField& Foam::HGSSolidThermo::sigmaS() const
{
    return sigmaS_;
}


const Foam::volScalarField& Foam::HGSSolidThermo::emissivity() const
{
    return emissivity_;
}


const Foam::volScalarField&  Foam::HGSSolidThermo::K() const
{
    notImplemented("HGSSolidThermo::K()");
    return volScalarField::null();
}


const Foam::volSymmTensorField& Foam::HGSSolidThermo::directionalK() const
{
    notImplemented("HGSSolidThermo::directionalK()");
    return const_cast<volSymmTensorField&>(volSymmTensorField::null());
}


Foam::basicSolidMixture& Foam::HGSSolidThermo::composition()
{
    notImplemented("HGSSolidThermo::composition()");
    return *reinterpret_cast<basicSolidMixture*>(0);
}


const Foam::basicSolidMixture& Foam::HGSSolidThermo::composition() const
{
    notImplemented("HGSSolidThermo::composition() const");
    return *reinterpret_cast<const basicSolidMixture*>(0);
}


Foam::tmp<Foam::volScalarField> Foam::HGSSolidThermo::hs() const
{
    notImplemented("HGSSolidThermo::hs()");
    return volScalarField::null();
}


Foam::tmp<Foam::scalarField> Foam::HGSSolidThermo::hs(const label patchI)
const
{
    notImplemented("HGSSolidThermo::hs(const label)");
    return scalarField::null();
}


Foam::tmp<Foam::scalarField> Foam::HGSSolidThermo::K
(
    const label patchI
)const
{
    notImplemented("HGSSolidThermo::K(const label)");
    return scalarField::null();
}


Foam::tmp<Foam::symmTensorField> Foam::HGSSolidThermo::directionalK
(
    const label
)const
{
    notImplemented("HGSSolidThermo::directionalK(const label)");
    return symmTensorField::null();
}


bool Foam::HGSSolidThermo::read()
{
    return regIOobject::read();
}


bool Foam::HGSSolidThermo::writeData(Ostream& os) const
{
    return true;
}

// ************************************************************************* //
