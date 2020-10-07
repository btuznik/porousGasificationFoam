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

#include "interpolatedHGSSolidThermo.H"
#include "interpolateXY.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interpolatedHGSSolidThermo::interpolatedHGSSolidThermo
(
    const fvMesh& mesh,
    const word dictName
 )
:
    HGSSolidThermo(mesh),
    interpolateSolid(subDict(dictName)),
    dict_(subDict(dictName))
{
    calculate();
}


Foam::interpolatedHGSSolidThermo::interpolatedHGSSolidThermo
(
    const fvMesh& mesh,
    const word dictName,
    const dictionary& dict
 )
:
    HGSSolidThermo(mesh, dict),
    interpolateSolid(subDict(dictName)),
    dict_(subDict(dictName))
{
    calculate();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interpolatedHGSSolidThermo::~interpolatedHGSSolidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::interpolatedHGSSolidThermo::calculate()
{
    // Correct rho
    rho_.primitiveFieldRef() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        rhoValues_
    );

    forAll(rho_.boundaryField(), patchI)
    {
        rho_.boundaryFieldRef()[patchI] == this->rho(patchI)();
    }

    // Correct emissivity
    emissivity_.primitiveFieldRef() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        emissivityValues_
    );

    forAll(emissivity_.boundaryField(), patchI)
    {
        emissivity_.boundaryFieldRef()[patchI] == this->emissivity(patchI)();
    }


    // Correct absorptivity
    kappa_.primitiveFieldRef() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        kappaValues_
    );

    forAll(kappa_.boundaryField(), patchI)
    {
        kappa_.boundaryFieldRef()[patchI] == this->kappa(patchI)();
    }


    // Correct scatter
    sigmaS_.primitiveFieldRef() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        sigmaSValues_
    );

    forAll(sigmaS_.boundaryField(), patchI)
    {
        sigmaS_.boundaryFieldRef()[patchI] == this->sigmaS(patchI)();
    }
}


Foam::tmp<Foam::volScalarField> Foam::interpolatedHGSSolidThermo::Cp() const
{
    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                "Cp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimEnergy/(dimMass*dimTemperature)
        )
    );
    volScalarField& Cp = tCp.ref();

    Cp.primitiveFieldRef() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        cpValues_
    );

    forAll(Cp.boundaryField(), patchI)
    {
        Cp.boundaryFieldRef()[patchI] == this->Cp(patchI)();
    }

    return tCp;
}


Foam::tmp<Foam::volScalarField> Foam::interpolatedHGSSolidThermo::Hf() const
{
    tmp<volScalarField> tHf
    (
        new volScalarField
        (
            IOobject
            (
                "Hf",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimEnergy/dimMass
        )
    );
    volScalarField& Hf = tHf.ref();

    Hf.primitiveFieldRef() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        HfValues_
    );

    forAll(Hf.boundaryField(), patchI)
    {
        Hf.boundaryFieldRef()[patchI] == this->Hf(patchI)();
    }

    return tHf;
}


Foam::tmp<Foam::scalarField> Foam::interpolatedHGSSolidThermo::rho
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            interpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                rhoValues_
            )
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::interpolatedHGSSolidThermo::Cp
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            interpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                cpValues_
            )
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::interpolatedHGSSolidThermo::Hf
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            interpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                HfValues_
            )
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::interpolatedHGSSolidThermo::emissivity
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            interpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                emissivityValues_
            )
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::interpolatedHGSSolidThermo::kappa
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            interpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                kappaValues_
            )
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::interpolatedHGSSolidThermo::sigmaS
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            interpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                sigmaSValues_
            )
        )
    );
}


bool Foam::interpolatedHGSSolidThermo::read()
{
    return read(dict_);
}


bool Foam::interpolatedHGSSolidThermo::read(const dictionary& dict)
{
    bool ok = interpolateSolid::read(dict);
    return ok;
}


bool Foam::interpolatedHGSSolidThermo::writeData(Ostream& os) const
{
    bool ok = HGSSolidThermo::writeData(os);
    ok = interpolateSolid::writeData(os);

    return ok && os.good();
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const interpolatedHGSSolidThermo& s)
{
    s.writeData(os);
    return os;
}


// ************************************************************************* //
