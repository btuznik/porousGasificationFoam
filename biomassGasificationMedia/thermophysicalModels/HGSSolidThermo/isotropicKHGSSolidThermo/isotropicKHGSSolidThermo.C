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

#include "isotropicKHGSSolidThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(isotropicKHGSSolidThermo, 0);
    addToRunTimeSelectionTable
    (
        HGSSolidThermo,
        isotropicKHGSSolidThermo,
        mesh
    );

    addToRunTimeSelectionTable
    (
        HGSSolidThermo,
        isotropicKHGSSolidThermo,
        dictionary
    );

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isotropicKHGSSolidThermo::isotropicKHGSSolidThermo
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    interpolatedHGSSolidThermo(mesh, typeName + "Coeffs", dict),
    K_
    (
        IOobject
        (
            "K",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimTime/(dimLength*dimTemperature)
    ),
    KValues_ (Field<scalar>(subDict(typeName + "Coeffs").lookup("KValues")))
{
    correct();
}


Foam::isotropicKHGSSolidThermo::isotropicKHGSSolidThermo(const fvMesh& mesh)
:
    interpolatedHGSSolidThermo(mesh, typeName + "Coeffs"),
    K_
    (
        IOobject
        (
            "K",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimTime/(dimLength*dimTemperature)
    ),
    KValues_ (Field<scalar>(subDict(typeName + "Coeffs").lookup("KValues")))
{
    correct();
}


void Foam::isotropicKHGSSolidThermo::correct()
{

    // Correct K
    K_.primitiveFieldRef() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        KValues_
    );

    forAll(K_.boundaryField(), patchI)
    {
        K_.boundaryFieldRef()[patchI] == this->K(patchI)();
    }

    interpolatedHGSSolidThermo::calculate();
}


Foam::tmp<Foam::scalarField> Foam::isotropicKHGSSolidThermo::K
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
                KValues_
            )
        )
    );
}


bool Foam::isotropicKHGSSolidThermo::read()
{
    KValues_  = Field<scalar>(subDict(typeName + "Coeffs").lookup("KValues"));
    return true;
}


bool Foam::isotropicKHGSSolidThermo::writeData(Ostream& os) const
{
    os.writeKeyword("KValues") << KValues_ << token::END_STATEMENT << nl;
    bool ok = interpolatedHGSSolidThermo::writeData(os);

    return ok && os.good();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::isotropicKHGSSolidThermo::~isotropicKHGSSolidThermo()
{}


// ************************************************************************* //
