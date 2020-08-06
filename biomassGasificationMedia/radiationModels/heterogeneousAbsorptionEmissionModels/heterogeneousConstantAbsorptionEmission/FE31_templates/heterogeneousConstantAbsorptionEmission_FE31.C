/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "heterogeneousConstantAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiationModels
    {
        defineTypeNameAndDebug(heterogeneousConstantAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            heterogeneousAbsorptionEmissionModel,
            heterogeneousConstantAbsorptionEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::heterogeneousConstantAbsorptionEmission::heterogeneousConstantAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    heterogeneousAbsorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    a_(coeffsDict_.lookup("a")),
    as_(coeffsDict_.lookup("as")),
    borderAs_(coeffsDict_.lookup("borderAs")),
    e_(coeffsDict_.lookup("e")),
    es_(coeffsDict_.lookup("es")),
    borderEs_(coeffsDict_.lookup("borderEs")),
    E_(coeffsDict_.lookup("E"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::heterogeneousConstantAbsorptionEmission::~heterogeneousConstantAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousConstantAbsorptionEmission::aCont(const label bandI) const
{
    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "a",
                mesh_.time().timeName(),
                mesh_,  // HR 101116 -- Todo: Need to bring object registry here!
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            a_
        )
    );

    return ta;
}

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousConstantAbsorptionEmission::asCont(const label bandI) const
{
    tmp<volScalarField> tas
    (
        new volScalarField
        (
            IOobject
            (
                "as",
                mesh_.time().timeName(),
                mesh_,  // HR 101116 -- Todo: Need to bring object registry here!
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            as_
        )
    );

    return tas;
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousConstantAbsorptionEmission::borderAsCont(const label bandI) const
{
    tmp<volScalarField> tBorderAs
    (
        new volScalarField
        (
            IOobject
            (
                "borderAs",
                mesh_.time().timeName(),
                mesh_,  // HR 101116 -- Todo: Need to bring object registry here!
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            borderAs_
        )
    );

    return tBorderAs;
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousConstantAbsorptionEmission::eCont(const label bandI) const
{
    tmp<volScalarField> te
    (
        new volScalarField
        (
            IOobject
            (
                "e",
                mesh_.time().timeName(),
                mesh_,  // HR 101116 -- Todo: Need to bring object registry here!
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            e_
        )
    );

    return te;
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousConstantAbsorptionEmission::esCont(const label bandI) const
{
    tmp<volScalarField> tes
    (
        new volScalarField
        (
            IOobject
            (
                "es",
                mesh_.time().timeName(),
                mesh_,  // HR 101116 -- Todo: Need to bring object registry here!
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            es_
        )
    );

    return tes;
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousConstantAbsorptionEmission::borderEsCont(const label bandI) const
{
    tmp<volScalarField> tBorderEs
    (
        new volScalarField
        (
            IOobject
            (
                "borderEs",
                mesh_.time().timeName(),
                mesh_,  // HR 101116 -- Todo: Need to bring object registry here!
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            borderEs_
        )
    );

    return tBorderEs;
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousConstantAbsorptionEmission::ECont(const label bandI) const
{
    tmp<volScalarField> tE
    (
        new volScalarField
        (
            IOobject
            (
                "E",
                mesh_.time().timeName(),
                mesh_,  // HR 101116 -- Todo: Need to bring object registry here!
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            E_
        )
    );

    return tE;
}


// ************************************************************************* //
