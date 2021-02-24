/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
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

#include "heterogeneousConstantAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace heterogeneousAbsorptionEmissionModels
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
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::heterogeneousAbsorptionEmissionModels::
heterogeneousConstantAbsorptionEmission::heterogeneousConstantAbsorptionEmission
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
    E_(coeffsDict_.lookup("E")),
    borderL_(coeffsDict_.lookup("borderL"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::heterogeneousAbsorptionEmissionModels
::heterogeneousConstantAbsorptionEmission::~heterogeneousConstantAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModels::
heterogeneousConstantAbsorptionEmission::aCont(const label bandI) const
{
    return volScalarField::New
    (
        "a",
        mesh_,
        a_
    );
}

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModels::
heterogeneousConstantAbsorptionEmission::asCont(const label bandI) const
{
    return volScalarField::New
    (
        "as",
        mesh_,
        as_
    );
}

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModels::
heterogeneousConstantAbsorptionEmission::borderAsCont(const label bandI) const
{
    return volScalarField::New
    (
        "borderAs",
        mesh_,
        borderAs_
    );
}

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModels::
heterogeneousConstantAbsorptionEmission::ECont(const label bandI) const
{
    return volScalarField::New
    (
        "E",
        mesh_,
        E_
    );
}

Foam::dimensionedScalar
Foam::radiationModels::heterogeneousAbsorptionEmissionModels::
heterogeneousConstantAbsorptionEmission::borderL(const label bandI) const
{
    return  dimensionedScalar("borderL", borderL_.dimensions(), borderL_.value());
}
// ************************************************************************* //
