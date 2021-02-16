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

#include "heterogeneousAbsorptionEmissionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiationModels
    {
        defineTypeNameAndDebug(heterogeneousAbsorptionEmissionModel, 0);
        defineRunTimeSelectionTable(heterogeneousAbsorptionEmissionModel, dictionary);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::heterogeneousAbsorptionEmissionModel::heterogeneousAbsorptionEmissionModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dict_(dict),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiationModels::heterogeneousAbsorptionEmissionModel::~heterogeneousAbsorptionEmissionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::a(const label bandI) const
{
    return aDisp(bandI) + aCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::aCont(const label bandI) const
{
    return volScalarField::New
    (
        "aCont",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::aDisp(const label bandI) const
{
    return volScalarField::New
    (
        "aDisp",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::as(const label bandI) const
{
    return asCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::asCont(const label bandI) const
{
    return volScalarField::New
    (
        "asCont",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::borderAs(const label bandI) const
{
    return borderAsCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::borderAsCont(const label bandI) const
{
    return volScalarField::New
    (
        "borderAsCont",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::e(const label bandI) const
{
    return eDisp(bandI) + eCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::eCont(const label bandI) const
{
    return volScalarField::New
    (
        "eCont",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::eDisp(const label bandI) const
{
    return volScalarField::New
    (
        "eDisp",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::es(const label bandI) const
{
    return esCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::esCont(const label bandI) const
{
    return volScalarField::New
    (
        "esCont",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::borderEs(const label bandI) const
{
    return borderEsCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::borderEsCont(const label bandI) const
{
    return volScalarField::New
    (
        "borderEsCont",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::E(const label bandI) const
{
    return EDisp(bandI) + ECont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::ECont(const label bandI) const
{
    return volScalarField::New
    (
        "ECont",
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::EDisp(const label bandI) const
{
    return volScalarField::New
    (
        "EDisp",
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
    );
}

Foam::dimensionedScalar
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::borderL(const label bandI) const
{
    return  dimensionedScalar("borderL", dimLength, 0.0);
}

Foam::label Foam::radiationModels::heterogeneousAbsorptionEmissionModel::nBands() const
{
    return pTraits<label>::one;
}


const Foam::Vector2D<Foam::scalar>&
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::bands(const label n) const
{
    return Vector2D<scalar>::one;
}


bool Foam::radiationModels::heterogeneousAbsorptionEmissionModel::isGrey() const
{
    return false;
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::addIntensity
(
    const label rayI,
    const volScalarField& ILambda
) const
{
    return ILambda;
}


void Foam::radiationModels::heterogeneousAbsorptionEmissionModel::correct
(
    volScalarField& a,
    PtrList<volScalarField>& aj
) const
{
    a = this->a();
    aj[0] =  a;
}


// ************************************************************************* //
