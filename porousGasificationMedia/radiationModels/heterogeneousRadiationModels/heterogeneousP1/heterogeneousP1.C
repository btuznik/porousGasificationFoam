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

#include "heterogeneousP1.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "heterogeneousAbsorptionEmissionModel.H"
#include "scatterModel.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
    defineTypeNameAndDebug(heterogeneousP1, 0);
    addToRunTimeSelectionTable
    (
        heterogeneousRadiationModel,
        heterogeneousP1,
        porosity
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::heterogeneousP1::heterogeneousP1
(
    const volScalarField& T,
    const volScalarField& porosityF,
    const List<label>& surfF,
    const volScalarField& Ts
)
:
    heterogeneousRadiationModel(typeName, T),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), 0)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    ),
    as_
    (
        IOobject
        (
            "as",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    ),
    borderAs_
    (
        IOobject
        (
            "borderAs",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    ),
    es_
    (
        IOobject
        (
            "es",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    ),
    borderEs_
    (
        IOobject
        (
            "borderEs",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
    ),
    surfToVol_
    (
        IOobject
        (
            "surfToVol",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, 1)
    ),
    porosityF_(porosityF),
    surfL_(surfF),
    surfF_
    (
        IOobject
        (
            "surfF",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        scalar(0.0)
    ),
    whereIs_
    (
        IOobject
        (
            "whereIs",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        scalar(1.0)
    ),
    whereIsNot_
    (
        IOobject
        (
            "whereIsNot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        scalar(0.0)
    ),
    solidSh_
    (
        IOobject
        (
            "radiationSolidSh",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
    )
{
    forAll(porosityF_,cellI)
    {
        if (porosityF_[cellI] > (1.0 - pow(10.0,-8.0)))
        {
	    whereIs_[cellI] = 0.0;
	    whereIsNot_[cellI] = 1.0;
        }
        else
        {
	    whereIs_[cellI] = 1.0;
	    whereIsNot_[cellI] = 0.0;
        }
    }
    surfToVol_.ref() = 1./(pow(mesh_.V(),1./3.)*6.);

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::heterogeneousP1::~heterogeneousP1()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiationModels::heterogeneousP1::read()
{
    if (heterogeneousRadiationModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiationModels::heterogeneousP1::calculate()
{
    a_ = heterogeneousAbsorptionEmission_->aCont();
    as_ = heterogeneousAbsorptionEmission_->asCont();
    borderAs_ = heterogeneousAbsorptionEmission_->borderAsCont();
    e_ = heterogeneousAbsorptionEmission_->eCont();
    es_ = heterogeneousAbsorptionEmission_->esCont();
    borderEs_ = heterogeneousAbsorptionEmission_->borderEsCont();
    E_ = heterogeneousAbsorptionEmission_->ECont();
    const volScalarField sigmaEff(scatter_->sigmaEff());

    const dimensionedScalar a0 ("a0", a_.dimensions(), rootVSmall);

    surfF_ = surfF_*0;
    forAll(surfL_,cellI)
    {
        surfF_[surfL_[cellI]] = 1.0;
    }

    forAll(porosityF_,cellI)
    {
        if (porosityF_[cellI] > (1.0 - pow(10.0,-8.0)))  //this is an ad hoc threshold
        {
	           whereIs_[cellI] = 0.0;
	           whereIsNot_[cellI] = 1.0;
        }
        else
        {
	           whereIs_[cellI] = 1.0;
	           whereIsNot_[cellI] = 0.0;
        }
    }

    // Construct diffusion
    const volScalarField gamma
    (
        IOobject
        (
            "gammaRad",
            G_.mesh().time().timeName(),
            G_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        1.0/(3.0*(a_ + as_*(whereIs_-surfF_) + borderAs_*surfF_) + sigmaEff + a0)
    );

    volScalarField solidRadiation = (es_*(whereIs_-surfF_) + borderEs_*surfF_ )*physicoChemical::sigma*pow4(Ts_);

    // Solve G transport equation
    solve
    (
        fvm::laplacian(gamma, G_)
      - fvm::Sp((a_ + as_*(whereIs_-surfF_) + borderAs_*surfF_ ), G_)
     ==
      - 4.0*(e_*physicoChemical::sigma*pow4(T_) + solidRadiation) - E_
    );

    volScalarField Gair=G_;

    forAll(G_,cellI)
    {
        if (surfF_[cellI]==0)
        {
            solidSh_[cellI] = (G_[cellI]*as_[cellI] - 4.0*solidRadiation[cellI])*whereIs_[cellI];
        }
        else
        {
            solidSh_[cellI] = (G_[cellI]*borderAs_[cellI] - 4.0*solidRadiation[cellI])*whereIs_[cellI];
        }
    }

    volScalarField::Boundary& qrBf = qr_.boundaryFieldRef();

    // Calculate radiative heat flux on boundaries.
    forAll(mesh_.boundaryMesh(), patchi)
    {
        if (!G_.boundaryField()[patchi].coupled())
        {
            qrBf[patchi] =
                -gamma.boundaryField()[patchi]
                *G_.boundaryField()[patchi].snGrad();
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::radiationModels::heterogeneousP1::Rp() const
{
    return volScalarField::New
    (
        "Rp",
        4.0*heterogeneousAbsorptionEmission().eCont()*physicoChemical::sigma
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiationModels::heterogeneousP1::Ru() const
{
    const volScalarField::Internal& G = G_();
    const volScalarField::Internal E = heterogeneousAbsorptionEmission_->ECont()()();
    const volScalarField::Internal a = heterogeneousAbsorptionEmission_->aCont()()();

    return a*G - E;
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousP1::solidSh() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "tSolidSh",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            solidSh_
        )
    );
}


// ************************************************************************* //