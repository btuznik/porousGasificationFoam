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

#include "heterogeneousMeanTemp.H"
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
    defineTypeNameAndDebug(heterogeneousMeanTemp, 0);
    addToRunTimeSelectionTable
    (
        heterogeneousRadiationModel,
        heterogeneousMeanTemp,
        porosity
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::heterogeneousMeanTemp::heterogeneousMeanTemp
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
    Qr_
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
    borderL_
    (
        dimensionedScalar("borderL", dimLength, 0.0)
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
        if (porosityF_[cellI] > (1.0 - pow(10.0, -8.0)))
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::heterogeneousMeanTemp::~heterogeneousMeanTemp()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiationModels::heterogeneousMeanTemp::read()
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


void Foam::radiationModels::heterogeneousMeanTemp::calculate()
{
    borderAs_ = heterogeneousAbsorptionEmission_->borderAsCont();
    borderL_  = heterogeneousAbsorptionEmission_->borderL();
    const volScalarField sigmaEff(scatter_->sigmaEff());
    surfF_ = surfF_ * 0;

    scalar totalSurf = 0;

//    const dimensionedScalar a0 ("a0", a_.dimensions(), rootVSmall);

    forAll(surfL_,cellI)
    {
        scalar volume = mesh_.V()[surfL_[cellI]];
        scalar surfArea = 0.;
        forAll(mesh_.cells()[surfL_[cellI]],faceI)
        {
            surfArea += mesh_.magSf()[mesh_.cells()[surfL_[cellI]][faceI]];
        }
        surfF_[surfL_[cellI]] = borderL_.value() * surfArea / 6. / volume;
        totalSurf += borderL_.value() * surfArea / 6.;
    }

    scalar totalVol = 0;
    forAll(porosityF_, cellI)
    {
        if (porosityF_[cellI] > (1.0 - pow(10.0, -8.0)))  //this is an ad hoc threshold
        {
            whereIs_[cellI] = 0.0;
            whereIsNot_[cellI] = 1.0;
        }
        else
        {
            whereIs_[cellI] = 1.0;
            whereIsNot_[cellI] = 0.0;
            totalVol += mesh_.V()[cellI];
        }
    }

    reduce(totalSurf, sumOp<scalar>());
    reduce(totalVol, sumOp<scalar>());

    Info << "Radiation active volume to porous media volume ratio: " << totalSurf/max(totalVol,SMALL) << endl;

    dimensionedScalar boundaryMeanTemp("boundaryMeanTemp",dimless,0.0);
    dimensionedScalar boundarySurface("boundarySurface",dimless,0.0);

    // Calculate radiative heat flux on boundaries.
    bool thereIsWall = false;

    forAll(mesh_.boundaryMesh(), patchI)
    {
        if (mesh_.boundaryMesh()[patchI].type() == "wall")
        {
            boundaryMeanTemp += gSum(T_.boundaryField()[patchI]*mesh_.magSf().boundaryField()[patchI]);
            boundarySurface += gSum(mesh_.magSf().boundaryField()[patchI]);
            thereIsWall = true;
        }
    }

    if (not thereIsWall)
    {
        FatalErrorIn("Foam::radiation::heterogeneousMeanTemp")
            << "there is no patch of type: wall "
            << nl << "cannot calculate radiation source."
            << exit(FatalError);
    }
    // Solve G transport equation


    // eqZx2uHGn016
    volScalarField solidRadiation = borderAs_ * surfF_ * physicoChemical::sigma * pow4(Ts_);

    scalar radiationEnergy = (physicoChemical::sigma * pow4(boundaryMeanTemp / boundarySurface)).value();

    forAll(G_,cellI)
    {
        G_[cellI] = radiationEnergy;
        if (surfF_[cellI]==1)
        {
            solidSh_[cellI] = 4.0*(G_[cellI]*borderAs_[cellI] - solidRadiation[cellI])*whereIs_[cellI];
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::radiationModels::heterogeneousMeanTemp::Rp() const
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
                physicoChemical::sigma.dimensions()/dimLength,
                0.0
            )
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiationModels::heterogeneousMeanTemp::Ru() const
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
                "Ru",
                dimMass/dimLength/pow3(dimTime),
                0.0
            )
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousMeanTemp::solidSh() const
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
