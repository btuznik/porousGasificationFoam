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
    addToRadiationRunTimeSelectionTables(heterogeneousMeanTemp);
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
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
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
    porosityF_(porosityF),
    surfF_(surfF),
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::heterogeneousMeanTemp::~heterogeneousMeanTemp()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiationModels::heterogeneousMeanTemp::read()
{
    if (heterogeneousMeanTemp::read())
    {
        // nothing to read

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiationModels::heterogeneousMeanTemp::calculate()
{
    const volScalarField sigmaEff(scatter_->sigmaEff());

    const dimensionedScalar a0 ("a0", a_.dimensions(), rootVSmall);

    forAll(porosityF_,cellI)
    {
        if (porosityF_[cellI] > (1.0 - pow(10.0,-3.0)))  // zero porosity threshold can be set oterwise if needed
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

    scalar radiationEnergy = (physicoChemical::sigma*pow(boundaryMeanTemp/boundarySurface,4)).value();

    forAll(G_,cellI)
    {
	G_[cellI] = radiationEnergy;
    }

    volScalarField Gair=G_;

    volScalarField Gsolid=G_*0.0 + whereIsNot_*Gair;

    forAll(surfF_,cellI)
    {
        solidSh_[surfF_[cellI]] = Gair[cellI];
    }
}


Foam::tmp<Foam::volScalarField> Foam::radiationModels::heterogeneousMeanTemp::Rp() const
{
    return volScalarField::New
    (
        "Rp",
        4.0*absorptionEmission_->eCont()*physicoChemical::sigma
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiationModels::heterogeneousMeanTemp::Ru() const
{
    const volScalarField::Internal& G =
        G_();
    const volScalarField::Internal E =
        absorptionEmission_->ECont()()();
    const volScalarField::Internal a =
        absorptionEmission_->aCont()()();

    return a*G - E;
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
