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
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"

#include "heterogeneousAbsorptionEmissionModel.H"
#include "scatterModel.H"
#include "mathematicalConstants.H"
#include "radiationConstants.H"

using namespace Foam::radiation;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
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

Foam::radiation::heterogeneousMeanTemp::heterogeneousMeanTemp(const volScalarField& T, const volScalarField& porosityF, const List<label>& surfF, const volScalarField& Ts)
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
            "Qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
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
        dimensionedScalar("a", dimless/dimLength, 0.0)
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
        dimensionedScalar("a", dimless/dimLength, 0.0)
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
        dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
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
        dimensionedScalar("radiationSolidSh", dimMass/dimLength/pow3(dimTime), 0.0)
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

Foam::radiation::heterogeneousMeanTemp::~heterogeneousMeanTemp()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::heterogeneousMeanTemp::read()
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


void Foam::radiation::heterogeneousMeanTemp::calculate()
{
    const volScalarField sigmaEff(scatter_->sigmaEff());

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

    scalar radiationEnergy = (radiation::sigmaSB*pow(boundaryMeanTemp/boundarySurface,4)).value();

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


Foam::tmp<Foam::volScalarField> Foam::radiation::heterogeneousMeanTemp::Rp() const
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
                IOobject::NO_WRITE,
                false
            ),
            4.0*absorptionEmission_->eCont()*radiation::sigmaSB*whereIsNot_
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::radiation::heterogeneousMeanTemp::Ru() const
{
    const DimensionedField<scalar, volMesh>& G =
        G_.dimensionedInternalField();
    const DimensionedField<scalar, volMesh> E =
        absorptionEmission_->ECont()().dimensionedInternalField();
    const DimensionedField<scalar, volMesh> a =
        absorptionEmission_->aCont()().dimensionedInternalField();

    return  (a*G - 4.0*E)*whereIsNot_;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::heterogeneousMeanTemp::solidSh() const
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
