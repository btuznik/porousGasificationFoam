/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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

#include "fieldPorosityModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fieldPorosityModel, 0);
    defineRunTimeSelectionTable(fieldPorosityModel, mesh);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldPorosityModel::fieldPorosityModel
(
    const fvMesh& mesh,
    volScalarField& porosityF
)
:
    regIOobject
    (
        IOobject
        (
            "NO_NAME",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    name_(),
    mesh_(mesh),
    f_(0.),
    porosityF_(porosityF)
{
    IOobject porosityPropertiesHeader
    (
        "porosityProperties",
        mesh_.time().constant(),
        mesh_,
        IOobject::MUST_READ
    );

    word modelType;
    if (porosityPropertiesHeader.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary dict
        (
            IOobject
            (
                "porosityProperties",
                mesh_.time().constant(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );
        f_ = dict.lookupOrDefault("forchheimerCoeff",0.);
        Info << "Forchheimer coefficient f = " << f_  << " specified" << nl
             << "Darcy resistance term will be modified by + f*rho*mag(u)*sqrt(3.)*Df/|Df|"
             << nl << endl;

    }
    else
    {
        Info << "no Forchheimer coefficient specified" << nl
             << "Darcy resistance term only. For the Forchheimer term create porosityProperties dictionary."
             << nl << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldPorosityModel::~fieldPorosityModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fieldPorosityModel::addResistance
(
        fvVectorMatrix& UEqn,
        volTensorField& Df
)
const
{
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();
    const vectorField& U = UEqn.psi();

    word muName = "thermo:mu";
    word nuName = "nu";
    word rhoName = "rho";

    label n=0;
    forAll(porosityF_,celli)
            if(porosityF_[celli] >= 0 && porosityF_[celli] < 1)
                    n++;

    labelList cells(n);
    label i = 0;

    forAll(porosityF_,celli)
    {
            if (porosityF_[celli] >= 0 && porosityF_[celli] < 1)
            {
                cells[i]=celli;
                i++;
            }
    }

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho = mesh_.lookupObject<volScalarField>(rhoName);

        if (mesh_.foundObject<volScalarField>(muName))
        {
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            addViscousInertialResistance
            (
                Udiag,
                Usource,
                cells,
                V,
                rho,
                mu,
                U,
                Df
            );
        }
        else
        {
            const volScalarField& nu = mesh_.lookupObject<volScalarField>(nuName);

            addViscousInertialResistance
            (
                Udiag,
                Usource,
                cells,
                V,
                rho,
                rho * nu,
                U,
                Df
            );
        }
    }
    else
    {
        if (mesh_.foundObject<volScalarField>(nuName))
        {
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);

            addViscousInertialResistance
            (
                Udiag,
                Usource,
                cells,
                V,
                geometricOneField(),
                nu,
                U,
                Df
            );
        }
        else
        {
            const volScalarField& rho =
                mesh_.lookupObject<volScalarField>(rhoName);
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            addViscousInertialResistance
            (
                Udiag,
                Usource,
                cells,
                V,
                geometricOneField(),
                mu / rho,
                U,
                Df
            );
        }
    }
}

bool Foam::fieldPorosityModel::writeData(Ostream& os) const
{
    return true;
}

bool Foam::fieldPorosityModel::read(const dictionary& dict)
{
   return true;
}

// ************************************************************************* //
