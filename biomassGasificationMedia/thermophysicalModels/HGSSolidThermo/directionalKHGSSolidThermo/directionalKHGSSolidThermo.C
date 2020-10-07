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

#include "directionalKHGSSolidThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "transform.H"
#include "transformField.H"
#include "biomassInterpolateXY.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(directionalKHGSSolidThermo, 0);
    addToRunTimeSelectionTable
    (
        HGSSolidThermo,
        directionalKHGSSolidThermo,
        mesh
    );

    addToRunTimeSelectionTable
    (
        HGSSolidThermo,
        directionalKHGSSolidThermo,
        dictionary
    );

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::directionalKHGSSolidThermo::directionalKHGSSolidThermo
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    interpolatedHGSSolidThermo(mesh, typeName + "Coeffs", dict),
    directionalK_
    (
        IOobject
        (
            "K",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimEnergy/dimTime/(dimLength*dimTemperature)
    ),
    ccTransforms_
    (
        IOobject
        (
            "ccTransforms",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimLength
    )
{
    init();
}


Foam::directionalKHGSSolidThermo::directionalKHGSSolidThermo(const fvMesh& mesh)
:
    interpolatedHGSSolidThermo(mesh, typeName + "Coeffs"),
    directionalK_
    (
        IOobject
        (
            "K",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimEnergy/dimTime/(dimLength*dimTemperature)
    ),
    ccTransforms_
    (
        IOobject
        (
            "ccTransforms",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimLength
    )
{
    init();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::directionalKHGSSolidThermo::~directionalKHGSSolidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::directionalKHGSSolidThermo::init()
{
    KValues_ = Field<vector>(subDict(typeName + "Coeffs").lookup("KValues"));

    const fvMesh& mesh = directionalK_.mesh();

    // Determine transforms for cell centres
    forAll(mesh.C(), cellI)
    {
        vector dir = mesh.C()[cellI] - coordSys_.origin();
        dir /= mag(dir);

        // Define local coordinate system with
        // - e1 : axis from cc to centre
        // - e3 : rotation axis
        coordinateSystem cs
        (
            "cc",
            coordSys_.origin(),
            coordSys_.R().e3(),     //z',e3
            dir                 //x',e1
        );
        
        ccTransforms_[cellI] = cs.R().R();
    }
        
    forAll(mesh.C().boundaryField(), patchI)
    {
        const fvPatchVectorField& patchC = mesh.C().boundaryField()[patchI];
        fvPatchTensorField& patchT = ccTransforms_.boundaryFieldRef()[patchI];

        tensorField tc(patchT.size());
        forAll(tc, i)
        {
            vector dir = patchC[i] - coordSys_.origin();
            dir /= mag(dir);

            coordinateSystem cs
            (
                "cc",
                coordSys_.origin(),
                coordSys_.R().e3(),     //z',e3
                dir                 //x',e1
            );

            tc[i] = cs.R().R();
        }
        patchT = tc;
    }

    if (debug)
    {
        Info<< "directionalKHGSSolidThermo : dumping converted Kxx, Kyy, Kzz"
            << endl;
        {
            volVectorField Kxx
            (
                IOobject
                (
                    "Kxx",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                mesh,
                dimless
            );
            Kxx.primitiveFieldRef() = transform
            (
                ccTransforms_.internalField(),
                vectorField
                (
                    ccTransforms_.internalField().size(),
                    point(1, 0, 0)
                )
            );
            forAll(Kxx.boundaryField(), patchI)
            {
                Kxx.boundaryFieldRef()[patchI] = transform
                (
                    ccTransforms_.boundaryField()[patchI],
                    vectorField
                    (
                        ccTransforms_.boundaryField()[patchI].size(),
                        point(1, 0, 0)
                    )
                );
            }
            Kxx.write();
        }
        {
            volVectorField Kyy
            (
                IOobject
                (
                    "Kyy",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                mesh,
                dimless
            );
            Kyy.primitiveFieldRef() = transform
            (
                ccTransforms_.internalField(),
                vectorField
                (
                    ccTransforms_.internalField().size(),
                    point(0, 1, 0)
                )
            );
            forAll(Kyy.boundaryField(), patchI)
            {
                Kyy.boundaryFieldRef()[patchI] = transform
                (
                    ccTransforms_.boundaryField()[patchI],
                    vectorField
                    (
                        ccTransforms_.boundaryField()[patchI].size(),
                        point(0, 1, 0)
                    )
                );
            }
            Kyy.write();
        }
        {
            volVectorField Kzz
            (
                IOobject
                (
                    "Kzz",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                mesh,
                dimless
            );
            Kzz.primitiveFieldRef() = transform
            (
                ccTransforms_.internalField(),
                vectorField
                (
                    ccTransforms_.internalField().size(),
                    point(0, 0, 1)
                )
            );
            forAll(Kzz.boundaryField(), patchI)
            {
                Kzz.boundaryFieldRef()[patchI] = transform
                (
                    ccTransforms_.boundaryField()[patchI],
                    vectorField
                    (
                        ccTransforms_.boundaryField()[patchI].size(),
                        point(0, 0, 1)
                    )
                );
            }
            Kzz.write();
        }
    }

    correct();
}


Foam::symmTensor Foam::directionalKHGSSolidThermo::transformPrincipal
(
    const tensor& tt,
    const vector& st
) const
{
    return symmTensor
    (
        tt.xx()*st.x()*tt.xx()
      + tt.xy()*st.y()*tt.xy()
      + tt.xz()*st.z()*tt.xz(),

        tt.xx()*st.x()*tt.yx()
      + tt.xy()*st.y()*tt.yy()
      + tt.xz()*st.z()*tt.yz(),

        tt.xx()*st.x()*tt.zx()
      + tt.xy()*st.y()*tt.zy()
      + tt.xz()*st.z()*tt.zz(),

        tt.yx()*st.x()*tt.yx()
      + tt.yy()*st.y()*tt.yy()
      + tt.yz()*st.z()*tt.yz(),

        tt.yx()*st.x()*tt.zx()
      + tt.yy()*st.y()*tt.zy()
      + tt.yz()*st.z()*tt.zz(),

        tt.zx()*st.x()*tt.zx()
      + tt.zy()*st.y()*tt.zy()
      + tt.zz()*st.z()*tt.zz()
    );
}


void Foam::directionalKHGSSolidThermo::transformField
(
    symmTensorField& fld,
    const tensorField& tt,
    const vectorField& st
) const
{
    fld.setSize(tt.size());
    forAll(fld, i)
    {
        fld[i] = transformPrincipal(tt[i], st[i]);
    }
}


void Foam::directionalKHGSSolidThermo::correct()
{
    calculate();
    interpolatedHGSSolidThermo::calculate();
}


const Foam::volSymmTensorField&
Foam::directionalKHGSSolidThermo::directionalK() const
{
    return directionalK_;
}


void Foam::directionalKHGSSolidThermo::calculate()
{
    // Correct directionalK
    Field<vector> localK
    (
        biomassInterpolateXY
        (
            T_.internalField(),
            TValues_,
            KValues_
        )
    );

    // Transform into global coordinate system
    transformField
    (
        directionalK_.primitiveFieldRef(),
        ccTransforms_.internalField(),
        localK
    );

    forAll(directionalK_.boundaryField(), patchI)
    {
        directionalK_.boundaryFieldRef()[patchI] == this->directionalK(patchI)();
    }
}


const Foam::volScalarField& Foam::directionalKHGSSolidThermo::K() const
{
    forAll(KValues_, i)
    {
        const vector& v = KValues_[i];
        if
        (
            v.x() != v.y()
         || v.x() != v.z()
         || v.y() != v.z()
        )
        {
            FatalErrorIn("directionalKHGSSolidThermo::K() const")
                << "Supplied K values " << KValues_
                << " are not isotropic." << exit(FatalError);
        }
    }

    // Get temperature interpolated properties (principal directions)
    Field<vector> localK
    (
        biomassInterpolateXY
        (
            T_.internalField(),
            TValues_,
            KValues_
        )
    );

    tmp<volScalarField> tK
    (
        new volScalarField
        (
            IOobject
            (
                "K",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimEnergy/dimTime/(dimLength*dimTemperature)
        )
    );
    volScalarField& K = tK.ref();

    K.primitiveFieldRef() = biomassInterpolateXY
    (
        T_.internalField(),
        TValues_,
        KValues_.component(0)()
    );

    forAll(K.boundaryField(), patchI)
    {
        K.boundaryFieldRef()[patchI] == this->K(patchI)();
    }

    return tK;
}


Foam::tmp<Foam::scalarField> Foam::directionalKHGSSolidThermo::K
(
    const label patchI
) const
{
    forAll(KValues_, i)
    {
        const vector& v = KValues_[i];
        if
        (
            v.x() != v.y()
         || v.x() != v.z()
         || v.y() != v.z()
        )
        {
            FatalErrorIn("directionalKHGSSolidThermo::K() const")
                << "Supplied K values " << KValues_
                << " are not isotropic." << exit(FatalError);
        }
    }

    return tmp<scalarField>
    (
        new scalarField
        (
            biomassInterpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                KValues_.component(0)()
            )
        )
    );
}


Foam::tmp<Foam::symmTensorField> Foam::directionalKHGSSolidThermo::directionalK
(
    const label patchI
) const
{
    const fvPatchScalarField& patchT = T_.boundaryField()[patchI];

    Field<vector> localK(biomassInterpolateXY(patchT, TValues_, KValues_));

    tmp<symmTensorField> tglobalK(new symmTensorField(localK.size()));
    transformField(tglobalK.ref(), ccTransforms_.boundaryField()[patchI], localK);

    return tglobalK;
}


bool Foam::directionalKHGSSolidThermo::read()
{
    return read(subDict(typeName + "Coeffs"));
}


bool Foam::directionalKHGSSolidThermo::read(const dictionary& dict)
{
    coordSys_ = coordinateSystem(mesh_, dict);
    KValues_  = Field<vector>(subDict(typeName + "Coeffs").lookup("KValues"));
    return true;
}


bool Foam::directionalKHGSSolidThermo::writeData(Ostream& os) const
{
    bool ok = interpolatedHGSSolidThermo::writeData(os);
    os.writeKeyword("KValues") << KValues_ << token::END_STATEMENT << nl;
    return ok && os.good();
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const directionalKHGSSolidThermo& s)
{
    s.writeData(os);
    return os;
}


// ************************************************************************* //
