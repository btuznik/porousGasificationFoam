/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "volFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SolidThermo, class SolidThermoType, class GasThermoType>
inline Foam::PtrList<Foam::scalarField>&
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::RRs()
{
    return RRs_;
}


template<class SolidThermo, class SolidThermoType, class GasThermoType>
inline Foam::PtrList<Foam::scalarField>&
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::RRg()
{
    return RRg_;
}

template<class SolidThermo, class SolidThermoType, class GasThermoType>
inline Foam::tmp<Foam::volScalarField>
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::RRg() const
{
    tmp<volScalarField> tRRg
    (
        new volScalarField
        (
            IOobject
            (
                "RRg",
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        for (label i=0; i < nGases_; i++)
        {
            tRRg.ref().primitiveFieldRef() += RRg_[i];
        }
        tRRg.ref().correctBoundaryConditions();
    }
    return tRRg;
}

template<class SolidThermo, class SolidThermoType, class GasThermoType>
inline Foam::tmp<Foam::volScalarField>
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::RRs() const
{
    tmp<volScalarField> tRRs
    (
        new volScalarField
        (
            IOobject
            (
                "RRs",
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        for (label i=0; i < nSolids_; i++)
        {
            tRRs.ref().primitiveFieldRef() += RRs_[i];
        }
        tRRs.ref().correctBoundaryConditions();
    }
    return tRRs;
}

template<class SolidThermo, class SolidThermoType, class GasThermoType>
inline Foam::tmp<Foam::volScalarField>
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::RRs
(
    const label i
) const
{
    tmp<volScalarField> tRRs
    (
        new volScalarField
        (
            IOobject
            (
                "RRs(" + Ys_[i].name() + ')',
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        tRRs.ref().primitiveFieldRef() = RRs_[i];
        tRRs.ref().correctBoundaryConditions();
    }
    return tRRs;
}

template<class SolidThermo, class SolidThermoType, class GasThermoType>
inline Foam::tmp<Foam::volScalarField>
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::RRg
(
    const label i
) const
{
    tmp<volScalarField> tRRg
    (
        new volScalarField
        (
            IOobject
            (
                "RRg(" + this->pyrolisisGases_[i] + ')',
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        tRRg.ref().primitiveFieldRef() = RRg_[i];
        tRRg.ref().correctBoundaryConditions();
    }
    return tRRg;
}

template<class SolidThermo, class SolidThermoType, class GasThermoType>
inline Foam::label
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::nEqns() const
{
    // nEqns = number of species + temperature + pressure
    return nSpecie_ + 2;
}


template<class SolidThermo, class SolidThermoType, class GasThermoType>
inline Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::RR()
{
    return RR_;
}


template<class SolidThermo, class SolidThermoType, class GasThermoType>
inline const Foam::PtrList<Foam::solidHeterogeneousReaction>&
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::reactions() const
{
    return reactions_;
}

template<class SolidThermo, class SolidThermoType, class GasThermoType>
inline const Foam::PtrList<GasThermoType>&
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::gasThermo() const
{
    return gasThermo_;
}

template<class SolidThermo, class SolidThermoType, class GasThermoType>
inline const Foam::speciesTable&
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::gasTable() const
{
    return pyrolisisGases_;
}


template<class SolidThermo, class SolidThermoType, class GasThermoType>
inline const Foam::PtrList<SolidThermoType>&
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::specieThermos() const
{
    return solidThermo_;
}


template<class SolidThermo, class SolidThermoType, class GasThermoType>
inline Foam::label
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::nSpecie() const
{
    return nSpecie_;
}

template<class SolidThermo, class SolidThermoType, class GasThermoType>
inline Foam::label
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::nGases() const
{
    return nGases_;
}

template<class SolidThermo, class SolidThermoType, class GasThermoType>
inline Foam::label
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::nReaction() const
{
    return nReaction_;
}


template<class SolidThermo, class SolidThermoType, class GasThermoType>
inline Foam::scalar
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::Treact() const
{
    return Treact_;
}


template<class SolidThermo, class SolidThermoType, class GasThermoType>
inline Foam::scalar&
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::Treact()
{
    return Treact_;
}


template<class SolidThermo, class SolidThermoType, class GasThermoType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::RR
(
    const label i
) const
{
    return RR_[i];
}

template<class SolidThermo, class SolidThermoType, class GasThermoType>
Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::RR
(
    const label i
)
{
    return RR_[i];
}

template<class SolidThermo, class SolidThermoType, class GasThermoType>
inline Foam::scalarField&
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::coeffs()
{
    return coeffs_;
}


template<class SolidThermo, class SolidThermoType, class GasThermoType>
inline const Foam::scalarField&
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::coeffs() const
{
    return coeffs_;
}
// ************************************************************************* //
