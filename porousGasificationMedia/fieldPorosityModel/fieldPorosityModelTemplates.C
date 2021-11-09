/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
#include "fvm.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fieldPorosityModel::addViscousInertialResistance
(
    scalarField& Udiag,
    vectorField& Usource,
    const labelList& cells,
    const scalarField& V,
    const RhoFieldType& rho,
    const scalarField& mu,
    const vectorField& U,
    tensorField& Df 
) const
{
    forAll (cells, i)
    {
        // This is Darcy level only.
        tensor dragCoeff = mu[cells[i]] * Df[cells[i]];

        // This is Darcy-Forcheimer level which can be put into use
        // tensor dragCoeff = mu[cells[i]]*Df[cells[i]] + (rho[cells[i]]*mag(U[cells[i]]))*Ff[cells[i]];

        // Isotropic part that goes into diagonal part of U matrix.
        scalar isoDragCoeff = tr(dragCoeff);
        
        Udiag[cells[i]] += V[cells[i]]*isoDragCoeff;

        // Non-isotropic part that goes into source part of U matrix
        Usource[cells[i]] -=
            V[cells[i]] * ((dragCoeff - I * isoDragCoeff) & U[cells[i]]);
    }
}

// ************************************************************************* //
