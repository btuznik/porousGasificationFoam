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

#include "heterogeneousConstRho.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::heterogeneousConstRho::heterogeneousConstRho
(
    scalar rho
)
:
    rho_(rho)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::heterogeneousConstRho::heterogeneousConstRho(const heterogeneousConstRho& pg)
:
    rho_(pg.rho_)
{}


inline Foam::autoPtr<Foam::heterogeneousConstRho> Foam::heterogeneousConstRho::clone() const
{
    return autoPtr<heterogeneousConstRho>(new heterogeneousConstRho(*this));
}


inline Foam::autoPtr<Foam::heterogeneousConstRho> Foam::heterogeneousConstRho::New(const dictionary& is)
{
    return autoPtr<heterogeneousConstRho>(new heterogeneousConstRho(is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::heterogeneousConstRho::rho(scalar) const
{
    return rho_;
}


// ************************************************************************* //
