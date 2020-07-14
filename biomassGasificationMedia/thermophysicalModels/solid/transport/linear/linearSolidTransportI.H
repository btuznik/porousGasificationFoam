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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class thermo>
inline Foam::linearSolidTransport<thermo>::linearSolidTransport
(
    const thermo& t,
    const scalar K0,
    const scalar n0
)
:
    thermo(t),
    K0_(K0),
    n0_(n0)
{}


template<class thermo>
inline Foam::linearSolidTransport<thermo>::linearSolidTransport
(
    const word& name,
    const linearSolidTransport& ct
)
:
    thermo(name, ct),
    K0_(ct.K0_),
    n0_(ct.n0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class thermo>
inline Foam::scalar Foam::linearSolidTransport<thermo>::K
(
    const scalar T
) const
{
    return K0_ + n0_*T;
}


template<class thermo>
inline Foam::scalar Foam::linearSolidTransport<thermo>::alpha
(
    const scalar T
) const
{
    scalar Cp = this->Cp(T);

    scalar rho = this->rho(T);

    return (K(T)/(rho*Cp));
}


// ************************************************************************* //
