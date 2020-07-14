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

inline Foam::solidConstRate::solidConstRate
(
    const scalar A,
    const scalar Ta,
    const scalar Tcrit
)
:
    A_(A),
    Ta_(Ta),
    Tcrit_(Tcrit)
{}


inline Foam::solidConstRate::solidConstRate
(
    const speciesTable&,
    Istream& is
)
:
    A_(readScalar(is.readBegin("solidConst(Istream&)"))),
    Ta_(readScalar(is)),
    Tcrit_(readScalar(is))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::solidConstRate::operator()
(
    const scalar T,
    const scalar,
    const scalarField&
) const
{
    scalar ak = A_;

    return ak;
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const solidConstRate& arr
)
{
    os  << token::BEGIN_LIST
        << arr.A_ << token::SPACE << arr.Ta_ << token::SPACE << arr.Tcrit_
        << token::END_LIST;
    return os;
}


// ************************************************************************* //
