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

#include "solidHeterogeneousReaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const List<label>& solidHeterogeneousReaction::slhs() const
{
    return slhs_;
}

inline const List<scalar>& solidHeterogeneousReaction::slhsSto() const
{
    return slhsSto_;
}

inline const List<label>& solidHeterogeneousReaction::glhs() const
{
    return glhs_;
}

inline const List<scalar>& solidHeterogeneousReaction::glhsSto() const
{
    return glhsSto_;
}

inline const List<label>& solidHeterogeneousReaction::srhs() const
{
    return srhs_;
}

inline const List<scalar>& solidHeterogeneousReaction::srhsSto() const
{
    return srhsSto_;
}

inline const List<label>& solidHeterogeneousReaction::grhs() const
{
    return grhs_;
}

inline const List<scalar>& solidHeterogeneousReaction::grhsSto() const
{
    return grhsSto_;
}

inline const speciesTable& solidHeterogeneousReaction::pyrolisisGases() const
{
    return pyrolisisGases_;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Ostream& operator<<
(
    Ostream& os,
    const solidHeterogeneousReaction& r
)
{
    r.write(os);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
