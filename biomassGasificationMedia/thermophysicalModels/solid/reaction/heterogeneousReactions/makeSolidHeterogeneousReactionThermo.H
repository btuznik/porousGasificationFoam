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

InClass
    Foam::makeSolidHeterogeneousReactionThermo

Description
    Macros for instantiating reactions on given solid thermo packages

\*---------------------------------------------------------------------------*/

#ifndef makeHeterogeneousReactionThermo_H
#define makeHeterogeneousReactionThermo_H

#include "solidHeterogeneousReaction.H"
#include "IrreversibleSolidHeterogeneousReaction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeReaction(ReactionType, ReactionRate)                              \
                                                                              \
    typedef solidHeterogeneousReaction heterogeneousReaction;                     \
                                                                              \
    typedef ReactionType<ReactionRate>                                        \
        ReactionType##ReactionRate;                                           \
                                                                              \
    template<>                                                                \
    const word ReactionType##ReactionRate::typeName                           \
    (                                                                         \
        ReactionType::typeName_()                                             \
      + ReactionRate::type()                                                  \
      + heterogeneousReaction::typeName_()                                      \
    );                                                                        \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        heterogeneousReaction,                                                             \
        ReactionType##ReactionRate,                                           \
        Istream                                                               \
    );


#define makeIRReactions(ReactionRate)                                         \
                                                                              \
    makeReaction(IrreversibleSolidHeterogeneousReaction, ReactionRate)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif


// ************************************************************************* //
