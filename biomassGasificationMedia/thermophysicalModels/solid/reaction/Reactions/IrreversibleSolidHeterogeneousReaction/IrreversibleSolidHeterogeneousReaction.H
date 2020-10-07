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

Class
    Foam::IrreversibleSolidHeterogeneousReaction

Description
    Simple extension of Reaction to handle reversible reactions

SourceFiles
    IrreversibleSolidHeterogeneousReaction.C

\*---------------------------------------------------------------------------*/

#ifndef IrreversibleSolidHeterogeneousReaction_H
#define IrreversibleSolidHeterogeneousReaction_H

#include "solidHeterogeneousReaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class IrreversibleSolidHeterogeneousReaction Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionRate>
class IrreversibleSolidHeterogeneousReaction
:
    public solidHeterogeneousReaction
{
    // Private data

        // Reaction rate
        ReactionRate k_;

        // Reaction heat
        scalar heatReact_;

        // Reaction order
        List<scalar> nReact_;

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=
        (
            const IrreversibleSolidHeterogeneousReaction<ReactionRate>&
        );


public:

    //- Runtime type information
    TypeName("irreversible");


    // Constructors

        //- Construct from components
        IrreversibleSolidHeterogeneousReaction
        (
            const solidHeterogeneousReaction& reaction,
            const ReactionRate& reactionRate,
            const scalar nReact
        );


        //- Construct from Istream
        IrreversibleSolidHeterogeneousReaction
        (
	        const PtrList<volScalarField>& gasPhaseGases,
            const speciesTable& components,
            Istream& is,
            const speciesTable& pyrolysisGases
        );


    //- Destructor
    virtual ~IrreversibleSolidHeterogeneousReaction()
    {}


    // Member Functions

        // IrreversibleSolidHeterogeneousReaction rate coefficients

            //- Forward rate constant
            virtual scalar kf
            (
                const scalar T,
                const scalar p,
                const scalarField& c
            ) const;


           //- Reaction order
            virtual scalar heatReact() const;

           //- Reaction order
           virtual List<scalar> nReact() const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "IrreversibleSolidHeterogeneousReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
