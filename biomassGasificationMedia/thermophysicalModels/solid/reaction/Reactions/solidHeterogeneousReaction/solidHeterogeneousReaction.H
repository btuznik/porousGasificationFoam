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
    Foam::solidHeterogeneousReaction

Description

    Read solid reactions of the type S1 + G1 = S2 + G1

SourceFiles
    solidHeterogeneousReactionI.H
    solidHeterogeneousReaction.C

\*---------------------------------------------------------------------------*/

#ifndef solidHeterogeneousReaction_H
#define solidHeterogeneousReaction_H

#include "speciesTable.H"
#include "scalarField.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

#include "volFieldsFwd.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class solidHeterogeneousReaction;

inline Ostream& operator<<(Ostream&, const solidHeterogeneousReaction&);

/*---------------------------------------------------------------------------*\
                        Class solidHeterogeneousReactionDeclaration
\*---------------------------------------------------------------------------*/

class solidHeterogeneousReaction
{

private:

    // Private data

	//- List of gas species from gas phase available for gasification reactions
	const PtrList<volScalarField>& gasPhaseGases_;

        //- List of solid names present in reaction system
        const speciesTable& components_;

        //- List of gas species present in reaction system
        speciesTable pyrolisisGases_;

        //- Solid components index for the left-hand-side of the reaction
        List<label> slhs_;

	//- Solid mass coefficients for left-hand-side of the reaction
	List<scalar> slhsSto_;
	
	//- Gas components index for the left-hand-side of the reaction
	List<label> glhs_;

	//- Gas stoichmietry coefficients for left-hand-side of the reaction
	List<scalar> glhsSto_;
	
        //- Solid components index for the right-hand-side of the reaction
        List<label> srhs_;

	//- Solid mass coefficients for left-hand-side of the reaction
	List<scalar> srhsSto_;
	
        //- Specie index for the right-hand-side of the reaction
        List<label> grhs_;

	//- Gas stoichmietry coefficients for right-hand-side of the reaction
	List<scalar> grhsSto_;
	
    // Private Member Functions

        //- Set rhs and lhs of the reaction
        void setLRhs(Istream&);

        //- Look for the component index in the reaction
        label componentIndex(bool& isGas, token t);

        //- Disallow default bitwise copy construct
        solidHeterogeneousReaction(const solidHeterogeneousReaction&);

        //- Disallow default bitwise assignment
        void operator=(const solidHeterogeneousReaction&);


public:

    //- Runtime type information
    TypeName("HeterogeneousReaction");

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            solidHeterogeneousReaction,
            Istream,
            (
		const PtrList<volScalarField>& gasPhaseGases,
                const speciesTable& components,
                Istream& is,
                const speciesTable& pyrolysisGases
            ),
            (gasPhaseGases, components, is, pyrolysisGases)
        );


    // Public classes

        //- Class used for the read-construction of PtrLists of reaction
        class iNew
        {
	    const PtrList<volScalarField>& gasPhaseGases_;
            const speciesTable& components_;
            speciesTable pyrolisisGases_;

        public:

            iNew
            (
                const PtrList<volScalarField>& gasPhaseGases,
		const speciesTable& components,
                Istream& pyrolisisGases
            )
            :
		gasPhaseGases_(gasPhaseGases),
                components_(components),
                pyrolisisGases_(pyrolisisGases)
            {}

            autoPtr<solidHeterogeneousReaction> operator()(Istream& is) const
            {
                return autoPtr<solidHeterogeneousReaction>
                (
                    solidHeterogeneousReaction::New
                    (
			gasPhaseGases_,
                        components_,
                        is,
                        pyrolisisGases_
                    )
                );
            }
        };


    // Constructors

        //- Construct from components
        solidHeterogeneousReaction
        (
	    const PtrList<volScalarField>& gasPhaseGases,
            const speciesTable& components,
            const speciesTable& pyrolisisGases,
            const List<label>& slhs,
	    const List<label>& glhs,
            const List<label>& srhs,
            const List<label>& grhs
        );

        //- Construct as copy given new speciesTable
        solidHeterogeneousReaction
        (
            const solidHeterogeneousReaction&,
	    const PtrList<volScalarField>& gasPhaseGases,
            const speciesTable& components,
            const speciesTable& pyrolisisGases
        );

        //- Construct from Istream
        solidHeterogeneousReaction
        (
	    const PtrList<volScalarField>& gasPhaseGases,
            const speciesTable& components,
            Istream& is,
            const speciesTable& pyrolisisGases
        );

        //- Construct and return a clone
        virtual autoPtr<solidHeterogeneousReaction> clone() const
        {
            return autoPtr<solidHeterogeneousReaction>
            (
                new solidHeterogeneousReaction(*this)
            );
        }


    // Selectors

        //- Return a pointer to a new patchField created on freestore from input
        static autoPtr<solidHeterogeneousReaction> New
        (
	    const PtrList<volScalarField>& gasPhaseGases,
            const speciesTable& components,
            Istream&,
            const speciesTable& pyrolisisGases
        );


    //- Destructor
    virtual ~solidHeterogeneousReaction()
    {}


    // Member Functions

        // Access

            inline const List<label>& slhs() const;
	    inline const List<scalar>& slhsSto() const;
	    inline const List<label>& glhs() const;
	    inline const List<scalar>& glhsSto() const;
            inline const List<label>& srhs() const;
	    inline const List<scalar>& srhsSto() const;
            inline const List<label>& grhs() const;
	    inline const List<scalar>& grhsSto() const;

            inline const speciesTable& pyrolisisGases() const;


        // solidHeterogeneousReactionrate coefficients

            virtual scalar kf
            (
                const scalar T,
                const scalar p,
                const scalarField& c
            ) const;

	    virtual scalar heatReact() const;

            virtual List<scalar> nReact() const;

        //- Write
        virtual void write(Ostream&) const;


    // Ostream Operator

        friend Ostream& operator<<
        (
            Ostream&,
            const solidHeterogeneousReaction&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidHeterogeneousReactionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
