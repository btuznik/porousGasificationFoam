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
#include "DynamicList.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(solidHeterogeneousReaction, 0);
    defineRunTimeSelectionTable(solidHeterogeneousReaction, Istream);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidHeterogeneousReaction::solidHeterogeneousReaction
(
    const PtrList<volScalarField>& gasPhaseGases,
    const speciesTable& componets,
    const speciesTable& pyrolisisGases,
    const List<label>& slhs,
    const List<label>& glhs,
    const List<label>& srhs,
    const List<label>& grhs
)
:
    gasPhaseGases_(gasPhaseGases),
    components_(componets),
    pyrolisisGases_(pyrolisisGases),
    slhs_(slhs),
    slhsSto_(),
    glhs_(glhs),
    glhsSto_(),
    srhs_(srhs),
    srhsSto_(),
    grhs_(grhs),
    grhsSto_()
{
}


Foam::solidHeterogeneousReaction::solidHeterogeneousReaction
(
    const solidHeterogeneousReaction& r,
    const PtrList<volScalarField>& gasPhaseGases,
    const speciesTable& componets,
    const speciesTable& pyrolisisGases
)
:
    gasPhaseGases_(gasPhaseGases),
    components_(componets),
    pyrolisisGases_(pyrolisisGases),
    slhs_(r.slhs_),
    slhsSto_(),
    glhs_(r.glhs_),
    glhsSto_(),
    srhs_(r.srhs_),
    srhsSto_(),
    grhs_(r.grhs_),
    grhsSto_()
{
}


Foam::solidHeterogeneousReaction::solidHeterogeneousReaction
(
    const PtrList<volScalarField>& gasPhaseGases,
    const speciesTable& components,
    Istream& is,
    const speciesTable& pyrolisisGases
)
:
    gasPhaseGases_(gasPhaseGases),
    components_(components),
    pyrolisisGases_(pyrolisisGases)
{
    setLRhs(is);
}


Foam::label Foam::solidHeterogeneousReaction::componentIndex
(
    bool& isGas,
    token t
)
{
    if (t.isWord())
    {
        word componentName = t.wordToken();

        size_t i = componentName.find('=');

        if (i != word::npos)
        {
            string exponentStr = componentName
            (
                i + 1,
                componentName.size() - i - 1
            );
            componentName = componentName(0, i);
        }
        if (components_.contains(componentName))
        {
            isGas = false;
            return (components_[componentName]);
        }
        else if (pyrolisisGases_.contains(componentName))
        {
            isGas = true;
            return (pyrolisisGases_[componentName]);
        }
        else
        {
            FatalIOErrorIn
            (
                "solidHeterogeneousReaction::componentIndex(bool&, Istream& is)",
                t.wordToken()
            )
                << "Cannot find component" << componentName
                << "in tables :" << pyrolisisGases_ << " or "
                << components_
                << exit(FatalIOError);
            return -1;
        }

    }
    else
    {
        FatalIOErrorIn("solidHeterogeneousReaction::componentIndex(bool&, Istream& is)", t.wordToken())
            << "Expected a word but found " << t.info()
            << exit(FatalIOError);
        return -1;
    }
}


void Foam::solidHeterogeneousReaction::setLRhs(Istream& is)
{
    label index = 0;

    bool LHS = true;
    bool number = false;
    bool isGas = false;
    bool proceed = true;
    bool isComponent = false;

    scalar stoich = 0.0;

    DynamicList<label> dslhs;
    DynamicList<label> dsrhs;
    DynamicList<label> dglhs;
    DynamicList<label> dgrhs;
    DynamicList<scalar> dslhsSto;
    DynamicList<scalar> dsrhsSto;
    DynamicList<scalar> dglhsSto;
    DynamicList<scalar> dgrhsSto;

    while (is && proceed)
    {
	token t(is);

	if (t.isPunctuation())
	{
            if (t == token::BEGIN_LIST)
            {
                proceed = false;
		is.putBack(t);
            }
	    if (t == token::ASSIGN)
	    {
		LHS = false;
	    }
	}
	else if (t.isNumber())
	{
	    number = true;
	    stoich = t.number();
	}
	else
	{
	    index = componentIndex(isGas,t);
	    isComponent = true;
	}

	if (isComponent)
	{
	    if (LHS)
	    {
		if(isGas)
		{
                    forAll(gasPhaseGases_,i)
                            {
                                if (pyrolisisGases_[index] == gasPhaseGases_[i].name())
                                {
                                        dglhs.append(i);
                                }
                            }
		    if (number)
		    {
			dglhsSto.append(stoich);
		    }
		    else
		    {
			dglhsSto.append(1.0);
		    }
		}
		else
		{
		    if(number)
		    {
			dslhsSto.append(stoich);
		    }
		    else
		    {
			dslhsSto.append(1.0);
		    }
		    dslhs.append(index);
		}
	    }
	    else
	    {
		if(isGas)
		{
		    dgrhs.append(index);
		    if (number)
		    {
			dgrhsSto.append(stoich);
		    }
		    else
		    {
			dgrhsSto.append(1.0);
		    }
		}
		else
		{
		    if(number)
		    {
			dsrhsSto.append(stoich);
		    }
		    else
		    {
			dsrhsSto.append(1.0);
		    }
		    dsrhs.append(index);
		}
	    }
            number = false;
    	    isGas = false;
    	    isComponent = false;
	}
    }

    slhs_ = dslhs.shrink();
    srhs_ = dsrhs.shrink();
    glhs_ = dglhs.shrink();
    grhs_ = dgrhs.shrink();
    slhsSto_ = dslhsSto.shrink();
    srhsSto_ = dsrhsSto.shrink();
    glhsSto_ = dglhsSto.shrink();
    grhsSto_ = dgrhsSto.shrink();

/*
        {
            FatalIOErrorIn("solidHeterogeneousReaction::lsrhs(Istream& is)", is)
                << "Cannot find component in tables"
                << exit(FatalIOError);
        }
    }

    FatalIOErrorIn("solidHeterogeneousReaction::lsrhs(Istream& is)", is)
        << "Cannot continue reading reaction data from stream"
        << exit(FatalIOError);
*/
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidHeterogeneousReaction> Foam::solidHeterogeneousReaction::New
(
    const PtrList<volScalarField>& gasPhaseGases,
    const speciesTable& species,
    Istream& is,
    const speciesTable& pyrolisisGases
)
{
    if (is.eof())
    {
        FatalIOErrorIn
        (
            "solidHeterogeneousReaction::New(const speciesTable& species,"
            " const HashPtrTable& thermoDatabase, Istream&)",
            is
        )   << "solidHeterogeneousReactiontype not specified" << nl << nl
            << "Valid solidHeterogeneousReactiontypes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word reactionTypeName(is);

    IstreamConstructorTable::iterator cstrIter
        = IstreamConstructorTablePtr_->find(reactionTypeName);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "solidHeterogeneousReaction::New(const speciesTable& species,"
            " const HashPtrTable& thermoDatabase, Istream&)",
            is
        )   << "Unknown reaction type "
            << reactionTypeName << nl << nl
            << "Valid reaction types are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<solidHeterogeneousReaction>
    (
        cstrIter()(gasPhaseGases, species, is, pyrolisisGases)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidHeterogeneousReaction::write(Ostream& os) const
{
    os << type() << nl << "    ";

    os << "slhs " << slhs_.size() << " srhs " << srhs_.size() << " glhs " << glhs_.size() << " grhs " << grhs_.size() << nl;

    forAll(slhs_, i)
    {
        os  <<  " + " << slhsSto_[i] << " " << components_[slhs_[i]];
    }

    forAll(glhs_, i)
    {
	os << " + " << glhsSto_[i]  << " " << gasPhaseGases_[glhs_[i]].name();
    }

    os << " = ";

    forAll(srhs_, i)
    {
        os <<  " + " << srhsSto_[i] << " " <<  components_[srhs_[i]];
    }

    forAll(grhs_, i)
    {
        os <<  " + " << grhsSto_[i]  << " " << pyrolisisGases_[grhs_[i]];
    }

    os  << endl << "   ";
}


Foam::scalar Foam::solidHeterogeneousReaction::kf
(
    const scalar T,
    const scalar p,
    const scalarField& c
) const
{
    return 0.0;
}

Foam::scalar Foam::solidHeterogeneousReaction::heatReact() const
{
    return 0.0;
}

Foam::List<scalar> Foam::solidHeterogeneousReaction::nReact() const
{
    return List<scalar>();
}


// ************************************************************************* //
