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

#include "reactingSolidHeterogeneousMixture.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoSolidType>
Foam::reactingSolidHeterogeneousMixture<ThermoSolidType>::reactingSolidHeterogeneousMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const PtrList<volScalarField>& gasPhaseGases
)
:
    multiComponentSolidMixture<ThermoSolidType>
    (
        thermoDict,
        mesh
    ),
    PtrList<solidHeterogeneousReaction>
    (
        mesh.lookupObject<dictionary> ("chemistryProperties").lookup("solidReactions"),
        solidHeterogeneousReaction::iNew
        (
            gasPhaseGases,
            this->components_,
            mesh.lookupObject<dictionary> ("chemistryProperties").lookup("species")
        )
    )
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoSolidType>
void Foam::reactingSolidHeterogeneousMixture<ThermoSolidType>::read
(
  const dictionary& thermoDict
)
{}


// ************************************************************************* //
