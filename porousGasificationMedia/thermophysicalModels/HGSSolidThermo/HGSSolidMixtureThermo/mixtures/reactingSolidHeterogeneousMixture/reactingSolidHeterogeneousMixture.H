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
    Foam::reactingSolidHeterogeneousMixture

Description
    Foam::reactingSolidHeterogeneousMixture

SourceFiles
    reactingSolidHeterogeneousMixture.C

\*---------------------------------------------------------------------------*/

#ifndef reactingSolidHeterogeneousMixture_H
#define reactingSolidHeterogeneousMixture_H

#include "multiComponentSolidMixture.H"
#include "solidHeterogeneousReaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class reactingSolidHeterogeneousMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoSolidType>
class reactingSolidHeterogeneousMixture
:
    public multiComponentSolidMixture<ThermoSolidType>,
    public PtrList<solidHeterogeneousReaction>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        reactingSolidHeterogeneousMixture(const reactingSolidHeterogeneousMixture&);

        //- Disallow default bitwise assignment
        void operator=(const reactingSolidHeterogeneousMixture&);


public:

    //- The type of thermo package this mixture is instantiated for
    typedef ThermoSolidType thermoType;


    // Constructors

        //- Construct from dictionary and mesh
        reactingSolidHeterogeneousMixture(const dictionary&, const fvMesh&, const PtrList<volScalarField>& gasPhaseGases);


    //- Destructor
    virtual ~reactingSolidHeterogeneousMixture()
    {}


    // Member functions

        //- Read dictionary
        void read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "reactingSolidHeterogeneousMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
