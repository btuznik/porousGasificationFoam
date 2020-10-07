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
    Foam::heterogeneousConstRho

Description
    rho constant

SourceFiles
    heterogeneousConstRhoI.H
    heterogeneousConstRho.C

\*---------------------------------------------------------------------------*/

#ifndef heterogeneousConstRho_H
#define heterogeneousConstRho_H

#include "dictionary.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declaration of friend functions and operators
class heterogeneousConstRho;

inline Ostream& operator<<(Ostream&, const heterogeneousConstRho&);

/*---------------------------------------------------------------------------*\
                           Class constRho Declaration
\*---------------------------------------------------------------------------*/

class heterogeneousConstRho
{
    // Private data

        //- Density [kg/m3]
        scalar rho_;


public:

    // Constructors

        //- Construct from components
        inline heterogeneousConstRho(scalar rho);

        //- Construct from dictionary
        //heterogeneousConstRho(Istream&);
        heterogeneousConstRho(const dictionary& dict);

        //- Construct as named copy
        inline heterogeneousConstRho(const heterogeneousConstRho&);

        //- Construct and return a clone
        inline autoPtr<heterogeneousConstRho> clone() const;

        // Selector from dictionary
        inline static autoPtr<heterogeneousConstRho> New(const dictionary& is);


    // Member functions
        //- Return the instantiated type name
        static word typeName()
        {
            return "heterogeneousConstRho";
        }

        //- Return density [kg/m3]
        inline scalar rho(scalar T) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const heterogeneousConstRho&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "heterogeneousConstRhoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
