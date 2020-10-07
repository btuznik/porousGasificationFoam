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
    Foam::solidConstRate

Description
    Const reaction rate for solids

SourceFiles
    solidConstRateI.H

\*---------------------------------------------------------------------------*/

#ifndef solidConstRate_H
#define solidConstRate_H

#include "scalarField.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declaration of friend functions and operators
class solidConstRate;
Ostream& operator<<(Ostream&, const solidConstRate&);

/*---------------------------------------------------------------------------*\
                       Class solidConstRate Declaration
\*---------------------------------------------------------------------------*/

class solidConstRate
{
    // Private data

        scalar A_;
        scalar Ta_;
        scalar Tcrit_;


public:

    // Constructors

        //- Construct from components
        inline solidConstRate
        (
            const scalar A,
            const scalar Ta,
            const scalar Tcrit
        );

        //- Construct from Istream
        inline solidConstRate
        (
            const speciesTable& species,
            Istream& is
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "SolidConst";
        }

        inline scalar operator()
        (
            const scalar T,
            const scalar p,
            const scalarField& c
        ) const;


    // Ostream Operator

        inline friend Ostream& operator<<
        (
            Ostream&,
            const solidConstRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidConstRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
