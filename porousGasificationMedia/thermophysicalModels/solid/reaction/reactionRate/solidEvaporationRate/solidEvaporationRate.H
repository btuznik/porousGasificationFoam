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
    Foam::solidEvaporationRate

Description
    Evaporation reaction rate for solids

SourceFiles
    solidEvaporationRateI.H

\*---------------------------------------------------------------------------*/

#ifndef solidEvaporationRate_H
#define solidEvaporationRate_H

#include "scalarField.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declaration of friend functions and operators
class solidEvaporationRate;
Ostream& operator<<(Ostream&, const solidEvaporationRate&);

/*---------------------------------------------------------------------------*\
                       Class solidEvaporationRate Declaration
\*---------------------------------------------------------------------------*/

class solidEvaporationRate
{
    // Private data

        scalar A_;
        scalar Ta_;
        scalar Tcrit_;


public:

    // Constructors

        //- Construct from components
        inline solidEvaporationRate
        (
            const scalar A,
            const scalar Ta,
            const scalar Tcrit
        );

        //- Construct from Istream
        inline solidEvaporationRate
        (
            const speciesTable& species,
            Istream& is
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "SolidEvaporation";
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
            const solidEvaporationRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidEvaporationRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
