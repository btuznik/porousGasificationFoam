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

Description
    costant specie transfer model

SourceFiles
    emptyST.C

\*---------------------------------------------------------------------------*/

#ifndef emptyST_H
#define emptyST_H

#include "specieTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class const Declaration
\*---------------------------------------------------------------------------*/

class emptyST
:
    public specieTransferModel
{
private:

public:

    //- Runtime type information
    TypeName("emptyST");


    // Constructors

        //- Construct from components
        emptyST
        (
            const volScalarField& por,
	        const volScalarField& por0
        );


    // Selectors

        //- Return a reference to the selected specieTransfer model
        static autoPtr<emptyST> New
        (
            const volScalarField& por,
	        const volScalarField& por0
        );


    //- Destructor
    virtual ~emptyST()
    {}


    // Member Functions

        //- Return the field for specie transfer between solid/gas
        virtual tmp<volScalarField> ST() const;

        //- Read specieTransferProperties dictionary
        virtual bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
