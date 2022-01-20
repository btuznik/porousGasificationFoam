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
    pipe specie transfer model

SourceFiles
    pipeST.C

\*---------------------------------------------------------------------------*/

#ifndef pipeST_H
#define pipeST_H

#include "specieTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class const Declaration
\*---------------------------------------------------------------------------*/

class pipeST
:
    public specieTransferModel
{
private:

    scalar pipeRadius_;
    const volVectorField& Up_;
    const volScalarField& rhop_;
    const volScalarField& mup_;
    const volScalarField& alphap_;

public:

    //- Runtime type information
    TypeName("pipeST");


    // Constructors

        //- Construct from components
        pipeST
        (
            const volScalarField& por,
	        const volScalarField& por0
        );

    // Selectors

        //- Return a reference to the selected specieTransfer model
        static autoPtr<pipeST> New
        (
            const volScalarField& por,
	        const volScalarField& por0
        );


    //- Destructor
    virtual ~pipeST()
    {}


    // Member Functions

        //- Return the field for specie transfer between solid/gas on surfaces
        virtual tmp<volScalarField> ST() const;

        //- Read specieTransferProperties dictionary
        virtual bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
