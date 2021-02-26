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
    pipe heat transfer model

SourceFiles
    pipe.C

\*---------------------------------------------------------------------------*/

#ifndef pipe_H
#define pipe_H

#include "heatTransferModel.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//namespace heatTransfer
//{

/*---------------------------------------------------------------------------*\
                           Class pipe Declaration
\*---------------------------------------------------------------------------*/

class pipeCONV
:
    public heatTransferModel
{
private:

    scalar hCoeff_;
    scalar pipeRadius_;
    const volVectorField& Up_;
    const volScalarField& rhop_;
    const volScalarField& alphap_;
    const volScalarField& mup_;
    const fluidThermo& thermop_;

public:

    //- Runtime type information
    TypeName("pipeCONV");


    // Constructors

        //- Construct from components
        pipeCONV
        (
            const volScalarField& porosity,
            const volScalarField& initialPorosity
        );


    // Selectors

        //- Return a reference to the selected heatTransfer model
        static autoPtr<pipeCONV> New
        (
            const volScalarField& porosity,
            const volScalarField& initialPorosity
        );


    //- Destructor
    virtual ~pipeCONV()
    {}


    // Member Functions

        //- Return the field for heat transfer between solid/gas.
        virtual tmp<volScalarField> CONV() const;

        //- Read heatTransferProperties dictionary
        virtual bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
