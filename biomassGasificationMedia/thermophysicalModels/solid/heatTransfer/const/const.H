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
    Foam::incompressible::const

Description
    Turbulence model for const incompressible flow.

SourceFiles
    const.C

\*---------------------------------------------------------------------------*/

#ifndef const_H
#define const_H

#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//namespace heatTransfer
//{

/*---------------------------------------------------------------------------*\
                           Class const Declaration
\*---------------------------------------------------------------------------*/

class constHTC
:
    public heatTransferModel
{
private:

    scalar HTCCoeff_;
    scalar borderHTCCoeff_;
    scalar surfaceCoeff_;

public:

    //- Runtime type information
    TypeName("constHTC");


    // Constructors

        //- Construct from components
        constHTC
        (
            const volScalarField& porosity,
            const volScalarField& initialPorosity
        );


    // Selectors

        //- Return a reference to the selected heatTransfer model
        static autoPtr<constHTC> New
        (
            const volScalarField& porosity,
            const volScalarField& initialPorosity
        );


    //- Destructor
    virtual ~constHTC()
    {}


    // Member Functions

        //- Return the heatTransfer viscosity, i.e. 0 for const flow
        virtual tmp<volScalarField> HTC() const;

        //- Return the field for heat transfer between solid/gas on surfaces
        virtual tmp<volScalarField> borderHTC() const;

        //- Correct the const viscosity
        virtual void correct();

        //- Read heatTransferProperties dictionary
        virtual bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//} // End namespace heatTransfer
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
