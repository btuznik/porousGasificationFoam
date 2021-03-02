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
    Foam::pipeCONV

Description
    Pipe heat transfer model

    Porosity is modelled as a collection of parallel pipes. This models assumes
    that the porous medium with initial porosity \f$ \epsilon_0 \f$ has a regular
    lattice of parallel pipes with initial radius \f$ a_{c0} \f$. Constant
    number of pipes is assumed. Its radii evolve in time while porosity changes.

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

/*---------------------------------------------------------------------------*\
                           Class pipe Declaration
\*---------------------------------------------------------------------------*/

class pipeCONV
:
    public heatTransferModel
{
private:

    //- Convective heat transfer coefficient
    scalar hCoeff_;

    scalar pipeRadius_;

    //- Velocity.
    const volVectorField& Up_;

    //- Gas density
    const volScalarField& rhop_;

    //- Thermal diffusivity [kg/m/s]
    const volScalarField& alphap_;

    //- Dynamic viscosity [kg/m/s]
    const volScalarField& mup_;

    //- Reference to the fundamental fluid thermodynamic properties
    const fluidThermo& thermop_;

    // Private member functions

        //- Gas thermal conductivity
        scalar kf(const label cellI) const;

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
