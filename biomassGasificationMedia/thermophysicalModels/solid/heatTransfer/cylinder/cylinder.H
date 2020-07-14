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
    cylinder.C

\*---------------------------------------------------------------------------*/

#ifndef cylinder_H
#define cylinder_H

#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//namespace heatTransfer
//{

/*---------------------------------------------------------------------------*\
                           Class const Declaration
\*---------------------------------------------------------------------------*/

/**
 * cylinderHTC class.
 * Prosity is modelled as a collection of parallel pipes. This models assumes
 * that the porous medium with initial void fraction \f$ \gamma_0 \f$ has a regular
 * lattice of parallel pipes with initial radius \f$ a_0 \f$.
 */

class cylinderHTC
:
    public heatTransferModel
{
private:

    scalar HTCCoeff_;
    scalar borderHTCCoeff_;
    scalar pipeRadius_;

public:

    //- Runtime type information
    TypeName("cylinderHTC");


    // Constructors

        /** @brief Construct from porosity and porosity_0
        */
        cylinderHTC
        (
            const volScalarField& por,
            const volScalarField& por0
        );

    //- Destructor
    virtual ~cylinderHTC()
    {}


    // Member Functions

        //- Return the heatTransfer viscosity, i.e. 0 for const flow
        /** @brief Calculates the surface per volume in the porous medium.
         * Creates temporary volScalarField field and calculates surface per volume
         * in the porous medium, according to the formula:
          \f[
           \Sigma = \frac{2}{a_0} \left( \frac{\gamma_0}{\gamma} \right)^{\frac{1}{2}}
          \f]
         * @return volScalarField with surface per volume in the porous medium.
         */
        virtual tmp<volScalarField> HTC();

        //- Return the field for heat transfer between solid/gas on surfaces
        virtual tmp<volScalarField> borderHTC();

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
