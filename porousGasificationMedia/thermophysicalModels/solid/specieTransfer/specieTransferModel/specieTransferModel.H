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

Namespace
    Foam::incompressible::specieTransferModels

Description
    Namespace for incompressible specieTransfer specieTransfer models.

Class
    Foam::incompressible::specieTransferModel

Description
    Abstract base class for incompressible specieTransfer models
    (RAS, LES and laminar).

SourceFiles
    specieTransferModel.C
    newTurbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef specieTransferModel_H
#define specieTransferModel_H

#include "primitiveFieldsFwd.H"
#include "Time.H"
#include "IOdictionary.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvMesh.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class specieTransferModel Declaration
\*---------------------------------------------------------------------------*/

class specieTransferModel
:
    public IOdictionary
{

protected:

    // Protected data

	    const Time& runTime_;
        const fvMesh& mesh_;
        const volScalarField& por_;
        const volScalarField& por0_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        specieTransferModel(const specieTransferModel&);

        //- Disallow default bitwise assignment
        void operator=(const specieTransferModel&);


public:

    //- Runtime type information
    TypeName("specieTransferModel");

    // Declare run-time New selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            specieTransferModel,
            porosity,
            (
                const volScalarField& por, const volScalarField& por0
            ),
            (por,por0)
        );


    // Constructors

    //- Construct from porosity
        specieTransferModel
        (
            const volScalarField& por,
	        const volScalarField& por0
        );


    // Selectors

        //- Return a reference to the selected specieTransfer model
	static autoPtr<specieTransferModel> New(const volScalarField&,const volScalarField&);


    //- Destructor
    virtual ~specieTransferModel()
    {}


    // Member Functions

        //- Access function to flux field
        inline const volScalarField& por() const
        {
            return por_;
        }

        inline const volScalarField& por0() const
        {
            return por0_;
        }
        //- Return the specieTransfer coefficient
        virtual tmp<volScalarField> ST() const = 0;

        //- Read specieTransferProperties dictionary
        virtual bool read() = 0;

        //- Default dummy write function
        virtual bool writeData(Ostream&) const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
