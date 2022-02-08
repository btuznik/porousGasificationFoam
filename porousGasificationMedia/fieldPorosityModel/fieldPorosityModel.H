/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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
    Foam::porosityModel

Description
    Porous zone definition based on cell zones.

    Porous zone definition based on cell zones and parameters obtained from a
    control dictionary constructed from the given stream. The orientation of
    the porous region is defined with the same notation as a coordinateSystem,
    but only a Cartesian coordinate system is valid.

    Implemented porosity models:

    Darcy (@e d parameter)
    @f[
        S = - (\mu \, d ) U
    @f]

    Since negative Darcy/Forchheimer parameters are invalid, they can be used
    to specify a multiplier (of the max component).

    The porousReactingZones method porousReactingZones::ddt() mirrors the normal fvm::ddt()
    method, but accounts for the effective volume of the cells.

SourceFiles
   fieldPorosityModel.C
   fieldPorosityModelTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fieldPorosityModel_H
#define fieldPorosityModel_H

#include "fvMesh.H"
#include "dictionary.H"
#include "fvMatricesFwd.H"
#include "runTimeSelectionTables.H"
#include "coordinateSystem.H"
#include "dimensionedVector.H"
#include "keyType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class porosityModel Declaration
\*---------------------------------------------------------------------------*/

class fieldPorosityModel
:
    public regIOobject
{
    // Private Member Functions

        //- Disallow default bitwise copy construction
        fieldPorosityModel(const fieldPorosityModel&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const fieldPorosityModel&) = delete;


protected:

    // Protected data

        //- Porosity name
        word name_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Forchheimer coefficient
        scalar f_;

        //- porosity Field
        volScalarField& porosityF_;

    // Protected Member Functions

        //- Df based viscous and intertial resistance
        template<class RhoFieldType>
        void addViscousInertialResistance
        (
            scalarField& Udiag,
            vectorField& Usource,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const scalarField& mu,
            const vectorField& U,
            tensorField& Df
        ) const;

        //- Return label index
        label fieldIndex(const label index) const;


public:

    //- Runtime type information
    TypeName("fieldPorosityModel");

    //- Selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        fieldPorosityModel,
        mesh,
        (
            const fvMesh& mesh,
    	    volScalarField& porosityF
        ),
        (mesh, porosityF)
    );

    //- Constructor
    fieldPorosityModel
    (
        const fvMesh& mesh,
        volScalarField& porosityF
    );


    //- Selector
    static autoPtr<fieldPorosityModel> New
    (
        const fvMesh& mesh,
        volScalarField& porosityF
    );

    //- Destructor
    virtual ~fieldPorosityModel();


    // Member Functions

        //- Return const access to the porosity model name
        inline const word& name() const;

        //- Return const access to the porosity active flag
        inline bool active() const;

        //- Transform the model data wrt mesh changes
        void transformModelData();

        //- Add resistance based on Df field
        void addResistance(fvVectorMatrix& UEqn, volTensorField& Df) const;

    // I-O

        //- Write
        bool writeData(Ostream& os) const;

        //- Inherit read from regIOobject
        using regIOobject::read;

        //- Read porosity dictionary
        bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#include "fieldPorosityModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
