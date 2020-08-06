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
    Foam::radiation

Description
    Namespace for radiation modelling

Class
    Foam::radiation::heterogeneousRadiationModel

Description
    Top level model for radiation modelling

SourceFiles
    heterogeneousRadiationModel.C

\*---------------------------------------------------------------------------*/

#ifndef heterogeneousRadiationModel_H
#define heterogeneousRadiationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "basicThermo.H"
#include "fvMatrices.H"
#include "blackBodyEmission.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

// Forward declaration of classes
class heterogeneousAbsorptionEmissionModel;
class scatterModel;

/*---------------------------------------------------------------------------*\
                       Class heterogeneousRadiationModel Declaration
\*---------------------------------------------------------------------------*/

class heterogeneousRadiationModel
:
    public IOdictionary
{
protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the time database
        const Time& time_;

        //- Reference to the gas temperature field
        const volScalarField& T_;

        //- Reference to the solid temperature field
        const volScalarField& Ts_;

        //- Model specific dictionary input parameters
        Switch radiation_;

        //- Radiation model dictionary
        dictionary coeffs_;

        //- Radiation solver frequency - number flow solver iterations per
        //  radiation solver iteration
        label solverFreq_;

        //- Flag to enable radiation model to be evaluated on first iteration
        bool firstIter_;


        // References to the radiation sub-models

            //- Absorption/emission model
            autoPtr<heterogeneousAbsorptionEmissionModel> absorptionEmission_;

            //- Scatter model
            autoPtr<scatterModel> scatter_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        heterogeneousRadiationModel(const heterogeneousRadiationModel&);

        //- Disallow default bitwise assignment
        void operator=(const heterogeneousRadiationModel&);


public:

    //- Runtime type information
    TypeName("heterogeneousRadiationModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             heterogeneousRadiationModel,
             dictionary,
             (
                 const volScalarField& T
             ),
             (T)
         );

         declareRunTimeSelectionTable
         (
             autoPtr,
             heterogeneousRadiationModel,
             porosity,
             ( 
                 const volScalarField& T,
                 const volScalarField& porosityF,
                 const List<label>& surfF,
                 const volScalarField& Ts
             ),
             (T,porosityF,surfF,Ts)
         );

    // Constructors

        //- Null constructor
        heterogeneousRadiationModel(const volScalarField& T);

        //- Construct from components
        heterogeneousRadiationModel(const word& type, const volScalarField& T);

        //- Construct from porosity, surface and solid temperature
        heterogeneousRadiationModel(const volScalarField& T, const volScalarField& porosityF, const List<label>& surfF, const volScalarField& Ts);

    // Selectors

         //- Return a reference to the selected radiation model
         static autoPtr<heterogeneousRadiationModel> New(const volScalarField& T);

         //- Return a reference to the selected radiation model
         static autoPtr<heterogeneousRadiationModel> New(const volScalarField& T, const volScalarField& porosityF, const List<label>& surfF, const volScalarField& Ts);


    //- Destructor
    virtual ~heterogeneousRadiationModel();


    // Member Functions

        // Edit

            //- Main update/correction routine
            virtual void correct();

            //- Solve radiation equation(s)
            virtual void calculate() = 0;

            //- Read radiationProperties dictionary
            virtual bool read() = 0;


        // Access

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const = 0;

            //- Source term component (constant)
            virtual tmp<DimensionedField<scalar, volMesh> > Ru() const = 0;

            //- Enthalpy source term
            virtual tmp<fvScalarMatrix> Sh(basicThermo& thermo) const;

            //- Sensible enthalpy source term
            virtual tmp<fvScalarMatrix> Shs(basicThermo& thermo) const;

            //- Enthalpy source term for solid
            virtual tmp<volScalarField> solidSh() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
