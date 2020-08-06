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
    Foam::volPyrolysis

Description
    Reacting and radiating volumetric pyrolysis model
    created by Kamil Kwiatkowski & Pawel Zuk (biomassgasification.eu project)

SourceFiles
    volPyrolysis.C

\*---------------------------------------------------------------------------*/

#ifndef volPyrolysis_H
#define volPyrolysis_H

#include "heterogeneousPyrolysisModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heterogeneousPyrolysisModels
{


/*---------------------------------------------------------------------------*\
                      Class volPyrolysis Declaration
\*---------------------------------------------------------------------------*/

class volPyrolysis
:
    public heterogeneousPyrolysisModel
{
private:

    // Private member functions

        /// Disallow default bitwise copy construct
        volPyrolysis(const volPyrolysis&);

        /// Disallow default bitwise assignment
        void operator=(const volPyrolysis&);

        /// Read model controls
        void readReactingOneDimControls();


protected:

    ///- thermophysical properties of gaseous species.
    psiReactionThermo& gasThermo_;

    /// Switch for equilibrium/nonequilibrium model
    Switch equilibrium_;


    // Reference to solid thermo properties

        /// Density of solid region[kg/m3]
        volScalarField rho_;
        /// Non-const access to solid's temperature
        volScalarField Ts_;

        // Solution parameters

            /// Number of non-orthogonal correctors
            label nNonOrthCorr_;
            /// Maximum diffussivity
            scalar maxDiff_;

    // Reference to main region fields

        /// porosity
        volScalarField& porosity_;
        /// initial porosity
        volScalarField initialPorosity_;
        volScalarField porosityArch_;
        /// void fraction
        volScalarField voidFraction_;
        
        //- radiation
	volScalarField& radiation_;

    // Fields

        //- masks for equations
        volScalarField whereIs_;
        volScalarField whereIsNot_;

        //- Reference to heat transport model
        autoPtr<heatTransferModel> HTmodel_;
        volScalarField HTC_;

        //- border cells inside solid body
        List<label> surfF_;

    // Checks

        //- scalar max dT in a cell
        scalar maxDT_;


    // Protected member functions

        //- Read control parameters
        bool read();

        //- Read control parameters from dict
        bool read(const dictionary& dict);

   tmp<volScalarField> heatUpGasCalc() const;

public:

    //- Runtime type information
    TypeName("volPyrolysis");


    // Constructors

        //- Construct from type name, mesh and main region gas list
        volPyrolysis
        (
            const word& modelType,
            const fvMesh& mesh,
            psiReactionThermo& gasThermo,
            volScalarField& whereIs
        );
        
        //- Construct from type name, mesh, main region gas list and radiation
        volPyrolysis
        (
            const word& modelType,
            const fvMesh& mesh,
            psiReactionThermo& gasThermo,
            volScalarField& whereIs,
            volScalarField& radiation
        );

    //- Destructor
    virtual ~volPyrolysis();

    // Member Functions

        // Access

            //- Fields

                /// Return const temperature of solid phase[K]
                virtual const volScalarField& T() const;

    // parametrs from main region

        //- Returm time modification factor
        virtual Switch equilibrium() const;

    // Solution parameters

        //- Return the number of non-orthogonal correctors
        inline label nNonOrthCorr() const;

       // Evolution

            //- Pre-evolve region
            virtual void preEvolveRegion();

            //- Evolve the pyrolysis equations
            virtual void evolveRegion();

            //- Evolve porosity field
            virtual void evolvePorosity();

        //Comunication with main region

            //- Return enthalpy source - Eulerian phase only
            virtual Foam::tmp<Foam::volScalarField> heatTransfer();

            //- Return HTC
            virtual Foam::tmp<Foam::volScalarField> HTC() const;

        // I-O

            //- Provide some feedback
            virtual void info() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heterogeneousPyrolysisModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volPyrolysisI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
