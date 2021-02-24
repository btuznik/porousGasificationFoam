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
    Foam::heterogeneousPyrolysisModel

Description
    base class for heterogeneusPyrolysisModel

SourceFiles
    heterogeneousPyrolysisModelI.H
    heterogeneousPyrolysisModel.C

\*---------------------------------------------------------------------------*/

/** @file
 * Interface of base class for pyrolysis models
 */

#ifndef heterogeneousPyrolysisModel_H
#define heterogeneousPyrolysisModel_H

#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "psiReactionThermo.H"
#include "HGSSolidThermo.H"
#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;
class Time;

namespace heterogeneousPyrolysisModels
{

/*---------------------------------------------------------------------------*\
                      Class heterogeneousPyrolysisModel Declaration
\*---------------------------------------------------------------------------*/

/**
 * The heterogeneousPyrolysisModel class.
 * It is an implementation of a base class for different pyrolysis models. This
 * class creates all necessary fields of solid properties (e.g., mass fields of
 * solid components) and integrates physical submodels for the solid state
 * (e.g., radiation, kinetics, conduction). For each cell containing the solid
 * phase (where γ < 1), the source terms for the gas- and solid-phase equations
 * are calculated.
 */
class heterogeneousPyrolysisModel
:
    public IOdictionary
{
protected:

    // Protected Data

        /// Reference to the mesh database
        const fvMesh& mesh_;

        /// Reference to the time database
        const Time& time_;

        /// Active flag
        Switch active_;

        /// Active information output
        Switch infoOutput_;

        /// Model coefficients dictionary
        dictionary coeffs_;

        //-  Numerical control parameter
        scalar reactionDeltaMin_;

    // Protected Member Functions

        ///- Read control parameters
        virtual bool read();

        ///- Read control parameters from dictionary
        virtual bool read(const dictionary& dict);

private:

    // Private Member Functions

        //- Read pyrolysis controls
        void readPyrolysisControls();

        ///- Disallow default bitwise copy construct
        heterogeneousPyrolysisModel(const heterogeneousPyrolysisModel&);

        ///- Disallow default bitwise assignment
        void operator=(const heterogeneousPyrolysisModel&);

public:

    //- Runtime type information
    TypeName("heterogeneousPyrolysisModel");

        declareRunTimeSelectionTable
        (
             autoPtr,
             heterogeneousPyrolysisModel,
             mesh,
             (
                const word& modelType,
                const fvMesh& mesh
             ),
             (modelType, mesh)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            heterogeneousPyrolysisModel,
            noRadiation,
            (
                const word& modelType,
                const fvMesh& mesh,
                HGSSolidThermo& solidThermo,
                psiReactionThermo& gasThermo,
                volScalarField& whereIs
             ),
             (modelType, mesh, solidThermo, gasThermo, whereIs)
        );
        
        declareRunTimeSelectionTable
        (
            autoPtr,
            heterogeneousPyrolysisModel,
            radiation,
            (
                const word& modelType,
                const fvMesh& mesh,
                HGSSolidThermo& solidThermo,
                psiReactionThermo& gasThermo,
                volScalarField& whereIs,
                volScalarField& radiation
             ),
             (modelType, mesh, solidThermo, gasThermo, whereIs, radiation)
        );

    // Constructors

        ///- Construct from type name and mesh
        heterogeneousPyrolysisModel
        (
            const word& modelType,
            const fvMesh& mesh
        );

        ///- Cunstruct from type name and mesh and gas species list
        heterogeneousPyrolysisModel
        (
            const word& modelType,
            const fvMesh& mesh,
            HGSSolidThermo& solidThermo,
            psiReactionThermo& gasThermo,
            volScalarField& whereIs
        );
        
        //- Constructor from type name, mesh, gas species list and radiation
        heterogeneousPyrolysisModel
        (
            const word& modelType,
            const fvMesh& mesh,
            HGSSolidThermo& solidThermo,
            psiReactionThermo& gasThermo,
            volScalarField& whereIs,
            volScalarField& radiation
        );

    // Selectors

        //- Return a reference to the selected pyrolysis model
        static autoPtr<heterogeneousPyrolysisModel> New
        (
          const fvMesh& mesh,
          HGSSolidThermo& solidThermo,
          psiReactionThermo& gasThermo,
          volScalarField& whereIs
        );
        
        static autoPtr<heterogeneousPyrolysisModel> New
        (
          const fvMesh& mesh,
          HGSSolidThermo& solidThermo,
          psiReactionThermo& gasThermo,
          volScalarField& whereIs,
          volScalarField& radiation
        );

    //- Destructor
    virtual ~heterogeneousPyrolysisModel();

    // Member Functions

        // Access

            // Fields

                //- Return density [kg/m3]
                virtual volScalarField& rhoConst() const = 0;

                //- Return const temperature [K]
                virtual const volScalarField& T() const = 0;

                //- Return specific heat capacity [J/kg/K]
                virtual const tmp<volScalarField> Cp() const = 0;

                //- Return the region absorptivity [1/m]
                virtual const volScalarField& kappa() const = 0;

                //- Return the region thermal conductivity [W/m/k]
                virtual const volScalarField& K() const = 0;

                //- Return list of cells in surface layer inside solid
                virtual const List<label>& surf() const = 0;


        // Evolution

            /// Pre-evolve region
            virtual void preEvolveRegion();

            /// Evolve the region
            virtual void evolveRegion();

            /// Evolve the film
            virtual void evolve();


        // Helper function

            //- Mean diffusion number of the solid region
            virtual scalar solidRegionDiffNo() const;

            //- Return max diffusivity allowed in the solid
            virtual scalar maxDiff() const;

            /// Returm time modification factor
            virtual scalar maxTime() const;

            /// Return equilbirium model switch
            virtual Switch equilibrium() const;

            ///- Return the solution dictionary
            inline const dictionary& solution() const;

            //- Return the model coefficients dictionary
            inline const dictionary& coeffs() const;

            //- Return the active flag
            inline const Switch& active() const;

        // Source fields

            //- Return total mass source - Eulerian phase only
            virtual Foam::tmp<Foam::volScalarField> Srho() const;

            //- Return mass source for specie i - Eulerian phase only
            virtual Foam::tmp<Foam::volScalarField> Srho
            (
                const label i
            ) const = 0;

            //- Return enthalpy source - Eulerian phase only
            virtual Foam::tmp<Foam::volScalarField> heatTransfer();

            //- Return enthalpy needed to heat up gas
            virtual Foam::tmp<Foam::volScalarField> heatUpGas() const;

            //- Return enthalpy needed to heat up gas
            virtual Foam::tmp<Foam::volScalarField> solidChemistrySh() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heterogeneousPyrolysisModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "heterogeneousPyrolysisModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
