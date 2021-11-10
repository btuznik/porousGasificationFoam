/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::ODESolidHeterogeneousChemistryModel

Description
    Extends base chemistry model by adding a thermo package, and ODE functions.
    Introduces chemistry equation system and evaluation of chemical source
    terms.

SourceFiles
    StandardChemistryModelI.H
    ODESolidHeterogeneousChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef ODESolidHeterogeneousChemistryModel_H
#define ODESolidHeterogeneousChemistryModel_H

#include "BasicSolidChemistryModel.H"
#include "solidHeterogeneousReaction.H"
#include "psiReactionThermo.H"
#include "ODESystem.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
            Class ODESolidHeterogeneousChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class SolidThermo, class SolidThermoType, class GasThermoType>
class ODESolidHeterogeneousChemistryModel
:
    public BasicSolidChemistryModel<SolidThermo>,
    public ODESystem
{
    // Private Member Functions

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);


protected:

    typedef SolidThermoType thermoType;


    // Protected data

        //- Mesh reference
        const fvMesh& mesh_;

        //- List of gas species present in gas phase of reaction system
        PtrList<volScalarField>& gasPhaseGases_;

        //- Reference to the field of specie mass fractions
        const PtrList<volScalarField>& Ys_;

        //- Thermodynamic data of solids
        const PtrList<SolidThermoType>& solidThermo_;

        //- Thermodynamic data of gases
       const PtrList<GasThermoType> gasThermo_;

        //- List of gas species present in reaction system
        speciesTable pyrolisisGases_;

        //- Reactions
        const PtrList<solidHeterogeneousReaction>& reactions_;

        //- Number of gas species
        label nGases_;

        //- Number of species
        label nSpecie_;

        //- Number of solid components
        label nSolids_;

        //- Number of reactions
        label nReaction_;

        //- Temperature below which the reaction rates are assumed 0
        scalar Treact_;

        //- List of reaction rate per solid [kg/m3/s]
        PtrList<scalarField> RRs_;

        //- List of reaction rate per gas [kg/m3/s]
        PtrList<scalarField> RRg_;

        //- List of reaction rate per specie [kg/m^3/s]
        PtrList<volScalarField::Internal> RR_;

        //- Temporary concentration field
        mutable scalarField c_;

        //- Temporary rate-of-change of concentration field
        mutable scalarField dcdt_;

        //- Sh based on reaction heat
        volScalarField shReactionHeat_;

        //- ODE coefficients
        scalarField coeffs_;


    // Protected Member Functions

        //- Write access to chemical source terms
        //  (e.g. for multi-chemistry model)
        inline PtrList<volScalarField::Internal>& RR();

        //- Write access to source terms for solids
        inline PtrList<scalarField>& RRs();

        //- Write access to source terms for gases
        inline PtrList<scalarField>& RRg();

private:

        //- List of accumulative solid concentrations
        mutable PtrList<volScalarField> Ys0_;

        //- Cell counter
        label cellCounter_;

        //- List of active reacting cells
        List<bool> reactingCells_;

        //- cellVolume
        volScalarField V_;

        //- Choose source of information about solid chemistry reaction energy
        Switch solidReactionEnergyFromEnthalpy_;

        //- Choose source of information about solid chemistry reaction energy
        Switch stoichiometricReactions_;

        Switch diffusionLimitedReactions_;

        //- Change of energy in time while heat of reaction approach is chosen
        scalar solidReactionDeltaEnergy_;

        //- Show the relative reaction ratios
        Switch showRRR_;

        //- Dictionary for gasName communication
        List<label> gasDictionary_;

        //- Gas density
        const volScalarField& rhoG_;

        //- Solid density
        const volScalarField& rho_;

        //- Porosity
        const volScalarField& porosityF_;

        //- Specie Transfer
        const volScalarField& ST_;

        //- Temporary concentration field
        mutable scalarField specieConcentration_;

        //- Chemical time scale
        scalar tauC_;

        //- Chemical time step
        scalar dt_;

        // Private members

            //- Set reacting status of cell, cellI
            void setCellReacting(const label cellI, const bool active);

            //- Resets the values of reaction rates and adjust size for the changing meshes.
            void resetReactionRates(const volScalarField& rho);

            void solveOneCell
            (
                const scalar t0,
                const scalar deltaT,
                scalar& deltaTMin,
                const label celli,
                const volScalarField& rho,
                volScalarField& newDeltaTMin
            );

            void calculateSourceTerms
            (
                const scalar t0,
                const scalar deltaT,
                scalar& deltaTMin,
                const label celli,
                const scalar solidRho,
                scalar& Ti,
                const scalarField& initialSpecieConcentration,
                const scalarField& omegaPreq
            );

            void updateReactionRates
            (
                const scalar deltaT,
                scalar& deltaTMin,
                const label celli,
                const scalar solidRho,
                const scalar gasRho,
                const scalarField& initialSpecieConcentration,
                volScalarField& newDeltaTMin
            );

public:

    //- Runtime type information
    TypeName("ODESolidHeterogeneousChemistryModel");

    // Constructors

        //- Construct from thermo
        ODESolidHeterogeneousChemistryModel(const SolidThermo& thermo, PtrList<volScalarField>& gasPhaseGases);

        //- Disallow default bitwise copy construction
        ODESolidHeterogeneousChemistryModel(const ODESolidHeterogeneousChemistryModel&);


    //- Destructor
    virtual ~ODESolidHeterogeneousChemistryModel();


    // Member Functions

        //- The reactions
        inline const PtrList<solidHeterogeneousReaction>& reactions() const;

        //- Thermodynamic data of gases
        inline const PtrList<GasThermoType>& gasThermo() const;

        //- Gases table
        inline const speciesTable& gasTable() const;

        //- Thermodynamic data of the species
        inline const PtrList<SolidThermoType>& specieThermos() const;

        //- The number of species
        virtual inline label nSpecie() const;

        //- The number of solids
        inline label nGases() const;

        //- The number of reactions
        virtual inline label nReaction() const;

        //- Temperature below which the reaction rates are assumed 0
        inline scalar Treact() const;

        //- Temperature below which the reaction rates are assumed 0
        inline scalar& Treact();

        //- dc/dt = omega, rate of change in concentration, for each species
        virtual scalarField omega
        (
            const scalarField& c,
            const scalar T,
            const scalar p,
            scalarField& rR,
            const bool updateC0 = false
        ) const;

        //- Return the reaction rate for reaction r and the reference
        //  species and charateristic times
        virtual scalar omega
        (
            const solidHeterogeneousReaction& r,
            const scalarField& c,
            const scalar T,
            const scalar p,
            scalar& pf,
            scalar& cf,
            label& lRef,
            scalar& pr,
            scalar& cr,
            label& rRef
        ) const;

        //- Calculates the reaction rates
        virtual void calculate();

        // Chemistry model functions (overriding abstract functions in
        // basicChemistryModel.H)

            //- Return const access to the chemical source terms for solids
            inline tmp<volScalarField> RRs(const label i) const;

            //- Return const access to the chemical source terms for gases
            inline tmp<volScalarField> RRg(const label i) const;

            //- Return total solid source term
            inline tmp<volScalarField> RRs() const;

            //- Return total gas source term
            inline tmp<volScalarField> RRg() const;

            //- Return const access to the chemical source terms for specie, i
            inline const volScalarField::Internal& RR
            (
                const label i
            ) const;

            //- Return non const access to chemical source terms [kg/m^3/s]
            virtual volScalarField::Internal& RR
            (
                const label i
            );

            //- Return reaction rate of the speciei in reactionI
            virtual tmp<volScalarField::Internal> calculateRR
            (
                const label reactionI,
                const label speciei
            ) const;

            //- Return sensible enthalpy for gas i [J/Kg]
            virtual tmp<volScalarField> gasHs
            (
                const volScalarField& p,
                const volScalarField& T,
                const label i
            ) const;

            //- Solve the reaction system for the given start time and time
            //  step and return the characteristic time
            virtual scalar solve(const scalar t0, const scalar deltaT);

            //- Return the chemical time scale
            virtual tmp<volScalarField> tc() const;

            //- Return source for enthalpy equation [kg/m/s3]
            virtual tmp<volScalarField> Sh() const;

            //- Return the heat release rate [kg/m/s^3]
            virtual tmp<volScalarField> Qdot() const;

            //- Return source for porosity equation [1/s]
            virtual tmp<volScalarField> RRpor(const volScalarField T) const;


        // ODE functions (overriding abstract functions in ODE.H)

            //- Number of ODE's to solve
            inline virtual label nEqns() const;

            //- Return access to coefficients
            virtual scalarField& coeffs();

            //- Return reference to coefficients
            virtual const scalarField& coeffs() const;

            virtual void derivatives
            (
                const scalar t,
                const scalarField& c,
                const label li,
                scalarField& dcdt
            ) const;

            virtual void jacobian
            (
                const scalar t,
                const scalarField& c,
                const label li,
                scalarField& dcdt,
                scalarSquareMatrix& J
            ) const;

            virtual scalar solve
            (
                scalarField &c,
                const scalar T,
                const scalar p,
                const label li,
                const scalar t0,
                const scalar dt
            ) const = 0;

            //- Update ODE after the solution, advancing by delta
            virtual void update(const scalar delta)
            {}



    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ODESolidHeterogeneousChemistryModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ODESolidHeterogeneousChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ODESolidHeterogeneousChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
