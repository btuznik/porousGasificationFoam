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

\*---------------------------------------------------------------------------*/

#include "ODESolidHeterogeneousChemistryModel.H"
#include "multiComponentMixture.H"
#include "UniformField.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "reactingSolidHeterogeneousMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SolidThermo, class SolidThermoType, class GasThermoType>
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::ODESolidHeterogeneousChemistryModel
(
    const SolidThermo& solidThermo,
    PtrList<volScalarField>& gasPhaseGases
)
:
    BasicSolidChemistryModel<SolidThermo>(solidThermo),
    ODESystem(),
    mesh_(solidThermo.mesh()),
    gasPhaseGases_(gasPhaseGases),
    Ys_(this->solidThermo().composition().Y()),
    solidThermo_
    (
        static_cast<const reactingSolidHeterogeneousMixture<SolidThermoType>&>
            (this->solidThermo().composition()).solidData()
    ),
    gasThermo_(
        dynamic_cast<const multiComponentMixture<GasThermoType>&>
        (mesh_.lookupObject <psiReactionThermo> ("thermophysicalProperties")).specieThermos()
    ),
    pyrolisisGases_
    (
        mesh_.lookupObject<dictionary>
            ("chemistryProperties").lookup("species")
    ),
    reactions_
    (
        static_cast<const reactingSolidHeterogeneousMixture<SolidThermoType>&>
            (this->solidThermo().composition())
    ),
    nGases_(pyrolisisGases_.size()),
    nSpecie_(Ys_.size()+ nGases_),
    nSolids_(Ys_.size()),
    nReaction_(reactions_.size()),
    Treact_
    (
        BasicSolidChemistryModel<SolidThermo>::template lookupOrDefault<scalar>
        (
            "Treact",
            0
        )
    ),
    RRs_(nSolids_),
    RRg_(nGases_),
    RR_(nSpecie_),
    c_(nSpecie_),
    dcdt_(nSpecie_),
    shReactionHeat_
    (
        IOobject
        (
            "shReactionHeat",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),
    coeffs_(nSpecie_ + 2),
    Ys0_(nSolids_),
    cellCounter_(0),
    reactingCells_(mesh_.nCells(), true),
    V_
    (
        IOobject
        (
            "V_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
    dimensionedScalar("zero",dimVolume,0.0)
    ),
    solidReactionEnergyFromEnthalpy_
    (
        mesh_.lookupObject<dictionary>
            ("chemistryProperties").lookupOrDefault("solidReactionEnergyFromEnthalpy",true)
    ),
    stoichiometricReactions_
    (
        mesh_.lookupObject<dictionary>
            ("chemistryProperties").lookupOrDefault("stoichiometricReactions",false)
    ),
    diffusionLimitedReactions_
    (
        mesh_.lookupObject<dictionary>
            ("chemistryProperties").lookupOrDefault("diffusionLimitedReactions",false)
    ),
    diffusionLimitedReactionsAlpha_
    (
        true
    ),
    solidReactionDeltaEnergy_(0.0),
    showRRR_
    (
        mesh_.lookupObject<dictionary>
            ("chemistryProperties").lookupOrDefault("showRelativeReactionRates",false)
    ),
    rhoG_
    (
        mesh_.lookupObject<volScalarField>("rho")
    ),
    rho_
    (
        solidThermo.rho()
    ),
    porosityF_
    (
      mesh_.lookupObject<volScalarField>("porosityF")
    ),
    ST_
    (
        mesh_.lookupObject<volScalarField>("STvol")
    ),
    specieConcentration_(nSpecie_, 0.0),
    tauC_(0.0),
    dt_(0.0)
{
    // create the fields for the chemistry sources
    Info << "gases in gas phase: " << gasPhaseGases_.size() << " \n" << endl;

    forAll(gasPhaseGases_,i)
    {
        Info << gasPhaseGases_[i].name() << " \n";
    }

    Info << "gases from pyrolysis: " << endl;

    forAll(pyrolisisGases_,i)
    {
        Info << pyrolisisGases_[i] << " \n";
    }

    forAll(RRs_, fieldI)
    {
        RRs_.set
        (
            fieldI,
            new scalarField(mesh_.nCells(), 0.0)
        );


        IOobject header
        (
            Ys_[fieldI].name() + "0",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ
        );

        // check if field exists and can be read
        if (header.typeHeaderOk<volScalarField>(true))
        {
            Ys0_.set
            (
                fieldI,
                new volScalarField
                (
                    IOobject
                    (
                        Ys_[fieldI].name() + "0",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_
                )
            );
        }
        else
        {
            volScalarField Y0Default
            (
                IOobject
                (
                    "Y0Default",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            );

            Ys0_.set
            (
                fieldI,
                new volScalarField
                (
                    IOobject
                    (
                        Ys_[fieldI].name() + "0",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    Y0Default
                )
            );

            // Calculate inital values of Ysi0 = rho*delta*Yi
            Ys0_[fieldI] = this->solidThermo().rho() * Ys_[fieldI] * dimensionedScalar("tmp",dimVol/dimMass,1.);
            Ys0_[fieldI].ref() *= mesh_.V();
        }
    }

    V_.ref() = mesh_.V();

    forAll(RRg_, fieldI)
    {
        RRg_.set(fieldI, new scalarField(mesh_.nCells(), 0.0));
    }

    Info<< "ODESolidHeterogeneousChemistryModel: Number of solids = " << nSolids_
        << " and reactions = " << nReaction_ << endl;

    Info<< "Number of gases from pyrolysis = " << nGases_ << endl;

    forAll(reactions_, i)
    {
        Info<< indent << "Reaction " << i << nl << reactions_[i] << nl;
    }

    Info<< "solidReactionEnergyFromEnthalpy " << solidReactionEnergyFromEnthalpy_ << nl;
    Info<< "stoichiometricReactions " << stoichiometricReactions_ << nl;

    gasDictionary_.resize(nGases_);
    forAll(gasDictionary_,gasI)
    {
        forAll(gasPhaseGases_,gasJ)
        {
            if (gasPhaseGases_[gasJ].name() == pyrolisisGases_[gasI])
            {
                gasDictionary_[gasI] = gasJ;
            }
        }
    }

    Info << "diffusionLimitedReactions " << diffusionLimitedReactions_ << nl;

    if ( diffusionLimitedReactions_ )
    {
        word STmodelName
        (
            IOdictionary
            (
                IOobject
                (
                    "specieTransferProperties",
                    mesh_.time().constant(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            ).lookup("specieTransferModel")
        );
        diffusionLimitedReactionsAlpha_ =  !( STmodelName == "constST" );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class SolidThermo, class SolidThermoType, class GasThermoType>
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::
~ODESolidHeterogeneousChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class SolidThermo, class SolidThermoType, class GasThermoType>
void Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::
setCellReacting(const label cellI, const bool active)
{
    reactingCells_[cellI] = active;
}

template<class SolidThermo, class SolidThermoType, class GasThermoType>
Foam::scalarField  Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::omega
(
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalarField& rR,
    const bool updateC0
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    const label cellI = cellCounter_;
    label nEq = nEqns();

    if (not solidReactionEnergyFromEnthalpy_)
    {
        nEq = nEqns()+1;
    }

    scalarField om(nEq,0.0);

    forAll(reactions_, i)
    {
        const solidHeterogeneousReaction& R = reactions_[i];

        scalar omegai = omega
        (
            R, c, T, 0.0, pf, cf, lRef, pr, cr, rRef
        );

    scalar massCoefficient = 1;

    scalar substrates = 0;
    scalar products = 0;

    scalar solidSubstrates = 0;
    scalar solidProducts = 0;

        scalar massStream = 0;

        if (stoichiometricReactions_)
        {
            forAll(R.grhs(), g)
            {
                label gi = gasDictionary_[R.grhs()[g]];
                products += gasThermo_[gi].W()*R.grhsSto()[g];
            }

            forAll(R.glhs(), g)
            {
                label gi = R.glhs()[g];
                substrates += gasThermo_[gi].W()*R.glhsSto()[g];
            }

            forAll(R.slhs(), s)
            {
                solidSubstrates += R.slhsSto()[s];
            }

            forAll(R.srhs(), s)
            {
                solidProducts += R.srhsSto()[s];
            }

            if ((substrates + solidSubstrates  == 0) || (products + solidProducts == 0))
            {
                FatalErrorIn("omega")
                    << "Reaction:\n" << R
                    << "\nis not really a reaction" << exit(FatalError);
            }

            if (solidSubstrates > solidProducts)
            {
                scalar sr = solidProducts/solidSubstrates;
                massCoefficient = 1./(products-substrates);

                forAll(R.slhs(), s)
                {
                    label si = R.slhs()[s];
                    om[si] -= omegai*R.slhsSto()[s]/solidSubstrates;
                    massStream -= omegai*R.slhsSto()[s]/solidSubstrates;
                    if (updateC0)
                    {
                        Ys0_[si][cellI] = this->solidThermo().rho()[cellI] *Ys_[si][cellI] * V_[cellI];
                    }
                }

                forAll(R.srhs(), s)
                {
                    label si = R.srhs()[s];
                    om[si] += sr*omegai*R.srhsSto()[s]/solidProducts;
                    if (updateC0)
                    {
                        Ys0_[si][cellI] = this->solidThermo().rho()[cellI] *Ys_[si][cellI] * V_[cellI];
                    }
                }

                forAll(R.grhs(), g)
                {
                    label gi = gasDictionary_[R.grhs()[g]];
                    om[gi + nSolids_] +=  (1.0 - sr) * omegai * massCoefficient * gasThermo_[gasDictionary_[R.grhs()[g]]].W() * R.grhsSto()[g];
                }

                forAll(R.glhs(), g)
                {
                    label gi = R.glhs()[g];
                    om[gi + nSolids_] -=  (1.0 - sr)*massCoefficient*omegai*gasThermo_[gi].W()*R.glhsSto()[g];
                    massStream -= (1.0 - sr)*massCoefficient*omegai*gasThermo_[gi].W()*R.glhsSto()[g];
                }

                if (not solidReactionEnergyFromEnthalpy_)
                {
                    om[nEqns()] -= omegai*R.heatReact();
                }
            }
            else if (solidSubstrates < solidProducts)
            {
                if (solidSubstrates > 0)
                {
                    scalar sr = solidProducts/solidSubstrates;
                    massCoefficient = 1./(products-substrates);

                    forAll(R.slhs(), s)
                    {
                        label si = R.slhs()[s];
                        om[si] -= omegai*R.slhsSto()[s]/solidSubstrates;
                        massStream -= omegai*R.slhsSto()[s]/solidSubstrates;
                        if (updateC0)
                        {
                            Ys0_[si][cellI] = this->solidThermo().rho()[cellI] *Ys_[si][cellI] * V_[cellI];
                        }
                    }

                    forAll(R.srhs(), s)
                    {
                        label si = R.srhs()[s];
                        om[si] += omegai*R.srhsSto()[s]/solidSubstrates;

                        if (updateC0)
                        {
                            Ys0_[si][cellI] = this->solidThermo().rho()[cellI] *Ys_[si][cellI] * V_[cellI];
                        }
                    }

                    forAll(R.grhs(), g)
                    {
                        label gi = gasDictionary_[R.grhs()[g]];
                        om[gi + nSolids_] +=  (1.0 - sr)*omegai*massCoefficient*gasThermo_[gasDictionary_[R.grhs()[g]]].W()*R.grhsSto()[g];
                    }

                    forAll(R.glhs(), g)
                    {
                        label gi = R.glhs()[g];
                        om[gi + nSolids_] -=  (1.0 - sr)*omegai*massCoefficient*gasThermo_[gi].W()*R.glhsSto()[g];
                        massStream -= (1.0 - sr)*omegai*massCoefficient*gasThermo_[gi].W()*R.glhsSto()[g];
                    }

                    if (not solidReactionEnergyFromEnthalpy_)
                    {
                        om[nEqns()] -= omegai*R.heatReact();
                    }
                }
                else
                {
                    scalar sr = products/substrates;
                    forAll(R.grhs(), g)
                    {
                        label gi = gasDictionary_[R.grhs()[g]];
                        om[gi + nSolids_] +=  sr*omegai*gasThermo_[gasDictionary_[R.grhs()[g]]].W()*R.grhsSto()[g]/products;
                    }

                    forAll(R.glhs(), g)
                    {
                        label gi = R.glhs()[g];
                        om[gi + nSolids_] -=  omegai*gasThermo_[gi].W()*R.glhsSto()[g]/substrates;
                        massStream -= omegai*gasThermo_[gi].W()*R.glhsSto()[g]/substrates;
                    }

                    forAll(R.srhs(), s)
                    {
                        label si = R.srhs()[s];
                        om[si] += omegai*(1-sr)*R.srhsSto()[s]/solidProducts;

                        if (updateC0)
                        {
                            Ys0_[si][cellI] = this->solidThermo().rho()[cellI] *Ys_[si][cellI] * V_[cellI];
                        }
                    }
                    if (not solidReactionEnergyFromEnthalpy_)
                    {
                        om[nEqns()] -= omegai*R.heatReact();
                    }
                }
            }
            else
            {
                if (products == substrates)
                {

                    forAll(R.slhs(), s)
                    {
                        label si = R.slhs()[s];
                        om[si] -= omegai*R.slhsSto()[s]/solidSubstrates;
                        massStream -= omegai*R.slhsSto()[s]/solidSubstrates;
                        if (updateC0)
                        {
                            Ys0_[si][cellI] = this->solidThermo().rho()[cellI] *Ys_[si][cellI] * V_[cellI];
                        }
                    }

                    forAll(R.srhs(), s)
                    {
                        label si = R.srhs()[s];
                        om[si] += omegai*R.srhsSto()[s]/solidProducts;
                        if (updateC0)
                        {
                            Ys0_[si][cellI] = this->solidThermo().rho()[cellI] *Ys_[si][cellI] * V_[cellI];
                        }
                    }
                    if (products > 0)
                    {
                        forAll(R.grhs(), g)
                        {
                            label gi = gasDictionary_[R.grhs()[g]];
                            om[gi + nSolids_] +=  omegai*massCoefficient*gasThermo_[gasDictionary_[R.grhs()[g]]].W()*R.grhsSto()[g]/products;
                        }

                        forAll(R.glhs(), g)
                        {
                            label gi = R.glhs()[g];
                            om[gi + nSolids_] -=  omegai*massCoefficient*gasThermo_[gi].W()*R.glhsSto()[g]/substrates;
                            massStream -= omegai*massCoefficient*gasThermo_[gi].W()*R.glhsSto()[g]/substrates;
                        }
                    }

                    if (not solidReactionEnergyFromEnthalpy_)
                    {
                        om[nEqns()] -= omegai*R.heatReact();
                    }

                }
                else
                {
                    FatalErrorIn("omega")
                        << "Reaction:\n" << R
                        << "\ntype is not implemented" << exit(FatalError);
                }
            }
        }
        else
        {
            forAll(R.grhs(), g)
            {
                products += R.grhsSto()[g];
            }

            forAll(R.glhs(), g)
            {
                substrates += R.glhsSto()[g];
            }

            forAll(R.slhs(), s)
            {
                solidSubstrates += R.slhsSto()[s];
            }

            forAll(R.srhs(), s)
            {
                solidProducts += R.srhsSto()[s];
            }

            if ((substrates + solidSubstrates  == 0) || (products + solidProducts == 0))
            {
                FatalErrorIn("omega")
                    << "Reaction:\n" << R
                    << "\nis not really a reaction" << exit(FatalError);
            }

            scalar totalSubstrates = substrates + solidSubstrates;

            if ((totalSubstrates > 0) and (mag(totalSubstrates - (products + solidProducts)) < SMALL))
            {
                forAll(R.slhs(), s)
                {
                    label si = R.slhs()[s];
                    om[si] -= omegai*R.slhsSto()[s]/totalSubstrates;
                    massStream -= omegai*R.slhsSto()[s]/totalSubstrates;
                    if (updateC0)
                    {
                        Ys0_[si][cellI] = this->solidThermo().rho()[cellI] *Ys_[si][cellI] * V_[cellI];

                    }
                }

                forAll(R.srhs(), s)
                {
                    label si = R.srhs()[s];
                    om[si] += omegai*R.srhsSto()[s]/totalSubstrates;
                    if (updateC0)
                    {
                        Ys0_[si][cellI] = this->solidThermo().rho()[cellI] *Ys_[si][cellI] * V_[cellI];
                    }
                }

                forAll(R.grhs(), g)
                {
                    label gi = gasDictionary_[R.grhs()[g]];
                    om[gi + nSolids_] +=  omegai*R.grhsSto()[g]/totalSubstrates;
                }

                forAll(R.glhs(), g)
                {
                    label gi = R.glhs()[g];
                    om[gi + nSolids_] -= omegai*R.glhsSto()[g]/totalSubstrates;
                    massStream -= omegai*R.glhsSto()[g]/totalSubstrates;
                }
                if (not solidReactionEnergyFromEnthalpy_)
                {
                    om[nEqns()] -= omegai*R.heatReact();
                }
            }
            else
            {
                    FatalErrorIn("omega")
                        << "Reaction:\n" << R
                        << "\nviolates mass conservation" << exit(FatalError);
            }
        }
        rR[i] = mag(massStream);
    }
    return om;
}

template<class SolidThermo, class SolidThermoType, class GasThermoType>
scalar Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::omega
(
    const solidHeterogeneousReaction& R,
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    scalarField c1(nSpecie_, 0.0);

    label cellI = cellCounter_;

    for (label i=0; i<nSpecie_; i++)
    {
        c1[i] = max(0.0, c[i]);
    }

    scalar kf = R.kf(T, 0.0, c1);

    const label Nl = R.slhs().size();
    if (Nl > 0)
    {
        if ( R.glhs().size() > 0 )
        {
            for (label s=0; s < Nl; s++)
            {
            label si = R.slhs()[s];
                kf *= pow(Ys_[si][cellI],R.nReact()[s]);
            }
            forAll(R.glhs(),i)
            {
                kf *= pow(gasPhaseGases_[R.glhs()[i]].internalField()[cellI],R.nReact()[Nl+i]);
            }
        }
        else
        {
            for (label s=0; s<Nl; s++)
            {
                label si = R.slhs()[s];
                kf *= pow(Ys_[si][cellI],R.nReact()[s]);
            }
        }
        kf *= this->solidThermo().rho()[cellI];
    }
    else
    {
        forAll(R.glhs(),i)
        {
            kf *= pow(gasPhaseGases_[R.glhs()[i]].internalField()[cellI],R.nReact()[i]);
        }
        kf *= rhoG_[cellI];
    }

    if (diffusionLimitedReactions_ and (R.glhs().size() > 0 ) and (kf != 0))
    {
        scalar avKf = 1./kf;
        forAll(R.glhs(),i)
        {
            scalar addAvKf = 0.;

            if (diffusionLimitedReactionsAlpha_)
            {
                addAvKf = (ST_[cellI]*gasThermo_[R.glhs()[i]].alphah(p, T)*gasPhaseGases_[R.glhs()[i]].internalField()[cellI]);
            }
            else
            {
                addAvKf = (ST_[cellI]*gasPhaseGases_[R.glhs()[i]].internalField()[cellI]);
            }

            if (addAvKf != 0)
            {
                avKf += 1./addAvKf;
            }
            else
            {
                kf = 0;
            }
        }
        if (kf != 0)
        {
            kf = 1./avKf;
        }
    }

    return kf;
}


template<class SolidThermo, class SolidThermoType, class GasThermoType>
void Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::derivatives
(
    const scalar time,
    const scalarField& c,
    const label li,
    scalarField& dcdt
) const
{
    scalar T = c[nSpecie_];

    dcdt = 0.0;
    scalarField rR(nReaction_,0.0);

    label celli = cellCounter_;

    scalarField omegaPreq(omega(c,T,0,rR) * (1. - porosityF_[celli]));
    if (solidReactionEnergyFromEnthalpy_)
    {
        dcdt = omegaPreq;
    }
    else
    {
        forAll(dcdt,i)
        {
            dcdt[i] = omegaPreq[i];
        }
    }

    //Total mass concentration
    scalar cTot = 0.0;
    for (label i=0; i<nSolids_; i++)
    {
        cTot += c[i];
    }

    scalar newCp = 0.0;
    scalar newhi = 0.0;

    if (solidReactionEnergyFromEnthalpy_)
    {
    for (label i=0; i<nSolids_; i++)
    {
            scalar dYidt = dcdt[i];
            scalar Yi = c[i];
            newCp += Yi*solidThermo_[i].Cp(T);
            newhi -= dYidt*solidThermo_[i].hf();
    }
    }
    else
    {
        for (label i=0; i<nSolids_; i++)
        {
                scalar Yi = c[i];
                newCp += Yi*solidThermo_[i].Cp(T);
        }
        newhi += omegaPreq[nEqns()];
    }


    scalar dTdt = newhi/newCp;
    scalar dtMag = min(500.0, mag(dTdt));
    dcdt[nSpecie_] = dTdt*dtMag/(mag(dTdt) + 1.0e-10);

    // dp/dt = ...
    dcdt[nSpecie_ + 1] = 0.0;
}


template<class SolidThermo, class SolidThermoType, class GasThermoType>
void Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::jacobian
(
    const scalar t,
    const scalarField& c,
    const label li,
    scalarField& dcdt,
    scalarSquareMatrix& dfdc
) const
{

    label cellI = cellCounter_;

    scalar T = c[nSpecie_];
    scalar p = c[nSpecie_ + 1];

    scalarField c2(nSpecie_, 0.0);
    scalarField rR(nReaction_, 0.0);

    scalar solidSubstrates = 0.;
    scalar solidProducts = 0.;
    scalar substrates = 0.;
    scalar products = 0.;
    scalarField stCoeffs(nSpecie_, 0.0);

    for (label i=0; i<nSolids_; i++)
    {
        c2[i] = max(c[i], 0.0);
    }

    for (label i=0; i<nEqns(); i++)
    {
        for (label j=0; j<nEqns(); j++)
        {
            dfdc[i][j] = 0.0;
        }
    }

    scalarField omegaPreq(omega(c2,T,0.0,rR)*(1.-porosityF_[cellI]));
    if (solidReactionEnergyFromEnthalpy_)
    {
        dcdt = omegaPreq;
    }
    else
    {
        forAll(dcdt,i)
        {
            dcdt[i] = omegaPreq[i];
        }
    }

    for (label ri=0; ri<reactions_.size(); ri++)
    {
        const solidHeterogeneousReaction& R = reactions_[ri];

        const label Ns = R.slhs().size();
        const label Ng = R.glhs().size();

        solidSubstrates = 0.;
        solidProducts = 0.;
        substrates = 0.;
        products = 0.;
        stCoeffs = 0.;

        if (stoichiometricReactions_)
        {
            forAll(R.grhs(), g)
            {
                label gi = gasDictionary_[R.grhs()[g]];
                products += gasThermo_[gi].W()*R.grhsSto()[g];
            }
            forAll(R.glhs(), g)
            {
                label gi = R.glhs()[g];
                substrates += gasThermo_[gi].W()*R.glhsSto()[g];
            }
            forAll(R.slhs(), s)
            {
                solidSubstrates += R.slhsSto()[s];
            }
            forAll(R.srhs(), s)
            {
                solidProducts += R.srhsSto()[s];
            }

            if (solidSubstrates > solidProducts)
            {
                scalar sr = solidProducts/solidSubstrates;
                scalar massCoefficient = 1./(products-substrates);

                forAll(R.slhs(), s)
                {
                    label si = R.slhs()[s];
                    stCoeffs[si] -= R.slhsSto()[s]/solidSubstrates;
                }
                forAll(R.srhs(), s)
                {
                    label si = R.srhs()[s];
                    stCoeffs[si] += sr*R.srhsSto()[s]/solidProducts;
                }
                forAll(R.grhs(), g)
                {
                    label gi = gasDictionary_[R.grhs()[g]];
                    stCoeffs[gi + nSolids_] +=  (1.0 - sr)*massCoefficient*gasThermo_[gasDictionary_[R.grhs()[g]]].W()*R.grhsSto()[g];
                }
                forAll(R.glhs(), g)
                {
                    label gi = R.glhs()[g];
                    stCoeffs[gi + nSolids_] -=  (1.0 - sr)*massCoefficient*gasThermo_[gi].W()*R.glhsSto()[g];
                }
            }
            else if (solidSubstrates < solidProducts)
            {
                if (solidSubstrates > 0)
                {
                    scalar sr = solidProducts/solidSubstrates;
                    scalar massCoefficient = 1./(products-substrates);
                    forAll(R.slhs(), s)
                    {
                        label si = R.slhs()[s];
                        stCoeffs[si] -= R.slhsSto()[s]/solidSubstrates;
                    }
                    forAll(R.srhs(), s)
                    {
                        label si = R.srhs()[s];
                        stCoeffs[si] += R.srhsSto()[s]/solidSubstrates;
                    }
                    forAll(R.grhs(), g)
                    {
                        label gi = gasDictionary_[R.grhs()[g]];
                        stCoeffs[gi + nSolids_] +=  (1.0 - sr)*massCoefficient*gasThermo_[gasDictionary_[R.grhs()[g]]].W()*R.grhsSto()[g];
                    }
                    forAll(R.glhs(), g)
                    {
                        label gi = R.glhs()[g];
                        stCoeffs[gi + nSolids_] -=  (1.0 - sr)*massCoefficient*gasThermo_[gi].W()*R.glhsSto()[g];
                    }
                }
                else
                {
                    scalar sr = products/substrates;
                    forAll(R.grhs(), g)
                    {
                        label gi = gasDictionary_[R.grhs()[g]];
                        stCoeffs[gi + nSolids_] +=  sr*gasThermo_[gasDictionary_[R.grhs()[g]]].W()*R.grhsSto()[g]/products;
                    }
                    forAll(R.glhs(), g)
                    {
                        label gi = R.glhs()[g];
                        stCoeffs[gi + nSolids_] -=  gasThermo_[gi].W()*R.glhsSto()[g]/substrates;
                    }
                    forAll(R.srhs(), s)
                    {
                        label si = R.srhs()[s];
                        stCoeffs[si] += (1-sr)*R.srhsSto()[s]/solidProducts;
                    }
                }
            }
            else if (products == substrates)
            {
                forAll(R.slhs(), s)
                {
                    label si = R.slhs()[s];
                    stCoeffs[si] -= R.slhsSto()[s]/solidSubstrates;
                }
                forAll(R.srhs(), s)
                {
                    label si = R.srhs()[s];
                    stCoeffs[si] += R.srhsSto()[s]/solidProducts;
                }
                if (products > 0)
                {
                    forAll(R.grhs(), g)
                    {
                        label gi = gasDictionary_[R.grhs()[g]];
                        stCoeffs[gi + nSolids_] +=  gasThermo_[gasDictionary_[R.grhs()[g]]].W()*R.grhsSto()[g]/products;
                    }
                    forAll(R.glhs(), g)
                    {
                        label gi = R.glhs()[g];
                        stCoeffs[gi + nSolids_] -=  gasThermo_[gi].W()*R.glhsSto()[g]/substrates;
                    }
                }
            }
        }
        else
        {
            forAll(R.glhs(), g)
            {
                substrates += R.glhsSto()[g];
                stCoeffs[R.glhs()[g]+nSolids_] = -R.glhsSto()[g];
            }
            forAll(R.grhs(), g)
            {
                label gi = gasDictionary_[R.grhs()[g]];
                products += R.grhsSto()[g];
                stCoeffs[gi+nSolids_] = R.grhsSto()[g];
            }
            forAll(R.slhs(), s)
            {
                solidSubstrates += R.slhsSto()[s];
                stCoeffs[R.slhs()[s]] = -R.slhsSto()[s];
            }
            forAll(R.srhs(), s)
            {
                solidProducts += R.srhsSto()[s];
                stCoeffs[R.srhs()[s]] = R.srhsSto()[s];
            }
                stCoeffs = stCoeffs/(substrates+solidSubstrates);
        }

        scalar kf0 = R.kf(T, 0.0, c2);
        if (Ns > 0)
        {
            kf0 *= this->solidThermo().rho()[cellI];
        }
        else
        {
            kf0 *= rhoG_[cellI];
        }
        kf0 *= (1.-porosityF_[cellI]);

        for (label rSj=0; rSj < Ns + Ng; rSj++)
        {
            label sj;
            if (rSj < Ns)
            {
                sj = R.slhs()[rSj];
            }
            else
            {
                sj = R.glhs()[rSj-Ns] +  nSolids_;
            }

            scalar kf = kf0;
            scalar kfTot = kf0;

            for (label rSi=0; rSi < Ns + Ng; rSi++)
            {
                label si;
                if (rSi < Ns)
                {
                    si = R.slhs()[rSi];
                }
                else
                {
                    si = R.glhs()[rSi-Ns] +  nSolids_;
                }

                scalar el = R.nReact()[rSi];
                if (rSi < Ns)
                {
                    kfTot *= pow(Ys_[si][cellI],el);
                }
                else
                {
                    kfTot *= pow(gasPhaseGases_[R.glhs()[rSi-Ns]].internalField()[cellI],el);
                }

                if (rSi == rSj)
                {
                    if (el < 1.0)
                    {
                        if (c2[si]>SMALL)
                        {
                            if (rSi < Ns)
                            {
                                kf *= el*pow(Ys_[si][cellI] + VSMALL, el - 1.0);
                            }
                            else
                            {
                                kf *= el*pow(gasPhaseGases_[R.glhs()[rSi-Ns]].internalField()[cellI] + VSMALL, el - 1.0);
                            }
                        }
                        else
                        {
                            kf = 0.0;
                        }
                    }
                    else
                    {
                        if (rSi < Ns)
                        {
                            kf *= el*pow(Ys_[si][cellI], el - 1.0);
                        }
                        else
                        {
                            kf *= el*pow(gasPhaseGases_[R.glhs()[rSi-Ns]].internalField()[cellI], el - 1.0);
                        }
                    }
                }
                else
                {
                     if (rSi < Ns)
                     {
                         kf *= pow(Ys_[si][cellI],el);
                     }
                     else
                     {
                         kf *= pow(gasPhaseGases_[R.glhs()[rSi-Ns]].internalField()[cellI],el);
                     }
                }
            }

            if (diffusionLimitedReactions_ and (Ng > 0 ) and (kf0 != 0))
            {
                scalar avKf = 1./kf0;
                scalar kf00 = kf0;
                scalar chosenKf = 0;
                scalar chosenKf0 = 0;
                scalar addAvKf = 0.;

                for (label rSi=Ns; rSi < Ns + Ng; rSi++)
                {
                    if (diffusionLimitedReactionsAlpha_)
                    {
                        addAvKf = (ST_[cellI]*gasThermo_[R.glhs()[rSi-Ns]].alphah(p, T)*gasPhaseGases_[R.glhs()[rSi-Ns]].internalField()[cellI]);
                    }
                    else
                    {
                        addAvKf = (ST_[cellI]*gasPhaseGases_[R.glhs()[rSi-Ns]].internalField()[cellI]);
                    }

                    if (addAvKf != 0)
                    {
                        avKf += 1./addAvKf;
                        if (diffusionLimitedReactionsAlpha_)
                        {
                            chosenKf = ST_[cellI]*gasThermo_[R.glhs()[rSi-Ns]].alphah(p, T);
                            chosenKf0 = (ST_[cellI]*gasThermo_[R.glhs()[rSi-Ns]].alphah(p, T)*gasPhaseGases_[R.glhs()[rSi-Ns]].internalField()[cellI]);
                        }
                        else
                        {
                            chosenKf = ST_[cellI];
                            chosenKf0 = (ST_[cellI]*gasPhaseGases_[R.glhs()[rSi-Ns]].internalField()[cellI]);
                        }
                    }
                    else
                    {
                        kf0 = 0.;
                    }
                }
                if (kf0 != 0)
                {
                    kf0 = 1./avKf;
                }
                if ( kf00*chosenKf0 != 0 )
                {
                    kf = kf0*kf0*(kf/(kf00*kf00)+chosenKf/(chosenKf0*chosenKf0));
                }
                else
                {
                    kf = 0.;
                }
            }

            for (label rSi=0; rSi < Ns + Ng; rSi++)
            {
                label si;
                if (rSi < Ns)
                {
                    si = R.slhs()[rSi];
                }
                else
                {
                    si = R.glhs()[rSi-Ns] +  nSolids_;
                }
                dfdc[si][sj] += kf * stCoeffs[si];
            }

            forAll(R.srhs(), i)
            {
                label si = R.srhs()[i];
                dfdc[si][sj] += kf * stCoeffs[si];
            }
            forAll(R.grhs(), g)
            {
                label gi = gasDictionary_[R.grhs()[g]];
                dfdc[gi+nSolids_][sj] += kf*stCoeffs[gi+nSolids_];
            }
        }
    }

    // calculate the dcdT elements numerically
    scalar delta = 1.0e-8;

    scalarField dcdT0(dcdt);
    omegaPreq = omega(c2,T - delta ,0,rR)*(1.-porosityF_[cellI]);
    if (solidReactionEnergyFromEnthalpy_)
    {
        dcdT0 = omegaPreq;
    }
    else
    {
        forAll(dcdT0,i)
        {
            dcdT0[i] = omegaPreq[i];
        }
    }

    scalarField dcdT1(dcdt);
    omegaPreq = omega(c2,T + delta,0,rR)*(1.-porosityF_[cellI]);
    if (solidReactionEnergyFromEnthalpy_)
    {
        dcdT1 = omegaPreq;
    }
    else
    {
        forAll(dcdT1,i)
        {
            dcdT1[i] = omegaPreq[i];
        }
    }

    for (label i=0; i<nEqns(); i++)
    {
        dfdc[i][nSpecie_] = 0.5*(dcdT1[i] - dcdT0[i]) / delta;
    }
}


template<class SolidThermo, class SolidThermoType, class GasThermoType>
Foam::tmp<Foam::volScalarField>
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::tc() const
{
    tmp<volScalarField> ttc
    (
        volScalarField::New
        (
            "tc",
            this->mesh(),
            dimensionedScalar(dimTime, small),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );


    return ttc;
}

template<class SolidThermo, class SolidThermoType, class GasThermoType>
Foam::tmp<Foam::volScalarField>
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::gasHs
(
    const volScalarField& p,
    const volScalarField& T,
    const label index
) const
{

    tmp<volScalarField> tHs
    (
        new volScalarField
        (
            IOobject
            (
                "Hs_" + pyrolisisGases_[index],
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy/dimMass, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& gasHs = tHs.ref();

    const GasThermoType& mixture = gasThermo_[index];

    forAll(gasHs.internalField(), cellI)
    {
        gasHs[cellI] = mixture.Hs(p[cellI], T[cellI]);
    }

    return tHs;
}

template<class SolidThermo, class SolidThermoType, class GasThermoType>
Foam::tmp<Foam::volScalarField>
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        scalarField& Sh = tSh.ref();

        if (solidReactionEnergyFromEnthalpy_)
        {
            forAll(Ys_, i)
            {
            forAll(Sh, cellI)
            {
                scalar hf = solidThermo_[i].hf();
                Sh[cellI] -= hf*RRs_[i][cellI];
            }
            }
            forAll(Sh, cellI)
            {
            forAll(pyrolisisGases_, i)
            {
                scalar Hc = gasThermo_[gasDictionary_[i]].Hf();
                Sh[cellI] -= Hc*RRg_[i][cellI];
            }
            }
        }
        else
        {
            forAll(Sh,cellI)
            {
            Sh[cellI] = shReactionHeat_[cellI];
            }
        }
    }
    return tSh;
}


template<class SolidThermo, class SolidThermoType, class GasThermoType>
Foam::tmp<Foam::volScalarField>
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::Qdot() const
{
    tmp<volScalarField> tQdot
    (
        volScalarField::New
        (
            "Qdot",
            this->mesh_,
            dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
        )
    );

    if (this->chemistry_)
    {
        scalarField& Qdot = tQdot.ref();

        forAll(Ys_, i)
        {
            forAll(Qdot, celli)
            {
                const scalar hi = 0;
                Qdot[celli] -= hi*RR_[i][celli];
            }
        }
    }

    return tQdot;
}

template<class SolidThermo, class SolidThermoType, class GasThermoType>
Foam::tmp<Foam::volScalarField>
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::RRpor(const volScalarField T) const
{
    tmp<volScalarField> tRRpor
    (
        new volScalarField
        (
            IOobject
            (
                "RRpor",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimless/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        scalarField& RRpor = tRRpor.ref();

        forAll(Ys_, i)
        {
            forAll(RRpor, cellI)
            {
                scalar rho = solidThermo_[i].rho(T[cellI]);
                RRpor[cellI] -= RRs_[i][cellI]/rho;
            }
        }

    }
    return tRRpor;
}

template<class SolidThermo, class SolidThermoType, class GasThermoType>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::calculateRR
(
    const label ri,
    const label si
) const
{
    tmp<volScalarField::Internal> tRR
    (
        volScalarField::Internal::New
        (
            "RR",
            this->mesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, 0)
        )
    );



    return tRR;
}


template<class SolidThermo, class SolidThermoType, class GasThermoType>
void Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::calculate()
{

    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->solidThermo().rho()
    );

    resetReactionRates(rho);

    if (this->chemistry_)
    {
        forAll(rho, celli)
        {
            cellCounter_ = celli;

            const scalar delta = this->mesh().V()[celli];

            if (reactingCells_[celli])
            {
                scalar rhoi = rho[celli];
                scalar Ti = this->solidThermo().T()[celli];

                scalarField c(nSpecie_, 0.0);
                for (label i=0; i<nSolids_; i++)
                {
                    c[i] = rhoi*Ys_[i][celli]*delta;
                }
                // it is not seen from outside the class so no need to wrap it
                scalarField rR(nReaction_, 0.0);
                const scalarField dcdt = omega(c, Ti, 0.0, rR, true);
                forAll(RRs_, i)
                {
                    RRs_[i][celli] = dcdt[i]/delta;
                }

                forAll(RRg_, i)
                {
                    RRg_[i][celli] = dcdt[nSolids_ + i]/delta;
                }
            }
        }
    }
}

template<class SolidThermo, class SolidThermoType, class GasThermoType>
void
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::resetReactionRates
(
    const volScalarField& rho
)
{
    if (this->mesh().changing())
    {
        forAll(RRs_, i)
        {
            RRs_[i].setSize(rho.size());
        }
        forAll(RRg_, i)
        {
            RRg_[i].setSize(rho.size());
        }
    }

    forAll(RRs_, i)
    {
        RRs_[i] = 0.0;
    }
    forAll(RRg_, i)
    {
        RRg_[i] = 0.0;
    }
}

template<class SolidThermo, class SolidThermoType, class GasThermoType>
Foam::scalar
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::solve
(
    const scalar t0,
    const scalar deltaT
)
{
    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->solidThermo().rho()
    );

    volScalarField newDeltaTMin
    (
        IOobject
        (
            "newDeltaTMin",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, GREAT)
    );

    resetReactionRates(rho);
    shReactionHeat_ *= 0.0;

    if (!this->chemistry_)
    {
        return GREAT;
    }
    else
    {
        scalar deltaTMin = GREAT;
        newDeltaTMin = GREAT;

        forAll(rho, celli)
        {
            solveOneCell(t0, deltaT, deltaTMin, celli, rho, newDeltaTMin);
        }

        // Doesn't allow the time-step to rise.
        deltaTMin = min(deltaTMin, deltaT);
        reduce(deltaTMin,minOp<scalar>());
        scalar output = gMin(newDeltaTMin);

        return output;
    }
}

template<class SolidThermo, class SolidThermoType, class GasThermoType>
void
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::solveOneCell
(
    const scalar t0,
    const scalar deltaT,
    scalar& deltaTMin,
    const label celli,
    const volScalarField& rho,
    volScalarField& newDeltaTMin
)
{
    if (reactingCells_[celli])
    {
        cellCounter_ = celli;

        // Initialises properties of the cell.
        scalar solidRho = rho[celli] * (1. - porosityF_[celli]);
        scalar gasRho = rhoG_[celli] * porosityF_[celli];
        scalar Ti = this->solidThermo().T()[celli];

        scalarField initialSpecieConcentration(nSpecie_, 0.0);

        scalarField rR(nReaction_, 0.0);

        scalarField omegaPreq(omega(initialSpecieConcentration, Ti, 0.0, rR) * (1. - porosityF_[celli]));

        if (showRRR_ && gSum(rR) > 0)
        {
            Info<< "relative reaction rates: " << rR / gSum(rR) << endl;
        }

        for (label i = 0; i < nSolids_; ++i)
        {
            specieConcentration_[i] = solidRho * Ys_[i][celli];
        }

        for (label i = 0; i < nGases_; ++i)
        {
            specieConcentration_[nSolids_ + i] = gasRho * gasPhaseGases_[i][celli];
        }

        initialSpecieConcentration = specieConcentration_;

        calculateSourceTerms(t0, deltaT, deltaTMin, celli, solidRho, Ti, initialSpecieConcentration, omegaPreq);

        updateReactionRates(deltaT, deltaTMin, celli, solidRho, gasRho, initialSpecieConcentration, newDeltaTMin);

        if (not solidReactionEnergyFromEnthalpy_)
        {
            shReactionHeat_[celli] = omegaPreq[nEqns()];
        }
    }
}

template<class SolidThermo, class SolidThermoType, class GasThermoType>
void
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::calculateSourceTerms
(
    const scalar t0,
    const scalar deltaT,
    scalar& deltaTMin,
    const label celli,
    const scalar solidRho,
    scalar& Ti,
    const scalarField& initialSpecieConcentration,
    const scalarField& omegaPreq
)
{
    scalar t = t0;

    tauC_ = this->deltaTChem_[celli];
    dt_ = min(deltaT, tauC_);

    scalar timeLeft = deltaT;

    // Calculate the source terms.
    while (timeLeft > SMALL)
    {
        tauC_ = this->solve(specieConcentration_, Ti, 0.0, celli, t, dt_);
        t += dt_;
        // Update the temperature.
        scalar cTot = 0.0;

        // Total mass density.
        for (label i=0; i<nSolids_; i++)
        {
            cTot += specieConcentration_[i];
        }

        scalar newCp = 0.0;
        scalar newhi = 0.0;
        scalar invRho = 0.0;
        scalarList dcdt = (specieConcentration_ - initialSpecieConcentration)/dt_;

        if (solidReactionEnergyFromEnthalpy_)
        {
            for (label i=0; i < nSolids_; i++)
            {
                scalar dYi = dcdt[i];
                scalar Yi = specieConcentration_[i];
                newhi -= dYi * solidThermo_[i].hf();
                newCp += Yi * solidThermo_[i].Cp(Ti);
                invRho += Yi / solidThermo_[i].rho(Ti);
            }
        }
        else
        {
            for (label i=0; i<nSolids_; i++)
            {
                scalar Yi = specieConcentration_[i];
                newCp += Yi * solidThermo_[i].Cp(Ti);
                invRho += Yi / solidThermo_[i].rho(Ti);
            }
            newhi += omegaPreq[nEqns()];
        }

        scalar dTi = (newhi/(newCp*solidRho))*dt_;

        Ti += dTi;

        timeLeft -= dt_;
        this->deltaTChem_[celli] = tauC_;
        dt_ = min(timeLeft, tauC_);
        dt_ = max(dt_, SMALL);
    }

    deltaTMin = min(tauC_, deltaTMin);
}

template<class SolidThermo, class SolidThermoType, class GasThermoType>
void
Foam::ODESolidHeterogeneousChemistryModel<SolidThermo, SolidThermoType, GasThermoType>::updateReactionRates
(
    const scalar deltaT,
    scalar& deltaTMin,
    const label celli,
    const scalar solidRho,
    const scalar gasRho,
    const scalarField& initialSpecieConcentration,
    volScalarField& newDeltaTMin
)
{
   const scalarField changeOfConcentration = specieConcentration_ - initialSpecieConcentration;

    // the integration inside time step tends to loose accuracy
    // for conserving the mass we correct for it
    // to normalize to calculated solid loss
    scalar sourceSolid = 0.0;
    scalar sourceGas = 0.0;
    scalar sourceCorrect = 1.0;
    for (label i = 0; i<nSolids_; i++)
    {
       sourceSolid += changeOfConcentration[i] / deltaT;
    }
    for (label i = 0; i < nGases_; i++)
    {
       sourceGas += changeOfConcentration[nSolids_ + i] / deltaT;
    }

    if (sourceGas != 0)
    {
       sourceCorrect = mag(sourceSolid / sourceGas);
    }

    (void) sourceCorrect;

    forAll(RRs_, i)
    {
        RRs_[i][celli] = changeOfConcentration[i] / deltaT ;

        if ((solidRho * RRs_[i][celli] != 0) && (mag(RRs_[i][celli]) > ROOTVSMALL))
        {
            if (1. <= specieConcentration_[i] / solidRho)
            {
                newDeltaTMin[celli] = min(newDeltaTMin[celli], (1.01 - Ys_[i][celli]) / RRs_[i][celli] * solidRho);
                deltaTMin = min(deltaTMin, newDeltaTMin[celli]);

                if (1.02 < specieConcentration_[i] / solidRho)
                {
                    Info << indent << indent << " too much   " << specieConcentration_[i]/solidRho << " of " << Ys_[i].name() << " in cell " << celli << " limits time step to [s] " << deltaTMin << nl;
                }
            }

            if (0. >= specieConcentration_[i] / solidRho)
            {
                newDeltaTMin[celli] = min(newDeltaTMin[celli], (-.01-Ys_[i][celli])/RRs_[i][celli]*solidRho);
                deltaTMin = min(deltaTMin,newDeltaTMin[celli]);
                if (-0.02 > specieConcentration_[i]/solidRho)
                {
                    Info << indent << indent << " too little " << specieConcentration_[i]/solidRho << " of " << Ys_[i].name() << " in cell " << celli << " limits time step to [s] " << deltaTMin << nl;
                }
            }
            if (deltaTMin < 0)
                Info<< dt_ << " " << deltaTMin << " " << tauC_ << " an error occured: negative deltaT from solid chemistry" << endl;
        }
    }

    forAll(RRg_, i)
    {
        RRg_[i][celli] = changeOfConcentration[nSolids_ + i] / deltaT * sourceCorrect;

        if ((changeOfConcentration[nSolids_ + i] * gasRho != 0) && (mag(RRg_[i][celli]) > ROOTVSMALL))
        {
            scalar dtm = deltaTMin;
            if (1. <= specieConcentration_[nSolids_ + i] / gasRho)
            {
                newDeltaTMin[celli] = min(newDeltaTMin[celli], (1.01-gasPhaseGases_[i][celli])/RRg_[i][celli]*gasRho);
                deltaTMin = min(deltaTMin,newDeltaTMin[celli]);
                if (1.02 < specieConcentration_[nSolids_ + i]/gasRho)
                {
                    Info << indent << indent << " too much   " << specieConcentration_[nSolids_ + i]/gasRho << " of " << gasPhaseGases_[i].name() << " in cell " << celli << " limits time step to [s] " << deltaTMin << " " << newDeltaTMin[celli] << nl;
                }
            }
            if (0. >= specieConcentration_[nSolids_ + i] / gasRho)
            {
                newDeltaTMin[celli] = min(newDeltaTMin[celli], (-.01-gasPhaseGases_[i][celli])/RRg_[i][celli]*gasRho);
                deltaTMin = min(deltaTMin,newDeltaTMin[celli]);
                if (-0.02 > specieConcentration_[nSolids_ + i]/gasRho)
                {
                    Info << indent << indent << " too little " << specieConcentration_[nSolids_ + i]/gasRho << " of " << gasPhaseGases_[i].name() << " in cell " << celli << " limits time step to [s] " << deltaTMin << nl;
                }
            }
            if (deltaTMin < 0)
            {
                Info<< dt_  << " " << dtm << " " << deltaTMin << " " << tauC_ << " "
                    << (1.-gasPhaseGases_[i][celli]) << " " << specieConcentration_[nSolids_ + i]/gasRho
                    << " " << gasPhaseGases_[i][celli]  << " an error occured: negative deltaT from solid chemistry" << endl;
            }
        }
    }
}
// ************************************************************************* //
