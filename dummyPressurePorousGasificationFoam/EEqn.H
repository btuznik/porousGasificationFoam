// eqZx2uHGn044

if (solveEnergy)
{
    if (pyrolysisZone.equilibrium())
    {}
    else
    {
        volScalarField& he = thermo.he();
        volScalarField heatUpGas = pyrolysisZone.heatUpGas()();
        heatTransfer = pyrolysisZone.heatTransfer()();

        fvScalarMatrix EEqn
        (
            fvm::ddt(porosityF * rho, he)
          + fvm::div(phi, he)
          + fvc::ddt(porosityF * rho, K)
          + fvc::div(phi, K)
          + (
                he.name() == "e"
              ? fvc::div
                (
                    fvc::absolute(phi/fvc::interpolate(porosityF * rho), U),
                    p,
                    "div(phiv,p)"
                )
              : -dpdt
            )
          + thermophysicalTransport->divq(he)
         ==
          + fvOptions(rho, he)
          + heatTransfer
          + heatUpGas
          + radiation->Sh(thermo, he)
        );

        EEqn.relax();

        fvOptions.constrain(EEqn);

        EEqn.solve();

        fvOptions.correct(he);

        thermo.correct();

        Info<< " T gas min/max = " << gMin(T) << ", "
            << gMax(T) << endl;
        Info<< " rho gas min/max = " << min(rho).value() << ", " << max(rho).value() << endl;
                Info<< " Sh gas min/max = " << min(reaction->Qdot()).value() << ", " << max(reaction->Qdot()).value() << endl;
    }
}
