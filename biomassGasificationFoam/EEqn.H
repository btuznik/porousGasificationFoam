if (solveEnergy)
{
    ddtKinEn = fvc::ddt(porosityF*rho, K);
    divPhiKinEn = fvc::div(phi, K);

    if (pyrolysisZone.equilibrium())
    {}
    else
    {
        volScalarField& he = thermo.he();
        volScalarField heatUpGas = pyrolysisZone.heatUpGas()();
        heatTransfer = pyrolysisZone.heatTransfer()();

        fvScalarMatrix EEqn
        (
            fvm::ddt(porosityF * rho, he) + mvConvection->fvmDiv(phi, he)
          + fvc::ddt(rho, K) + fvc::div(phi, K)
          + (
                he.name() == "e"
              ? fvc::div
                (
                    fvc::absolute(phi/fvc::interpolate(rho), U),
                    p,
                    "div(phiv,p)"
                )
              : -dpdt
            )
            + thermophysicalTransport->divq(he)
         ==
            reaction->Qdot()*(1.-porosityF)
            - (ddtKinEn + divPhiKinEn)
            + fvOptions(rho, he)
            + heatTransfer
            + heatUpGas
            //+ radiation->Shs(thermo)
        );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    Info<< "min/max(T gas) = "
        << min(T).value() << ", " << max(T).value() << endl;
    Info<< " rho gas min/max = " << min(rho).value() << ", " << max(rho).value() << endl;
            Info<< " Sh gas min/max = " << min(reaction->Qdot()).value() << ", " << max(reaction->Qdot()).value() << endl;
   }
}
