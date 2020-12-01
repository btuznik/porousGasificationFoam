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

        Info<< " rho gas min/max = " << gMin(rho) << ", " << gMax(rho) << endl;
        Info<< " Sh gas min/max = " << gMin(chemistrySh) << ", " << gMax(chemistrySh) << endl;
        Info<< " porosity min/max = " << gMin(porosityF) << ", " << gMax(porosityF) << endl;
        Info<< " heatUpGas min/max = " << gMin(heatUpGas) << ", " << gMax(heatUpGas) << endl;
        Info<< " heatTransfer min/max = " << gMin(heatTransfer) << ", " << gMax(heatTransfer) << endl;

        fvScalarMatrix EEqn
        (
            fvm::ddt(porosityF * rho, he) + mvConvection->fvmDiv(phi, he)
          + fvc::ddt(porosityF * rho, K) + fvc::div(phi, K)
          + (
                he.name() == "e"
              ? fvc::div
                (
                    fvc::absolute(phi/fvc::interpolate(porosityF * rho), U),
                    pKin,
                    "div(phiv,pKin)"
                )
              : -dpdt
            )
            + thermophysicalTransport->divq(he)
         ==
            chemistrySh*(1.-porosityF)
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

    Info<< " hs min/max = " << gMin(he) << ", " << gMax(he) << endl;
    Info<< " T gas min/max   = " << gMin(T) << ", "
        << gMax(T) << endl;
   }
}
