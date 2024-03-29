// eqZx2uHGn041

tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi)")
    )
);

{
    reaction->correct();
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (Y[i].name() != inertSpecie)
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(porosityF * rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence->muEff(), Yi)
             ==
                reaction->R(Yi)
              + fvOptions(rho, Yi)
              + pyrolysisZone.Srho(i)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve("Yi");

            fvOptions.correct(Yi);

            Yi.max(0.0);
            if (gMax(Yi) > 2.0)
                Info << "warning: mass fraction exceeded 2" << endl;
            Yt += Yi;
        }
        else
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(porosityF*rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence->muEff(), Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve("Yi");

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
    }

    forAll(Y, i)
    {
        if (Y[i].name() != inertSpecie)
        {
            Info<< " gas  " << Y[i].name() << " Srho min/max   = "
                << min(pyrolysisZone.Srho(i)()).value() << ", "
                << max(pyrolysisZone.Srho(i)()).value()
                << " min Y = " << min(Y[i]).value() <<" max Y = " << max(Y[i]).value() << endl;

            Y[i] = Y[i] / Yt;

            Info<< " min Ynorm = " << min(Y[i]).value() <<" max Ynorm = " << max(Y[i]).value() << endl;
        }
        else
        {
            Info<< " gas  " <<Y[i].name() << " is inert specie, min Y = "
                << min(Y[i]).value() <<" max Y = " << max(Y[i]).value();

            Y[i] = Y[i]/Yt;

            Info<< " min Ynorm = " << min(Y[i]).value() <<" max Ynorm = " << max(Y[i]).value() << endl;
        }
    }
}
