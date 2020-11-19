rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, pKin));

if (pimple.nCorrPiso() <= 1)
{
    tUEqn.clear();
}

if (pimple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *(
            fvc::flux(HbyA)
          + MRF.zeroFilter
            (
                rhorAUf*fvc::ddtCorr(rho, U, phi)/fvc::interpolate(rho)
            )
        )
    );

    MRF.makeRelative(fvc::interpolate(psi), phid);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, pKin)
          + fvm::div(phid, pKin)
          - fvm::laplacian(rhorAUf, pKin)
         ==
            fvOptions(psi, pKin, rho.name())
        );

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi == pEqn.flux();
        }
    }
}
else
{
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::flux(rho*HbyA)
          + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi))
        )
    );

    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(pKin, rho, U, phiHbyA, rhorAUf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi * porosityF, pKin)
          + fvc::div(phiHbyA)
          - fvm::laplacian(rhorAUf, pKin)
         ==
            fvOptions(psi * porosityF, pKin, rho.name())
          + Srho * (1.0 - porosityF)
        );

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
pKin.relax();

U = HbyA - rAU*fvc::grad(pKin);
U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);

if (pressureControl.limit(pKin))
{
    pKin.correctBoundaryConditions();
    rho = thermo.rho();
}

if (thermo.dpdt())
{
    dpdt = fvc::ddt(pKin);
}
