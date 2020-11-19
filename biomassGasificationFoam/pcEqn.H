rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
volScalarField rAtU(1.0/(1.0/rAU - UEqn.H1()));
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
                fvc::interpolate(rho*rAU)*fvc::ddtCorr(rho, U, phi)
               /fvc::interpolate(rho)
            )
        )
    );

    MRF.makeRelative(fvc::interpolate(psi), phid);

    surfaceScalarField phic
    (
        "phic",
        fvc::interpolate(rho*(rAtU - rAU))*fvc::snGrad(pKin)*mesh.magSf()
    );

    HbyA -= (rAU - rAtU)*fvc::grad(pKin);

    volScalarField rhorAtU("rhorAtU", rho*rAtU);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, pKin)
          + fvm::div(phid, pKin)
          + fvc::div(phic)
          - fvm::laplacian(rhorAtU, pKin)
         ==
            fvOptions(psi, pKin, rho.name())
        );

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi == phic + pEqn.flux();
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
          + MRF.zeroFilter(fvc::interpolate(rho*rAU)*fvc::ddtCorr(rho, U, phi))
        )
    );

    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    phiHbyA += fvc::interpolate(rho*(rAtU - rAU))*fvc::snGrad(pKin)*mesh.magSf();
    HbyA -= (rAU - rAtU)*fvc::grad(pKin);

    volScalarField rhorAtU("rhorAtU", rho*rAtU);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(pKin, rho, U, phiHbyA, rhorAtU, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, pKin)
          + fvc::div(phiHbyA)
          - fvm::laplacian(rhorAtU, pKin)
         ==
            fvOptions(psi, pKin, rho.name())
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

U = HbyA - rAtU*fvc::grad(pKin);
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
