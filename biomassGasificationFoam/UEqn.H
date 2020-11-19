// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U) + fvm::div(phi, U)
  + MRF.DDt(rho, U)
  + pyrolysisZone.viscosityDrop()() * turbulence->divDevTau(U)
 ==
    fvOptions(rho, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

tmp<volScalarField> trAU;
tmp<volTensorField> trTU;

pZones.addResistance(UEqn,Df);
    
fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(pKin));

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}

trAU = 1.0/UEqn.A();
trAU.ref().rename("rAU");
