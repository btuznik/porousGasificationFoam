// Solve the Momentum equation
// eqZx2uHGn043


MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U) + fvm::div(phi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevTau(U)
 ==
    fvOptions(rho, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

tmp<volScalarField> trAU;
tmp<volTensorField> trTU;

pZones.addResistance(UEqn,Df); // eqZx2uHGn015

volScalarField rUA = 1.0 / UEqn.A();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn == -fvc::grad(p)
    );

    fvOptions.correct(U);
}
else
{
    U = rUA * (UEqn.H() - fvc::grad(p));
    U.correctBoundaryConditions();
}

trAU = 1.0/UEqn.A();
trAU.ref().rename("rAU");
