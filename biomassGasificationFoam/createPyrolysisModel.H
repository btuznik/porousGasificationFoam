Info<< "\nCreating pyrolysis model" << endl;

typedef typename heterogeneousPyrolysisModels::heterogeneousPyrolysisModel pyrolysisModelType;

autoPtr<pyrolysisModelType> pyrolysis
(
    pyrolysisModelType::New(mesh, thermo, porosityF)
);

pyrolysisModelType& pyrolysisZone = pyrolysis();

// TODO
/*
volScalarField Srho
(
    IOobject
    (
        "Srho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pyrolysisZone.Srho()()
);
*/

volScalarField heatTransfer
(
    IOobject
    (
        "heatTransfer",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pyrolysisZone.heatTransfer()
);

