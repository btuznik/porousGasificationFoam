Info<< "Reading solid termophysical properties" << endl;

IOdictionary chemProperties
(
    IOobject
    (
        "chemistryProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

autoPtr<HGSSolidThermo> sThermo
(
    HGSSolidThermo::New
    (
         mesh,
         IOdictionary
             (
                 IOobject
                 (
                     "solidThermophysicalProperties",
                     mesh.time().constant(),
                     mesh,
                     IOobject::MUST_READ,
                     IOobject::NO_WRITE,
                     false
                 )
             ),
        Y
    )
);

HGSSolidThermo& solidThermo = sThermo();


Info<< "\nCreating pyrolysis model" << endl;

typedef typename heterogeneousPyrolysisModels::heterogeneousPyrolysisModel pyrolysisModelType;

autoPtr<pyrolysisModelType> pyrolysis
(
    pyrolysisModelType::New(mesh, solidThermo, thermo, porosityF, radiationF)
);

pyrolysisModelType& pyrolysisZone = pyrolysis();

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
    pyrolysisZone.Srho()
);

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
