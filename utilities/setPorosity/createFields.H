    Info<< "Reading field porosityF\n" << endl;
    volScalarField porosityF
    (
        IOobject
        (
            "porosityF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("porosityF", dimensionSet(0,0,0,0,0,0,0), 0.0)
    );

    Info<< "Reading field Df\n" << endl;
    volTensorField Df
    (
        IOobject
        (
            "Df",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimless/dimLength/dimLength
    );
   
    volTensorField anisotropyK
    (
        IOobject
        (
            "anisotropyK",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimless
    ); 
