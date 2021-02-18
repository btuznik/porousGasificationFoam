    IOobject rhoSHeader
    (
        "rhos",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (!rhoSHeader.typeHeaderOk<volScalarField>(true))
    {
        Info<< "    No rhoS field" << endl;
        continue;
    }

    volScalarField rhoS
    (
        rhoSHeader,
        mesh
    );

    IOobject porosityFHeader
    (
        "porosityF",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (!porosityFHeader.typeHeaderOk<volScalarField>(true))
    {
        Info<< "    No porosityF field" << endl;
        continue;
    }

    volScalarField porosityF
    (
        porosityFHeader,
        mesh
    );

