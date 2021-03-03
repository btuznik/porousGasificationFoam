Info<< "Creating porous zones" << nl << endl;

fieldPorosityModel pZones(mesh, porosityF);
Switch pressureImplicitPorosity(false);

