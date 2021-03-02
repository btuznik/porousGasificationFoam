/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    postChannel

Description
    Post-processes data from channel flow calculations.

    For each time: calculate: txx, txy,tyy, txy,
    eps, prod, vorticity, enstrophy and helicity. Assuming that the mesh
    is periodic in the x and z directions, collapse Umeanx, Umeany, txx,
    txy and tyy to a line and print them as standard output.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "volFields.H"
//#include "makeGraph.H"
 #include "graph.H"
#include "writeFile.H"

#include "fvcVolumeIntegrate.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    argList::noParallel();

    // Disable reading from constant/ and 0/ directories.
    timeSelector::addOptions(false, true);

    scalarField time;
    scalarField totalMass;

    // Get times list.
    instantList timeDirs = timeSelector::select0(runTime, args);

    // For each time step reads fields and caluclates the total mass.
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        #include "readFields.H"

        dimensionedScalar mass = fvc::domainIntegrate(rhoS * (1. - porosityF));

        time.append(runTime.value());
        totalMass.append(mass.value());
    }

    // Creates graph.
    const fileName path(runTime.rootPath()/runTime.caseName());
    mkDir(path);

    const word fileName("totalMass");
    const word format("raw");

    graph("Total mass loss", "time", "mass", time, totalMass).write(path/fileName, format);

    Info<< "\nEnd\n" << endl;

    return 0;
}
