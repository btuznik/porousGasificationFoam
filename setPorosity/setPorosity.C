/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    setPorosity

Description
    Apllication for creation of porosity field and Darcy resistance tensor field
    created by Kamil Kwiatkowski & Pawel Zuk (biomassgasification.eu project)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OSspecific.H"
#include "fixedValueFvPatchFields.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#include "setRootCase.H"

#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"

#include "medium.H"

Info<<"Writing modified porosityF, Df and anisotropyK \n"<<endl;
	
porosityF.write();
Df.write();
//anisotropyK.write();
    
Info<< "End\n" << endl;

return 0;
}



//*************************************************************************//

