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
    medium.H

Description
    template for porous medium shape
\*---------------------------------------------------------------------------*/


// volVectorField containig coordinates of cell centers
const volVectorField& Ccells = porosityF.mesh().C();

// loop over cells containing logical conditions for porosity shape and value
forAll(Ccells, celli)
{

// example: const value in whole space
//	porosityF[celli]=0.46;

// example: a half plane x < 0
//<<<<<<< HEAD
	if ((Ccells[celli].x() < 0.1) && (Ccells[celli].x() > -0.1) && (Ccells[celli].y() < 0.1) 
		&& (Ccells[celli].y() > -0.1) && (Ccells[celli].z() < 0.0) && (Ccells[celli].z() > -0.1))
	{
		porosityF[celli]=0.6;
//=======
//	if ((Ccells[celli].x() < 0.0016) && (Ccells[celli].x() > -0.0016) && (Ccells[celli].y() < 0.0016) 
//		&& (Ccells[celli].y() > -0.0016) && (Ccells[celli].z() < 0.0016) && (Ccells[celli].z() > -0.0016))
//	{
//       		porosityF[celli]=0.1;
//>>>>>>> bd1366d589e0f6b9f14c68b48f4e94f7c43166f8
	}
	else
	{
		porosityF[celli]=1.0;
	}

}


// base vectors of the local coordinate system in which porosity tensor is diagonal
// v1 - x, v2 - y, v3 = v1xv2
vector v1;
vector v2;
vector v3;

//v1
v1.x()=1;
v1.y()=0;
v1.z()=0;

//v2
v2.x()=0;
v2.y()=1;
v2.z()=0;

//normalization
v1=v1/mag(v1);
v2=v2/mag(v2);

//v3
v3=v1^v2;
	
//matrix of transformation to local coordinate system
tensor T_loc(v1,v2,v3);	

//porosity tensor in local coordinates (chosen for porosity tensor to be diagonal)
tensor D_loc(0,0,0,0,0,0,0,0,0);
tensor anisotropyK_loc(0,0,0,0,0,0,0,0,0);

forAll(Ccells,celli){
//values of porosity in the local coordinates 
	D_loc.xx()=16e8;
	D_loc.yy()=16e8;
	D_loc.zz()=16e8;

	anisotropyK_loc.xx()=1;
	anisotropyK_loc.yy()=1;
	anisotropyK_loc.zz()=1;
//transformation of porosity tensor to normal base in cartesian coordinates	
	Df[celli] = (T_loc & D_loc & T_loc.T());
	anisotropyK[celli] = (T_loc & anisotropyK_loc & T_loc.T());
}

Df.correctBoundaryConditions();
