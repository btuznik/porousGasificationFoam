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

\*---------------------------------------------------------------------------*/

#include "heterogeneousRadiationModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//Tutaj nic nie zmieniałem ~ TADEK
Foam::autoPtr<Foam::heterogeneousRadiationModel> 
	Foam::heterogeneousRadiationModel::New
(
    const volScalarField& T
)
{
    IOobject radIO
    (
        "radiationProperties",
        T.time().constant(),
        T.mesh(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        false
    );

    word modelType("none");
    if (radIO.typeHeaderOk<IOdictionary>(false))
    {
        IOdictionary(radIO).lookup("radiationModel") >> modelType;
    }
    else
    {
        Info<< "Radiation model not active: radiationProperties not found"
            << endl;
    }

    Info<< "Selecting radiationModel " << modelType << endl;

    TConstructorTable::iterator cstrIter =
        TConstructorTablePtr_->find(modelType);

    if (cstrIter == TConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown radiationModel type "
            << modelType << nl << nl
            << "Valid radiationModel types are:" << nl
            << TConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<heterogeneousRadiationModel>(cstrIter()(T));
}

//Przekopiowałem czytanie słownika z bgf dla fe31
Foam::autoPtr<Foam::heterogeneousRadiationModel> 
	Foam::heterogeneousRadiationModel::New
(
    const volScalarField& T,
	const volScalarField& porosityF,
	const List<label>& surfF, 
	const volScalarField& Ts
)
{
    const word modelType
    (
        IOdictionary
        (
            IOobject
            (
                "radiationProperties",
                T.time().constant(),
                T.mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("heterogeneousRadiationModel")
    );

    Info<< "Selecting heterogeneousRadiationModel " << modelType << endl;

    porosityConstructorTable::iterator cstrIter =
        porosityConstructorTablePtr_->find(modelType);

    if (cstrIter == porosityConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown radiationModel type "
            << modelType << nl << nl
            << "Valid radiationModel types are:" << nl
            << porosityConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<heterogeneousRadiationModel>(cstrIter()(T,porosityF,surfF,Ts));
}


// ************************************************************************* //
