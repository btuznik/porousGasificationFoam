/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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

#include "basicSolidChemistryModel.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::autoPtr<ChemistryModel> Foam::basicSolidChemistryModel::New
(
    const typename ChemistryModel::reactionThermo& thermo,
    PtrList<volScalarField>& gasPhaseGases,
    const word thermoName
)
{
    IOdictionary chemistryDict
    (
        IOobject
        (
            "chemistryProperties",
            thermo.db().time().constant(),
            thermo.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    if (!chemistryDict.isDict("chemistryType"))
    {
        FatalErrorInFunction
            << "Template parameter based chemistry solver selection is no "
            << "longer supported. Please create a chemistryType dictionary"
            << "instead." << endl << endl << "For example, the entry:" << endl
            << "    chemistrySolver ode<StandardChemistryModel<"
            << "rhoChemistryModel,sutherland<specie<janaf<perfectGas>,"
            << "sensibleInternalEnergy>>>>" << endl << endl << "becomes:"
            << endl << "    chemistryType" << endl << "    {" << endl
            << "        solver ode;" << endl << "        method standard;"
            << endl << "    }" << exit(FatalError);
    }

    const dictionary& chemistryTypeDict =
        chemistryDict.subDict("solidChemistryType");

    const word& solverName
    (
        chemistryTypeDict.found("solver")
      ? chemistryTypeDict.lookup("solver")
      : chemistryTypeDict.found("chemistrySolver")
      ? chemistryTypeDict.lookup("chemistrySolver")
      : chemistryTypeDict.lookup("solver") // error if neither entry is found
    );

    const word& methodName(chemistryTypeDict.lookup("method"));
    const word& solidThermoType(chemistryTypeDict.lookup("solidThermoType"));

    dictionary chemistryTypeDictNew;
    chemistryTypeDictNew.add("solver", solverName);
    chemistryTypeDictNew.add("method", methodName);
    chemistryTypeDictNew.add("solidThermoType", solidThermoType);

    Info<< "Selecting chemistry solver for solid phase " << chemistryTypeDictNew << endl;

    typedef typename ChemistryModel::thermoConstructorTable cstrTableType;
    cstrTableType* cstrTable = ChemistryModel::thermoConstructorTablePtr_;

    const word chemSolverCompThermoName =
        solverName + '<' + methodName + '<'
      + HGSSolidThermo::typeName + ','
      + solidThermoType + ',' + thermoName + ">>";

    typename cstrTableType::iterator cstrIter =
        cstrTable->find(chemSolverCompThermoName);
    if (cstrIter == cstrTable->end())
    {
        Info << "Chemistry model not found. " << endl;
        FatalErrorInFunction
            << "Unknown " << typeName_() << " type " << solverName << endl
            << endl
            << "Valid models are:" << endl
            << cstrTable->toc() << endl;

        FatalErrorInFunction << exit(FatalError);
    }
    Info << "basicSolidNew" << endl;
    return autoPtr<ChemistryModel>(cstrIter()(thermo, gasPhaseGases));
}

// ************************************************************************* //
