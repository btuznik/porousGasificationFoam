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

\*---------------------------------------------------------------------------*/

#include "pipeST.H"
#include "Time.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
      

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pipeST, 0);
addToRunTimeSelectionTable(specieTransferModel, pipeST, porosity);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pipeST::pipeST
(
    const volScalarField& por,
    const volScalarField& por0
)
:
    specieTransferModel(por,por0),
    pipeRadius_(1.0),
    Up_(db().lookupObject<volVectorField>("U")),
    rhop_(db().lookupObject<volScalarField>("rho")),
    mup_(db().lookupObject<volScalarField>("mu")),
    alphap_(db().lookupObject<volScalarField>("alpha"))
{
   read();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<pipeST> pipeST::New
(
    const volScalarField& por,
    const volScalarField& por0
)
{
    return autoPtr<pipeST>
    (
        new pipeST( por,por0)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> pipeST::ST() const
{
// eqZx2uHGn007
    Foam::tmp<Foam::volScalarField> STloc_ = Foam::tmp<Foam::volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "STloc",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "zero", dimless/dimTime, 0.0
            )
        )
    );
    
    forAll (STloc_(),cellI)
    {
        STloc_.ref()[cellI] = pow(por()[cellI], 0.5) * pow(por0()[cellI], 0.5) * 2.0 / pipeRadius_ *
            (1.83) / pipeRadius_ * alphap_[cellI];  //eqZx2uHGn019 eqZx2uHGn020
    }

    return STloc_;
}

bool pipeST::read()
{

	IOdictionary dict
        (
            IOobject
            (
                "specieTransferProperties",
                mesh_.time().constant(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

    const dictionary& params = dict.subDict("Parameters");

    params.lookup("pipeRadius") >> pipeRadius_;

    return true;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
