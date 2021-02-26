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

#include "makeHGSSolidMixtureThermo.H"

#include "constRho.H"

#include "constSolidThermo.H"
#include "exponentialSolidThermo.H"

#include "constSolidTransport.H"
#include "exponentialSolidTransport.H"
#include "linearSolidTransport.H"

#include "constSolidRad.H"

#include "HGSSolidThermo.H"

#include "multiComponentSolidMixture.H"

#include "reactingSolidHeterogeneousMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

makeSolidHeterogeneousMixtureThermoName
(
    HGSSolidThermo,
    HGSSolidMixtureThermo,
    reactingSolidHeterogeneousMixture,
    exponentialSolidTransport,
    constSolidRad,
    exponentialSolidThermo,
    constRho,
    expoHeterogeneous
);


makeSolidHeterogeneousMixtureThermoName
(
    HGSSolidThermo,
    HGSSolidMixtureThermo,
    reactingSolidHeterogeneousMixture,
    constSolidTransport,
    constSolidRad,
    constSolidThermo,
    constRho,
    constHeterogeneous
);


makeSolidHeterogeneousMixtureThermoName
(
    HGSSolidThermo,
    HGSSolidMixtureThermo,
    reactingSolidHeterogeneousMixture,
    linearSolidTransport,
    constSolidRad,
    constSolidThermo,
    constRho,
    linearHeterogeneous
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
