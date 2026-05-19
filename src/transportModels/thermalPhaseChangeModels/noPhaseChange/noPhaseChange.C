/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 Alex Rattner
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

#include "noPhaseChange.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{
    defineTypeNameAndDebug(noPhaseChange, 0);
    addToRunTimeSelectionTable
    (
        thermalPhaseChangeModel, 
        noPhaseChange, 
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPhaseChangeModels::noPhaseChange::noPhaseChange
(
        const word& name,
        const dictionary& dict,
        const immiscibleIncompressibleTwoPhaseMixture& mixture,
        const volScalarField& T,
        const volScalarField& alpha1
)
:
    thermalPhaseChangeModel
    (
        name, 
        dict, 
        mixture, 
        T, 
        alpha1
    ),
    Q_pc_
    (
        IOobject
        (
            "PhaseChangeHeat",
            T_.time().timeName(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        T.mesh(),
        dimensionedScalar( "Q_pc", dimensionSet(1,-1,-3,0,0,0,0), Zero )
    )
{

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField>
Foam::thermalPhaseChangeModels::noPhaseChange::Q_pc() const
{
    return Q_pc_;
}


Foam::tmp<Foam::volScalarField>
Foam::thermalPhaseChangeModels::noPhaseChange::PCV() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            "PCV",
            T_.time().timeName(),
            T_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T_.mesh(),
        dimensionedScalar("PCV", dimless/dimTime, Zero)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::thermalPhaseChangeModels::noPhaseChange::alpha1Gen() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            "alpha1Gen",
            T_.time().timeName(),
            T_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T_.mesh(),
        dimensionedScalar("alpha1Gen", dimless/dimTime, Zero)
    );
}


void Foam::thermalPhaseChangeModels::noPhaseChange::correct()
{}   // no-op


bool Foam::thermalPhaseChangeModels::noPhaseChange::read(const dictionary& dict)
{
    return thermalPhaseChangeModel::read(dict);
}


// ************************************************************************* //
