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

#include "thermalPhaseChangeModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalPhaseChangeModel, 0);
    defineRunTimeSelectionTable(thermalPhaseChangeModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPhaseChangeModel::thermalPhaseChangeModel
(
    const word& name,
    const dictionary& dict,
    const immiscibleIncompressibleTwoPhaseMixture& mixture,
    const volScalarField& T,
    const volScalarField& alpha1
)
:
    name_(name),
    dict_(dict),
    mixture_(mixture),
    T_(T),
    alpha1_(alpha1),
    T_sat_("T_sat", dimTemperature, dict),
    h_lv_("h_lv", dimEnergy/dimMass, dict),
    sw_PCV_(dict.lookupOrDefault<Switch>("DilatationSource", Switch(true))),
    sw_alpha1Gen_(dict.lookupOrDefault<Switch>("PhaseFractionSource", Switch(true)))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::thermalPhaseChangeModel::Qcorr() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            "Qcorr",
            T_.time().timeName(),
            T_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T_.mesh(),
        dimensionedScalar("Qcorr", dimensionSet(1,-1,-3,0,0,0,0), Zero)
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::thermalPhaseChangeModel::kappaf() const
{
    return mixture_.kappaf();
}


// PCV = (Q_pc / h_lv) * (1/rho_v - 1/rho_l)    [1/s]
// Positive for evaporation: vapour volume is created.
Foam::tmp<Foam::volScalarField>
Foam::thermalPhaseChangeModel::PCV() const
{
    if (sw_PCV_)
    {
        return (Q_pc()/h_lv_)
              *( scalar(1)/mixture_.rho2() - scalar(1)/mixture_.rho1() );
    }

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


// alpha1Gen = -Q_pc / (rho_l * h_lv)    [1/s]
// Negative for evaporation (liquid destroyed), positive for condensation.
Foam::tmp<Foam::volScalarField>
Foam::thermalPhaseChangeModel::alpha1Gen() const
{
    if (sw_alpha1Gen_)
    {
        return Q_pc() / (mixture_.rho1() * h_lv_);
    }

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

Foam::tmp<Foam::surfaceScalarField>
Foam::thermalPhaseChangeModel::phiStefan() const
{
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phiStefanZero",
                alpha1_.time().timeName(),
                alpha1_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            alpha1_.mesh(),
            dimensionedScalar
            (
                "zero",
                dimVolume/dimTime,
                Zero
            )
        )
    );
}


bool Foam::thermalPhaseChangeModel::read(const dictionary& dict)
{
    dict_ = dict;
    T_sat_ = dimensionedScalar("T_sat", dimTemperature, dict);
    h_lv_  = dimensionedScalar("h_lv",  dimEnergy/dimMass, dict);
    sw_PCV_       = dict.lookupOrDefault<Switch>("DilatationSource",    Switch(true));
    sw_alpha1Gen_ = dict.lookupOrDefault<Switch>("PhaseFractionSource", Switch(true));
    return true;
}


Foam::tmp<Foam::volScalarField>
Foam::thermalPhaseChangeModel::Q_pc_thermal() const
{
    return Q_pc();
}


// ************************************************************************* //
