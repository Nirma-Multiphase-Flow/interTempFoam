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
#include "surfaceInterpolate.H"
#include "fvmSup.H"

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
    sw_alpha1Gen_(dict.lookupOrDefault<Switch>("PhaseFractionSource", Switch(true))),
    cp1_
    (
        "cp1", dimEnergy/dimMass/dimTemperature,
        mixture.subDict(mixture.phase1Name()).get<scalar>("cp")
    ),
    cp2_
    (
        "cp2", dimEnergy/dimMass/dimTemperature,
        mixture.subDict(mixture.phase2Name()).get<scalar>("cp")
    ),
    Pr1_("Pr1", dimless, mixture.subDict(mixture.phase1Name()).get<scalar>("Pr")),
    Pr2_("Pr2", dimless, mixture.subDict(mixture.phase2Name()).get<scalar>("Pr"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::thermalPhaseChangeModel::kappaf() const
{
    // Arithmetic mixture conductivity k_i = rho_i nu_i cp_i / Pr_i
    const surfaceScalarField alpha1f
    (
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );

    return surfaceScalarField::New
    (
        "kappaf",
        alpha1f*mixture_.rho1()*cp1_/Pr1_
       *fvc::interpolate(mixture_.nuModel1().nu())
      + (scalar(1) - alpha1f)*mixture_.rho2()*cp2_/Pr2_
       *fvc::interpolate(mixture_.nuModel2().nu())
    );
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


Foam::tmp<Foam::fvScalarMatrix>
Foam::thermalPhaseChangeModel::TSource(const volScalarField& T) const
{
    // Legacy latent-heat handling, previously inline in TEqn.H: implicit
    // penalty that removes |Q_pc_thermal| and pins T towards Tsat where
    // Q_pc_thermal is non-zero (see audit D5/T2 for its limitations)
    const volScalarField Q(Q_pc_thermal());
    const dimensionedScalar eps_T("eps_T", dimTemperature, 0.25);
    const volScalarField Acoeff(mag(Q)/max(mag(T - T_sat_), eps_T));

    return fvm::Sp(Acoeff, T) - Acoeff*T_sat_;
}


// ************************************************************************* //
