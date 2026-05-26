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

Class
    Foam::thermalPhaseChangeModel

Description
    Parent class for thermal phase change models.

SourceFiles
    thermalPhaseChangeModel.C
    newThermalPhaseChangeModel.C

\*---------------------------------------------------------------------------*/

#ifndef thermalPhaseChangeModel_H
#define thermalPhaseChangeModel_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "tmp.H"
#include "autoPtr.H"
#include "Switch.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class thermalPhaseChangeModel Declaration
\*---------------------------------------------------------------------------*/

class thermalPhaseChangeModel
{

protected:

    word       name_;
    dictionary dict_;

    const immiscibleIncompressibleTwoPhaseMixture& mixture_;
    const volScalarField& T_;
    const volScalarField& alpha1_;

    // Read from phaseChangeProperties dictionary
    dimensionedScalar T_sat_;   // saturation temperature [K]
    dimensionedScalar h_lv_;    // latent heat per unit mass [J/kg]

    Switch sw_PCV_;             // toggle: inject PCV into pEqn
    Switch sw_alpha1Gen_;       // toggle: inject alpha1Gen into alphaEqn

    // Disallow copy
    thermalPhaseChangeModel(const thermalPhaseChangeModel&) = delete;
    void operator=(const thermalPhaseChangeModel&) = delete;


public:

    TypeName("thermalPhaseChangeModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        thermalPhaseChangeModel,
        dictionary,
        (
            const word& name,
            const dictionary& dict,
            const immiscibleIncompressibleTwoPhaseMixture& mixture,
            const volScalarField& T,
            const volScalarField& alpha1
        ),
        (name, dict, mixture, T, alpha1)
    );

    static autoPtr<thermalPhaseChangeModel> New
    (
        const word& name,
        const dictionary& dict,
        const immiscibleIncompressibleTwoPhaseMixture& mixture,
        const volScalarField& T,
        const volScalarField& alpha1
    );

    thermalPhaseChangeModel
    (
        const word& name,
        const dictionary& dict,
        const immiscibleIncompressibleTwoPhaseMixture& mixture,
        const volScalarField& T,
        const volScalarField& alpha1
    );

    virtual ~thermalPhaseChangeModel() = default;


    // Pure-virtual: derived classes implement this
    virtual tmp<volScalarField> Q_pc() const = 0;   // W/m^3

    // Thermally-corrected latent source for TEqn.
    // Default returns Q_pc() — models that don't implement betaThermal blending
    // are unaffected.
    virtual tmp<volScalarField> Q_pc_thermal() const;

    // Enthalpy correction term [W/m³].
    // Non-zero only for models with Helmholtz mass redistribution (HardtWondra).
    // Default implementation returns zero so existing TEqn works unchanged.
    virtual tmp<volScalarField> Qcorr() const;

    //- Stefan transport flux
    virtual tmp<surfaceScalarField> phiStefan() const;

    // Harmonic-mean face conductivity [W/m/K].
    // Default: arithmetic mixture.kappaf().
    // Override in models that store individual phase conductivities (HardtWondra).
    virtual tmp<surfaceScalarField> kappaf() const;

    // Non-virtual: derived algebraically from Q_pc()
    virtual tmp<volScalarField> PCV()       const;  // 1/s, > 0 for evap
    virtual tmp<volScalarField> alpha1Gen() const;  // 1/s, > 0 for cond

    // Saturation temperature accessor (shared by all derived models)
    const dimensionedScalar& T_sat() const { return T_sat_; }

    virtual void correct() = 0;
    virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //