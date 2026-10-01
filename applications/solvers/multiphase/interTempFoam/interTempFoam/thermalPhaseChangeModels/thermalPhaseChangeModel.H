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
#include "fvMatrices.H"

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

    // Per-phase specific heat and Prandtl number, read from the phase
    // sub-dictionaries of transportProperties (formerly members of a
    // locally modified incompressibleTwoPhaseMixture, see audit A1)
    dimensionedScalar cp1_;
    dimensionedScalar cp2_;
    dimensionedScalar Pr1_;
    dimensionedScalar Pr2_;

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



    // Harmonic-mean face conductivity [W/m/K].
    // Default: arithmetic mixture of k_i = rho_i nu_i cp_i / Pr_i.
    // Override in models that store individual phase conductivities (HardtWondra).
    virtual tmp<surfaceScalarField> kappaf() const;

    // Non-virtual: derived algebraically from Q_pc()
    virtual tmp<volScalarField> PCV()       const;  // 1/s, > 0 for evap
    virtual tmp<volScalarField> alpha1Gen() const;  // 1/s, > 0 for cond

    // Saturation temperature accessor (shared by all derived models)
    const dimensionedScalar& T_sat() const { return T_sat_; }

    // Specific heats of phase 1 and 2
    const dimensionedScalar& cp1() const { return cp1_; }
    const dimensionedScalar& cp2() const { return cp2_; }

    virtual void correct() = 0;

    //- Latent-heat / interface terms of the temperature equation, added to
    //  its left-hand side. Default: the legacy implicit Tsat pin
    //  Acoeff (T - Tsat), Acoeff = |Q_pc_thermal|/max(|T - Tsat|, 0.25 K).
    virtual tmp<fvScalarMatrix> TSource(const volScalarField& T) const;

    //- Update after the temperature of the last outer iteration is
    //  converged (sharp models evaluate the interfacial mass transfer here
    //  for the next time step). Default: nothing.
    virtual void correctAfterT() {}

    //- Maximum interface normal speed |mdot''|/rho1 [m/s] for the
    //  phase-change time-step limit. Default: 0 (no limit).
    virtual scalar maxInterfaceSpeed() const { return 0; }
    virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //