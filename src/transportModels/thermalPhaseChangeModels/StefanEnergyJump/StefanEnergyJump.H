/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
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
    Foam::thermalPhaseChangeModels::StefanEnergyJump

Description
    Stefan energy-jump phase-change model for interTempFoam.

    Computes the volumetric latent-heat source Q_pc [W/m^3] from the
    conductive heat flux at the liquid-vapour interface:

        Q_pc = sign(T - Tsat) * kEff * |∇T · n̂| * δ_Γ

    where:
      - sign(T - Tsat)  selects evaporation (+) or condensation (-)
      - kEff            = k_liq  when T > Tsat  (liquid drives evaporation)
                        = k_vap  when T < Tsat  (vapour drives condensation)
      - |∇T · n̂|       = magnitude of temperature gradient normal to interface
      - δ_Γ             = 6 α (1-α) |∇α_smooth| regularised interface delta [1/m]

    The regularised delta suppresses spurious sources in single-phase cells
    and localises phase change strictly to the interface band.

    Alpha is pre-smoothed with pseudo-diffusion scaled to the mesh cell size
    so that VOF-staircase noise in |∇α| is filtered before use.

    Derived quantities (consistent with Q_pc):
      - mdot     = Q_pc / h_fg               [kg/m^3/s]   (diagnostic)
      - alpha1Gen = -Q_pc / (rho_l * h_fg)   [1/s]   → MULES source
      - PCV       = (Q_pc/h_fg)*(1/ρ_v-1/ρ_l) [1/s]  → pressure-Poisson RHS

Usage
    In constant/phaseChangeProperties:
    \verbatim
    model      StefanEnergyJump;

    T_sat      373.15;       // saturation temperature [K]
    h_lv       2.26e6;       // latent heat [J/kg]

    k_liq      0.68;         // liquid thermal conductivity [W/m/K]
    k_vap      0.025;        // vapour thermal conductivity [W/m/K]

    nSmoothIter     3;       // alpha-smoothing iterations (default 3)
    alphaSmoothWidth 1.0;    // smoothing width in cell units (default 1.0)
    deltaT_act      0.01;    // activation threshold |T-Tsat| [K] (default 0.01)
    RelaxFac        0.05;    // under-relaxation [0,1] (default 0.05)

    DilatationSource    true;
    PhaseFractionSource true;
    \endverbatim

SourceFiles
    StefanEnergyJump.C

\*---------------------------------------------------------------------------*/

#ifndef StefanEnergyJump_H
#define StefanEnergyJump_H

#include "thermalPhaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{

/*---------------------------------------------------------------------------*\
                    Class StefanEnergyJump Declaration
\*---------------------------------------------------------------------------*/

class StefanEnergyJump
:
    public thermalPhaseChangeModel
{
    // Private data

        const fvMesh& mesh_;

        //- Volumetric latent-heat source  [W/m^3]
        volScalarField Q_pc_;

        //- Mass transfer rate (diagnostic output)  [kg/m^3/s]
        volScalarField mdot_;

        //- Regularised interface delta function (diagnostic output)  [1/m]
        volScalarField deltaInterface_;

        //- Number of alpha pseudo-diffusion smoothing iterations
        label nSmoothIter_;

        //- Smoothing width expressed as cell-width multiples
        scalar alphaSmoothWidth_;

        //- Minimum |T - Tsat| required to activate phase change  [K]
        scalar deltaT_act_;

        //- Under-relaxation factor applied to Q_pc  [0,1]
        scalar RelaxFac_;

        //- Liquid thermal conductivity  [W/m/K]
        dimensionedScalar k_liq_;

        //- Vapour thermal conductivity  [W/m/K]
        dimensionedScalar k_vap_;


    // Private member functions

        //- Recompute Q_pc_ from current T_ and alpha1_
        void calcQ_pc();

        //- Disallow copy and assignment
        StefanEnergyJump(const StefanEnergyJump&) = delete;
        void operator=(const StefanEnergyJump&) = delete;


public:

    TypeName("StefanEnergyJump");


    // Constructor

        StefanEnergyJump
        (
            const word& name,
            const dictionary& dict,
            const immiscibleIncompressibleTwoPhaseMixture& mixture,
            const volScalarField& T,
            const volScalarField& alpha1
        );


    //- Destructor
    virtual ~StefanEnergyJump() = default;


    // Member functions

        //- Return the volumetric phase-change heat source  [W/m^3]
        virtual tmp<volScalarField> Q_pc() const override;

        //- Recompute Q_pc from current fields
        virtual void correct() override;

        //- Re-read model coefficients
        virtual bool read(const dictionary& dict) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalPhaseChangeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
