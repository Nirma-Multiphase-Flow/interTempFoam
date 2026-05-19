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
    Foam::thermalPhaseChangeModels::StefanEquilibriumMerged

Description
    Merged Stefan + Equilibrium phase-change model for interTempFoam.

    Eliminates the multiplicative gradient noise of StefanEnergyJump by
    replacing both noisy gradient fields (|∇T·n̂| and |∇α_smooth|) with a
    single scalar driving force (T - Tsat)/dt, localised to the interface
    by 6α(1-α), using the VOF-mixture thermal capacity:

        Q_merged = 6·αs·(1-αs) · rhoCp_mix · (T - Tsat) / Δt    [W/m³]

        rhoCp_mix = α·ρ₁·cp₁ + (1-α)·ρ₂·cp₂

    Using the mixture rhoCp (rather than a one-sided pos/neg selector) avoids
    the ~3350× asymmetry between the evaporation and condensation arms and
    removes the step discontinuity at T = Tsat that caused checkerboard Q_pc.

    Sign convention (same as StefanEnergyJump):
      Q_pc > 0  heat sink   (evaporation)
      Q_pc < 0  heat source (condensation)

    (T-Tsat) carries the sign naturally; no sign_pc multiplier needed.

Usage
    In constant/phaseChangeProperties:
    \verbatim
    model      StefanEquilibriumMerged;

    T_sat      373.15;       // saturation temperature [K]
    h_lv       2.26e6;       // latent heat [J/kg]

    nSmoothIter     5;       // alpha-smoothing iterations (default 5)
    alphaSmoothWidth 1.0;    // smoothing width in cell units (default 1.0)
    deltaT_act      0.01;    // activation threshold |T-Tsat| [K] (default 0.01)
    RelaxFac        0.2;     // under-relaxation [0,1] (default 0.2)

    DilatationSource    true;
    PhaseFractionSource true;
    \endverbatim

SourceFiles
    StefanEquilibriumMerged.C

\*---------------------------------------------------------------------------*/

#ifndef StefanEquilibriumMerged_H
#define StefanEquilibriumMerged_H

#include "thermalPhaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{

/*---------------------------------------------------------------------------*\
                Class StefanEquilibriumMerged Declaration
\*---------------------------------------------------------------------------*/

class StefanEquilibriumMerged
:
    public thermalPhaseChangeModel
{
    // Private data

        const fvMesh& mesh_;

        //- Volumetric latent-heat source  [W/m^3]
        volScalarField Q_pc_;

        //- Mass transfer rate (diagnostic output)  [kg/m^3/s]
        volScalarField mdot_;

        //- Smooth interface weight 6αs(1-αs) [dimensionless], max ≤ 1.5
        volScalarField deltaInterface_;

        //- Number of alpha pseudo-diffusion smoothing iterations
        label nSmoothIter_;

        //- Smoothing width expressed as cell-width multiples
        scalar alphaSmoothWidth_;

        //- Minimum |T - Tsat| required to activate phase change  [K]
        scalar deltaT_act_;

        //- Under-relaxation factor applied to Q_pc  [0,1]
        scalar RelaxFac_;


    // Private member functions

        //- Recompute Q_pc_ from current T_ and alpha1_
        void calcQ_pc();

        //- Disallow copy and assignment
        StefanEquilibriumMerged(const StefanEquilibriumMerged&) = delete;
        void operator=(const StefanEquilibriumMerged&) = delete;


public:

    TypeName("StefanEquilibriumMerged");


    // Constructor

        StefanEquilibriumMerged
        (
            const word& name,
            const dictionary& dict,
            const immiscibleIncompressibleTwoPhaseMixture& mixture,
            const volScalarField& T,
            const volScalarField& alpha1
        );


    //- Destructor
    virtual ~StefanEquilibriumMerged() = default;


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
