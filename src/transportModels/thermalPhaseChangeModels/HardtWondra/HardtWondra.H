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
    Foam::thermalPhaseChangeModels::HardtWondra

Description
    Diffuse-interface phase-change model following Hardt & Wondra (2008).

    Physical causality chain:
      heat diffusion  →  conductive flux at interface
      →  raw Stefan mass flux  →  Helmholtz redistribution
      →  smooth volumetric source  →  TEqn / alphaEqn / pEqn

    Core formula:

        mdot_raw = |∇α_s| · (k_eff ∇T · n̂) / h_lv          [kg/m³/s]

        mdot - λ² ∇²mdot = mdot_raw                (Helmholtz)

        Q_pc  = mdot · h_lv                         [W/m³]

    Key differences from temperature-relaxation models:
      - No (T-Tsat)/Δt forcing: Q_pc does not scale with 1/Δt
      - Interface area density |∇α_s| [1/m] is gradient-based and
        mesh-convergent, not the ad-hoc 6α(1-α) indicator
      - Helmholtz equation smooths gradient noise while conserving total
        mass transfer over the diffuse interface band

    The enthalpy correction Qcorr_ (Hardt & Wondra eq. 42) compensates for
    artificial enthalpy changes in cells where the Helmholtz redistribution
    places mass sources outside the physical phase-change zone.

Usage
    In constant/phaseChangeProperties:
    \verbatim
    model      HardtWondra;

    T_sat      373.15;       // saturation temperature [K]
    h_lv       2.26e6;       // latent heat [J/kg]

    k_liq      0.68;         // liquid thermal conductivity [W/m/K]
    k_vap      0.025;        // vapour thermal conductivity [W/m/K]

    nSmoothIter      5;      // alpha pseudo-diffusion iterations
    alphaSmoothWidth 1.0;    // smoothing width in cell units

    lambdaSmearCells 1.5;    // Helmholtz length: λ = n × h_cell
    liquidBiasCoeff  1.0;    // kEff bias: 1.0 = full mixture, 0.95 = liquid-biased

    mdotMax          100;    // hard clip on |mdot| [kg/m³/s]
    useEnthalpyCorrection true;
    RelaxFac         1.0;

    DilatationSource    true;
    PhaseFractionSource true;
    \endverbatim

SourceFiles
    HardtWondra.C

\*---------------------------------------------------------------------------*/

#ifndef HardtWondra_H
#define HardtWondra_H

#include "thermalPhaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{

/*---------------------------------------------------------------------------*\
                       Class HardtWondra Declaration
\*---------------------------------------------------------------------------*/

class HardtWondra
:
    public thermalPhaseChangeModel
{
    // Private data

        const fvMesh& mesh_;

        //- Volumetric latent-heat source  [W/m³]
        volScalarField Q_pc_;

        //- Enthalpy correction (Hardt & Wondra eq. 42)  [W/m³]
        volScalarField Qcorr_;

        //- Helmholtz-smoothed mass flux  [kg/m³/s]
        volScalarField mdot_;

        surfaceScalarField phiStefan_;

        //- Raw Stefan mass flux before Helmholtz smoothing  [kg/m³/s]
        volScalarField mdotRaw_;

        //- Interface area density |∇α_s|  [1/m]
        volScalarField interfaceArea_;

        //- Interface-normal conductive heat flux  [W/m³] = [W/m²]·[1/m]
        volScalarField qn_;

        //- Alpha pseudo-diffusion smoothing iterations
        label nSmoothIter_;

        //- Smoothing width in cell-width units
        scalar alphaSmoothWidth_;

        //- Helmholtz spreading length in cell units: λ = n × h_cell
        scalar lambdaSmearCells_;

        //- Conductivity bias toward liquid side [0,1]; 1.0 = full mixture
        scalar liquidBiasCoeff_;

        //- Hard clip on |mdot_|  [kg/m³/s]
        scalar mdotMax_;

        //- Under-relaxation factor applied to mdot_
        scalar RelaxFac_;

        //- Toggle for enthalpy correction term
        bool useEnthalpyCorrection_;

        //- Liquid thermal conductivity  [W/m/K]
        dimensionedScalar k_liq_;

        //- Vapour thermal conductivity  [W/m/K]
        dimensionedScalar k_vap_;


    // Private member functions

        void calcQ_pc();

        HardtWondra(const HardtWondra&) = delete;
        void operator=(const HardtWondra&) = delete;


public:

    TypeName("HardtWondra");


    // Constructor

        HardtWondra
        (
            const word& name,
            const dictionary& dict,
            const immiscibleIncompressibleTwoPhaseMixture& mixture,
            const volScalarField& T,
            const volScalarField& alpha1
        );


    //- Destructor
    virtual ~HardtWondra() = default;


    // Member functions

        //- Return volumetric latent-heat source  [W/m³]
        virtual tmp<volScalarField> Q_pc()   const override;

        //- Return enthalpy correction term  [W/m³]
        virtual tmp<volScalarField> Qcorr()  const override;

        virtual tmp<surfaceScalarField> phiStefan() const override;

        //- Return harmonic-mean face conductivity  [W/m/K]
        virtual tmp<surfaceScalarField> kappaf() const override;

        //- Recompute phase-change fields from current T and alpha1
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
