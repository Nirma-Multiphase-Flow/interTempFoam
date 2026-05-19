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

\*---------------------------------------------------------------------------*/

#include "StefanEquilibriumMerged.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcLaplacian.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{
    defineTypeNameAndDebug(StefanEquilibriumMerged, 0);

    addToRunTimeSelectionTable
    (
        thermalPhaseChangeModel,
        StefanEquilibriumMerged,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPhaseChangeModels::StefanEquilibriumMerged::StefanEquilibriumMerged
(
    const word& name,
    const dictionary& dict,
    const immiscibleIncompressibleTwoPhaseMixture& mixture,
    const volScalarField& T,
    const volScalarField& alpha1
)
:
    thermalPhaseChangeModel(name, dict, mixture, T, alpha1),

    mesh_(T.mesh()),

    Q_pc_
    (
        IOobject
        (
            "PhaseChangeHeat",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Q_pc", dimensionSet(1,-1,-3,0,0,0,0), Zero)
    ),

    mdot_
    (
        IOobject
        (
            "mdot",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("mdot", dimDensity/dimTime, Zero)
    ),

    deltaInterface_
    (
        IOobject
        (
            "deltaInterface",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("deltaI", dimless, Zero)
    ),

    nSmoothIter_(dict.lookupOrDefault<label>("nSmoothIter", 5)),
    alphaSmoothWidth_(dict.lookupOrDefault<scalar>("alphaSmoothWidth", 1.0)),
    deltaT_act_(dict.lookupOrDefault<scalar>("deltaT_act", 0.01)),
    RelaxFac_(dict.lookupOrDefault<scalar>("RelaxFac", 0.2))
{
    correct();
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::thermalPhaseChangeModels::StefanEquilibriumMerged::calcQ_pc()
{
    // =========================================================================
    // STEP 1 – Smooth alpha to de-noise VOF staircase
    //
    // Identical to StefanEnergyJump Step 1: N iterations of explicit
    // pseudo-diffusion scaled by Fourier number.  Working copy only —
    // alpha1_ is never modified.
    // =========================================================================

    const scalar dt = max(mesh_.time().deltaTValue(), SMALL);

    const scalar h_ref =
        pow(gMax(mesh_.V().field()), scalar(1)/scalar(3));

    const scalar Fo_per_iter = min(
        alphaSmoothWidth_ * alphaSmoothWidth_
        / max(scalar(nSmoothIter_), scalar(1)),
        scalar(0.25)
    );

    const dimensionedScalar D_smooth
    (
        "D_smooth",
        dimArea/dimTime,
        Fo_per_iter * h_ref * h_ref / dt
    );
    const dimensionedScalar dt_dim("dt_dim", dimTime, dt);

    volScalarField alphaSmooth
    (
        IOobject
        (
            "alphaSmooth_merged",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1_
    );

    for (label i = 0; i < nSmoothIter_; ++i)
    {
        alphaSmooth += fvc::laplacian(D_smooth, alphaSmooth) * dt_dim;

        alphaSmooth.primitiveFieldRef() =
            max(min(alphaSmooth.primitiveField(), scalar(1)), scalar(0));
        alphaSmooth.correctBoundaryConditions();
    }


    // =========================================================================
    // STEP 2 – Smooth interface weight  [dimensionless]
    //
    // KEY DIFFERENCE from StefanEnergyJump:
    //   Stefan:  deltaInterface_ = 6α(1-α)|∇α_smooth|   [1/m]  — gradient noise
    //   Merged:  deltaInterface_ = 6α(1-α)               [-]    — no gradients
    //
    // Peak = 1.5 at α=0.5 (interface centre).
    // Exactly zero in bulk phases (α=0 or α=1).
    // Any value > 1.5 in diagnostics indicates a bug in smoothing.
    // =========================================================================

    const volScalarField interfaceWeight =
        scalar(6) * alphaSmooth * (scalar(1) - alphaSmooth);

    deltaInterface_ = interfaceWeight;


    // =========================================================================
    // STEP 3 – VOF-mixture ρCp weighted by local volume fraction
    //
    // rhoCp_mix = α·ρ₁cp₁ + (1-α)·ρ₂cp₂
    //
    // Interface cells contain both phases; their thermal inertia must reflect
    // both contributions.  The previous one-sided selector
    //   pos(T-Tsat)*ρ₁cp₁ + neg(T-Tsat)*ρ₂cp₂
    // imposed pure-phase ρcp on mixture cells, creating a ~3350× source
    // asymmetry between the evaporation and condensation arms and a step
    // discontinuity at T = Tsat that caused checkerboard Q_pc patterns.
    // rhoCp_mix is continuous, α-weighted, and bounded by the actual
    // thermal content of each cell.  Sign of Q_pc is carried by (T-Tsat).
    // =========================================================================

    const dimensionedScalar& rho1 = mixture_.rho1();
    const dimensionedScalar& rho2 = mixture_.rho2();
    const dimensionedScalar& cp1  = mixture_.cp1();
    const dimensionedScalar& cp2  = mixture_.cp2();

    const volScalarField superheat = T_ - T_sat_;

    const volScalarField rhoCp_mix =
        alpha1_ * rho1 * cp1 + (scalar(1) - alpha1_) * rho2 * cp2;


    // =========================================================================
    // STEP 4 – Raw merged source  [W/m^3]
    //
    // Q_merged = 6αs(1-αs) · rhoCp_mix · (T-Tsat) / Δt
    //
    // Sign is carried by (T-Tsat):
    //   T > Tsat → Q > 0 → heat sink (evaporation)   ✓
    //   T < Tsat → Q < 0 → heat source (condensation) ✓
    //
    // No sign_pc multiplier: that would double-negate the condensation arm.
    // No k_liq/k_vap: thermal diffusion handled by fvm::laplacian in TEqn.H.
    // =========================================================================

    Q_pc_ = interfaceWeight * rhoCp_mix * superheat / dt_dim;


    // =========================================================================
    // STEP 5 – Activation threshold
    //
    // Zero out cells where |T-Tsat| < deltaT_act_ to prevent parasitic
    // phase change driven by numerical noise near equilibrium.
    // =========================================================================

    {
        scalarField& qFld = Q_pc_.primitiveFieldRef();
        const scalarField& tFld = T_.primitiveField();
        const scalar Tsat = T_sat_.value();

        forAll(qFld, celli)
        {
            if (Foam::mag(tFld[celli] - Tsat) < deltaT_act_)
            {
                qFld[celli] = 0;
            }
        }
    }


    // =========================================================================
    // STEP 6 – Fluid availability limiter
    //
    // Cannot evaporate more liquid than present, nor condense more vapour.
    //   LimEvap = α1 · ρ₁ · h_lv / Δt    [W/m^3]
    //   LimCond = (1-α1) · ρ₂ · h_lv / Δt [W/m^3]
    //
    // No WallField_ term: 6α(1-α) ≈ 0 in pure-liquid wall cells (α≈1),
    // so wall artefacts that required WallField_ in InterfaceEquilibrium
    // do not arise here.
    // =========================================================================

    const dimensionedScalar dt_lim("dt_lim", dimTime, max(dt, SMALL));

    const volScalarField LimEvap
    (
        IOobject
        (
            "LimEvap_merged",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1_ * rho1 * h_lv_ / dt_lim
    );

    const volScalarField LimCond
    (
        IOobject
        (
            "LimCond_merged",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (scalar(1) - alpha1_) * rho2 * h_lv_ / dt_lim
    );

    Q_pc_ = max(-LimCond, min(Q_pc_, LimEvap));


    // =========================================================================
    // STEP 7 – Volume-dilatation CFL limiter
    //
    // Scale down Q_pc in cells where |PCV|·Δt > 0.5 to prevent alpha
    // shooting through [0,1] bounds.
    // =========================================================================

    {
        const scalar alpha_CFL = 0.5;
        const scalar vol_fac =
            (1.0/rho2.value() - 1.0/rho1.value()) / h_lv_.value();

        scalarField& qFld = Q_pc_.primitiveFieldRef();

        forAll(qFld, celli)
        {
            const scalar pcv_abs = Foam::mag(qFld[celli]) * vol_fac * dt;
            if (pcv_abs > alpha_CFL)
            {
                qFld[celli] *= alpha_CFL / pcv_abs;
            }
        }
    }


    // =========================================================================
    // STEP 8 – Under-relaxation
    //
    // Default 0.2 (vs Stefan's 0.05): safe because there is no gradient noise
    // to amplify through Q_pc.
    // =========================================================================

    Q_pc_ *= RelaxFac_;


    // =========================================================================
    // STEP 9 – Derived mass-transfer rate (diagnostic)
    // =========================================================================

    mdot_ = Q_pc_ / h_lv_;


    // =========================================================================
    // STEP 10 – Diagnostics
    //
    // max(6α(1-α)) ≤ 1.5 is a runtime sanity check.
    // Any value > 1.5 indicates a bug in the smoothing step.
    // =========================================================================

    Info<< "StefanEquilibriumMerged:" << nl
        << "  ∫Q_pc dV        = "
        << gSum(Q_pc_.primitiveField() * mesh_.V().field())  << " W"    << nl
        << "  ∫mdot dV        = "
        << gSum(mdot_.primitiveField() * mesh_.V().field())  << " kg/s" << nl
        << "  max(6α(1-α))    = "
        << gMax(deltaInterface_.primitiveField())             << " [-]"  << nl
        << "  min/max(T)      = "
        << gMin(T_.primitiveField()) << " / "
        << gMax(T_.primitiveField())                          << " K"    << nl
        << "  max/min(Q_pc)   = "
        << gMax(Q_pc_.primitiveField()) << " / "
        << gMin(Q_pc_.primitiveField())                       << " W/m3" << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::thermalPhaseChangeModels::StefanEquilibriumMerged::Q_pc() const
{
    return Q_pc_;
}


void Foam::thermalPhaseChangeModels::StefanEquilibriumMerged::correct()
{
    calcQ_pc();
}


bool Foam::thermalPhaseChangeModels::StefanEquilibriumMerged::read
(
    const dictionary& dict
)
{
    thermalPhaseChangeModel::read(dict);

    nSmoothIter_      = dict.lookupOrDefault<label>("nSmoothIter",      5);
    alphaSmoothWidth_ = dict.lookupOrDefault<scalar>("alphaSmoothWidth", 1.0);
    deltaT_act_       = dict.lookupOrDefault<scalar>("deltaT_act",       0.01);
    RelaxFac_         = dict.lookupOrDefault<scalar>("RelaxFac",         0.2);

    return true;
}


// ************************************************************************* //
