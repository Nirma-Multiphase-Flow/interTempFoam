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

#include "StefanEnergyJump.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "fvcLaplacian.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{
    defineTypeNameAndDebug(StefanEnergyJump, 0);

    addToRunTimeSelectionTable
    (
        thermalPhaseChangeModel,
        StefanEnergyJump,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPhaseChangeModels::StefanEnergyJump::StefanEnergyJump
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
        dimensionedScalar("deltaI", dimless/dimLength, Zero)
    ),

    nSmoothIter_(dict.lookupOrDefault<label>("nSmoothIter", 3)),
    alphaSmoothWidth_(dict.lookupOrDefault<scalar>("alphaSmoothWidth", 1.0)),
    deltaT_act_(dict.lookupOrDefault<scalar>("deltaT_act", 0.01)),
    RelaxFac_(dict.lookupOrDefault<scalar>("RelaxFac", 0.05)),

    k_liq_("k_liq", dimPower/dimLength/dimTemperature, dict),
    k_vap_("k_vap", dimPower/dimLength/dimTemperature, dict)
{
    correct();
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::thermalPhaseChangeModels::StefanEnergyJump::calcQ_pc()
{
    // =========================================================================
    // STEP 1 – Smooth alpha to de-noise VOF staircase gradients
    //
    // PROBLEM: fvc::grad(alpha1) on a VOF field is highly noisy because the
    // compressed VOF interface is essentially a staircase.  |∇α| is
    // concentrated on single faces, creating point-like source spikes.
    //
    // FIX: N iterations of explicit pseudo-diffusion:
    //         α_{n+1} = α_n + D_eff · ∆t · ∇²α_n
    //
    // SCALING: For explicit diffusion stability, Fo = D·∆t/h² ≤ 0.25.
    //   Target Fo per iteration = alphaSmoothWidth²/nSmoothIter, capped at 0.25.
    //   Total smoothed area ≈ nSmoothIter · Fo · h² = alphaSmoothWidth² · h²
    //
    // With alphaSmoothWidth=1, nSmoothIter=3: Fo≈0.33 (slightly above 0.25).
    // Since we bound alpha to [0,1] every iteration this minor over-diffusion
    // is harmless; it does not contaminate the governing equations.
    //
    // SIZING h_ref: use the global maximum cell volume to get the largest cell
    // dimension (conservative – ensures enough smoothing in all cells).
    // =========================================================================

    const scalar dt = max(mesh_.time().deltaTValue(), SMALL);

    // Characteristic cell size (cube-root of maximum cell volume)
    const scalar h_ref =
        pow(gMax(mesh_.V().field()), scalar(1)/scalar(3));

    // Clamp Fo for strict stability (adjust alphaSmoothWidth/nSmoothIter freely)
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

    // Working copy – does NOT modify alpha1_
    volScalarField alphaSmooth
    (
        IOobject
        (
            "alphaSmooth_stefan",
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

        // Bound strictly to [0,1] – prevents diffusion instability from
        // propagating source spikes outside the interface zone.
        alphaSmooth.primitiveFieldRef() =
            max(min(alphaSmooth.primitiveField(), scalar(1)), scalar(0));
        alphaSmooth.correctBoundaryConditions();
    }


    // =========================================================================
    // STEP 2 – Regularised interface delta function  δ_Γ  [1/m]
    //
    // PROBLEM: raw |∇α| is noisy (concentrated on single faces, staircase
    // artefacts, anisotropic with mesh aspect ratio).  Thermal fingers and
    // checkerboard sources directly trace to this noise.
    //
    // FIX: δ_Γ = 6·α·(1-α)·|∇α_smooth|
    //
    // PROPERTIES:
    //   • Exactly zero in single-phase cells (α=0 or α=1): suppresses all
    //     spurious bulk-phase sources regardless of gradient noise there.
    //   • Maximum at α=0.5 (6·0.5·0.5 = 1.5): localises source to interface.
    //   • Smooth across the diffuse interface: no sharp source boundaries.
    //   • Dimensions: [1/m] (same as |∇α|).
    //
    // The 6·α·(1-α) factor is computed from the SMOOTH alpha to prevent
    // the raw VOF α (which can be exactly 0 or 1 in adjacent cells) from
    // zeroing the delta function in cells where gradients are largest.
    // =========================================================================

    const dimensionedScalar eps_grad
    (
        "eps_grad",
        dimless/dimLength,
        SMALL
    );

    const volVectorField gradAlphaSmooth
    (
        IOobject
        (
            "gradAlphaSmooth_stefan",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(alphaSmooth)
    );

    const volScalarField magGradAlpha
    (
        IOobject
        (
            "magGradAlpha_stefan",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mag(gradAlphaSmooth)
    );

    // Regularised delta: localised, smooth, zero in single-phase regions.
    deltaInterface_ =
        scalar(6)*alphaSmooth*(scalar(1) - alphaSmooth)*magGradAlpha;

    // Unit interface normal (points from vapour into liquid for α1=liquid)
    const volVectorField nHat
    (
        IOobject
        (
            "nHat_stefan",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        gradAlphaSmooth / (magGradAlpha + eps_grad)
    );


    // =========================================================================
    // STEP 3 – Temperature gradient in interface-normal direction
    //
    // dTdn = ∇T · n̂  where n̂ = ∇α/|∇α| (vapour→liquid direction).
    //
    // For the canonical geometry (liquid below, vapour above, heated from wall):
    //   • ∇T points downward (hot wall below → T decreases upward)
    //   • n̂ points downward (α1=1 below, α1=0 above → gradient downward)
    //   • dTdn > 0 in the liquid near the interface (evaporation scenario)
    //
    // We use |dTdn| rather than dTdn so the amplitude is sign-agnostic about
    // the geometric orientation of the interface.  The PHYSICAL sign of Q_pc
    // is determined exclusively by sign(T - Tsat) in the next step.
    //
    // This makes the model correct for:
    //   • Horizontal interfaces (liquid below or above)
    //   • Vertical interfaces (liquid left or right)
    //   • Arbitrary interface orientation
    // =========================================================================

    const volVectorField gradT
    (
        IOobject
        (
            "gradT_stefan",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(T_)
    );

    // Normal-direction temperature gradient magnitude [K/m]
    const volScalarField magDTdn
    (
        IOobject
        (
            "magDTdn_stefan",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mag(gradT & nHat)
    );


    // =========================================================================
    // STEP 4 – Effective conductivity and sign selection
    //
    // Physics of the Stefan condition:
    //
    //   EVAPORATION (T > Tsat):
    //     Superheated liquid conducts heat toward the interface.
    //     Liquid-side conduction drives evaporation: use k_liq.
    //     Q_pc > 0 → heat removed from fluid (latent heat consumed).
    //
    //   CONDENSATION (T < Tsat):
    //     Subcooled vapour draws heat away from the interface.
    //     Vapour-side conduction drives condensation: use k_vap.
    //     Q_pc < 0 → heat added to fluid (latent heat released).
    //
    //   EQUILIBRIUM (T = Tsat):
    //     Both pos(superheat) and neg(superheat) are zero → kEff = 0 → Q_pc = 0.
    //     No driving force, no phase change. Exactly thermodynamically consistent.
    //
    // Combined formula:
    //   Q_pc = sign_pc * kEff * |dTdn| * δ_Γ
    //
    //   sign_pc = pos(T-Tsat) - neg(T-Tsat)    {+1, 0, -1}
    //   kEff    = pos(T-Tsat)*k_liq + neg(T-Tsat)*k_vap
    //
    //   => Q_pc = pos*k_liq*|dTdn|*δ - neg*k_vap*|dTdn|*δ
    //
    // NOTE: pos(x)*neg(x) = 0 for all x (never both nonzero simultaneously),
    // so the two terms cannot interfere.
    //
    // NOTE on dimension of pos() / neg():
    //   OpenFOAM's pos(volScalarField) returns a dimensionless volScalarField.
    //   T_ - T_sat_ has dimensions [K]; pos/neg strip that dimension. ✓
    // =========================================================================

    const volScalarField superheat = T_ - T_sat_;

    // kEff: k_liq for evap, k_vap for cond, 0 at Tsat
    const volScalarField kEff
    (
        IOobject
        (
            "kEff_stefan",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pos(superheat)*k_liq_ + neg(superheat)*k_vap_
    );

    // sign_pc: +1 evap, -1 cond, 0 equilibrium
    const volScalarField sign_pc
    (
        IOobject
        (
            "sign_pc_stefan",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pos(superheat) - neg(superheat)
    );


    // =========================================================================
    // STEP 5 – Raw Stefan source  [W/m^3]
    //
    //   Q_raw = sign_pc · kEff · |∇T·n̂| · δ_Γ
    //
    // Dimension check:
    //   [−] · [W/m/K] · [K/m] · [1/m] = [W/m^3] ✓
    //
    // For TEqn:  fvm::ddt(rhoCp,T) + ... - fvm::laplacian(κ,T) + Q_pc = 0
    //   Q_pc > 0 → sink  (evaporation cools the liquid)
    //   Q_pc < 0 → source (condensation heats the fluid)
    // =========================================================================

    Q_pc_ = sign_pc * kEff * magDTdn * deltaInterface_;


    // =========================================================================
    // STEP 6 – Activation threshold
    //
    // PURPOSE: Prevent parasitic phase change when T ≈ Tsat.
    //   Near equilibrium, ∇T and |∇α| are small, but numerical gradients
    //   can be nonzero.  The threshold removes sub-threshold sources
    //   that lack physical driving force.
    //
    // Default deltaT_act = 0.01 K (essentially zero for practical purposes,
    // can be raised to 0.1–1 K to damp oscillations near equilibrium).
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
    // STEP 7 – Fluid availability limiter
    //
    // PHYSICS: Cannot evaporate more liquid than is present, nor condense
    // more vapour than is present.
    //
    //   Max evaporation rate:  LimEvap = α1 · ρ_l · h_fg / dt   [W/m^3]
    //     → entire liquid fraction converts in one timestep
    //
    //   Max condensation rate: LimCond = (1-α1) · ρ_v · h_fg / dt [W/m^3]
    //     → entire vapour fraction converts in one timestep
    //
    // Constraint: -LimCond ≤ Q_pc ≤ LimEvap
    //
    // This is the PRIMARY stability limiter for large Q_pc spikes caused by
    // residual noise in |∇α| that the smoothing did not fully remove.
    // =========================================================================

    const dimensionedScalar& rho1 = mixture_.rho1();
    const dimensionedScalar& rho2 = mixture_.rho2();
    const dimensionedScalar dt_lim("dt_lim", dimTime, max(dt, SMALL));

    const volScalarField LimEvap
    (
        IOobject
        (
            "LimEvap_stefan",
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
            "LimCond_stefan",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (scalar(1) - alpha1_) * rho2 * h_lv_ / dt_lim
    );

    Q_pc_ = max(-LimCond, min(Q_pc_, LimEvap));


    // =========================================================================
    // STEP 8 – Volume-dilatation CFL limiter
    //
    // PHYSICS: Phase change produces/destroys volume via:
    //   PCV = Q_pc/h_fg · (1/ρ_v − 1/ρ_l)  [1/s]
    //
    // Numerical stability requires |PCV|·∆t ≤ α_CFL (default 0.5) per cell.
    // Violating this causes alpha to shoot through [0,1] bounds even with MULES,
    // producing pressure spikes in pEqn and timestep collapse.
    //
    // Implementation: directly scale down Q_pc in cells where PCV would
    // exceed the CFL limit.  This is a last-resort safety backstop after
    // the fluid availability limiter.
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
    // STEP 9 – Under-relaxation
    //
    // RelaxFac damps the coupling between the phase-change source and the
    // thermal/VOF equations.  For strong or stiff evaporation, large RelaxFac
    // causes oscillatory feedback:
    //   high Q → rapid alpha change → shifted interface → new Q spike → ...
    //
    // Recommended values:
    //   RelaxFac = 0.01–0.05  for initial stability testing
    //   RelaxFac = 0.1–0.3    for well-resolved, isotropic meshes
    //   RelaxFac = 1.0        only with very fine mesh + small dt + low cAlpha
    // =========================================================================

    Q_pc_ *= RelaxFac_;


    // =========================================================================
    // STEP 10 – Derived mass-transfer rate (diagnostic)
    // =========================================================================

    mdot_ = Q_pc_ / h_lv_;


    // =========================================================================
    // STEP 11 – Diagnostics
    // =========================================================================

    const scalar Q_int =
        gSum(Q_pc_.primitiveField() * mesh_.V().field());

    const scalar mdot_int =
        gSum(mdot_.primitiveField() * mesh_.V().field());

    Info<< "Stefan phase-change:" << nl
        << "  ∫Q_pc dV     = " << Q_int        << " W"      << nl
        << "  ∫mdot dV     = " << mdot_int      << " kg/s"   << nl
        << "  max(deltaI)  = "
        << gMax(deltaInterface_.primitiveField()) << " /m"   << nl
        << "  max(|∇αs|)   = "
        << gMax(magGradAlpha.primitiveField())    << " /m"   << nl
        << "  min(T)       = "
        << gMin(T_.primitiveField())              << " K"    << nl
        << "  max(T)       = "
        << gMax(T_.primitiveField())              << " K"    << nl
        << "  max(Q_pc)    = "
        << gMax(Q_pc_.primitiveField())           << " W/m3" << nl
        << "  min(Q_pc)    = "
        << gMin(Q_pc_.primitiveField())           << " W/m3" << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::thermalPhaseChangeModels::StefanEnergyJump::Q_pc() const
{
    return Q_pc_;
}


void Foam::thermalPhaseChangeModels::StefanEnergyJump::correct()
{
    calcQ_pc();
}


bool Foam::thermalPhaseChangeModels::StefanEnergyJump::read
(
    const dictionary& dict
)
{
    thermalPhaseChangeModel::read(dict);

    nSmoothIter_      = dict.lookupOrDefault<label>("nSmoothIter",      3);
    alphaSmoothWidth_ = dict.lookupOrDefault<scalar>("alphaSmoothWidth", 1.0);
    deltaT_act_       = dict.lookupOrDefault<scalar>("deltaT_act",       0.01);
    RelaxFac_         = dict.lookupOrDefault<scalar>("RelaxFac",         0.05);

    k_liq_ = dimensionedScalar("k_liq", dimPower/dimLength/dimTemperature, dict);
    k_vap_ = dimensionedScalar("k_vap", dimPower/dimLength/dimTemperature, dict);

    return true;
}


// ************************************************************************* //
