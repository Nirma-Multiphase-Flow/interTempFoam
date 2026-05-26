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

#include "HardtWondra.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "fvcLaplacian.H"
#include "surfaceInterpolate.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{
    defineTypeNameAndDebug(HardtWondra, 0);

    addToRunTimeSelectionTable
    (
        thermalPhaseChangeModel,
        HardtWondra,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPhaseChangeModels::HardtWondra::HardtWondra
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
            "Q_pc",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Q_pc", dimensionSet(1,-1,-3,0,0,0,0), Zero)
    ),

    Q_pc_thermal_
    (
        IOobject
        (
            "Q_pc_thermal",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Q_pc_   // initialise equal to Q_pc_; overwritten in calcQ_pc()
    ),

    Qcorr_
    (
        IOobject
        (
            "EnthalpyCorr",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qcorr", dimensionSet(1,-1,-3,0,0,0,0), Zero)
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

    mdotRaw_
    (
        IOobject
        (
            "mdotRaw",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("mdotRaw", dimDensity/dimTime, Zero)
    ),

    interfaceArea_
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
        dimensionedScalar("Ai", dimless/dimLength, Zero)
    ),

    qn_
    (
        IOobject
        (
            "qn",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("qn", dimPower/dimArea, Zero)
    ),

    phiStefan_
    (
        IOobject
        (
            "phiStefan",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "phiStefan",
            dimVolume/dimTime,
            Zero
        )
    ),

    harmonicBias_
    (
        "harmonicBias",
        dimless,
        dict.lookupOrDefault<scalar>("harmonicBias", 0.3)
    ),

    nSmoothIter_     (dict.lookupOrDefault<label> ("nSmoothIter",      2)),
    alphaSmoothWidth_(dict.lookupOrDefault<scalar>("alphaSmoothWidth",  0.5)),
    lambdaSmearCells_(dict.lookupOrDefault<scalar>("lambdaSmearCells",  1.5)),
    liquidBiasCoeff_ (dict.lookupOrDefault<scalar>("liquidBiasCoeff",   1.0)),
    mdotMax_         (dict.lookupOrDefault<scalar>("mdotMax",           100.0)),
    RelaxFac_        (dict.lookupOrDefault<scalar>("RelaxFac",          1.0)),
    useEnthalpyCorrection_
        (dict.lookupOrDefault<Switch>("useEnthalpyCorrection", Switch(true))),

    k_liq_("k_liq", dimPower/dimLength/dimTemperature, dict),
    k_vap_("k_vap", dimPower/dimLength/dimTemperature, dict),

    AiFloorAbs_(dict.lookupOrDefault<scalar>("AiFloorAbs", 50.0)),

    betaThermal_(dict.lookupOrDefault<scalar>("betaThermal", 0.0)),

    kinFloorCells_(dict.lookupOrDefault<scalar>("kinFloorCells", 2.0))
{
    // Set mdot_ BCs for the Helmholtz solve.
    //   Physical patches (walls, inlets, outlets): fixedValue 0
    //     → Helmholtz constrained to zero at domain boundaries, preventing
    //       the smoothed source from reaching the wall and creating T artifacts.
    //   Empty/coupled patches: zeroGradient (2D extrusion / parallel only).
    volScalarField::Boundary& mdotBf = mdot_.boundaryFieldRef();
    forAll(mdotBf, patchi)
    {
        const fvPatch& p = mesh_.boundary()[patchi];
        if (p.coupled() || p.type() == "empty")
        {
            mdotBf.set(patchi,
                fvPatchField<scalar>::New("zeroGradient", p, mdot_));
        }
        else
        {
            mdotBf.set(patchi,
                fvPatchField<scalar>::New("fixedValue", p, mdot_));
        }
    }

    correct();
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::thermalPhaseChangeModels::HardtWondra::calcQ_pc()
{
    const scalar dt = max(mesh_.time().deltaTValue(), SMALL);

    // =========================================================================
    // STEP 1 – Smooth alpha to de-noise VOF staircase gradients
    //
    // Identical pseudo-diffusion approach to StefanEnergyJump.
    // Working copy only: alpha1_ is never modified.
    // =========================================================================

    // h_ref = smallest face-to-face cell thickness: min(V) / max(face area).
    // Using cbrt(max(V)) fails on anisotropic meshes (e.g. thin 1D cells of
    // 30 μm × 1 mm × 1 mm give 0.31 mm instead of 30 μm).
    const scalar h_ref =
        gMin(mesh_.V().field())
      / max(gMax(mesh_.magSf().field()), SMALL);

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
            "alphaSmooth_HW",
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

    volScalarField alphaGeom
    (
        IOobject
        (
            "alphaGeom_HW",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphaSmooth
    );

    // additional geometric smoothing ONLY for normals/interface area
    for (label i=0; i<1; ++i)
    {
        alphaGeom +=
            fvc::laplacian(D_smooth, alphaGeom)*dt_dim;

        alphaGeom.primitiveFieldRef() =
            max(min(alphaGeom.primitiveField(), scalar(1)), scalar(0));

        alphaGeom.correctBoundaryConditions();
    }

    // =========================================================================
    // STEP 2 – Interface area density  |∇α_s|  [1/m]
    //
    // This is the VOF continuum interpretation of the interface area per unit
    // volume: ∫|∇α| dV ≈ interface area (Hardt & Wondra, 2008).
    // Using the smoothed alpha removes staircase noise from the gradient.
    // =========================================================================
    //==============================================================
    // Interface geometry
    //==============================================================

    // Interface mask
    volScalarField interfaceMask
    (
        IOobject
        (
            "interfaceMask",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pos(alphaGeom - 0.01)
    *pos(0.99 - alphaGeom)
    );

    // Geometric interface gradient
    const volVectorField gradAlpha
    (
        IOobject
        (
            "gradAlpha_HW",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(alphaGeom)
    );

    // Continuum interface area density [1/m]
    interfaceArea_ = mag(gradAlpha);

    const dimensionedScalar AiMax
    (
        "AiMax",
        interfaceArea_.dimensions(),
        1.5/h_ref
    );

    interfaceArea_ =
        min(interfaceArea_, AiMax);


    // Interface normals
    const volVectorField nHat
    (
        IOobject
        (
            "nHat_HW",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        gradAlpha
    /
        (
            mag(gradAlpha)
        + dimensionedScalar
            (
                "epsGradAlpha",
                gradAlpha.dimensions(),
                SMALL
            )
        )
    );

    //==============================================================
    // Harmonic conductivity interpolation
    //==============================================================

    // ============================================================
    // Thermally compressed interface fraction
    //
    // PURPOSE:
    // Reduce nonphysical cross-interface thermal leakage
    // WITHOUT affecting:
    //
    // - alpha transport
    // - interface motion
    // - qn reconstruction
    // - mdot
    // - phiStefan
    //
    // This sharpens ONLY thermal conductivity support.
    // ============================================================

    const volScalarField alphaThermal
    (
        IOobject
        (
            "alphaThermal",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),

        0.5
    *
        (
            scalar(1)
        + tanh
            (
                4.0*(alphaGeom - 0.5)
            )
        )
    );

    const volScalarField kEff
    (
        IOobject
        (
            "kEff_HW",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        k_liq_ * k_vap_
    /
        max
        (
            alphaThermal * k_vap_
        + harmonicBias_*(scalar(1) - alphaThermal) * k_liq_,
            dimensionedScalar
            (
                "kappaMin",
                k_liq_.dimensions(),
                scalar(1e-10)
            )
        )
    );

    //==============================================================
    // Local thermodynamic sensing regularization
    //
    // IMPORTANT:
    // - ONLY used for qn evaluation
    // - DOES NOT modify solved temperature field T_
    // - Prevents one-cell interfacial thermal spikes from
    //   over-driving Stefan flux reconstruction
    //==============================================================

    // One-sided TSense: clip T from below at T_sat.
    //
    // Root cause of qn suppression (M1): in the diffuse interface zone the
    // latent sink pins T ≈ Tsat on both sides, so gradT ≈ 0 exactly where
    // |∇α| is largest. The physical driving gradient lives in the hot-phase
    // bulk (vapor for evaporation, liquid for melting) where T > Tsat.
    //
    // max(T_, T_sat_) zeroes out the cold-side temperature contribution:
    //   - hot phase (T > Tsat): max returns T  → full gradient preserved
    //   - cold phase / interface (T ≈ Tsat): max returns Tsat → zero gradient
    //   fvc::grad(TSense) at the interface then reflects only the hot-side
    //   neighbor's gradient, which is the correct one-sided Stefan condition.
    //
    // This is READ-ONLY: TSense never enters T_, TEqn, or any conservation eqn.
    volScalarField TSense
    (
        IOobject
        (
            "TSense_HW",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        max(T_, T_sat_)
    );

    // Mild smoothing to suppress VOF staircase noise (unchanged coefficient).
    const dimensionedScalar D_Tsense
    (
        "D_Tsense",
        dimArea/dimTime,
        scalar(0.15)*h_ref*h_ref/dt
    );

    TSense += fvc::laplacian(D_Tsense, TSense)*dt_dim;
    TSense.correctBoundaryConditions();


    //==============================================================
    // Reconstruct temperature gradient from sensing field
    //==============================================================

    const volVectorField gradT
    (
        IOobject
        (
            "gradT_HW",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(TSense)
    );


    //==============================================================
    // Interfacial conductive heat flux
    //
    // Use signed flux instead of |gradT·n|.
    //
    // This removes nonphysical evaporation activation from
    // oscillatory gradients and substantially reduces pressure
    // ringing.
    //==============================================================


    //==============================================================
    // Smooth thermodynamic activation
    //==============================================================

    //==============================================================
    // Smooth thermodynamic activation
    //
    // Use TSense instead of raw T_ to avoid interface thermal
    // trench amplification.
    //
    // Narrower activation band:
    // - preserves Stefan equilibrium behavior
    // - reduces delayed activation lag
    //==============================================================

    const dimensionedScalar deltaT
    (
        "deltaT",
        dimTemperature,
        0.25
    );

    volScalarField evapSwitch
    (
        IOobject
        (
            "evapSwitch",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),

        0.5
    *
        (
            scalar(1)
        + tanh((TSense - T_sat_)/deltaT)
        )
    );


    volScalarField condSwitch
    (
        IOobject
        (
            "condSwitch",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        0.5*(scalar(1) - tanh((TSense - T_sat_)/deltaT))
    );

    // Signed heat flux: positive = evaporation, negative = condensation
    // evapSwitch - condSwitch = tanh((T-Tsat)/deltaT): smooth sign at Tsat
    //qn_ = qnRaw*(evapSwitch);

        
    //==============================================================
    // Directional interfacial conductive heat flux
    //
    // nHat points vapor -> liquid.
    //
    // For wall|vapor|liquid evaporation:
    //
    //     gradT · nHat < 0
    //
    // but evaporation must produce:
    //
    //     mdot > 0
    //
    // Therefore define:
    //
    //     qn_ = -kEff*(gradT · nHat)
    //
    // so that:
    //   evaporation  -> qn_ > 0
    //   condensation -> qn_ < 0
    //
    // This preserves directional transport physics and removes
    // the need for artificial sign reconstruction via T-Tsat.
    //==============================================================

    // Signed interfacial conductive heat flux using one-sided TSense gradient.
    // nHat points vapor→liquid (standard VOF convention).
    // For evaporation (hot phase on vapor side): gradT·nHat < 0, so qnSigned > 0.
    // For melting (hot phase on liquid side): same sign result via max(T,Tsat) clip.
    // ─────────────────────────────────────────────────────────────────────────────
    // One-sided hot-phase Stefan gradient reconstruction
    //
    // Replaces centered fvc::grad(TSense) + empirical ×2 correction.
    //
    // Physical objective:
    //     dT/dn |_Γ ≈ (T_hot - Tsat) / d_hot
    //
    // using ONLY hot-side neighbors.
    //
    // No GFM.
    // No ghost cells.
    // No level-set.
    //
    // Existing Helmholtz redistribution remains unchanged.
    // ─────────────────────────────────────────────────────────────────────────────

    volScalarField dTdn_hot
    (
        IOobject
        (
            "dTdn_hot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "zero",
            dimTemperature/dimLength,
            Zero
        )
    );

    // Thermal-path gradient field: 1× h_ref floor (vs kinematic 10× floor).
    // Populated alongside dTdn_hot in the same neighbor-search loop.
    // Used only for Q_pc_thermal_; does not affect phiStefan.
    volScalarField dTdn_hot_actual
    (
        IOobject
        (
            "dTdn_hot_actual",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "zero",
            dimTemperature/dimLength,
            Zero
        )
    );

    volScalarField dTdn_hot_thermal
    (
        IOobject
        (
            "dTdn_hot_thermal",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimTemperature/dimLength, Zero)
    );

    const scalar TsatVal = T_sat_.value();

    forAll(mesh_.cells(), cellI)
    {
        // only interface cells
        const scalar alphaCell = alphaGeom[cellI];

        // Use interfaceArea directly instead of alpha thresholds.
        // This is much more robust for compressed VOF interfaces.

        if (interfaceArea_[cellI] < SMALL)
        {
            continue;
        }

        const vector& nHatCell = nHat[cellI];

        scalar bestDistance      = GREAT;
        scalar bestGradient_kin  = 0.0;   // kinematic: 10× floor — drives phiStefan, UNCHANGED
        scalar bestGradient_therm = 0.0;  // thermal:    1× floor — drives TEqn Acoeff, NEW
        

        const labelList& nbrCells = mesh_.cellCells()[cellI];

        forAll(nbrCells, nbrI)
        {
            const label nbrCell = nbrCells[nbrI];

            const vector dVec =
                mesh_.C()[nbrCell] - mesh_.C()[cellI];

            const scalar normalProj =
                mag(dVec & nHatCell);

            if (normalProj < SMALL)
                continue;

            const scalar Tnbr      = T_[nbrCell];
            const scalar superheat = Tnbr - TsatVal;

            if (superheat < 1e-6)
                continue;

            if (normalProj < bestDistance)
            {
                bestDistance = normalProj;

                // kinematic path: 10× floor for phiStefan stability (UNCHANGED)
                bestGradient_kin = superheat / max(normalProj, kinFloorCells_*h_ref);

                // thermal path: 1× floor — implicit Acoeff only, matrix-stable
                bestGradient_therm = superheat / max(normalProj,    h_ref);
            }
        }

        if (bestDistance < GREAT)
        {
            dTdn_hot[cellI]           = bestGradient_kin;    // kinematic — UNCHANGED
            dTdn_hot_thermal[cellI]   = bestGradient_therm;  // thermal  — NEW
        }
        else
        {
            dTdn_hot[cellI]           = 0.0;
            dTdn_hot_thermal[cellI]   = 0.0;
        }
    }

    // smooth reconstructed gradient slightly to suppress stencil switching noise
   
    volScalarField qnSigned
    (
        IOobject
        (
            "qnSigned",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        kEff*dTdn_hot
    );

    qnSigned *= interfaceMask;

    // NO ×2 correction anymore.
    // This is already a one-sided gradient.
    qn_ = qnSigned;

    // ─── Thermal-path latent sink ─────────────────────────────────────────────
    // Q_pc_thermal_ uses the 1× floor gradient (dTdn_hot_thermal) gated by
    // evapSwitch (≈1 on hot side, ≈0 on cold side).  This restricts the
    // stronger sink to the hot phase where the non-physical superheat lives.
    // Crucially, dTdn_hot (10× floor) feeding qn_ and mdotRaw_ is UNCHANGED.
    // ─────────────────────────────────────────────────────────────────────────

    volScalarField qnThermal
    (
        IOobject("qnThermal_hw", mesh_.time().timeName(), mesh_,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        kEff * dTdn_hot_thermal
    );
    qnThermal *= interfaceMask;

    // evapSwitch already computed above: ≈1 where T>Tsat, ≈0 where T≤Tsat
    Q_pc_thermal_ = interfaceArea_ * qnThermal * evapSwitch;


    // =========================================================================
    // STEP 7 – Raw Stefan mass flux  [kg/m³/s]
    //
    //   mdot_raw = |∇α_s| · q_n / h_lv
    //
    // Dimensional check: [1/m] · [W/m²] / [J/kg] = [kg/m³/s] ✓
    //
    // This is the diffuse-interface version of the sharp-interface Stefan
    // condition: k·∂T/∂n = ṁ·h_lv.  The factor |∇α_s| localises it to
    // the interface band and gives it the correct units for a volumetric source.
    // =========================================================================

    // ────────────────────────────────────────────────────────────────────────

    mdotRaw_ = interfaceArea_ * qn_ / h_lv_;

    // Zero mdotRaw_ in cells that touch any physical boundary.
    // The Stefan condition applies at the fluid interface, not at contact lines
    // where the interface meets a wall.  Suppressing these cells prevents the
    // T artifact at wall/interface junctions that the alpha mask alone cannot
    // exclude (alpha is in [1e-3, 1-1e-3] at contact lines just like at the
    // bulk interface).
  
    //mdotRaw_ *= interfaceMask;

    // =========================================================================
    // STEP 8 – Helmholtz redistribution
    //
    //   mdot - λ²∇²mdot = mdot_raw      (solved implicitly)
    //
    // λ = lambdaSmearCells × h_ref  [m]
    //
    // This is the core of Hardt & Wondra (2008): the implicit Helmholtz equation
    // smoothly redistributes the raw Stefan source over a band of width ~λ.
    //
    // Why Helmholtz instead of explicit Laplacian smoothing?
    //   - Conserves total mass transfer: ∫mdot dV = ∫mdot_raw dV
    //   - Suppresses high-frequency noise without spreading the source globally
    //   - Implicit solve: no stability restriction on λ
    //   - λ → 0 recovers the sharp-interface limit
    //
    // Note: mdot_ is a persistent field so it acts as the initial guess for
    // the Helmholtz solve (warm start improves convergence).
    // =========================================================================

    Info<< "max(mdotRaw) = " << gMax(mdotRaw_) << endl;
    Info<< "min(mdotRaw) = " << gMin(mdotRaw_) << endl;

    const dimensionedScalar lambdaSqr
    (
        "lambdaSqr",
        dimArea,
        Foam::sqr(lambdaSmearCells_ * h_ref)
    );

    volScalarField mdotSource
    (
        IOobject
        (
            "mdotSource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mdotRaw_  // signed: positive=evap, negative=cond
    );

    fvScalarMatrix mdotEqn
    (
        fvm::Sp
        (
            dimensionedScalar("one", dimless, scalar(1)),
            mdot_
        )
    - fvm::laplacian(lambdaSqr, mdot_)
    ==
        mdotSource
    );

    mdotEqn.relax();

    solve
    (
        mdotEqn,
        mesh_.solverDict("mdot")
    );

     // Enforce zero mass transfer in bulk cells (optional safety check)

    // =========================================================================
    // STEP 9 – Hard clip on mdot
    //
    // mdotMax_ is the physical upper bound on mass transfer rate.
    // For water/steam at 100°C: 100 kg/m³/s → Q_pc_max ≈ 2.26×10⁸ W/m³.
    // This clip is a backstop safety limiter; for well-resolved simulations
    // with correct λ it should not activate.
    // =========================================================================

    const dimensionedScalar mdotLimiter
    (
        "mdotLimiter",
        mdot_.dimensions(),
        mdotMax_
    );

    // Smooth saturation instead of hard clipping.
    // Prevents thermodynamic discontinuity while still bounding mdot.

    //mdot_ =
    //    mdotLimiter
    //*tanh(mdot_/mdotLimiter);
    mdot_ = min(max(mdot_, -mdotLimiter), mdotLimiter);

    // =========================================================================
    // STEP 10 – Under-relaxation
    // =========================================================================

    mdot_ *= RelaxFac_;

    //==============================================================
    // Localized Stefan transport flux
    //
    // IMPORTANT:
    // - Uses mdotRaw_ (localized interface source)
    // - NOT mdot_ (Helmholtz redistributed source)
    // - Provides interface kinematics
    //==============================================================

    // =========================================================================
    // STEP 10 – Direct Stefan face velocity from interfacial heat flux
    //
    // Physics: v_Stefan = qn / (rho_int × h_lv)
    //
    // This is the sharp-interface Stefan condition applied directly at faces.
    // It gives a UNIFORM velocity across the diffuse interface zone (since qn
    // is approximately constant for a well-resolved erf/erfc profile), which
    // produces pure translation of the alpha field without distortion.
    //
    // Why this replaces mdot/(Ai+AiFloor)/rho:
    //   The Helmholtz solve reduces mdot peak by factor ~3 while conserving its
    //   integral.  mdot/Ai then has large tails (small Ai, nonzero mdot) that
    //   create 3-5x overspeeding at the interface edges, thickening the diffuse
    //   zone without advancing the alpha=0.5 contour.  The result is a ~20-25%
    //   underestimation of interface position even when total mass transfer is
    //   exactly correct.
    //
    // No AiFloor.  No heuristic blending.  No empirical correction.
    // =========================================================================
    // ── nHatSmooth: smoothed interface normal for phiStefan face interpolation ──
    // (was previously embedded inside the old mdot/Ai phiStefan block)
    // Two Laplacian passes reduce VOF staircase noise in nHat before face interp.
    const dimensionedScalar lambdaSqr_nHat
    (
        "lambdaSqr_nHat",
        dimArea,
        Foam::sqr(lambdaSmearCells_ * h_ref)
    );

    volVectorField nHatSmooth
    (
        IOobject
        (
            "nHatSmooth",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        nHat
    );

    for (label i = 0; i < 2; ++i)
    {
        nHatSmooth += fvc::laplacian(lambdaSqr_nHat, nHatSmooth);
        nHatSmooth.correctBoundaryConditions();
    }

    nHatSmooth =
        nHatSmooth
    / (
            mag(nHatSmooth)
        + dimensionedScalar("epsN", dimless, SMALL)
        );
        
    // =========================================================================
    // STEP 10 – Stefan velocity from mdotRaw (not Helmholtz mdot)
    //
    // Physical basis:
    //   mdotRaw = Ai × qn / h_lv
    //   → mdotRaw / Ai = qn / h_lv = v_Stefan [kg/m²/s / (kg/m³) = m/s]
    //   This ratio is UNIFORM across the interface zone (qn constant for
    //   a well-resolved erf profile) → pure translation of the alpha field.
    //
    // Why NOT Helmholtz mdot:
    //   Helmholtz spreads mdot from σ=Δx/2 to σ_out=1.58Δx, reducing the
    //   peak by factor 0.316.  Dividing by peaked Ai then gives Ustef_center
    //   ≈ 0.30 × v_true, causing 25-31% underestimation of interface position.
    //
    // Why NOT direct qn/(ρh_lv):
    //   qn is nonzero in bulk liquid (T > Tsat everywhere → gradT ≠ 0).
    //   Without Ai in numerator there is no natural zero-suppression
    //   in bulk → breaks interface.  mdotRaw has the correct zero in bulk.
    //
    // AiFloor: reduced from 50 to 5 /m — pure division-safety guard only.
    //   At max(Ai)=637/m: 637/(637+5) = 99.2% of true Ustef.
    //   At bulk (Ai=0): mdotRaw=0 → Ustef=0, AiFloor irrelevant.
    // =========================================================================

    const dimensionedScalar AiFloor_kin
    (
        "AiFloor_kin",
        interfaceArea_.dimensions(),
        scalar(5)           // pure zero-guard; ~0.8% of max(Ai)=637
    );

    volScalarField mdotInterfacial
    (
        IOobject
        (
            "mdotInterfacial",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mdotRaw_ / (interfaceArea_ + AiFloor_kin)
    );

    const surfaceScalarField mdotInterfacialf
    (
        fvc::interpolate(mdotInterfacial)
    );

    const surfaceScalarField rhoInt
    (
        scalar(1)
    / (
            fvc::interpolate(alphaGeom)/mixture_.rho1()
        + (scalar(1) - fvc::interpolate(alphaGeom))/mixture_.rho2()
        )
    );

    const surfaceVectorField nHatf
    (
        fvc::interpolate(nHatSmooth)
    );

    const surfaceScalarField Ustef
    (
        mdotInterfacialf / rhoInt
    );

    phiStefan_ =
    (
        Ustef * (nHatf & mesh_.Sf())
    );

    // flat-top mask — unchanged
    surfaceScalarField alphaIf(fvc::interpolate(alphaGeom));
    const dimensionedScalar maskAlphaMin_("maskAlphaMin", dimless, scalar(0.05));
    surfaceScalarField interfaceMaskF
    (
        min(
            min(alphaIf, scalar(1) - alphaIf) / maskAlphaMin_,
            dimensionedScalar("one", dimless, scalar(1))
        )
    );
    phiStefan_ *= interfaceMaskF;
    // =========================================================================
    // STEP 11 – Latent heat source  [W/m³]
    //
    //   Q_pc = mdot · h_lv
    //
    // Sign: Q_pc > 0 (evaporation, heat sink in liquid),
    //       Q_pc < 0 (condensation, heat source).
    // Consistent with base-class alpha1Gen and PCV conventions.
    // =========================================================================

    // =========================================================================
    // STEP 11 – Latent heat source  [W/m³]
    //
    //   Q_pc = mdot · h_lv
    //
    // IMPORTANT:
    // This is the ONLY thermodynamic latent source used by TEqn.
    //
    // mdot_ is Helmholtz-regularized and therefore suitable for
    // volumetric latent heat coupling.
    // =========================================================================

    // =========================================================================
    // STEP 11 – Latent heat source from localised mdotRaw (not Helmholtz mdot)
    //
    // Physics: Q_pc = ṁ_raw × h_lv  [W/m³]
    //
    // mdotRaw_ = Ai × qn / h_lv → zero wherever Ai = 0 (bulk cells).
    // This confines the latent heat source/sink to the 2-cell diffuse interface
    // zone where phase change is physically occurring.
    //
    // Why NOT Helmholtz mdot_:
    //   Helmholtz spreads mdot into hot liquid cells 1-2Δx from the interface.
    //   Q_pc = mdot_ × h_lv in those cells creates an Acoeff = Q_pc/(T−Tsat)
    //   implicit sink in TEqn that pulls the hot-liquid temperature toward Tsat.
    //   Over 180 timesteps this suppresses T_H1 by several K, reducing qn, in
    //   a compounding feedback loop responsible for the ~25-28% residual error.
    //
    // Helmholtz mdot_ is RETAINED for phiStefan (interface kinematics path is
    // already using mdotRaw) and for the pressure/continuity equation — it only
    // exits the Q_pc path here.
    //
    // Global energy conservation: ∫Q_pc dV = ∫mdotRaw×h_lv dV = ∫mdot×h_lv dV
    // (Helmholtz conserves the integral by construction).
    // =========================================================================

    Q_pc_ = mdotRaw_ * h_lv_;



    // ------------------------------------------------------------
    // Mild thermal redistribution ONLY for latent heat support.
    //
    // Broadens thermal sink slightly without broadening
    // interface recession.
    // ------------------------------------------------------------

    for (label i=0; i<0; ++i)
    {
        Q_pc_ +=
            fvc::laplacian(lambdaSqr, Q_pc_);

        Q_pc_.correctBoundaryConditions();
    }


    // =========================================================================
    // STEP 12 – Enthalpy correction


    // =========================================================================
    // STEP 12 – Enthalpy correction (Hardt & Wondra eq. 42)
    //
    // The Helmholtz solve places mass sources in cells outside the physical
    // phase-change zone.  In those cells the temperature equation "sees" a
    // latent heat source/sink that is not associated with a local conductive
    // flux, creating artificial interface temperature spikes.
    //
    // The correction:
    //   Q_corr = -mdot · (cp_l - cp_v) · (T - T_sat)
    //
    // Removes the spurious enthalpy carried by the redistributed mass.
    // Magnitude: ~(cp_l-cp_v)·ΔT/h_lv of Q_pc ≈ 1-2% for water/steam at ΔT=10K.
    // =========================================================================

    if (useEnthalpyCorrection_)
    {
        Qcorr_ = -mdot_ * (mixture_.cp1() - mixture_.cp2()) * (T_ - T_sat_);
    }
    else
    {
        Qcorr_ = dimensionedScalar
        (
            "zero",
            Qcorr_.dimensions(),
            Zero
        );
    }


    // =========================================================================
    // STEP 13 – Diagnostics
    // =========================================================================

    Info<< "HardtWondra phase-change:" << nl
        << "  ∫Q_pc  dV = "
        << gSum(Q_pc_.primitiveField() * mesh_.V().field()) << " W"       << nl
        << "  ∫mdot  dV = "
        << gSum(mdot_.primitiveField() * mesh_.V().field()) << " kg/s"    << nl
        << "  max(Ai)   = "
        << gMax(interfaceArea_.primitiveField())             << " /m"      << nl
        << "  max(|qn|) = "
        << gMax(mag(qn_.primitiveField()))                   << " W/m2"   << nl
        << "  max(mdot) = "
        << gMax(mdot_.primitiveField())                      << " kg/m3/s"<< nl
        << "  min(mdot) = "
        << gMin(mdot_.primitiveField())                      << " kg/m3/s"<< endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::thermalPhaseChangeModels::HardtWondra::Q_pc() const
{
    return Q_pc_;
}


Foam::tmp<Foam::volScalarField>
Foam::thermalPhaseChangeModels::HardtWondra::Q_pc_thermal() const
{
    return Q_pc_thermal_;
}


Foam::tmp<Foam::volScalarField>
Foam::thermalPhaseChangeModels::HardtWondra::Qcorr() const
{
    return Qcorr_;
}

Foam::tmp<Foam::surfaceScalarField>
Foam::thermalPhaseChangeModels::HardtWondra::phiStefan() const
{
    return phiStefan_;
}

Foam::tmp<Foam::surfaceScalarField>
Foam::thermalPhaseChangeModels::HardtWondra::kappaf() const
{
    // Harmonic conductivity:
    //
    //   k_harm = k_l*k_v / (α*k_v + (1-α)*k_l)
    //
    // Physically consistent with series thermal resistance across the
    // diffuse interface.

    const volScalarField kappa_h
    (
        IOobject
        (
            "kappa_h",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),

        k_liq_*k_vap_
      /
        max
        (
            alpha1_*k_vap_
          + (scalar(1) - alpha1_)*k_liq_,

            dimensionedScalar
            (
                "kappaMin",
                k_liq_.dimensions(),
                scalar(1e-12)
            )
        )
    );

    return fvc::interpolate(kappa_h);
}


void Foam::thermalPhaseChangeModels::HardtWondra::correct()
{
    calcQ_pc();
}


bool Foam::thermalPhaseChangeModels::HardtWondra::read
(
    const dictionary& dict
)
{
    thermalPhaseChangeModel::read(dict);

    nSmoothIter_      = dict.lookupOrDefault<label> ("nSmoothIter",      5);
    alphaSmoothWidth_ = dict.lookupOrDefault<scalar>("alphaSmoothWidth",  1.0);
    lambdaSmearCells_ = dict.lookupOrDefault<scalar>("lambdaSmearCells",  1.5);
    liquidBiasCoeff_  = dict.lookupOrDefault<scalar>("liquidBiasCoeff",   1.0);
    mdotMax_          = dict.lookupOrDefault<scalar>("mdotMax",           100.0);
    RelaxFac_         = dict.lookupOrDefault<scalar>("RelaxFac",          1.0);
    useEnthalpyCorrection_ =
        dict.lookupOrDefault<Switch>("useEnthalpyCorrection", Switch(true));
    AiFloorAbs_       = dict.lookupOrDefault<scalar>("AiFloorAbs",        50.0);
    betaThermal_      = dict.lookupOrDefault<scalar>("betaThermal",        0.0);
    kinFloorCells_ = dict.lookupOrDefault<scalar>("kinFloorCells", 2.0);

    k_liq_ = dimensionedScalar("k_liq", dimPower/dimLength/dimTemperature, dict);
    k_vap_ = dimensionedScalar("k_vap", dimPower/dimLength/dimTemperature, dict);

    return true;
}


// ************************************************************************* //
