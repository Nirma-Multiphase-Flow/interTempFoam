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

    nSmoothIter_     (dict.lookupOrDefault<label> ("nSmoothIter",      2)),
    alphaSmoothWidth_(dict.lookupOrDefault<scalar>("alphaSmoothWidth",  0.5)),
    lambdaSmearCells_(dict.lookupOrDefault<scalar>("lambdaSmearCells",  1.5)),
    liquidBiasCoeff_ (dict.lookupOrDefault<scalar>("liquidBiasCoeff",   1.0)),
    mdotMax_         (dict.lookupOrDefault<scalar>("mdotMax",           100.0)),
    RelaxFac_        (dict.lookupOrDefault<scalar>("RelaxFac",          1.0)),
    useEnthalpyCorrection_
        (dict.lookupOrDefault<Switch>("useEnthalpyCorrection", Switch(true))),

    k_liq_("k_liq", dimPower/dimLength/dimTemperature, dict),
    k_vap_("k_vap", dimPower/dimLength/dimTemperature, dict)
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
            alphaGeom * k_vap_
        + (scalar(1) - alphaGeom) * k_liq_,
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
        T_
    );


    // Mild localized smoothing ONLY for thermodynamic sensing
    //
    // Keep extremely weak:
    // - preserves physical gradients
    // - suppresses only grid-scale oscillations

    const scalar tempSmoothCoeff = 0.15;

    const dimensionedScalar D_Tsense
    (
        "D_Tsense",
        dimArea/dimTime,
        tempSmoothCoeff*h_ref*h_ref/dt
    );


    // Single pseudo-diffusion iteration ONLY
    TSense +=
        fvc::laplacian(D_Tsense, TSense)*dt_dim;

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

    volScalarField qnRaw
    (
        IOobject
        (
            "qnRaw",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),

        max
        (
            kEff*(gradT & nHat),

            dimensionedScalar
            (
                "zeroQn",
                dimPower/dimArea,
                0
            )
        )
    );
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


    qn_ = qnRaw*evapSwitch;

        

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

    mdotRaw_ = interfaceArea_ * qn_ / h_lv_;

    // Zero mdotRaw_ in cells that touch any physical boundary.
    // The Stefan condition applies at the fluid interface, not at contact lines
    // where the interface meets a wall.  Suppressing these cells prevents the
    // T artifact at wall/interface junctions that the alpha mask alone cannot
    // exclude (alpha is in [1e-3, 1-1e-3] at contact lines just like at the
    // bulk interface).
  
    mdotRaw_ *= interfaceMask;

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
        max
        (
            mdotRaw_,
            dimensionedScalar
            (
                "zero",
                mdotRaw_.dimensions(),
                0
            )
        )
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

    mdot_ =
        mdotLimiter
    *tanh(mdot_/mdotLimiter);


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

    //==============================================================
    // Convert volumetric source -> interfacial mass flux
    //==============================================================

    volScalarField mdotStefan
    (
        IOobject
        (
            "mdotStefan",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mdotRaw_
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
        mdotStefan
    /
        max
        (
            interfaceArea_,
            dimensionedScalar
            (
                "AiMin",
                interfaceArea_.dimensions(),
                SMALL
            )
        )
    );

    // Face interpolation
    const surfaceScalarField mdotInterfacialf
    (
        fvc::interpolate(mdotInterfacial)
    );

    const surfaceVectorField nHatf
    (
        fvc::interpolate(nHat)
    );

    const surfaceScalarField rhoInt
    (
        1.0
    /
        (
            fvc::interpolate(alphaGeom)/mixture_.rho1()
        + (1.0 - fvc::interpolate(alphaGeom))/mixture_.rho2()
        )
    );

    const surfaceScalarField Ustef
    (
        mdotInterfacialf / rhoInt
    );

    // Stefan volumetric face flux [m3/s]
    phiStefan_ =
    (
        Ustef
    *(nHatf & mesh_.Sf())
    );

    // Restrict to interface
    surfaceScalarField alphaIf
    (
        fvc::interpolate(alphaGeom)
    );

    // Compact support centered around alpha = 0.5
    surfaceScalarField interfaceMaskF
    (
        16.0*alphaIf*(scalar(1.0) - alphaIf)
    );

    // Sharpen strongly
    interfaceMaskF = sqr(interfaceMaskF);

    // Apply compact localization
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

    Q_pc_ = mdot_ * h_lv_;



    // ------------------------------------------------------------
    // Mild thermal redistribution ONLY for latent heat support.
    //
    // Broadens thermal sink slightly without broadening
    // interface recession.
    // ------------------------------------------------------------

    for (label i=0; i<1; ++i)
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
        Qcorr_ =
            dimensionedScalar("zero", dimensionSet(1,-1,-3,0,0,0,0), Zero);
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

    k_liq_ = dimensionedScalar("k_liq", dimPower/dimLength/dimTemperature, dict);
    k_vap_ = dimensionedScalar("k_vap", dimPower/dimLength/dimTemperature, dict);

    return true;
}


// ************************************************************************* //
