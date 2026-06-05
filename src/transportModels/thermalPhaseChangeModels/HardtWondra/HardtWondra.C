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

Description
    Hardt & Wondra phase-change model, gradient path driven by the genuine
    reconstructed distance function (RDF). The "RDF" field is built upstream
    by a plicRDF reconstruction + reconstructedDistanceFunction in the solver
    and looked up from the object registry here. The mass / phiStefan machinery
    is unchanged.

\*---------------------------------------------------------------------------*/

#include "HardtWondra.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "fvcLaplacian.H"
#include "surfaceInterpolate.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "fvcSurfaceIntegrate.H"
#include "upwind.H"
#include "fvCFD.H"

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
        dimensionedScalar("phiStefan", dimVolume/dimTime, Zero)
    ),

    UStefan_
    (
        IOobject
        (
            "UStefan",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("UStefan", dimVelocity, Zero)
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

    kinFloorCells_(dict.lookupOrDefault<scalar>("kinFloorCells", 2.0)),

    nExtrapIter_(dict.lookupOrDefault<label>("nExtrapIter", 5)),

    extrapCFL_(dict.lookupOrDefault<scalar>("extrapCFL", 0.2)),

    alphaPureLiquid_(dict.lookupOrDefault<scalar>("alphaPureLiquid", 0.99)),

    alphaPureVapor_(dict.lookupOrDefault<scalar>("alphaPureVapor", 0.01)),

    maxGradT_(dict.lookupOrDefault<scalar>("maxGradT", 1e8)),

    // Retained for .H compatibility; UNUSED on the genuine-RDF path (the
    // distance now comes from the reconstructedDistanceFunction library, not
    // from an analytic alpha inversion). Safe to delete from .H if you wish.
    interfaceWidthCells_(dict.lookupOrDefault<scalar>("interfaceWidthCells", 1.0)),

    coldPhaseIsHighAlpha1_
        (dict.lookupOrDefault<bool>("coldPhaseIsHighAlpha1", false))
{
    // Set mdot_ BCs for the Helmholtz solve (unchanged).
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
    // -------------------------------------------------------------------------
    // RDF availability guard.
    //
    // "RDF" is registered & filled upstream in the solver:
    //     surf.reconstruct();                         // plicRDF segment data
    //     RDF.markCellsNearSurf(interfaceCells, 4);   // ring level >= band
    //     RDF.constructRDF(RDF.nextToInterface(),
    //                      surf.centre(), surf.normal(), exchangeFields);
    // called AFTER alphaEqnSubCycle/mixture.correct() and BEFORE this correct().
    //
    // On the very first (construction-time) correct() the field may not be
    // registered yet; rather than abort, emit a warning and produce zero phase
    // change for that single call (it is recomputed with a live RDF every
    // outer iteration thereafter).
    // -------------------------------------------------------------------------
    if (!mesh_.foundObject<volScalarField>("RDF"))
    {
        WarningInFunction
            << "Field 'RDF' not in registry. Build the "
            << "reconstructedDistanceFunction and call constructRDF() before "
            << "phaseChangePtr->correct(). Producing zero phase change this call."
            << endl;

        Q_pc_         = dimensionedScalar(Q_pc_.dimensions(),         Zero);
        Q_pc_thermal_ = dimensionedScalar(Q_pc_thermal_.dimensions(), Zero);
        mdot_         = dimensionedScalar(mdot_.dimensions(),         Zero);
        mdotRaw_      = dimensionedScalar(mdotRaw_.dimensions(),      Zero);
        qn_           = dimensionedScalar(qn_.dimensions(),           Zero);
        phiStefan_    = dimensionedScalar(phiStefan_.dimensions(),    Zero);
        return;
    }

    const scalar dt = max(mesh_.time().deltaTValue(), SMALL);

    // =========================================================================
    // STEP 1 – Smooth alpha to de-noise VOF staircase gradients (UNCHANGED)
    // =========================================================================
    const scalar h_ref =
        gMin(mesh_.V().field())
      / max(gMax(mesh_.magSf().field()), SMALL);

    const scalar Fo_per_iter = min
    (
        alphaSmoothWidth_*alphaSmoothWidth_/max(scalar(nSmoothIter_), scalar(1)),
        scalar(0.25)
    );

    const dimensionedScalar D_smooth
    (
        "D_smooth", dimArea/dimTime, Fo_per_iter*h_ref*h_ref/dt
    );
    const dimensionedScalar dt_dim("dt_dim", dimTime, dt);

    volScalarField alphaSmooth
    (
        IOobject
        (
            "alphaSmooth_HW", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        alpha1_
    );

    for (label i = 0; i < nSmoothIter_; ++i)
    {
        alphaSmooth += fvc::laplacian(D_smooth, alphaSmooth)*dt_dim;
        alphaSmooth.correctBoundaryConditions();
    }
    alphaSmooth.primitiveFieldRef() =
        max(min(alphaSmooth.primitiveField(), scalar(1)), scalar(0));
    alphaSmooth.correctBoundaryConditions();


    // =========================================================================
    // STEP 2 – Interface normal nHat & area density Ai (UNCHANGED)
    //          Still required by mdotRaw_ (Step 8) and phiStefan_ (Step 9).
    // =========================================================================
    const volVectorField gradAlpha
    (
        IOobject
        (
            "gradAlpha_HW", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, false
        ),
        fvc::grad(alphaSmooth)
    );

    volVectorField gradAlphaSmooth
    (
        IOobject
        (
            "nHatSmooth_HW", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, false
        ),
        mesh_,
        dimensionedVector("zero", gradAlpha.dimensions(), Zero),
        "zeroGradient"
    );
    gradAlphaSmooth.primitiveFieldRef() = gradAlpha.primitiveField();
    gradAlphaSmooth.correctBoundaryConditions();

    {
        const dimensionedScalar lambdaSqr_nHat
        (
            "lambdaSqr_nHat", dimArea, Foam::sqr(lambdaSmearCells_*h_ref)
        );
        fvVectorMatrix gaEqn
        (
            fvm::Sp(scalar(1), gradAlphaSmooth)
          - fvm::laplacian(lambdaSqr_nHat, gradAlphaSmooth)
         == gradAlpha
        );
        gaEqn.solve();
    }

    interfaceArea_ = mag(gradAlphaSmooth);
    const dimensionedScalar AiMax("AiMax", interfaceArea_.dimensions(), 1.5/h_ref);
    interfaceArea_ = min(interfaceArea_, AiMax);

    const volVectorField nHat
    (
        IOobject
        (
            "nHat_HW", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, false
        ),
        gradAlphaSmooth
      / (interfaceArea_ + dimensionedScalar("epsN", dimless/dimLength, SMALL))
    );


    // =========================================================================
    // STEP 3 – Interface band (UNCHANGED) — used by Step 8 localisation.
    // =========================================================================
    volScalarField interfaceBand
    (
        IOobject
        (
            "interfaceBand", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, false
        ),
        pos(interfaceArea_ - dimensionedScalar("epsArea", dimless/dimLength, 1e-4))
    );

    const label nBandExpand = 1;
    for (label i = 0; i < nBandExpand; ++i)
    {
        interfaceBand =
            max(interfaceBand, fvc::average(fvc::interpolate(interfaceBand)));
    }
    interfaceBand = pos(interfaceBand - 0.01);


    // =========================================================================
    // STEP 4 – Signed normal distance from the genuine RDF library
    //
    // "RDF" carries the magnitude of the reconstructed distance to the PLIC
    // segment (>0 only inside the markCellsNearSurf ring band; 0 outside).
    // Re-sign it by phase so it is +ve in liquid (alpha>=0.5) and -ve in
    // vapor/solid (alpha<0.5). Using alpha for the sign makes the gradient
    // immune to the orientation convention of the reconstruction normal.
    // =========================================================================
    const volScalarField& RDFfield = mesh_.lookupObject<volScalarField>("RDF");

    volScalarField psi
    (
        IOobject
        (
            "RDF_HW", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("psi", dimLength, Zero),
        "zeroGradient"
    );
    {
        scalarField&       p = psi.primitiveFieldRef();
        const scalarField& r = RDFfield.primitiveField();
        const scalarField& a = alpha1_.primitiveField();
        forAll(p, celli)
        {
            const scalar sgn = (a[celli] >= 0.5 ? scalar(1) : scalar(-1));
            p[celli] = sgn*mag(r[celli]);
        }
        psi.correctBoundaryConditions();
    }


    // =========================================================================
    // STEP 5 – Local RDF normal gradient & heat flux qn_
    //
    //   dT/dn ~= (T - Tsat)/psi          (Scheufler Eqn 14 with exact distance)
    //   qDot_i = -k_i * dT/dn_i
    //   qn_    =  qDot_liq - qDot_vap    (sign IDENTICAL to the previous code,
    //                                     so Steps 8-11 / phiStefan / TEqn
    //                                     Acoeff are untouched)
    //
    // Gate on the RDF narrow band (|RDF| > SMALL): outside the band RDF is
    // exactly 0, where (T-Tsat)/psi would otherwise hit the floor and fabricate
    // a spurious gradient. psi-floor guards a cell centre lying on the segment.
    // =========================================================================
    const scalar psiFloor = 0.5*h_ref;

    volScalarField dTdn_liq
    (
        IOobject
        (
            "dTdn_liq", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("g", dimTemperature/dimLength, Zero),
        "zeroGradient"
    );
    volScalarField dTdn_vap
    (
        IOobject
        (
            "dTdn_vap", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("g", dimTemperature/dimLength, Zero),
        "zeroGradient"
    );

    {
        scalarField&       gl  = dTdn_liq.primitiveFieldRef();
        scalarField&       gv  = dTdn_vap.primitiveFieldRef();
        const scalarField& a   = alpha1_.primitiveField();
        const scalarField& Tc  = T_.primitiveField();
        const scalarField& p   = psi.primitiveField();
        const scalar       Ts  = T_sat_.value();

        forAll(gl, celli)
        {
            gl[celli] = 0.0;
            gv[celli] = 0.0;

            // RDF narrow band only
            if (mag(p[celli]) <= SMALL) continue;

            // One-sided: only the hot phase drives phase change.
            // coldPhaseIsHighAlpha1 = true  (melting):    hot = melt,   a < 0.5
            // coldPhaseIsHighAlpha1 = false (evaporation): hot = liquid, a >= 0.5
            // Skipping the cold-phase cells prevents double-counting when the
            // cold phase (solid / vapour) warms above Tsat over time.
            const bool cellIsHot = coldPhaseIsHighAlpha1_
                ? (a[celli] < 0.5)
                : (a[celli] >= 0.5);

            if (!cellIsHot) continue;

            const scalar sgn     = (p[celli] >= 0 ? scalar(1) : scalar(-1));
            const scalar psiSafe = sgn*max(mag(p[celli]), psiFloor);
            const scalar g = Foam::clamp
            (
                (Tc[celli] - Ts)/psiSafe,
                scalar(-maxGradT_),
                scalar(maxGradT_)
            );

            if (a[celli] >= 0.5) { gl[celli] = g; }   // liquid side, psi > 0
            else                 { gv[celli] = g; }   // vapor/solid side, psi < 0
        }
        dTdn_liq.correctBoundaryConditions();
        dTdn_vap.correctBoundaryConditions();
    }

    // qn_ = qDot_liq - qDot_vap = -k_liq*dTdn_liq + k_vap*dTdn_vap
    qn_ = (-k_liq_*dTdn_liq) - (-k_vap_*dTdn_vap);

    // Mask qn_ to the interface band BEFORE averaging.
    //
    // Why: qn_ is computed for every cell where |psi| > SMALL (the RDF band).
    // Even with the alpha-threshold fix in the solver, a few marginal cells
    // outside the true interface band can carry non-zero psi and a wrong
    // (T-Tsat)/psi value (e.g., cells near the hot wall with T >> Tsat and
    // a small stale psi).  The averaging pass then spreads these wrong values
    // to the true interface cells, diluting the actual heat flux.
    // Applying the interfaceBand mask first ensures that only cells with a
    // genuine |∇alpha| signal contribute to the post-averaged qn_.
    //
    // Physical justification: the interface heat flux qn_ has meaning only
    // in cells that lie on the diffuse interface; zeroing it elsewhere is
    // exact for a sharp interface and physically consistent for a diffuse one.
    // This fix is general (valid for evaporation, condensation, melting,
    // solidification) because interfaceBand is always derived from |∇alpha|.
    qn_ *= interfaceBand;

    // Single averaging pass — spreads qn_ by ~1 cell toward phase bulk.
    // With the band mask above, only genuine interface values are spread.
    for (label i = 0; i < 1; ++i)
    {
        qn_ = 0.5*qn_ + 0.5*fvc::average(fvc::interpolate(qn_));
    }


    // =========================================================================
    // STEP 8 – Mass flux & Helmholtz redistribution (UNCHANGED)
    // =========================================================================
    mdotRaw_ = interfaceArea_ * interfaceBand * qn_ / h_lv_;
    mdotRaw_.primitiveFieldRef() =
        max(min(mdotRaw_.primitiveField(), scalar(mdotMax_)), scalar(-mdotMax_));

    const dimensionedScalar lambdaSqr
    (
        "lambdaSqr", dimArea, Foam::sqr(lambdaSmearCells_*h_ref)
    );

    fvScalarMatrix mdotEqn
    (
        fvm::Sp(scalar(1), mdot_) - fvm::laplacian(lambdaSqr, mdot_) == mdotRaw_
    );
    mdotEqn.solve();


    // =========================================================================
    // STEP 9 – Stefan Velocity / phiStefan  (UNCHANGED — PRESERVED VERBATIM)
    // =========================================================================
    const dimensionedScalar AiFloor
    (
        "AiFloor", interfaceArea_.dimensions(), AiFloorAbs_
    );

    volScalarField mdotStefanVol
    (
        IOobject
        (
            "mdotStefanVol_HW", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, false
        ),
        // Use mdotRaw_ (sharp, band-localised) not mdot_ (Helmholtz-smoothed).
        // Helmholtz spreads mdot_ beyond interfaceBand; re-applying band then
        // discards that leaked portion, so phiStefan carries only ~27% of the
        // generated flux while TEqn (via Q_pc_thermal_=mdotRaw_*h_lv_) removes
        // 100% — a ~3.6x energy leak that stalls the front at ~0.49x analytical.
        // mdotRaw_*band/max(Ai,AiFloor) = qn_/h_lv_ at the interface (Ai>>AiFloor,
        // band=1), giving Ustef = qn/(h_lv*rho) — the exact Stefan velocity, and
        // now energy-consistent with the latent-heat sink in TEqn.H.
        (mdotRaw_ * interfaceBand) / max(interfaceArea_, AiFloor)
    );

    surfaceScalarField alphaIf = fvc::interpolate(alphaSmooth);
    surfaceScalarField rhoInt_f
    (
        alphaIf*mixture_.rho1() + (scalar(1) - alphaIf)*mixture_.rho2()
    );

    // phiN = interpolate(nHat) & Sf
    surfaceScalarField phiN(fvc::interpolate(nHat) & mesh_.Sf());

    surfaceScalarField Ustef = fvc::interpolate(mdotStefanVol) / rhoInt_f;
    phiStefan_ = Ustef * phiN;

    const dimensionedScalar maskAlphaMin("maskAlphaMin", dimless, scalar(0.05));
    surfaceScalarField interfaceMaskF
    (
        min
        (
            min(alphaIf, scalar(1) - alphaIf)/maskAlphaMin,
            dimensionedScalar("one", dimless, scalar(1))
        )
    );
    phiStefan_ *= interfaceMaskF;


    // =========================================================================
    // STEP 10 – Source Field Assignments (UNCHANGED)
    // =========================================================================
    Q_pc_         = mdot_    * h_lv_;   // SMOOTH: PCV fallback & pEqn
    Q_pc_thermal_ = mdotRaw_ * h_lv_;   // SHARP : TEqn Acoeff

    // Reconstruct cell-centred Stefan velocity vector from the face flux.
    // fvc::reconstruct divides phiStefan_ [m³/s] by face areas to recover
    // velocity [m/s], written to every time directory for ParaView.
    UStefan_ = fvc::reconstruct(phiStefan_);


    // =========================================================================
    // STEP 11 – Diagnostics
    // =========================================================================
    label nBand = 0;
    {
        const scalarField& p = psi.primitiveField();
        forAll(p, celli) { if (mag(p[celli]) > SMALL) ++nBand; }
    }

    Info<< "HardtWondra phase-change (RDF library gradient):" << nl
        << "  RDF band cells        = "
        << returnReduce(nBand, sumOp<label>()) << nl
        << "  Sum(Q_pc smooth)*dV   = "
        << gSum(Q_pc_.primitiveField()*mesh_.V().field()) << " W" << nl
        << "  Sum(Q_pc_thermal)*dV  = "
        << gSum(Q_pc_thermal_.primitiveField()*mesh_.V().field()) << " W" << nl
        << "  Sum(mdotRaw)*dV       = "
        << gSum(mdotRaw_.primitiveField()*mesh_.V().field()) << " kg/s" << nl
        << "  Sum(mdot smooth)*dV   = "
        << gSum(mdot_.primitiveField()*mesh_.V().field()) << " kg/s" << nl
        << "  max(Ai)               = "
        << gMax(interfaceArea_.primitiveField()) << " /m" << nl
        << "  max(|RDF psi|)        = "
        << gMax(mag(psi.primitiveField())) << " m" << nl
        << "  max(|qn|)             = "
        << gMax(mag(qn_.primitiveField())) << " W/m2" << nl
        << "  max(|dTdn_liq|)       = "
        << gMax(mag(dTdn_liq.primitiveField())) << " K/m" << nl
        << "  max(|dTdn_vap|)       = "
        << gMax(mag(dTdn_vap.primitiveField())) << " K/m" << nl
        << "  max(|phiStefan|)      = "
        << gMax(mag(phiStefan_.primitiveField())) << nl
        << "  sum(|phiStefan|)      = "
        << gSum(mag(phiStefan_.primitiveField())) << nl;
    Info<< "min/max(qn) = "
    << gMin(qn_)
    << " "
    << gMax(qn_)
    << endl;
    Info<< "min/max(mdotRaw) = "
    << gMin(mdotRaw_)
    << " "
    << gMax(mdotRaw_)
    << endl;
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
    // Harmonic conductivity (UNCHANGED):
    //   k_harm = k_l*k_v / (alpha*k_v + (1-alpha)*k_l)
    const volScalarField kappa_h
    (
        IOobject
        (
            "kappa_h", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE
        ),
        k_liq_*k_vap_
      / max
        (
            alpha1_*k_vap_ + (scalar(1) - alpha1_)*k_liq_,
            dimensionedScalar("kappaMin", k_liq_.dimensions(), scalar(1e-12))
        )
    );

    return fvc::interpolate(kappa_h);
}


Foam::tmp<Foam::volScalarField>
Foam::thermalPhaseChangeModels::HardtWondra::PCV() const
{
    // UNCHANGED: smeared mdot_ keeps pressure-velocity coupling stable.
    return tmp<volScalarField>::New
    (
        IOobject
        (
            "PCV_HW", mesh_.time().timeName(), mesh_,
            IOobject::NO_READ, IOobject::NO_WRITE, false
        ),
        mdot_ * (1.0/mixture_.rho2() - 1.0/mixture_.rho1())
    );
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

    kinFloorCells_ = dict.lookupOrDefault<scalar>("kinFloorCells", 2.0);

    k_liq_ = dimensionedScalar("k_liq", dimPower/dimLength/dimTemperature, dict);
    k_vap_ = dimensionedScalar("k_vap", dimPower/dimLength/dimTemperature, dict);

    interfaceWidthCells_ =
        dict.lookupOrDefault<scalar>("interfaceWidthCells", 1.0);

    coldPhaseIsHighAlpha1_ =
        dict.lookupOrDefault<bool>("coldPhaseIsHighAlpha1", false);

    return true;
}


// ************************************************************************* //