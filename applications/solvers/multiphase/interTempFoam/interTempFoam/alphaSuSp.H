// Phase-change sources for MULES — bounded implicit-sink / explicit-source split
// [1/s] = volumetric rate [kg/m³/s] / density [kg/m³]
//
// Sign convention (verified against Sheikh 1D melting, general for all modes):
//   alpha1Gen < 0  →  phase-1 is consumed  (melting: solid decreases)
//   alpha1Gen > 0  →  phase-1 is produced  (solidification: solid increases)
//
// Boundedness (THE FIX):
//   Previously alpha1Gen was applied as pure explicit Su with Sp=0.  In cells where
//   alpha1 ≈ 0 (hot-melt side) the predictor drove alpha.solid negative and
//   alpha1.max(0) discarded ~99% of the intended sink.  Fix: move the consumed-phase
//   sink into the implicit Sp channel so the contribution scales as Sp×alpha1 → 0
//   as alpha1 → 0, preventing the negative excursion entirely.
//
//   Consumed-phase (alpha1Gen < 0):
//     Sp = alpha1Gen / max(alpha1, alphaFloor)   [1/s, negative]
//     fvm::Sp(Sp, alpha1) → rate ∝ alpha1, vanishes as cell empties, no clip.
//
//   Produced-phase (alpha1Gen > 0):
//     Su = alpha1Gen   [1/s, explicit source]
//     Receiving cell has alpha1 ≈ 0, so it has headroom; no clip occurs.
//
//   alphaFloor = 0.01 caps |Sp| ≤ 100×|alpha1Gen|, preventing stiffness.
//
// Generality:
//   - Equal density (phiStefan=0, PCV=0): all interface motion rides on Sp/Su.
//   - Unequal density: phiStefan advection + PCV dilatation paths are unchanged;
//     Sp/Su handle the residual fraction conversion.
//   - Mode-agnostic: works for melting, solidification, evaporation, condensation;
//     assignment follows sign(alpha1Gen), not a mode-specific branch.

// divU is retained in the predictor's implicit matrix (fvm::Sp(Sp+divU,.))
// purely for numerical regularization -- it matches the canonical interFoam
// VoF predictor structure and keeps the U-p-alpha coupling well-conditioned
// for this solver's extreme density ratio (rho_l/rho_v ~ 1626).  Its bias
// (divU = PCV is ~1626x larger in magnitude than alpha1Gen = mdot/rho_l) is
// removed via deferred correction in alphaEqn.H: the explicit RHS subtracts
// divU*alpha1.oldTime() to cancel the implicit divU*alpha1 contribution, so
// the net source scales with divU*(alpha1-alpha1.oldTime()) [O(dt), small]
// instead of divU*alpha1 [O(1), the ~35-38x growth bug from an earlier
// attempt that added divU with no compensation at all]. A second earlier
// attempt removed divU entirely (pure Sp, no regularization) -- that is
// exactly mass-conservative in continuous form but destabilized the
// U-p-alpha coupling from the very first timestep (spurious 5-15 m/s
// velocities), sweeping the droplet away via convection. Deferred
// correction keeps both properties: implicit divU for stability, near-zero
// net bias for accuracy.
tmp<volScalarField> tDivU(fvc::div(phi));
const volScalarField& divU = tDivU();

tmp<volScalarField> tAlpha1Gen(phaseChangePtr->alpha1Gen());

const scalarField& alpha1GenI = tAlpha1Gen().primitiveField();
const scalarField& alpha1I    = alpha1.primitiveField();

const scalar alphaFloor = 0.01;  // prevents Sp stiffness when alpha1 → 0

volScalarField::Internal Su
(
    IOobject
    (
        "Su",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Su", dimless/dimTime, Zero)
);

volScalarField::Internal Sp
(
    IOobject
    (
        "Sp",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Sp", dimless/dimTime, Zero)
);

forAll(Su, celli)
{
    const scalar src = alpha1GenI[celli];
    if (src < 0.0)
    {
        // Consumed-phase sink: implicit to prevent alpha1 going negative
        Sp[celli] = src / max(alpha1I[celli], alphaFloor);
    }
    else if (src > 0.0)
    {
        // Produced-phase source: explicit, receiving cell has headroom
        Su[celli] = src;
    }
}

// --- Boundedness-loss compensation -----------------------------------------
// The redistributed alpha1Gen spreads symmetrically across the diffuse
// interface. The portion landing on cells already depleted of the consumed
// phase (alpha.solid -> 0 for a melt sink) cannot be applied: the bounded
// update Sp*alpha1 -> 0 there, and that mass-conversion credit evaporates.
// Consequence: APPLIED phase change is only ~0.48x the heat-flux rate,
// stalling the front at ~0.65x analytical despite the thermal path being
// correct. This block rescales Su and Sp so the APPLIED integral matches
// the REQUESTED integral, restoring the analytical front speed.
// Sign-agnostic => general for melting, solidification, evaporation,
// condensation and both equal / unequal density configurations.
{
    const scalarField& Vc = mesh.V().field();

    // Requested: heat-flux-determined total source
    const scalar reqInt = gSum(alpha1GenI*Vc);

    // Applied: what Su/Sp actually deliver given the current alpha1 field.
    // Explicit Su contributes directly; implicit Sp contributes Sp*alpha1.
    scalar appLoc = 0.0;
    forAll(Su, celli)
    {
        appLoc += (Sp[celli]*alpha1I[celli] + Su[celli])*Vc[celli];
    }
    const scalar appInt = returnReduce(appLoc, sumOp<scalar>());

    // Geff >= 1: boost to compensate MULES boundedness loss.
    // Clamp to [1, 2] — never inverts the source; caps boost so numerical noise
    // in reqInt/appInt (observed as high as 40x ratio) is not amplified catastrophically.
    // Gmax lowered 5→2: sufficient for the Stefan-validated planar front (Geff ≈ 1–1.5
    // in practice) while preventing ×5 amplification of near-zero noise sources.
    const scalar Gmax = 2.0;
    scalar Geff = 1.0;
    if (mag(appInt) > VSMALL && mag(reqInt) > VSMALL)
    {
        Geff = min(max(reqInt/appInt, scalar(1)), Gmax);
    }

    forAll(Su, celli)
    {
        Su[celli] *= Geff;
        Sp[celli] *= Geff;
    }

    Info<< "alphaSuSp: Geff=" << Geff
        << "  reqInt=" << reqInt << "  appInt=" << appInt << endl;
}

