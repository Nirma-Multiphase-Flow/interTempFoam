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
//
// Sign convention (verified for Sheikh 1D melting case):
//   psi > 0 in solid (alpha1=alpha.solid >= 0.5), psi < 0 in liquid.
//   Melting hot cell is liquid (psi<0) → qn < 0 → mdot < 0 → alpha1Gen < 0.
//   alpha.solid must decrease → Sp < 0 is correct. NO sign flip needed.

tmp<volScalarField> tDivU(phaseChangePtr->PCV());
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

volScalarField::Internal divU
(
    IOobject
    (
        "divU",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimless/dimTime,
    tDivU().primitiveField()
);
