// Face conductivity (harmonic mean — unchanged)
surfaceScalarField kappaf = phaseChangePtr->kappaf();

// ── Latent source fields ──────────────────────────────────────────────────
// Q_pc_now     : kinematic path (10× floor) — used for Qcorr, Q_cond
// Q_pc_th_now  : thermal path  ( 1× floor) — used for Acoeff only
const volScalarField Q_pc_now    = phaseChangePtr->Q_pc();
const volScalarField Q_pc_th_now = phaseChangePtr->Q_pc_thermal();

const dimensionedScalar Tsat  = phaseChangePtr->T_sat();
const dimensionedScalar eps_T("eps_T", dimTemperature, 2.5);
const scalar implicitFactor = 0.5;

// ── Main latent Acoeff (thermal path — ~10× stronger than before) ─────────
const volScalarField Acoeff
(
    implicitFactor
  * max
    (
        Q_pc_th_now / max(T - Tsat, eps_T),
        dimensionedScalar("zero", dimensionSet(1,-1,-3,-1,0,0,0), 0)
    )
);

// Condensation explicit source (still uses kinematic Q_pc — unchanged)
const volScalarField Q_cond
(
    min(Q_pc_now, dimensionedScalar("zero", dimensionSet(1,-1,-3,0,0,0,0), 0))
);

// ── Diffuse-interface indicator ───────────────────────────────────────────
const volScalarField indicator_TEqn = 4.0*alpha1*(scalar(1.0) - alpha1);

// ── Anti-trench floor (cold side) — UNCHANGED ────────────────────────────
// Prevents T dropping below Tsat in diffuse zone.
const dimensionedScalar C_floor("C_floor", dimless, scalar(5.0));
const volScalarField Acoeff_floor
(
    C_floor
  * (rhoCp / mesh.time().deltaT())
  * indicator_TEqn
  * pos(Tsat - T)
);

// ── Anti-superheat floor (hot side) — COLD-PHASE GATE ────────────────────
// Prevents non-physical T > Tsat in the solid side of the diffuse zone.
// Gate: alpha1 > 0.5 = solid side (for melting, coldPhaseIsHighAlpha1=true).
//       Excludes liquid-side cells (alpha1 < 0.5) where T > Tsat is the
//       physical driving gradient for qn → no suppression of interface speed.
// coldPhaseIsHighAlpha1 read from phaseChangeProperties in createFields.H.
// REPLACE:
// const volScalarField coldSideGate
// (
//     coldPhaseIsHighAlpha1
//     ? pos(alpha1 - dimensionedScalar("half", dimless, 0.5))
//     : pos(dimensionedScalar("half", dimless, 0.5) - alpha1)
// );

// WITH:
volScalarField coldSideGate
(
    IOobject("coldSideGate", mesh.time().timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("one", dimless, scalar(1.0))
);
if (coldPhaseIsHighAlpha1)
{
    coldSideGate = pos(alpha1 - dimensionedScalar("half", dimless, scalar(0.5)));
}
else
{
    coldSideGate = pos(dimensionedScalar("half", dimless, scalar(0.5)) - alpha1);
}

const dimensionedScalar C_hot("C_hot", dimless, scalar(5.0));
const volScalarField Acoeff_hot_floor
(
    C_hot
  * (rhoCp / mesh.time().deltaT())
  * indicator_TEqn
  * pos(T - Tsat)
  * coldSideGate           // ← restricts to cold-phase side only
);

// ── TEqn ─────────────────────────────────────────────────────────────────
fvScalarMatrix TEqn
(
    rhoCp * fvm::ddt(T)
  + fvm::div(rhoCpPhi, T)
  - fvm::laplacian(kappaf, T)
  + fvm::Sp(Acoeff + Acoeff_floor + Acoeff_hot_floor, T)
 ==
    (Acoeff*Tsat + Acoeff_floor*Tsat + Acoeff_hot_floor*Tsat)
  - Q_cond
  + phaseChangePtr->Qcorr()
);

TEqn.relax();
TEqn.solve();
T.correctBoundaryConditions();

Info<< "min/max(T) = " << min(T).value() << " " << max(T).value() << endl;