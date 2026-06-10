// Face conductivity (harmonic mean — unchanged)
surfaceScalarField kappaf = phaseChangePtr->kappaf();

const volScalarField    Q_pc_now = phaseChangePtr->Q_pc_thermal();  // sharp, interface-localized
const dimensionedScalar Tsat     = phaseChangePtr->T_sat();
const dimensionedScalar eps_T("eps_T", dimTemperature, 0.25);
const scalar implicitFactor = 1.0;   // =1 => implicit sink magnitude == |Q_pc| (exact latent heat)

// ── Implicit latent sink + interface Tsat-pin ──────────────────────────────
//   Q_pc_thermal is linear in (T-Tsat): qn = k*(T-Tsat)/psi, so the explicit
//   source -Q_pc_now = +(Ai*k/|psi|)*(T-Tsat) is a POSITIVE-feedback heating
//   term that runs T away from Tsat (the 526 K blow-up).
//
//   Acoeff = |Q_pc|/max(|T-Tsat|,eps_T) ≈ Ai*k/|psi| (the (T-Tsat) cancels),
//   Acoeff > 0 ⇒ unconditionally stable on the diagonal.  Net source
//       Acoeff*(Tsat - T) = -|Q_pc|*sign(T-Tsat)
//   removes exactly the latent heat |Q_pc| = |mdotRaw*h_lv| AND pins T->Tsat
//   only where Q_pc_thermal is sharp (the interface band), leaving the bulk
//   gradient free to drive qn.  See alphaEqn.H:278 — latent heat is handled
//   here in TEqn via Q_pc, not by the sensible-enthalpy flux.
const volScalarField Acoeff
(
    implicitFactor * mag(Q_pc_now) / max(mag(T - Tsat), eps_T)
);

// ── TEqn ───────────────────────────────────────────────────────────────────
// NO explicit -Q_pc_now: the implicit penalty IS the latent-heat handling.
// Adding -Q_pc_now alongside it either cancels it (net-zero source -> interface
// under-melts, "degraded") or, alone, gives the positive-feedback runaway.
fvScalarMatrix TEqn
(
    rhoCp*fvm::ddt(T)
  + fvm::div(rhoCpPhi, T)
  - fvm::laplacian(kappaf, T)
  + fvm::Sp(Acoeff, T)
 ==
    Acoeff*Tsat
);

TEqn.relax();
TEqn.solve();

// -------------------------------------------------------------------------
// Hardt & Wondra benchmark:
// Force pure liquid cells to saturation temperature while preserving
// interface cells.
// -------------------------------------------------------------------------

//const scalar TsatValue = Tsat.value();

//forAll(T, celli)
//{
//    if (alpha1[celli] > 0.99)
//    {
//        T[celli] = TsatValue;
//    }
//}

T.correctBoundaryConditions();

T.correctBoundaryConditions();

Info<< "min/max(T) = " << min(T).value() << " " << max(T).value() << endl;

Info<< "min/max(Q_pc) = "
    << gMin(Q_pc_now)
    << " "
    << gMax(Q_pc_now)
    << nl;
