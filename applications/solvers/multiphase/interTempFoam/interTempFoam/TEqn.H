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




// ── TEqn ─────────────────────────────────────────────────────────────────
fvScalarMatrix TEqn
(
    rhoCp * fvm::ddt(T)
  + fvm::div(rhoCpPhi, T)
  - fvm::laplacian(kappaf, T)
  + fvm::Sp(Acoeff, T)
 ==
    (Acoeff*Tsat)

);

TEqn.relax();
TEqn.solve();
T.correctBoundaryConditions();

Info<< "min/max(T) = " << min(T).value() << " " << max(T).value() << endl;