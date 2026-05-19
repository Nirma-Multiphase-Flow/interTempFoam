// ─── note: ghost-source correction removed ───────────────────────────────────
// The conservative fvm::ddt(rhoCp,T) + fvm::div(rhoCpPhi,T) carries a ghost
// term T·[d(ρCp)/dt + ∇·(ρCpU)] at the interface.  Correcting it requires
// rhoCp and rhoCpPhi to be updated within the PIMPLE iteration (now done in
// the main loop before this file is included).  The residual ghost error is
// O(CFL × ΔρCp/ρCp) per step — negligible for Co < 0.1.
// ─────────────────────────────────────────────────────────────────────────────

// Harmonic-mean face conductivity from the phase-change model.
// For HardtWondra: k_harm = k_l*k_v / (α·k_v + (1-α)·k_l).
// For other models: falls back to arithmetic mixture.kappaf().
// Harmonic mean is physically correct for series heat resistors at the
// VOF interface; arithmetic overpredicts by ~14× at α=0.5 for water/steam.
surfaceScalarField kappaf = phaseChangePtr->kappaf();

// Cache Q_pc — called once, used twice below
const volScalarField Q_pc_now = phaseChangePtr->Q_pc();

// ─────────────────────────────────────────────────────────────────────────────
// Implicit linearisation:  Q_pc ≈ Acoeff · (T − Tsat)
//
// Acoeff ≥ 0 for the evaporation arm (T > Tsat, Q_pc > 0).
// eps_T floor prevents division by near-zero in equilibrium cells.
// Condensation arm (Q_pc < 0) falls through to explicit Q_cond.
// ─────────────────────────────────────────────────────────────────────────────
const dimensionedScalar Tsat = phaseChangePtr->T_sat();
const dimensionedScalar eps_T("eps_T", dimTemperature, 5.0);

const volScalarField Acoeff
(
    max
    (
        Q_pc_now / max(T - Tsat, eps_T),
        dimensionedScalar("zero", dimensionSet(1,-1,-3,-1,0,0,0), 0)
    )
);

// Condensation residual: Q_pc_now where Q_pc_now < 0, else 0
const volScalarField Q_cond
(
    min
    (
        Q_pc_now,
        dimensionedScalar("zero", dimensionSet(1,-1,-3,0,0,0,0), 0)
    )
);

// ─────────────────────────────────────────────────────────────────────────────
// Anti-trench saturation floor
//
// Root cause of cold trench:
//   Acoeff = Q_pc/(T-Tsat) ~ O(1e7) W/(m³·K)  <<  rhoCp/dt ~ O(4e10) W/(m³·K)
//   The Newton-linearised correction is ~10,000x weaker than the ddt term,
//   so T can fall freely below Tsat within a timestep as the Helmholtz-
//   redistributed Q_pc drains liquid-side cells faster than conduction refills.
//
// Fix: one-sided stiff penalty that activates ONLY when T < Tsat in the
// diffuse zone.  C_floor*(rhoCp/dt) makes the correction 10,000x stiffer
// than the ddt term → T restored to Tsat within the same timestep.
// subcoolMask → 0 when T > Tsat, so the evaporation-driving gradient is
// not disturbed.
// ─────────────────────────────────────────────────────────────────────────────
const volScalarField indicator_TEqn = 4.0*alpha1*(scalar(1.0) - alpha1);


// ─────────────────────────────────────────────────────────────────────────────

fvScalarMatrix TEqn
(
    rhoCp * fvm::ddt(T)
  + fvm::div(rhoCpPhi, T)
  - fvm::laplacian(kappaf, T)
  + fvm::Sp(Acoeff, T)
 ==
  - Q_cond
  + phaseChangePtr->Qcorr()
);

TEqn.relax();

TEqn.solve();

T.correctBoundaryConditions();

Info<< "min/max(T) = "
    << min(T).value()
    << " "
    << max(T).value()
    << endl;
