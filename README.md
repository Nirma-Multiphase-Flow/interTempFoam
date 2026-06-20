# interTempFoam

A non-isothermal, phase-change-capable two-phase VOF solver for **OpenFOAM v2412**, derived from
`interFoam`. It adds a coupled energy (temperature) equation and a runtime-selectable thermal
phase-change model (evaporation / condensation / melting / solidification) built on geometric-VOF
(PLIC) interface reconstruction and a reconstructed distance function (RDF).

> **Status: active research code.** Core two-phase + energy transport is solid and validated against
> 1D analytical benchmarks. The phase-change source terms (Hardt–Wondra model) work for planar Stefan
> melting but the droplet d²-law evaporation case is **unfinished, ongoing work** as of the current
> HEAD. See [Status: Working vs. Broken / In-Progress](#status-working-vs-broken--in-progress) before
> relying on this for production results.

---

## Overview

`interTempFoam` solves two incompressible, immiscible fluids with VOF interface capturing, an energy
equation, and (optionally) a thermal phase-change source that converts mass between phases at the
interface based on the local conductive heat flux. Application areas:

- Stefan-problem melting / solidification
- Thermal stratification in two-phase systems
- Conjugate heat transfer across fluid interfaces
- Thermocapillary (Marangoni) flow, via temperature-dependent surface tension
- Droplet evaporation (d²-law) — **in progress, not yet validated**

---

## Mathematical Formulation

All source terms below are *as implemented in code*, not idealized textbook forms — file/line anchors
are given so behavior can be checked against the source directly.

### 1. VOF phase-fraction advection — `alphaEqn.H`, `alphaSuSp.H`

```
∂α₁/∂t + ∇·(α₁U) + ∇·[α₁(1−α₁)Uᵣ] = Su + Sp·α₁
```

Solved with MULES (optionally sub-cycled, `alphaEqnSubCycle.H`). `Uᵣ = φc·n̂` is the interface
compression flux. `Su`/`Sp` are the phase-change mass-conversion source, built in `alphaSuSp.H` from
`phaseChangePtr->alpha1Gen()` (units 1/s):

- **Consumed phase** (`alpha1Gen < 0`): applied **implicitly**, `Sp = alpha1Gen / max(α₁, 0.01)`, so the
  sink vanishes as the cell empties (prevents driving α₁ negative).
- **Produced phase** (`alpha1Gen > 0`): applied **explicitly**, `Su = alpha1Gen`.
- A **boundedness-loss compensation** factor `Geff = clamp(reqInt/appInt, 1, Gmax=2)` rescales `Su`/`Sp`
  so the volume-integrated *applied* source matches the heat-flux-determined *requested* source (MULES
  boundedness silently drops part of the source otherwise).
- The predictor matrix in `alphaEqn.H` carries `fvm::Sp(Sp + divU, α₁) − divU·α₁.oldTime()` — `divU` is
  kept in the implicit channel for numerical stability under the solver's large density ratios, with its
  bias removed via deferred correction (net effect O(Δt), not O(1)).

### 2. Continuity / pressure — `pEqn.H`

```
∇·U = PCV
```

```cpp
fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA) - phaseChangePtr->PCV()
```

`PCV = −ṅ·(1/ρ_vapor − 1/ρ_liquid)`, the volumetric dilatation source from mass crossing the interface
at different specific volumes. Zero for equal-density phases or when phase change is disabled.

### 3. Momentum — `UEqn.H`

```
∂(ρU)/∂t + ∇·(ρ_φU U) = −∇p_rgh − (g·x)∇ρ + σκn̂δ + ∇·τ
```
Standard CSF surface-tension body force (`mixture.surfaceTensionForce()`); no explicit recoil/Stefan
momentum source — phase change only enters momentum indirectly through the pressure dilatation term.

### 4. Energy / temperature — `TEqn.H`

```
∂(ρCp·T)/∂t + ∇·(ρCpφ·T) − ∇·(k_f∇T) = A·(T_sat − T)
```

```cpp
fvm::ddt(rhoCp, T)
+ fvm::div(rhoCpPhi, T)
- fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoCpPhi), T)   // boundedness compensation, variable rhoCp
- fvm::laplacian(kappaf, T)
+ fvm::Sp(Acoeff, T)
==
Acoeff*Tsat
```

`Acoeff = |Q_pc|/max(|T − Tsat|, 0.25 K)` is an **implicit latent-heat sink** that simultaneously (a)
removes exactly the latent heat `|Q_pc|` and (b) pins T → T_sat in the interface band. This deliberately
*replaces* an earlier explicit `−Q_pc` source, which was a positive-feedback term that diverged (the
code comments call this the "526 K blow-up"). **Latent heat is not carried by the sensible-enthalpy flux
`rhoCpPhi`** — see next section.

### 5. Property closures

- `ρ = α₁ρ₁ + (1−α₁)ρ₂`, `ρCp = α₁ρ₁Cp₁ + (1−α₁)ρ₂Cp₂` (volume-weighted, rebuilt every PIMPLE outer
  iteration after `pEqn.H`).
- `rhoCpPhi = fvc::interpolate(ρCp)·phi` — **uses the fluid flux `phi` only**, deliberately *not*
  `phiTotal` (which includes the Stefan flux). Mixing the Stefan flux into the sensible-enthalpy flux
  would inject `ρCp_liquid·T` into low-`ρCp` vapor cells and spike T; latent heat travels via the
  `Acoeff`/`Q_pc` path in `TEqn.H` instead.
- `rhoPhi = α₁ρ₁Φ_interp + φ·ρ₂` (momentum mass flux) **does** include the phase-change contribution.
- `κ = ρ·ν·cp/Pr` per phase (computed internally — `cp`/`Pr` are user inputs, conductivity is derived,
  not specified directly), `kappaf` = harmonic-mean face conductivity.
- `rhoCpPhi` and `rhoPhi` are both time-averaged across alpha sub-cycles in `alphaEqnSubCycle.H`
  (`rhoPhiSum`, `rhoCpPhiSum`), keeping the energy flux thermodynamically consistent with VOF advection.

### 6. Phase-change source — Hardt–Wondra diffuse-interface model

Implemented in `src/transportModels/thermalPhaseChangeModels/HardtWondra/`. Pipeline (per
`HardtWondra::calcQ_pc()`):

```
α (smoothed) → ∇α (Helmholtz-smoothed) → interfaceArea, interfaceBand
RDF → signed normal distance ψ (re-signed by phase)
dT/dn ≈ (T − T_sat)/ψ  (one-sided, hot phase only)
qₙ = −k_liq·dT/dn|_liq + k_vap·dT/dn|_vap
ṅ_raw = interfaceArea · band · qₙ / h_lv
Helmholtz redistribution:  ṅ − λ²∇²ṅ = ṅ_raw     (two λ: lambdaSmearCells for Q_pc/PCV,
                                                    lambdaAlphaCells for alpha1Gen)
Q_pc = ṅ · h_lv            → drives Acoeff in TEqn.H, and PCV in pEqn.H
alpha1Gen = Geff·(ṅ_alpha·band)/ρ₁  → drives Su/Sp in alphaSuSp.H
phiStefan = ṅ''·(1/ρ_vapor − 1/ρ_liquid)·n̂   (Shaikh 2016 Eq. 7, n̂ from ∇ψ)
```

### 7. PIMPLE coupling order (per outer iteration, `interTempFoam.C`)

```
1. alphaEqnSubCycle.H        — VOF advection (+ rhoPhiSum/rhoCpPhiSum accumulation)
2. surf.reconstruct()        — PLIC interface reconstruction
3. RDF.constructRDF(...)     — reconstructed distance function (+ bootstrap fallback if empty)
4. phaseChangePtr->correct() — ONLY on pimple.firstIter() (see Limitations)
5. phiTotal = phi ± phiStefan
6. UEqn.H, pEqn.H            — momentum + pressure (PCV source)
7. rebuild rho, rhoCp from current alpha1
8. TEqn.H                    — energy (Acoeff latent-heat sink)
9. turbulence->correct()
```
The `coldPhaseIsHighAlpha1` flag selects sign convention: `true` = melting (solid = α₁=1 = cold phase),
`false` = evaporation (vapor = α₁=0 = cold phase).

---

## Repository Structure

```
interTempFoam/
├── applications/solvers/multiphase/interTempFoam/interTempFoam/
│   ├── interTempFoam.C          # Main solver — PIMPLE loop, RDF, phase-change call
│   ├── createFields.H           # Field setup: U, p_rgh, T, mixture, RDF, phaseChangeProperties
│   ├── TEqn.H                   # Energy equation (Acoeff latent-heat sink)
│   ├── UEqn.H                   # Momentum predictor
│   ├── pEqn.H                   # Pressure (p_rgh) equation, PCV dilatation source
│   ├── alphaEqn.H               # VOF MULES solve, rhoPhi/rhoCpPhi rebuild
│   ├── alphaEqnSubCycle.H       # Sub-cycled alpha solve, rhoPhiSum/rhoCpPhiSum time-averaging
│   ├── alphaSuSp.H              # Phase-change Su/Sp source split for MULES + Geff rescaling
│   ├── alphaCourantNo.H         # Alpha-based CFL check
│   ├── createAlphaFluxes.H      # Alpha flux field setup
│   ├── correctPhi.H / initCorrectPhi.H   # Flux correction (moving mesh)
│   ├── rhofs.H                  # Face-interpolated densities
│   ├── setDeltaT.H / setRDeltaT.H        # Adaptive time-stepping / LTS
│   ├── docs/                    # Design notes (e.g. thermal-leakage-fix.md)
│   ├── report.md / report.pdf   # Fuller derivation & validation write-up
│   ├── Make/{files,options}
│   ├── interMixingFoam/         # Three-phase variant — NOT the primary solver
│   └── overInterDyMFoam/        # Overset-mesh variant — NOT the primary solver
│
└── src/transportModels/
    ├── Allwmake                                       # Builds all 8 libraries, in order
    ├── twoPhaseMixture/                                # Base two-phase mixture (stock OF)
    ├── interfaceProperties/                            # Surface tension + contact angle
    │   └── surfaceTensionModels/temperatureDependent/  # σ(T) — Marangoni flows
    ├── twoPhaseProperties/                             # alphaContactAngle BCs
    ├── incompressible/                                 # ★ extended with cp, Pr, kappa
    │   └── incompressibleTwoPhaseMixture/
    ├── compressible/                                   # Base (stock OF, unmodified)
    ├── immiscibleIncompressibleTwoPhaseMixture/        # VOF mixture used by the solver
    ├── geometricVoF/                                   # PLIC reconstruction, RDF, zoneDistribute
    └── thermalPhaseChangeModels/                        # ★ custom phase-change library
        ├── thermalPhaseChangeModel/                     # Abstract base + RTS factory
        ├── noPhaseChange/                                # Default / disabled
        ├── InterfaceEquilibrium/                         # Forces interface cells to equilibrium
        ├── StefanEnergyJump/                              # Q_pc from interfacial conductive flux
        ├── StefanEquilibriumMerged/                        # (T-Tsat)/dt driving force, α(1-α)-localized
        └── HardtWondra/                                    # ★ active model — see Math Formulation §6
```

---

## Installation & Compilation

### Prerequisites

- **OpenFOAM v2412** (openfoam.com release). Other versions will require porting.
- `WM_PROJECT_DIR`, `FOAM_USER_APPBIN`, `FOAM_USER_LIBBIN` set (i.e. OpenFOAM environment sourced).
- C++14-compatible compiler.

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
```

> **Portability caveat:** `applications/.../interTempFoam/Make/options` hardcodes
> `-I/usr/lib/openfoam/openfoam2412/src/surfMesh/lnInclude`. If your OpenFOAM install lives elsewhere,
> edit this line before building.

### Step 1 — Build the transport model libraries

```bash
cd src/transportModels
./Allwmake
```

This runs, **in this exact order** (per the actual `Allwmake` script — 8 libraries):
1. `twoPhaseMixture`
2. `interfaceProperties`
3. `twoPhaseProperties`
4. `incompressible` (adds `cp`, `Pr`, `kappa` accessors)
5. `compressible`
6. `immiscibleIncompressibleTwoPhaseMixture`
7. `geometricVoF`
8. `thermalPhaseChangeModels`

Installed to `$FOAM_USER_LIBBIN`.

### Step 2 — Build the solver

```bash
cd ../../applications/solvers/multiphase/interTempFoam/interTempFoam
wmake
```

Installed to `$FOAM_USER_APPBIN` as `interTempFoam`.

### Clean rebuild

```bash
cd src/transportModels && ./Allwmake -clean && ./Allwmake
cd ../../applications/solvers/multiphase/interTempFoam/interTempFoam && wclean && wmake
```

### Verify

```bash
which interTempFoam
interTempFoam -help
```

---

## Case Structure (User-Side Files)

> **No tutorial/example cases are bundled in this repository.** The layout below is what the solver
> *requires* at runtime (derived from `createFields.H`), not a ready-to-run tutorial.

### `0/` — initial fields

| Field | Type | Read mode | Notes |
|---|---|---|---|
| `U` | `volVectorField` | `MUST_READ` | Velocity [m/s] |
| `p_rgh` | `volScalarField` | `MUST_READ` | `p − ρgh` [Pa] |
| `T` | `volScalarField` | `MUST_READ` | Temperature [K] |
| `alpha.<phase1>` | `volScalarField` | via mixture | Phase fraction [0,1] |

### `constant/transportProperties`

```cpp
phases (water air);

water
{
    transportModel  Newtonian;
    nu              1e-06;     // kinematic viscosity [m^2/s]
    rho             1000;      // density [kg/m^3]
    cp              4182;      // specific heat [J/(kg.K)]
    Pr              7.0;       // Prandtl number [-]
}

air
{
    transportModel  Newtonian;
    nu              1.48e-05;
    rho             1.0;
    cp              1005;
    Pr              0.71;
}

sigma           0.07;          // surface tension [N/m]
```

Thermal conductivity is derived internally: `κ = ρ·ν·cp/Pr` (see Math Formulation §5).

### `constant/phaseChangeProperties` (optional)

If absent, the solver falls back to `noPhaseChange` (no mass/energy source). To enable Hardt–Wondra:

```cpp
model           HardtWondra;

T_sat           373.15;        // [K] required
h_lv            2.26e6;        // [J/kg] required
k_liq           0.679;         // [W/(m.K)] required
k_vap           0.026;         // [W/(m.K)] required

DilatationSource        true;  // enable PCV in pEqn
PhaseFractionSource     true;  // enable alpha1Gen in alphaSuSp
coldPhaseIsHighAlpha1   false; // false=evaporation (vapor=alpha1=0 is cold), true=melting

// Smoothing / numerics (see "Vestigial parameters" in Status section — several
// of these have no effect on the current code path; tune lambdaSmearCells,
// lambdaAlphaCells, nSmoothIter, mdotMax, maxGradT first):
nSmoothIter         2;
lambdaSmearCells    1.5;
lambdaAlphaCells    1.5;
mdotMax             100.0;
maxGradT            1e8;
```

### `system/`

- **`reconstructionDict`** — **required** by the geometricVoF library (`reconstructionSchemes::New`);
  selects the PLIC reconstruction scheme. Without it, `createFields.H` will fail to construct.
- `controlDict` — `application interTempFoam;`
- `fvSchemes` — needs an entry for the enthalpy flux divergence, e.g.
  `div(rhoCpPhi,T)  Gauss linearUpwind grad(T);`, alongside standard `interFoam` `div(rhoPhi,U)` etc.
- `fvSolution` — standard `alpha.<phase>` MULES block (+ `nAlphaSubCycles` if sub-cycling), PIMPLE
  controls, and a solver entry for `T`.

### Running

```bash
interTempFoam
# parallel:
decomposePar
mpirun -np 4 interTempFoam -parallel
reconstructPar
```

---

## Status: Working vs. Broken / In-Progress

### Working / validated

- Core VOF advection, momentum, and pressure-velocity coupling (inherited `interFoam` machinery).
- Energy equation with variable `ρCp`, consistent `rhoCpPhi` enthalpy flux (fluid-flux-only, by design),
  and sub-cycle accumulation of `rhoPhiSum`/`rhoCpPhiSum`.
- Internally derived conductivity `κ = ρν·cp/Pr` and harmonic-mean face conductivity.
- Temperature-dependent surface tension (Marangoni) model.
- 1D Stefan problem, interface-motion, and 1D melting — validated against Shaikh et al. (see git history;
  commits `801cc9d`, `b16dc0e`).
- Hardt–Wondra source chain (RDF → ψ → qₙ → ṅ → Helmholtz redistribution → `Q_pc`/`phiStefan`) for the
  planar melting/Stefan case.

### Broken / incomplete / in-progress

- **Droplet d²-law evaporation is unfinished.** HEAD commit is literally
  `"d-square law changes - ongoing"`; in-code comments reference unresolved interface lag (~0.49× of
  analytical), a stall at one test condition, and 35–38× spurious volume growth from an earlier attempt
  at the `divU` regularization in `alphaSuSp.H`. Do not treat evaporation results as validated.
- **`Qcorr` (Hardt & Wondra Eq. 42 enthalpy correction) is inert.** It's declared and exposed via the
  base-class virtual interface, but `calcQ_pc()` never assigns it (stays at zero), and it's never read in
  `TEqn.H` or anywhere in the solver.
- **A temporary diagnostic block is still live** in `HardtWondra.C`, explicitly commented
  `"DIAGNOSTIC Stage 0 — remove after PCV sign is confirmed and fixed"`.
- **The Hardt & Wondra Tsat-forcing benchmark block in `TEqn.H` is commented out** (the
  `if (alpha1>0.75) T=Tsat` loop) — left in place but disabled, not deleted.
- **Inconsistent defaults between `HardtWondra`'s constructor and its `read()` method**: e.g.
  `nSmoothIter` defaults to `2` in the constructor but `5` in `read()`; `alphaSmoothWidth` defaults to
  `0.5` vs `1.0`. Behavior can differ depending on whether `read()` is ever invoked.
- **One-timestep phase-change coupling lag by design**: `phaseChangePtr->correct()` runs only on
  `pimple.firstIter()`, not every outer correction, to avoid a `phiStefan` standing-wave instability —
  acceptable for now but a documented approximation, not a bug fix.
- **Several `phaseChangeProperties` dictionary entries are read but have no effect** on the current
  RDF-based code path: `interfaceWidthCells` (explicitly commented "UNUSED on the genuine-RDF path"),
  `betaThermal`, `kinFloorCells`, `liquidBiasCoeff`, `harmonicBias`, `RelaxFac` (no under-relaxation
  applied to `mdot_`), `useEnthalpyCorrection` (ties to the inert `Qcorr`), `nExtrapIter`, `extrapCFL`
  (Hamilton–Jacobi extrapolation controls, not used in the shown path). Setting these does nothing
  currently — don't spend time tuning them.
- **`Gmax` in `alphaSuSp.H` was recently lowered 5→2**, an uncalibrated trade-off between source-rescaling
  accuracy and noise amplification, called out as part of the ongoing d-square work.
- Minor hygiene: a stray editor backup `thermalPhaseChangeModel.H~` is committed; the solver's
  `Make/options` hardcodes an absolute OpenFOAM v2412 install path (see Installation caveat above).

---

## Differences from `interFoam`

| Feature | `interFoam` | `interTempFoam` |
|---|---|---|
| Energy equation | No | Yes — `TEqn.H`, implicit latent-heat/Tsat-pin sink |
| Phase change | No | Yes — runtime-selectable model, default `HardtWondra` |
| `rhoCpPhi` field | No | Enthalpy flux, fluid-flux-only, sub-cycle averaged |
| Interface representation | Algebraic VOF only | + PLIC reconstruction + RDF (geometricVoF) |
| Pressure equation source | `∇·U = 0` | `∇·U = PCV` (mass-transfer dilatation) |
| `transportProperties` | `nu`, `rho` | + `cp`, `Pr` per phase |
| Surface tension | Constant `sigma` | + optional `temperatureDependent` (Marangoni) |

---

## References

- Hardt, S. & Wondra, F. (2008). *Evaporation model for interfacial flows based on a continuum-field
  representation of the source terms.* J. Comput. Phys. — basis for the active phase-change model.
- Shaikh, J.N. et al. — 1D melting validation reference (see `report.md`).
- Natali, A. et al. (2014). *Non-isothermal two-phase flow simulation using VOF method.* — enthalpy flux
  formulation basis.
- Liu, J. et al. (2011). *Heat transfer in two-phase VOF simulations.* — `rhoCpPhi` consistency reference.
- Roenby, J., Bredmose, H., Jasak, H. — geometric VOF / `reconstructionSchemes` /
  `reconstructedDistanceFunction` (isoAdvector lineage).
- Rusche, H. (2002). *Computational Fluid Dynamics of Dispersed Two-Phase Flows at High Phase Fractions.*
  PhD thesis, Imperial College London.
- OpenFOAM `interFoam` solver — base implementation this solver extends.
- See in-repo `applications/solvers/multiphase/interTempFoam/interTempFoam/docs/` and `report.md` /
  `report.pdf` for the fuller derivation and validation history.

---

## License

Distributed under the GNU General Public License v3 — the same license as OpenFOAM. See individual
source file headers for copyright notices.
