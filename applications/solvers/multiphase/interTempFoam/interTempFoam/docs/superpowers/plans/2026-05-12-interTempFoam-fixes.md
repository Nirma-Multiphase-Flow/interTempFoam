# interTempFoam Numerical & Syntactical Fix Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix all numerical, syntactical, and physical correctness bugs in the `interTempFoam` VOF + energy solver (OpenFOAM v2412), validated against the Natali (2014) and Liu (2011) reference implementations.

**Architecture:** `interTempFoam` solves ρ(∂u/∂t + u·∇u) = −∇p + ∇(μ∇u), ∇·u=0, ∂α/∂t + u·∇α=0, and ρCp(∂T/∂t + u·∇T) = ∇·(k∇T). The energy equation is coupled to the VOF alpha field through the mixture thermal flux `rhoCpPhi` (≡ ρCp·φ on faces) and thermal conductivity `kappaf`. These fields must be updated consistently with `rhoPhi` after every alpha solve.

**Tech Stack:** OpenFOAM v2412 C++, wmake, custom `incompressibleTwoPhaseMixture` transport library at `$WM_PROJECT_USER_DIR/src/transportModels/incompressible/`

---

## Reference Validation Summary

Cross-checking current code against:
- **Natali (2014)** — "Adding heat transfer to InterFoam, OpenFOAM 2.3.0" (interTempFoam.pdf)
- **Liu (2011)** — "Implementation of myinterFoamDiabatic Solver" (Project_QingmingLIU-final.pdf)

### Verified Correct (no changes needed)
- `TEqn.H` — formulation `fvm::ddt(rhoCp,T) + fvm::div(rhoCpPhi,T) - fvm::laplacian(kappaf,T)` is correct per both references. `TEqn.relax()` and `T.correctBoundaryConditions()` are valid v2412 additions. ✓
- `alphaEqnSubCycle.H` lines 36-38 — `rhoCp == alpha1*rho1*cp1 + (1-alpha1)*rho2*cp2` is correct. ✓
- `alphaEqn.H` lines 260-263 — non-Euler (CrankNicolson) branch correctly uses `alphaPhi10` for both `rhoPhi` and `rhoCpPhi`. ✓
- `incompressibleTwoPhaseMixture.C` — `kappaf()` formula `α*ρ1*cp1*(1/Pr1)*ν1 + (1-α)*ρ2*cp2*(1/Pr2)*ν2` is correct (k = ρνcp/Pr = μcp/Pr). ✓
- `incompressibleTwoPhaseMixture.H` — `cp1_`, `cp2_`, `Pr1_`, `Pr2_` members and accessors `cp1()`, `cp2()`, `Pr1()`, `Pr2()`, `kappaf()` exist and are properly declared. ✓
- `createFields.H` — `rhoCp` field and `rhoCp.oldTime()` call correct. ✓
- `Make/options` — library linkage for custom transport model library is correct. ✓

### Bugs Found

| # | File | Lines | Severity | Description |
|---|------|--------|----------|-------------|
| 1 | `alphaEqn.H` | 247-249 | **CRITICAL** | Euler ddt branch updates `rhoPhi` but NOT `rhoCpPhi` → energy convection uses stale flux for all Euler-scheme runs (the default) |
| 2 | `createFields.H` | 120 | **IMPORTANT** | `rhoCpPhi` initialized as `rhoPhi*cp1` instead of `fvc::interpolate(rhoCp)*phi` → physically wrong initial flux (equals correct only if cp1==cp2) |
| 3 | `alphaEqnSubCycle.H` | 19-28 | **IMPORTANT** | Sub-cycle loop accumulates `rhoPhiSum` but not `rhoCpPhiSum` → when `nAlphaSubCycles > 1`, `rhoCpPhi` is not time-averaged across sub-cycles |
| 4 | `interTempFoam.C` | 28, 34-38, 60-66 | cosmetic | Application name/description still reads "interFoam" / "isothermal immiscible fluids" |

---

## File Map

| File | Action | What changes |
|------|--------|--------------|
| `alphaEqn.H` | Modify | Add `rhoCpPhi` update in Euler branch (line 249) |
| `createFields.H` | Modify | Fix `rhoCpPhi` initializer (line 120) |
| `alphaEqnSubCycle.H` | Modify | Add `rhoCpPhiSum` sub-cycle accumulation |
| `interTempFoam.C` | Modify | Fix application name/description strings |

---

## Task 1: Fix `rhoCpPhi` initialization in `createFields.H`

**Files:**
- Modify: `createFields.H:107-121`

The root path is `/home/param/OpenFOAM/param-v2412/applications/solvers/multiphase/interTempFoam/`.

**Problem:** `rhoCpPhi` is initialized as `rhoPhi*cp1`, which equals `fvc::interpolate(rho)*phi*cp1`. This is incorrect because it applies phase-1 specific heat to the entire mixture density flux. The correct value is `fvc::interpolate(rhoCp)*phi`, where `rhoCp = alpha1*rho1*cp1 + alpha2*rho2*cp2`.

**Validation (Natali tutorial, createFields slide):**
```
surfaceScalarField rhoCpPhi
(
    IOobject("rhoCpPhi", ..., NO_READ, NO_WRITE),
    fvc::interpolate(rhoCp)*phi
);
```

- [ ] **Step 1: Read the file to confirm current state**

  Open `createFields.H` and verify lines 107-121 read:
  ```cpp
  Info<< "Creating rhoCpPhi\n" << endl;
  surfaceScalarField rhoCpPhi
  (
      IOobject
      (
          "rhoCpPhi",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::NO_WRITE
      ),
      rhoPhi*cp1
  );
  ```

- [ ] **Step 2: Replace the incorrect initializer**

  In `createFields.H`, replace the initializer of `rhoCpPhi`:

  Old text (exact match):
  ```cpp
      rhoPhi*cp1
  ```

  New text:
  ```cpp
      fvc::interpolate(rhoCp)*phi
  ```

  Result: `rhoCpPhi` will be initialized as the face-interpolated `ρCp·φ` mixture flux, consistent with:
  `rhoCp = alpha1*rho1*cp1 + alpha2*rho2*cp2` (already computed on lines 91-105 of `createFields.H`).

- [ ] **Step 3: Verify the edit looks correct**

  Check lines 107-121 of `createFields.H` now read:
  ```cpp
  Info<< "Creating rhoCpPhi\n" << endl;

  surfaceScalarField rhoCpPhi
  (
      IOobject
      (
          "rhoCpPhi",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::NO_WRITE
      ),

      fvc::interpolate(rhoCp)*phi
  );
  ```

---

## Task 2: Add `rhoCpPhi` update in Euler branch of `alphaEqn.H`

**Files:**
- Modify: `alphaEqn.H:247-249`

**Problem:** The if-else block (lines 240-264) updating `rhoPhi` after the MULES alpha solve has two branches:
- **Euler branch** (lines 247-249): updates only `rhoPhi`, leaving `rhoCpPhi` stale.
- **Non-Euler branch** (lines 260-263): correctly updates both `rhoPhi` and `rhoCpPhi`.

Since Euler is the default ddt scheme for `ddt(rho,U)`, this means `rhoCpPhi` is **never updated** in standard runs. The energy convection term `fvm::div(rhoCpPhi, T)` in `TEqn.H` therefore uses the (wrong) initial value throughout the entire simulation.

**Reference (Natali, alphaEqn.H slide):**
```
rhoPhi = tphiAlpha()*(rho1 - rho2) + phi*rho2;
rhoCpPhi = tphiAlpha()*(rho1*cp1 - rho2*cp2) + phi*rho2*cp2;
```
The pattern `rhoCpPhi = alphaPhi*(rho1*cp1 - rho2*cp2) + phi*rho2*cp2` mirrors exactly the `rhoPhi` update pattern `rhoPhi = alphaPhi*(rho1 - rho2) + phi*rho2`.

- [ ] **Step 1: Confirm current state of the Euler branch**

  Lines 240-264 of `alphaEqn.H` should read:
  ```cpp
      if
      (
          word(mesh.ddtScheme("ddt(rho,U)"))
       == fv::EulerDdtScheme<vector>::typeName
       || word(mesh.ddtScheme("ddt(rho,U)"))
       == fv::localEulerDdtScheme<vector>::typeName
      )
      {
          rhoPhi = alphaPhi10*(rho1f - rho2f) + phiCN*rho2f;
      }
      else
      {
          ...
          rhoPhi = alphaPhi10*(rho1f - rho2f) + phi*rho2f;
          rhoCpPhi =
              alphaPhi10*(rho1*cp1 - rho2*cp2)
              + phi*rho2*cp2;
      }
  ```

- [ ] **Step 2: Add the missing `rhoCpPhi` update line in the Euler branch**

  In `alphaEqn.H`, replace the Euler branch body:

  Old text (exact):
  ```cpp
      {
          rhoPhi = alphaPhi10*(rho1f - rho2f) + phiCN*rho2f;
      }
      else
  ```

  New text:
  ```cpp
      {
          rhoPhi = alphaPhi10*(rho1f - rho2f) + phiCN*rho2f;
          rhoCpPhi = alphaPhi10*(rho1*cp1 - rho2*cp2) + phiCN*rho2*cp2;
      }
      else
  ```

  **Dimensional check:**
  - `alphaPhi10` has dimensions [m³/s] (volume flux)
  - `rho1*cp1` has dimensions [kg/m³ · J/(kg·K)] = [J/(m³·K)]
  - `alphaPhi10*(rho1*cp1 - rho2*cp2)` → [m³/s · J/(m³·K)] = [J/(s·K)] = [W/K]
  - `phiCN*rho2*cp2` → [m³/s · kg/m³ · J/(kg·K)] = [J/(s·K)] = [W/K]
  - `rhoCpPhi` has dimensions [W/K] ✓ (same as `rhoPhi` dimensioned as [kg/s] scaled by cp)

  Note: `rho1`, `rho2` are `dimensionedScalar` (not `rho1f`, `rho2f` which are aliases for the same thing defined in `rhofs.H`). Using `rho1*cp1` is consistent with the non-Euler branch on line 262.

- [ ] **Step 3: Verify the complete if-else block**

  After the edit, lines 240-264 should read:
  ```cpp
      if
      (
          word(mesh.ddtScheme("ddt(rho,U)"))
       == fv::EulerDdtScheme<vector>::typeName
       || word(mesh.ddtScheme("ddt(rho,U)"))
       == fv::localEulerDdtScheme<vector>::typeName
      )
      {
          rhoPhi = alphaPhi10*(rho1f - rho2f) + phiCN*rho2f;
          rhoCpPhi = alphaPhi10*(rho1*cp1 - rho2*cp2) + phiCN*rho2*cp2;
      }
      else
      {
          if (ocCoeff > 0)
          {
              alphaPhi10 =
                  (alphaPhi10 - (1.0 - cnCoeff)*alphaPhi10.oldTime())/cnCoeff;
          }

          rhoPhi = alphaPhi10*(rho1f - rho2f) + phi*rho2f;
          rhoCpPhi =
              alphaPhi10*(rho1*cp1 - rho2*cp2)
              + phi*rho2*cp2;
      }
  ```

---

## Task 3: Add `rhoCpPhi` sub-cycle accumulation in `alphaEqnSubCycle.H`

**Files:**
- Modify: `alphaEqnSubCycle.H:1-39`

**Problem:** When `nAlphaSubCycles > 1`, the sub-cycle loop accumulates `rhoPhiSum` and assigns it to `rhoPhi` to get the time-averaged mass flux. But there is no equivalent `rhoCpPhiSum` for `rhoCpPhi`. After the sub-cycle loop, `rhoCpPhi` holds only the value from the last sub-cycle iteration (and for the Euler scheme, it holds the wrong stale value until Task 2's fix is applied). The energy convection term therefore uses an inconsistent flux that does not match `rhoPhi`.

**Reference pattern (Liu, Natali):** The sub-cycle time-averaging for `rhoCpPhi` must mirror what is done for `rhoPhi`. In standard interFoam, only `rhoPhiSum` is needed. For `interTempFoam`, an additional `rhoCpPhiSum` must be accumulated.

- [ ] **Step 1: Read current `alphaEqnSubCycle.H`**

  Confirm the full file reads:
  ```cpp
  if (nAlphaSubCycles > 1)
  {
      dimensionedScalar totalDeltaT = runTime.deltaT();
      surfaceScalarField rhoPhiSum
      (
          mesh.newIOobject("rhoPhiSum"),
          mesh,
          dimensionedScalar(rhoPhi.dimensions(), Zero)
      );

      tmp<volScalarField> trSubDeltaT;

      if (LTS)
      {
          trSubDeltaT =
              fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
      }

      for
      (
          subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
          !(++alphaSubCycle).end();
      )
      {
          #include "alphaEqn.H"
          rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
      }

      rhoPhi = rhoPhiSum;
  }
  else
  {
      #include "alphaEqn.H"
  }

  rho == alpha1*rho1 + alpha2*rho2;
  rhoCp ==
      alpha1*rho1*cp1
    + (scalar(1) - alpha1)*rho2*cp2;
  ```

- [ ] **Step 2: Add `rhoCpPhiSum` alongside `rhoPhiSum`**

  Replace the full sub-cycle block:

  Old text (exact):
  ```cpp
  if (nAlphaSubCycles > 1)
  {
      dimensionedScalar totalDeltaT = runTime.deltaT();
      surfaceScalarField rhoPhiSum
      (
          mesh.newIOobject("rhoPhiSum"),
          mesh,
          dimensionedScalar(rhoPhi.dimensions(), Zero)
      );

      tmp<volScalarField> trSubDeltaT;

      if (LTS)
      {
          trSubDeltaT =
              fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
      }

      for
      (
          subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
          !(++alphaSubCycle).end();
      )
      {
          #include "alphaEqn.H"
          rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
      }

      rhoPhi = rhoPhiSum;
  }
  ```

  New text:
  ```cpp
  if (nAlphaSubCycles > 1)
  {
      dimensionedScalar totalDeltaT = runTime.deltaT();
      surfaceScalarField rhoPhiSum
      (
          mesh.newIOobject("rhoPhiSum"),
          mesh,
          dimensionedScalar(rhoPhi.dimensions(), Zero)
      );
      surfaceScalarField rhoCpPhiSum
      (
          mesh.newIOobject("rhoCpPhiSum"),
          mesh,
          dimensionedScalar(rhoCpPhi.dimensions(), Zero)
      );

      tmp<volScalarField> trSubDeltaT;

      if (LTS)
      {
          trSubDeltaT =
              fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
      }

      for
      (
          subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
          !(++alphaSubCycle).end();
      )
      {
          #include "alphaEqn.H"
          rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
          rhoCpPhiSum += (runTime.deltaT()/totalDeltaT)*rhoCpPhi;
      }

      rhoPhi = rhoPhiSum;
      rhoCpPhi = rhoCpPhiSum;
  }
  ```

- [ ] **Step 3: Verify the full updated file**

  The complete file should now read:
  ```cpp
  if (nAlphaSubCycles > 1)
  {
      dimensionedScalar totalDeltaT = runTime.deltaT();
      surfaceScalarField rhoPhiSum
      (
          mesh.newIOobject("rhoPhiSum"),
          mesh,
          dimensionedScalar(rhoPhi.dimensions(), Zero)
      );
      surfaceScalarField rhoCpPhiSum
      (
          mesh.newIOobject("rhoCpPhiSum"),
          mesh,
          dimensionedScalar(rhoCpPhi.dimensions(), Zero)
      );

      tmp<volScalarField> trSubDeltaT;

      if (LTS)
      {
          trSubDeltaT =
              fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
      }

      for
      (
          subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
          !(++alphaSubCycle).end();
      )
      {
          #include "alphaEqn.H"
          rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
          rhoCpPhiSum += (runTime.deltaT()/totalDeltaT)*rhoCpPhi;
      }

      rhoPhi = rhoPhiSum;
      rhoCpPhi = rhoCpPhiSum;
  }
  else
  {
      #include "alphaEqn.H"
  }

  rho == alpha1*rho1 + alpha2*rho2;
  rhoCp ==
      alpha1*rho1*cp1
    + (scalar(1) - alpha1)*rho2*cp2;
  ```

---

## Task 4: Fix application description in `interTempFoam.C`

**Files:**
- Modify: `interTempFoam.C:28-38`, `interTempFoam.C:60-66`

**Problem:** The header block and runtime note still describe the old interFoam solver as "isothermal." This is cosmetic but incorrect for a heat-transfer solver.

- [ ] **Step 1: Fix the Application doxygen tag (line 28)**

  Old text:
  ```cpp
  Application
      interFoam
  ```

  New text:
  ```cpp
  Application
      interTempFoam
  ```

- [ ] **Step 2: Fix the Description doxygen block (lines 34-38)**

  Old text:
  ```cpp
  Description
      Solver for two incompressible, isothermal immiscible fluids using a VOF
      (volume of fluid) phase-fraction based interface capturing approach,
      with optional mesh motion and mesh topology changes including adaptive
      re-meshing.
  ```

  New text:
  ```cpp
  Description
      Solver for two incompressible, non-isothermal immiscible fluids using a VOF
      (volume of fluid) phase-fraction based interface capturing approach,
      with energy equation for heat transfer. Optional mesh motion and mesh
      topology changes including adaptive re-meshing.
  ```

- [ ] **Step 3: Fix the argList::addNote string (lines 60-66)**

  Old text:
  ```cpp
      argList::addNote
      (
          "Solver for two incompressible, isothermal immiscible fluids"
          " using VOF phase-fraction based interface capturing.\n"
          "With optional mesh motion and mesh topology changes including"
          " adaptive re-meshing."
      );
  ```

  New text:
  ```cpp
      argList::addNote
      (
          "Solver for two incompressible, non-isothermal immiscible fluids"
          " using VOF phase-fraction based interface capturing with energy equation.\n"
          "With optional mesh motion and mesh topology changes including"
          " adaptive re-meshing."
      );
  ```

---

## Task 5: Rebuild and verify compilation

**Files:**
- No file edits — build verification only.

- [ ] **Step 1: Source OpenFOAM environment**

  Run:
  ```bash
  source /opt/openfoam/openfoam-v2412/etc/bashrc
  ```
  (adjust path to the actual OpenFOAM v2412 installation if different)

  Or check:
  ```bash
  echo $WM_PROJECT_VERSION
  ```
  Expected output: `v2412`

- [ ] **Step 2: Clean previous build artifacts**

  From the solver directory:
  ```bash
  cd /home/param/OpenFOAM/param-v2412/applications/solvers/multiphase/interTempFoam
  wclean
  ```
  Expected: messages like `Cleaning .../Make/linux64GccDPInt32Opt`

- [ ] **Step 3: Rebuild the custom transport library first**

  The transport library at `$WM_PROJECT_USER_DIR/src/transportModels/incompressible/` provides `cp1()`, `cp2()`, `kappaf()`. It must be compiled before the solver.

  ```bash
  cd $WM_PROJECT_USER_DIR/src/transportModels/incompressible
  wmake libso
  ```
  Expected: Compiles and links to `$FOAM_USER_LIBBIN/libincompressibleTransportModels.so` (or the name in its Make/files).

  Also rebuild immiscibleIncompressibleTwoPhaseMixture if it has a separate Make/:
  ```bash
  cd $WM_PROJECT_USER_DIR/src/transportModels/immiscibleIncompressibleTwoPhaseMixture
  wmake libso
  ```

- [ ] **Step 4: Build the solver**

  ```bash
  cd /home/param/OpenFOAM/param-v2412/applications/solvers/multiphase/interTempFoam
  wmake
  ```
  Expected: Compiles `interTempFoam.C` and links to `$FOAM_USER_APPBIN/interTempFoam`.

  Zero errors expected. Zero warnings about undefined variables.

- [ ] **Step 5: Confirm binary exists**

  ```bash
  ls -lh $FOAM_USER_APPBIN/interTempFoam
  ```
  Expected: Recent timestamp (today's date).

---

## Self-Review Against References

### Spec coverage check

| Reference requirement | Task covering it |
|----------------------|-----------------|
| `rhoCpPhi` initialized correctly | Task 1 |
| `rhoCpPhi` updated in Euler scheme alpha loop | Task 2 |
| `rhoCpPhi` updated in CrankNicolson alpha loop | Already correct — no task needed |
| `rhoCpPhi` sub-cycle time-averaged | Task 3 |
| `rhoCp` updated in sub-cycle | Already correct (lines 36-38) — no task needed |
| `TEqn.H` formulation correct | Already correct — no task needed |
| `kappaf()` formula correct | Already correct — no task needed |
| Transport model reads `cp`, `Pr` from transportProperties | Already correct — no task needed |
| Application description accurate | Task 4 |
| Solver compiles cleanly | Task 5 |

### Placeholder scan
No TBDs, TODOs, or "implement later" in any task above. All code blocks are complete.

### Type consistency
- `rhoCpPhi` initialized as `fvc::interpolate(rhoCp)*phi` in Task 1. Type: `surfaceScalarField`. ✓
- `rhoCpPhi` assigned in Task 2 as `alphaPhi10*(rho1*cp1 - rho2*cp2) + phiCN*rho2*cp2`. Type: `surfaceScalarField`. ✓
- `rhoCpPhiSum` declared as `surfaceScalarField` in Task 3 with `rhoCpPhi.dimensions()`. ✓
- `rhoCpPhi = rhoCpPhiSum` assignment in Task 3 is type-consistent. ✓

### Dimensional consistency of Task 2 fix
- `alphaPhi10`: [m³/s] (face volume flux)
- `rho1`: [kg/m³], `cp1`: [J/(kg·K)] → `rho1*cp1`: [J/(m³·K)]
- `alphaPhi10*(rho1*cp1 - rho2*cp2)`: [m³/s · J/(m³·K)] = [J/(s·K)] = [W/K]
- `phiCN`: [m³/s], `rho2`: [kg/m³], `cp2`: [J/(kg·K)]
- `phiCN*rho2*cp2`: [m³/s · kg/m³ · J/(kg·K)] = [W/K]
- `rhoCpPhi` total: [W/K] ✓ — This matches the expected dimension for the convective heat flux in the energy equation `fvm::div(rhoCpPhi, T)` where T is [K].
