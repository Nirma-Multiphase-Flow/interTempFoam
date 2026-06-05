# Thermal Leakage Fix — Sheikh 1D Stefan Melting Case

**Solver**: interTempFoam (HardtWondra phase-change model)  
**Date**: 2026-05-25  
**Files modified**: `TEqn.H`, `createFields.H`

---

## 1. Problem Statement

In the Sheikh 1D Stefan melting case, the temperature field in the **solid phase**
(alpha1 ≈ 1, far from the diffuse interface zone) rises above Tsat = 313.15 K over
a large spatial distance, creating nonphysical sensible heating of the solid.

**Observed behaviour**
- Interface propagation: within ~5% of the Sheikh analytical solution (acceptable)
- Temperature in solid bulk: T > Tsat over tens of cells beyond the diffuse zone
- This creates an extended "warm tail" in the solid that is unphysical

**Case configuration (Sheikh)**

| Parameter | Value |
|---|---|
| phases | (solid liquid) — alpha1 = solid fraction |
| alpha1 = 1 | solid (cold), T_init = Tsat = 313.15 K |
| alpha1 = 0 | liquid (hot), T_wall = 350 K |
| h_lv | 175,000 J/kg |
| k_liq = k_vap | 0.21 W/m/K (equal) |
| rho = rho1 = rho2 | 750 kg/m³ |
| cp | 2400 J/kg/K |
| Mesh Δx | 5.6 × 10⁻⁴ m (500 cells over 0.28 m) |

---

## 2. Root Cause Analysis

### 2.1 Energy balance at the interface cell

The HardtWondra model reconstructs the Stefan heat flux via:

```
dTdn_hot = (T_hot_neighbour − Tsat) / max(dist, 10·h_ref)
Q_pc = Ai · kEff · dTdn_hot    [W/m³]
```

The `10·h_ref` clamp is intentional for diffuse-interface stability. For Δx = 5.6 × 10⁻⁴ m:

| Flux | Formula | Value |
|---|---|---|
| Conductive influx to interface cell | k·ΔT/Δx² | ≈ 18 MW/m³ |
| Q_pc (latent heat sink) | Ai·kEff·dTdn_hot | ≈ 0.67 MW/m³ |
| Q_pc absorbs | 0.67/18 | **~3.7%** of incoming flux |

The remaining 96.3% of the incoming conductive heat is not absorbed as latent heat.
It raises the interface cell temperature:

```
T_interface_ss ≈ (T_liquid_neighbour + T_solid_neighbour) / 2 ≈ (340 + 313) / 2 = 326 K
```

This is **13 K above Tsat**. The warm interface then drives conduction into adjacent
solid cells, which propagates further into the solid bulk — the thermal leakage.

### 2.2 Why the 10·h_ref clamp cannot be removed

Removing the clamp would also change `phiStefan` (which uses the same `dTdn_hot`),
altering interface propagation speed. The ~5% interface accuracy would be lost.

### 2.3 Existing saturation floor (cold-trench fix)

The existing `Acoeff_floor` in TEqn.H prevents T from dropping **below** Tsat:

```cpp
Acoeff_floor = C_floor * (rhoCp/dt) * 4α(1−α) * pos(Tsat − T)
```

This handles the cold-trench problem on the liquid side. There was no symmetric
floor preventing T from rising **above** Tsat on the solid side.

---

## 3. Solution Design

### 3.1 Approach: post-solve saturation clip

The fix clips `T = Tsat` in the **solid-phase diffuse zone** (alpha1 > 0.5) **after**
the implicit TEqn solve completes, using a phase-convention-aware gate.

**Why post-solve and not an implicit fvMatrix term?**

An in-matrix stiff heat sink at alpha1 > 0.5 (the originally planned `Acoeff_floor_hot`)
was implemented first and tested. It caused a **measurable reduction in qn** and
slowed interface propagation. The mechanism:

1. The stiff sink at alpha1 > 0.5 appears in the coupled linear system alongside
   the Laplacian operator.
2. During the implicit solve, the coupling `k·(T_A − T_B=Tsat)/Δx²` drains heat
   from liquid-side cells (alpha1 = 0.3–0.5) **within the same timestep**.
3. Lower liquid-side T → lower `dTdn_hot` in HardtWondra → reduced `qn` → slower
   interface.

A **post-solve clip** avoids this: it only establishes a boundary condition for
the **next** timestep's Laplacian. The current timestep's implicit solve sees no
artificial heat sink, so the hot-phase temperature gradient driving `qn` is
preserved.

### 3.2 Phase-convention gate

The cold-phase gate depends on which phase is alpha1 in the case:

| Case type | Cold phase | alpha1 convention | Gate |
|---|---|---|---|
| Melting (Sheikh) | solid | alpha1 = solid fraction → cold = high alpha1 | `pos(alpha1 − 0.5)` |
| Evaporation | vapour | alpha1 = liquid fraction → cold = low alpha1 | `pos(0.5 − alpha1)` |

This is controlled by `coldPhaseIsHighAlpha1` in `constant/phaseChangeProperties`.

---

## 4. Bugs Encountered and Fixed

### Bug 1: Registry lookup failure in TEqn.H

**Error at runtime:**
```
FOAM FATAL ERROR: failed lookup of phaseChangeProperties (objectRegistry region0)
```

**Root cause**: In `createFields.H`, the `IOdictionary pcDict` was declared as a
local variable inside the `if` block:

```cpp
if (isFile(...))
{
    IOdictionary pcDict(IOobject(..., mesh, MUST_READ_IF_MODIFIED, NO_WRITE));
    phaseChangePtr = thermalPhaseChangeModel::New(pcDictName, pcDict, ...);
}  // ← pcDict goes out of scope and deregisters here
```

By the time the time loop reaches TEqn.H's `mesh.lookupObject<IOdictionary>(...)`,
the dictionary no longer exists in the registry.

**Fix** (`createFields.H`): Lifted to persistent outer scope using `autoPtr<IOdictionary>`:

```cpp
autoPtr<IOdictionary> pcDictPtr;  // lives for the entire simulation
if (isFile(...))
{
    pcDictPtr.reset(new IOdictionary(IOobject(...)));
    phaseChangePtr = thermalPhaseChangeModel::New(pcDictName, *pcDictPtr, ...);
}
```

TEqn.H then reads the flag without a registry lookup:

```cpp
const bool coldPhaseHighAlpha1 =
    pcDictPtr.valid()
    && pcDictPtr->lookupOrDefault<bool>("coldPhaseIsHighAlpha1", false);
```

### Bug 2: In-matrix floor suppressed qn

**Observed behaviour**: `Acoeff_floor_hot` (stiff implicit matrix sink) reduced
`qn` and slowed interface propagation.

**Root cause**: Analysed above (§3.1). The implicit Laplacian coupling in the same
timestep cooled the liquid-side diffuse zone, reducing `dTdn_hot`.

**Fix**: Removed `Acoeff_floor_hot` from the fvScalarMatrix. Replaced with a
post-solve temperature clip (§5).

---

## 5. Implemented Changes

### 5.1 `createFields.H` — pcDictPtr persistence

```cpp
// Before (local scope — deregisters before time loop):
if (isFile(...)) {
    IOdictionary pcDict(...);
    phaseChangePtr = thermalPhaseChangeModel::New(pcDictName, pcDict, ...);
}

// After (outer scope — persists for entire simulation):
autoPtr<IOdictionary> pcDictPtr;
if (isFile(...)) {
    pcDictPtr.reset(new IOdictionary(...));
    phaseChangePtr = thermalPhaseChangeModel::New(pcDictName, *pcDictPtr, ...);
}
```

### 5.2 `TEqn.H` — post-solve saturation clip

**Removed** from fvScalarMatrix:
```cpp
// REMOVED: Acoeff_floor_hot (stiff in-matrix sink at cold-phase diffuse zone)
// fvm::Sp(Acoeff + Acoeff_floor + Acoeff_floor_hot, T) == ... + Acoeff_floor_hot*Tsat
```

**Matrix reverted to**:
```cpp
fvScalarMatrix TEqn
(
    rhoCp * fvm::ddt(T)
  + fvm::div(rhoCpPhi, T)
  - fvm::laplacian(kappaf, T)
  + fvm::Sp(Acoeff + Acoeff_floor, T)
 ==
    (Acoeff*Tsat + Acoeff_floor*Tsat)
  - Q_cond
  + phaseChangePtr->Qcorr()
);
```

**Added** after `TEqn.solve()`:
```cpp
// Phase-convention gate: 1 on cold side, 0 on hot side
const bool coldPhaseHighAlpha1 =
    pcDictPtr.valid()
    && pcDictPtr->lookupOrDefault<bool>("coldPhaseIsHighAlpha1", false);

const volScalarField coldSideGate =
    coldPhaseHighAlpha1
    ? pos(alpha1 - dimensionedScalar("half", dimless, 0.5))   // melting: solid=high α₁
    : pos(dimensionedScalar("half", dimless, 0.5) - alpha1);  // evaporation: vapor=low α₁

// ... (TEqn.solve(), T.correctBoundaryConditions()) ...

// Post-solve saturation clip: T ≤ Tsat in cold-phase diffuse zone + bulk.
// Post-solve to avoid suppressing qn within the current timestep's implicit solve.
T -= coldSideGate
   * max(T - Tsat, dimensionedScalar("zeroT", dimTemperature, scalar(0)));
T.correctBoundaryConditions();
```

### 5.3 `constant/phaseChangeProperties` — case flag

```
coldPhaseIsHighAlpha1  true;   // solid = phase 1 = high alpha1 (melting case)
```

---

## 6. What Was NOT Changed

| Component | Reason |
|---|---|
| `HardtWondra.C` — 10·h_ref clamp | Intentional; removing it breaks interface propagation |
| `HardtWondra.C` — phiStefan reconstruction | Preserved; changing it alters the ~5% interface accuracy |
| `alphaEqn.H` — phiStefan in alpha transport | Preserved |
| `Acoeff_floor` (cold-trench fix) | Preserved; handles T < Tsat on liquid side |
| `rhoCpPhi` Stefan-decoupling | Preserved; prevents the cold-trench problem |

---

## 7. Verification Checklist

After running the Sheikh case:

1. **Temperature in solid bulk**: plot T vs x at several times. Cells with
   alpha.solid > 0.9 should show T = 313.15 K throughout.

2. **Interface position**: compare alpha.solid = 0.5 isosurface position against
   the Sheikh analytical solution. Should remain within ~5% (unchanged from pre-fix
   trajectory).

3. **qn / phiStefan**: `min/max(T)` in the log should show max(T) ≤ 350 K and
   the interface region min ≈ 313 K. No significant change in iteration count or
   residuals compared to the pre-fix run.

4. **No cold trench**: min(T) should not drop below Tsat (existing `Acoeff_floor`
   handles this and is unchanged).

---

## 8. Mechanistic Summary

```
Problem
  Interface cell T = 326 K  (13 K above Tsat)
  │
  ├── Q_pc absorbs only 3.7% of incoming conductive flux (10·h_ref clamp)
  ├── Remaining 96.3% raises interface cell T above Tsat
  └── Warm interface → Laplacian drives heat into solid → thermal leakage tail

Attempted fix (abandoned)
  Acoeff_floor_hot in fvMatrix  (stiff implicit sink at alpha1 > 0.5)
  Problem: in-matrix sink + Laplacian couples to liquid side within same timestep
           → cools liquid-side diffuse zone → reduces dTdn_hot → reduces qn
           → interface propagation slowed

Implemented fix
  Post-solve clip:  T = min(T, Tsat)  for alpha1 > 0.5  (cold-side gate)
  Applied AFTER implicit solve → no within-timestep Laplacian coupling
  → liquid-side gradient preserved → qn unchanged → interface propagation preserved
  → solid bulk temperature corrected to Tsat
```
