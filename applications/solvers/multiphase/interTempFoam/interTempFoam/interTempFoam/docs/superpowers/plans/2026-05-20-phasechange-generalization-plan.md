# HardtWondra Phase-Change Framework: Technical Evaluation & Generalization Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generalize the HardtWondra phase-change model to handle both wall|vapor|liquid and wall|liquid|vapor configurations with correct evaporation and condensation physics.

**Architecture:** Surgical corrections to `HardtWondra.C` only — the qn sign convention, evapSwitch condensation suppression, and Helmholtz source one-sidedness. `TEqn.H`, `alphaEqn.H`, `pEqn.H` already handle signed Q_pc correctly and need no changes.

**Tech Stack:** OpenFOAM v2412 C++, wmake, `thermalPhaseChangeModels` library at `src/transportModels/thermalPhaseChangeModels/`

---

## Part I — Structured Technical Evaluation

### 1. Phase Convention

| Field | Meaning | Phase 1 (water) | Phase 2 (steam) |
|---|---|---|---|
| `alpha1` | liquid volume fraction | 1 | 0 |
| `rho1/rho2` | density | 1000 kg/m³ | 1 kg/m³ |
| `cp1/cp2` | specific heat | 4216 J/kg/K | 2000 J/kg/K |
| `Q_pc > 0` | evaporation, heat sink | — | — |
| `Q_pc < 0` | condensation, heat source | — | — |
| `mdot > 0` | evaporation mass rate | [kg/m³/s] | — |

**nHat convention:** `nHat = grad(alpha)/|grad(alpha)|` always points **from vapor (alpha=0) toward liquid (alpha=1)**, regardless of which side the hot wall is on.

---

### 2. 1D Stefan Case Geometry

From `0/alpha.water` (1000 cells, leftWall T=383.15 K):

```
leftWall | cells 0-1: alpha=0 (steam) | cells 2-499: alpha=1 | cells 500-501: alpha=0 (steam) | cells 502-999: alpha=1 | rightOutlet
```

Two interface zones at t=0:

| Zone | Config | nHat direction | gradT direction |
|---|---|---|---|
| A (cells 1-2 boundary) | wall\|vapor\|liquid | +x | -x |
| B-left (cells 499-500) | liquid\|vapor | -x | -x |
| B-right (cells 501-502) | vapor\|liquid | +x | -x |

Zone A is the primary Stefan test geometry. nHat and gradT are **anti-aligned** in Zone A.

---

### 3. Source-Term Dependency Chain

```
T (prev step) + alpha1 (prev step)
      |
      v  [pimple.firstIter() only]
phaseChangePtr->correct() -> calcQ_pc()
  STEP 1: alpha pseudo-diffusion -> alphaSmooth, alphaGeom
  STEP 2: interfaceArea_ = |grad(alphaGeom)|        [1/m]
          nHat = grad(alphaGeom) / (|grad(alphaGeom)| + eps)
  STEP 3: kEff = k_liq*k_vap / (alpha*k_vap + bias*(1-alpha)*k_liq)
  STEP 4: TSense = mild-smoothed T
          gradT = fvc::grad(TSense)
  STEP 5: qnRaw = max(kEff*(gradT & nHat), 0)   <- ISSUE 1: SIGN
          evapSwitch = 0.5*(1+tanh((TSense-Tsat)/0.25))
          qn_ = qnRaw * evapSwitch                <- ISSUE 2: BLOCKS CONDENSATION
  STEP 6: mdotRaw_ = interfaceArea_ * qn_ / h_lv
  STEP 7: Helmholtz: mdot - lambda^2*laplacian(mdot) = max(mdotRaw_,0)  <- ISSUE 3
          mdot_ = mdotLimiter * tanh(mdot/mdotLimiter)
          mdot_ *= RelaxFac_
  STEP 8: phiStefan_ = Ustef*(nHatf & Sf) * interfaceMaskF   <- ISSUE 4: unbounded mask
  STEP 9: Q_pc_ = mdot_ * h_lv_
  STEP 10: Qcorr_ = -mdot_*(cp1-cp2)*(T-Tsat)

      |
      v  alphaSuSp.H
divU = PCV = Q_pc/h_lv*(1/rho2 - 1/rho1)   [used in alphaEqn + pEqn]
Su   = 0  <- ISSUE 5: alpha1Gen() never called despite PhaseFractionSource switch

      |
      v  alphaEqn
dAlpha/dt + div((phi+phiStefan)*alpha) = (0 + PCV*alpha)

      |
      v  pEqn
laplacian(rAUf, p_rgh) = div(phiHbyA) - PCV

      |
      v  TEqn
rhoCp*dT/dt + div(rhoCpPhi,T) - laplacian(kappaf,T)
  + Sp(Acoeff,T) == Acoeff*Tsat - Q_cond + Qcorr
where:
  Acoeff  = max(Q_pc/max(T-Tsat,5K), 0)   [implicit evap linearization]
  Q_cond  = min(Q_pc, 0)                   [explicit cond source]
  kappaf  = HardtWondra::kappaf()          [harmonic mean from alpha1_]
  rhoCpPhi = phi * interpolate(rhoCp)      [hydrodynamic flux only - correct]
```

---

### 4. Identified Issues

#### Issue 1 — CRITICAL: qnRaw sign kills evaporation for wall|vapor|liquid

**File:** `HardtWondra.C:488`

```cpp
qnRaw = max(kEff*(gradT & nHat), 0)
```

Sign analysis:

| Config | nHat | gradT (hot left) | gradT & nHat | qnRaw |
|---|---|---|---|---|
| wall\|liquid\|vapor | -x | -x | **positive** | > 0 ✓ |
| wall\|vapor\|liquid | +x | -x | **negative** | **= 0 ✗** |

For the current 1D Stefan test (Zone A = wall|vapor|liquid): `qnRaw = 0` everywhere. Consequence: `mdotRaw=0`, `mdot=0`, `Q_pc=0`, `phiStefan=0`. **No evaporation occurs.** The temperature in the vapor layer rises freely above Tsat with no latent heat sink.

---

#### Issue 2 — CRITICAL: evapSwitch blocks all condensation

**File:** `HardtWondra.C:522-543`

```cpp
evapSwitch = 0.5*(1 + tanh((TSense - T_sat_)/deltaT))
qn_ = qnRaw * evapSwitch
```

For T < Tsat: `evapSwitch -> 0`, so `qn_ = 0`. Combined with the Helmholtz `max(mdotRaw_,0)` clamp (Issue 3), `mdot_ >= 0` always. The `Q_cond = min(Q_pc,0)` branch in `TEqn.H` is permanently dead code for HardtWondra.

---

#### Issue 3 — MEDIUM: Helmholtz source is one-sided

**File:** `HardtWondra.C:609`

```cpp
mdotSource = max(mdotRaw_, 0)
```

Reinforces Issue 2: even if qn_ were allowed to go negative (condensation), the Helmholtz solve would zero out the negative source. Must be changed together with Issue 2.

---

#### Issue 4 — MINOR: interfaceMaskF for phiStefan is unbounded

**File:** `HardtWondra.C:764-771`

```cpp
surfaceScalarField interfaceMaskF(16.0*alphaIf*(1.0 - alphaIf));
interfaceMaskF = sqr(interfaceMaskF);   // max = 16 at alpha=0.5
```

`(16*alpha*(1-alpha))^2` reaches 16 at alpha=0.5. This is not a mask in [0,1]; it amplifies phiStefan by up to 16x. A normalized version: `sqr(4.0*alphaIf*(1.0 - alphaIf))` has max=1.

---

#### Issue 5 — MINOR: PhaseFractionSource switch has no effect

**File:** `alphaSuSp.H:8-16`

`thermalPhaseChangeModel::alpha1Gen()` returns `-Q_pc/(rho1*h_lv)` when toggled. But `Su` in `alphaSuSp.H` is hardwired to zero — `alpha1Gen()` is never called. The `PhaseFractionSource true/false` entry in `phaseChangeProperties` does nothing.

Note: this may be intentional. phiStefan + PCV together handle alpha transport; adding Su=alpha1Gen would triple-count. But the dead switch is misleading.

---

#### Issue 6 — MINOR: harmonicBias inconsistency

| Location | Uses harmonicBias? | Uses alpha smoothed? |
|---|---|---|
| `calcQ_pc()` kEff | yes (`harmonicBias_` param) | alphaGeom (2x smooth) |
| `kappaf()` for TEqn | no (hardcoded 1.0) | alpha1_ (unsmoothed) |

Minor inconsistency between the Stefan flux kEff and the TEqn thermal diffusion kappaf.

---

#### Issue 7 — MINOR: StefanEnergyJump kEff geometry-dependent (separate model)

**File:** `StefanEnergyJump.C:351`

```cpp
kEff = pos(superheat)*k_liq + neg(superheat)*k_vap
```

Selects by temperature sign, not by which phase conducts. For wall|vapor|liquid evaporation: T>Tsat in vapor -> uses k_liq (should use k_vap). StefanEnergyJump does not have Issue 1 (uses `mag(gradT & nHat)`) but has this kEff issue.

---

### 5. Mutual Consistency Table

| Component | Evaporation | Condensation | Geometry-independent? |
|---|---|---|---|
| `qnRaw = max(kEff*gradT·nHat, 0)` | wall\|liq\|vap only | blocked | **No** |
| `evapSwitch` | T>Tsat only | blocked | **No** |
| `mdotSource = max(mdotRaw_,0)` | positive only | blocked | **No** |
| `Q_pc = mdot*h_lv` | > 0 | = 0 (dead) | — |
| `Acoeff` in TEqn | correct | 0 (correct) | if Q_pc correct |
| `Q_cond = min(Q_pc,0)` | = 0 | = 0 (dead) | if Q_pc correct |
| `PCV` in pEqn | correct sign | = 0 (dead) | if Q_pc correct |
| `phiStefan` direction | correct | = 0 (dead) | **No** |

**TEqn, pEqn, alphaSuSp, alphaEqn downstream logic are all correct** — they just need a properly signed Q_pc delivered by the model. The fix is entirely within `HardtWondra.C`.

---

## Part II — Staged Refactoring Plan

> **STOP: Each phase requires explicit approval before execution.**

### File Map

| File | Phase | What changes |
|---|---|---|
| `HardtWondra/HardtWondra.C:488` | 1 | Replace `max(kEff*(gradT&nHat),0)` with `mag(...)` |
| `HardtWondra/HardtWondra.C:522-543` | 1+2 | Add condSwitch; sign qn_ by (evapSwitch-condSwitch) |
| `HardtWondra/HardtWondra.C:609` | 2 | Remove `max(mdotRaw_,0)` -> allow negative |
| `HardtWondra/HardtWondra.C:764-771` | 5 (opt) | Normalize interfaceMaskF to [0,1] |
| `alphaSuSp.H` | 4 | Document or remove PhaseFractionSource dead code |

---

### Phase 1: Fix qn for Wall|Vapor|Liquid (Evaporation Only)

**Scope:** Replace `max(kEff*(gradT&nHat),0)` with `mag(kEff*(gradT&nHat))` and drive sign via evapSwitch. No condensation yet.

---

#### Task 1.1: Confirm current zero-qnRaw symptom

**Files:** `log.interTempFoam`

- [ ] **Step 1: Check max(mdotRaw) in existing log**

```bash
grep "max(mdotRaw)" /home/param/OpenFOAM/param-v2412/run/1dstefan/log.interTempFoam | head -10
```

Expected output: lines like `max(mdotRaw) = 0` or very small values confirming Issue 1.

- [ ] **Step 2: Verify interface orientation**

```bash
grep "max(Ai)" /home/param/OpenFOAM/param-v2412/run/1dstefan/log.interTempFoam | head -5
```

Nonzero `max(Ai)` confirms the interface is detected but mdotRaw is still zero — consistent with qnRaw=0 due to sign mismatch.

- [ ] **Step 3: Commit analysis note (no code changes)**

```bash
echo "Analysis confirmed: qnRaw=0 for Zone A (wall|vapor|liquid) due to anti-aligned gradT and nHat."
```

---

#### Task 1.2: Replace qnRaw with magnitude formulation

**Files:**
- Modify: `/home/param/interTempFoam/src/transportModels/thermalPhaseChangeModels/HardtWondra/HardtWondra.C`

- [ ] **Step 1: Read existing qnRaw block**

Read lines 477–545 of `HardtWondra.C`:

```bash
sed -n '477,545p' /home/param/interTempFoam/src/transportModels/thermalPhaseChangeModels/HardtWondra/HardtWondra.C
```

- [ ] **Step 2: Implement the corrected qn block**

Replace the block from `volScalarField qnRaw(` through `qn_ = qnRaw*evapSwitch;` with:

```cpp
    //==============================================================
    // Signed normal heat flux: kEff*(gradT & nHat)
    //   > 0: heat flows in nHat direction (wall|liquid|vapor geometry)
    //   < 0: heat flows against nHat (wall|vapor|liquid geometry)
    // Use magnitude for geometry-independent Stefan flux amplitude.
    // Sign of mdot (evap vs cond) is determined by superheat sign.
    //==============================================================

    const volScalarField qnSigned
    (
        IOobject
        (
            "qnSigned_HW",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        kEff*(gradT & nHat)
    );

    const volScalarField qnMag
    (
        IOobject
        (
            "qnMag_HW",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mag(qnSigned)
    );

    const dimensionedScalar deltaT_sw("deltaT_sw", dimTemperature, 0.25);

    const volScalarField evapSwitch
    (
        IOobject
        (
            "evapSwitch_HW",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        0.5*(scalar(1) + tanh((TSense - T_sat_)/deltaT_sw))
    );

    // Phase 1: evaporation-only (condensation arm added in Phase 2)
    qn_ = qnMag * evapSwitch;
```

- [ ] **Step 3: Compile**

```bash
cd /home/param/interTempFoam/src/transportModels/thermalPhaseChangeModels
wmake 2>&1 | tail -15
```

Expected: clean compile ending with `.so` path.

- [ ] **Step 4: Run 1D Stefan case and verify mdotRaw nonzero**

```bash
cd /home/param/OpenFOAM/param-v2412/run/1dstefan
cp -r 0.orig 0
interTempFoam > log.phase1test 2>&1 &
grep "max(mdotRaw)" log.phase1test | head -10
```

Expected: `max(mdotRaw)` now > 0 for Zone A evaporation.

- [ ] **Step 5: Commit**

```bash
cd /home/param/interTempFoam
git add src/transportModels/thermalPhaseChangeModels/HardtWondra/HardtWondra.C
git commit -m "fix(HardtWondra): use mag(kEff*gradT.nHat) for geometry-independent Stefan flux

The previous max(kEff*(gradT & nHat), 0) produced zero for wall|vapor|liquid
because nHat and gradT are anti-aligned in that configuration.
Using the magnitude makes the flux amplitude geometry-independent.
Sign of evaporation vs condensation is driven by evapSwitch (T vs Tsat).

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
```

---

### Phase 2: Add Condensation Arm

**Scope:** Enable condensation (T < Tsat -> mdot < 0 -> Q_pc < 0). Requires signed qn, signed mdotSource in Helmholtz, and sign propagation through to Q_pc.

**Prerequisite:** Phase 1 validated (mdotRaw nonzero for Zone A, solver stable).

---

#### Task 2.1: Add condSwitch and signed mdotSource

**Files:**
- Modify: `HardtWondra.C` — qn_ signing and Helmholtz source

- [ ] **Step 1: Read Phase 1 qn_ block**

```bash
grep -n "evapSwitch\|qn_\|condSwitch" /home/param/interTempFoam/src/transportModels/thermalPhaseChangeModels/HardtWondra/HardtWondra.C
```

- [ ] **Step 2: Add condSwitch and sign qn_**

After the `evapSwitch` definition, add:

```cpp
    const volScalarField condSwitch
    (
        IOobject
        (
            "condSwitch_HW",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        0.5*(scalar(1) - tanh((TSense - T_sat_)/deltaT_sw))
    );

    // Signed qn_: positive for evaporation, negative for condensation
    // evapSwitch + condSwitch = 1 everywhere, so this = qnMag * sign(T-Tsat)
    qn_ = qnMag*(evapSwitch - condSwitch);
```

- [ ] **Step 3: Allow negative mdotSource in Helmholtz**

Find and change the line:

```cpp
    // BEFORE (line ~609):
    mdotSource = max(mdotRaw_, dimensionedScalar("zero", mdotRaw_.dimensions(), 0));

    // AFTER:
    // Allow negative source for condensation; Helmholtz solve propagates sign
    volScalarField mdotSource
    (
        IOobject("mdotSource", mesh_.time().timeName(), mesh_,
                  IOobject::NO_READ, IOobject::NO_WRITE),
        mdotRaw_   // signed: positive=evap, negative=cond
    );
```

- [ ] **Step 4: Verify Q_pc_ = mdot_ * h_lv_ (no change needed, sign propagates)**

```bash
grep -n "Q_pc_ = " /home/param/interTempFoam/src/transportModels/thermalPhaseChangeModels/HardtWondra/HardtWondra.C
```

Expected: `Q_pc_ = mdot_ * h_lv_` — already correct, sign carried through.

- [ ] **Step 5: Verify TEqn condensation path**

```bash
grep -n "Q_cond\|Acoeff" /home/param/interTempFoam/applications/solvers/multiphase/interTempFoam/interTempFoam/TEqn.H
```

Confirm:
- `Q_cond = min(Q_pc_now, 0)` — will be negative when Q_pc < 0 ✓
- `-Q_cond` enters TEqn RHS as positive heat source for condensation ✓
- `Acoeff = max(Q_pc/max(T-Tsat, eps_T), 0)` — for condensation: Q_pc<0, T<Tsat so T-Tsat<0, max(T-Tsat,eps_T)=eps_T>0, Q_pc/eps_T<0, max(...)=0 ✓

No changes needed to TEqn.H.

- [ ] **Step 6: Compile**

```bash
cd /home/param/interTempFoam/src/transportModels/thermalPhaseChangeModels
wmake 2>&1 | tail -10
```

- [ ] **Step 7: Run and check condensation diagnostic**

```bash
cd /home/param/OpenFOAM/param-v2412/run/1dstefan
cp -r 0.orig 0
interTempFoam > log.phase2test 2>&1 &
grep "min(mdot)" log.phase2test | head -10
```

Expected: `min(mdot)` shows negative values in condensation regions (right face of Zone B bubble).

- [ ] **Step 8: Commit**

```bash
cd /home/param/interTempFoam
git add src/transportModels/thermalPhaseChangeModels/HardtWondra/HardtWondra.C
git commit -m "feat(HardtWondra): add condensation arm with signed mdot/Q_pc

Enable condensation by:
- Adding condSwitch = 0.5*(1-tanh((T-Tsat)/0.25))
- Signing qn_ = qnMag*(evapSwitch - condSwitch)
- Removing max(mdotRaw_,0) clamp from Helmholtz source

Q_pc < 0 now activates the existing Q_cond = min(Q_pc,0) path
in TEqn.H (previously dead code for HardtWondra).

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
```

---

### Phase 3: kEff Phase Selection (Conductivity Correction)

**Scope:** Use k_vap for wall|vapor|liquid, k_liq for wall|liquid|vapor. Currently the harmonic mean overestimates by ~kHarm/k_vap ≈ 2x for wall|vapor|liquid.

**Prerequisite:** Phase 2 validated stable.

---

#### Task 3.1: Phase-selective kEff based on heat flow direction

**File:** `HardtWondra.C` — kEff computation block

- [ ] **Step 1: Understand the impact**

For current case: `k_liq=0.6`, `k_vap=0.05`, `harmonicBias=1`.
At alpha=0.5: `kHarm = 2*0.6*0.05/(0.6+0.05) = 0.092 W/mK`
For wall|vapor|liquid: k_vap should be used → qnMag reduced by `0.05/0.092 = 0.54x`
Phase-change rate reduced ~2x for Zone A. Adjust `RelaxFac` accordingly.

- [ ] **Step 2: Implement phase-selective kEff**

After computing `qnSigned` (from Phase 1 result), replace the kEff computation:

```cpp
    // Phase-selective conductivity: based on direction of heat flow at interface
    //   qnSigned > 0: heat flows in nHat direction (toward liquid) -> liquid drives evap -> k_liq
    //   qnSigned < 0: heat flows against nHat (toward vapor) -> vapor drives evap -> k_vap
    const volScalarField kEffPhysical
    (
        IOobject
        (
            "kEffPhysical_HW",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pos(qnSigned)*k_liq_ + (scalar(1) - pos(qnSigned))*k_vap_
    );

    // Recompute qnMag with physically correct conductivity
    const volScalarField qnMag
    (
        IOobject("qnMag_HW", mesh_.time().timeName(), mesh_,
                  IOobject::NO_READ, IOobject::NO_WRITE),
        kEffPhysical * mag(gradT & nHat)
    );
```

Note: `pos(qnSigned)` returns 1 when qnSigned > 0 (wall|liq|vap), 0 otherwise.

- [ ] **Step 3: Compile, run, compare Zone A mdot with Phase 2 baseline**

```bash
cd /home/param/interTempFoam/src/transportModels/thermalPhaseChangeModels
wmake 2>&1 | tail -5

cd /home/param/OpenFOAM/param-v2412/run/1dstefan
cp -r 0.orig 0
interTempFoam > log.phase3test 2>&1 &
grep "max(mdot)" log.phase3test | head -10
```

Expected: `max(mdot)` reduced ~0.5x compared to Phase 2 for Zone A.

- [ ] **Step 4: Commit**

```bash
cd /home/param/interTempFoam
git add src/transportModels/thermalPhaseChangeModels/HardtWondra/HardtWondra.C
git commit -m "fix(HardtWondra): select kEff based on heat flow direction at interface

Use k_liq when heat flows toward liquid (wall|liq|vap, qnSigned>0)
and k_vap when heat flows toward vapor (wall|vap|liq, qnSigned<0).
Harmonic mean was used previously, overestimating Stefan flux by ~2x
for wall|vapor|liquid with k_vap=0.05 vs k_liq=0.6.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
```

---

### Phase 4: Remove PhaseFractionSource Dead Code

**Scope:** Document that alpha1Gen() is intentionally unused (phiStefan+PCV handle alpha transport) and remove or comment the misleading switch.

---

#### Task 4.1: Verify alpha conservation and clean up

- [ ] **Step 1: Check alpha conservation in Phase 2 run**

```bash
grep "Phase-1 volume fraction" /home/param/OpenFOAM/param-v2412/run/1dstefan/log.phase2test | awk '{print $NF}' | head -20
```

If the weighted average drifts significantly, phiStefan+PCV may have a conservation gap and Su=alpha1Gen could be needed. If stable, proceed with documentation-only fix.

- [ ] **Step 2: Add comment to alphaSuSp.H clarifying the design**

```cpp
// alphaSuSp.H — Phase-change sources for MULES alpha equation
//
// Design note: alpha1Gen() is NOT used here (Su=0).
// Interface kinematics are handled by phiStefan (in phiCN = phi+phiStefan).
// Volume conservation is enforced by divU=PCV in fvm::Sp(divU,alpha1).
// Adding Su=alpha1Gen on top would triple-count the phase-change effect.
// The PhaseFractionSource dictionary switch is therefore a no-op and
// can be safely removed from phaseChangeProperties.
```

- [ ] **Step 3: Commit**

```bash
cd /home/param/interTempFoam
git add applications/solvers/multiphase/interTempFoam/interTempFoam/alphaSuSp.H
git commit -m "docs(alphaSuSp): clarify PhaseFractionSource is intentionally unused

alpha1Gen() is not called because phiStefan+PCV already handles
the full alpha transport. Adding Su=alpha1Gen would double/triple-count.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
```

---

### Phase 5: Normalize interfaceMaskF (Optional)

**Scope:** Fix `(16*alpha*(1-alpha))^2` (max=16) to `(4*alpha*(1-alpha))^2` (max=1).

**Prerequisite:** Phase 3 complete. Requires RelaxFac retuning after normalization since Stefan flux amplitude changes.

- [ ] **Step 1: Change formula in HardtWondra.C lines 764-771**

```cpp
// BEFORE:
surfaceScalarField interfaceMaskF(16.0*alphaIf*(scalar(1.0) - alphaIf));
interfaceMaskF = sqr(interfaceMaskF);   // max = 16 at alpha=0.5

// AFTER (normalized to [0,1]):
surfaceScalarField interfaceMaskF(sqr(4.0*alphaIf*(scalar(1.0) - alphaIf)));
// max = (4*0.25)^2 = 1 at alpha=0.5
```

- [ ] **Step 2: Retune RelaxFac by factor of ~4 to maintain Stefan flow rate**

In `constant/phaseChangeProperties`:
```
RelaxFac   3.2;  // was 0.8, increased 4x to compensate 16->1 amplitude reduction
```

- [ ] **Step 3: Compile, run, validate phiStefan magnitude similar to baseline**

```bash
wmake && interTempFoam > log.phase5test 2>&1 &
grep "max(mdot)" log.phase5test | head -5
```

- [ ] **Step 4: Commit**

```bash
git commit -m "fix(HardtWondra): normalize interfaceMaskF for phiStefan to [0,1]

Previous formula (16*alpha*(1-alpha))^2 peaked at 16 at alpha=0.5.
Corrected to (4*alpha*(1-alpha))^2 which peaks at 1.
RelaxFac adjusted 4x to maintain equivalent Stefan flow rate.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
```

---

## Part III — Validation Against 1D Stefan Analytical Solution

**Stefan-Neumann analytical interface position:**
```
s(t) = 2*beta*sqrt(kappa_l * t)
kappa_l = k_liq / (rho_liq * cp_liq) = 0.6 / (1000*4216) = 1.42e-7 m2/s
Stefan number St = cp_l * DeltaT / h_lv = 4216*10 / 1e6 = 0.04216
beta ~ sqrt(St/pi) ~ 0.116   [for St << 1]
```

| Phase | Validation check | Pass criterion |
|---|---|---|
| 1 | `max(mdotRaw)` > 0 for Zone A | Any nonzero value |
| 1 | T in vapor stays within 2 K of Tsat | min(T) > 371 K |
| 2 | `min(mdot)` < 0 in Zone B condensation region | Negative value appears |
| 2 | Zone B vapor bubble shrinks from right side | alpha integral decreasing for Zone B |
| 3 | Zone A phase-change rate ~2x lower than Phase 1 | max(mdot) reduced |
| 1+2+3 | Interface position vs sqrt(t) | Approximately linear slope |

---

## Appendix: Sign Convention Reference

```
Evaporation (liquid -> vapor, T > Tsat):
  Q_pc > 0, mdot > 0, PCV > 0, alpha1Gen < 0
  qn_ > 0, evapSwitch ~ 1

Condensation (vapor -> liquid, T < Tsat):
  Q_pc < 0, mdot < 0, PCV < 0, alpha1Gen > 0
  qn_ < 0, condSwitch ~ 1

nHat = grad(alpha)/|grad(alpha)|: always points vapor -> liquid
  wall|vapor|liquid: nHat = +x  (away from wall)
  wall|liquid|vapor: nHat = -x  (toward wall)

gradT (hot wall on left): always ~ -x direction

gradT & nHat:
  wall|vapor|liquid: negative  <- was the bug
  wall|liquid|vapor: positive  <- worked before
```
