# interTempFoam -- Publication-Grade Technical Audit

Diffuse-interface Stefan VOF solver with Hardt-Wondra (HW) Helmholtz
redistribution. Audit covers `interTempFoam.C`, `createFields.H`,
`alphaEqn{,SubCycle}.H`, `alphaSuSp.H`, `UEqn.H`, `pEqn.H`, `TEqn.H`,
`thermalPhaseChangeModel{.H,.C}`, `HardtWondra{.H,.C}`.

Audit date: 2026-05-27. Audited against publication scrutiny.

---

## 1. Architecture overview

### 1.1 PIMPLE coupling order (interTempFoam.C)

```
while pimple.loop():
  (firstIter) mesh.update();
  (firstIter) phaseChangePtr->correct();          // builds qn, mdotRaw, mdot, Q_pc, Q_pc_thermal, phiStefan, Qcorr
  phiTotal = phi + phiStefan;                     // shared by alphaEqn
  alphaEqnSubCycle  -> alphaEqn (MULES, phiTotal)  // VOF transport
  mixture.correct();                              // rho, mu from new alpha
  UEqn (no Stefan body force; phiStefan enters only via alphaEqn and pEqn)
  while pimple.correct():
    pEqn (Laplacian == div(phiHbyA) - PCV)        // PCV = volume-dilatation source
  rho, rhoCp recomputed
  TEqn (rhoCpPhi=phi*rhoCp_f, NOT phi+phiStefan)        // energy
```

### 1.2 The three HW physics pathways

| Pathway | Field used | Where injected | Purpose |
|---|---|---|---|
| **Kinematic** | `mdotRaw_` (then `Ustef = mdotRaw/(Ai+AiFloor_kin)/rho_int`) | `phiStefan -> alphaEqn.phiTotal` | Pure interface translation |
| **Thermodynamic (latent)** | `Q_pc_ = mdotRaw_ * h_lv` (and `Q_pc_thermal_` for TEqn implicit penalty) | TEqn `Acoeff` (implicit) + `Qcorr` (explicit) | Latent heat coupling |
| **Dynamic (continuity dilatation)** | `PCV = Q_pc/h_lv * (1/rho_2-1/rho_1)` | pEqn RHS | Mass-conserved volume source |

The architecture is **partially Helmholtz-redistributed** -- Helmholtz `mdot_`
exists, is solved, but is *no longer load-bearing*. See Sec.3.

---

## 2. Coupling / source-term dependency map

```
T ---+
    +-->  TSense=max(T,Tsat)  ->  gradT  (DEAD, see Sec.4.1)
    |
alpha --+->  alphaSmooth  ->  alphaGeom  ->  gradalpha  ->  Ai = |gradalpha|
    |                                 ->  nHat  -> nHatSmooth
    |                              ->  alphaThermal = 1/2(1+tanh 4(alpha-1/2))
    |                                 -> kEff (harmonic, biased)
    |
T,alpha -+->  cellCells one-sided search  ->  dTdn_hot      (kinematic, floor=kinFloorCells*h)
                                     ->  dTdn_hot_thermal (thermal, floor=h_ref)

dTdn_hot       ->  qn       =  kEff*dTdn_hot * interfaceMask
                  mdotRaw  =  Ai*qn / h_lv
                  +---> Helmholtz solve ->  mdot_   (used only by Qcorr & diagnostics)
                  +---> mdotInterfacial = mdotRaw/(Ai+AiFloor_kin)
                  |     -> Ustef = mdotInt_f / rho_int   -> phiStefan = Ustef*(nHat_f*Sf)*maskF
                  +---> Q_pc = mdotRaw*h_lv            -> TEqn Acoeff (via Q_pc_thermal), PCV, alpha1Gen (unused)
                  +---> (Q_pc_thermal feeds TEqn Acoeff implicit penalty only)

mdot_ (Helmholtz)-> Qcorr = -mdot*(cp_1-cp_2)*(T-Tsat) -> TEqn RHS
phiStefan        -> alphaEqn.phiTotal
PCV              -> pEqn RHS
```

---

## 3. Critical inconsistency -- Helmholtz redistribution is partially dead

Hardt & Wondra (2008) constructs `mdot = mdotRaw - lambda^2grad^2mdot` so that the
*volumetric latent source* `Q_pc = mdot*h_lv` is smooth and spreads into bulk
cells. Eq. (42) of HW (the `Qcorr` enthalpy correction) is **mathematically
derived assuming Q_pc uses the redistributed `mdot`**, because Qcorr removes the
artificial enthalpy carried by mass deposited outside the physical phase-change
band.

Current code (`HardtWondra.C:1122`):

```cpp
Q_pc_ = mdotRaw_ * h_lv_;                                          // localized
...
Qcorr_ = -mdot_ * (mixture_.cp1() - mixture_.cp2()) * (T_ - T_sat_);   // redistributed
```

This is **internally inconsistent**:

- If `Q_pc` is localized (uses `mdotRaw_`), then by HW's derivation `Qcorr == 0`.
  The current `Qcorr_` non-trivially modifies the energy budget by ~1-2 % using
  a mass field (`mdot_`) that does not appear elsewhere in the energy equation.
- The Helmholtz solve (`fvm::laplacian(lambda^2)`) is now load-bearing only through
  `Qcorr_`. If `useEnthalpyCorrection=false`, the entire Helmholtz solve is
  dead computational cost.

### 3.1 Recommended resolution

**Pick one consistent variant. The cleanest defensible choice is "localized
HW"** -- drop Helmholtz entirely and remove `Qcorr` (set to zero or remove).
Helmholtz redistribution was historically introduced to broaden the latent
sink for stability; the current code already broadens via the diffuse VoF
band itself, the `alphaThermal` tanh sharpening, and the Acoeff hot/cold
floors in `TEqn.H`. Removing `mdot_`/`Qcorr_` simplifies the model, eliminates
the inconsistency, removes one elliptic solve per outer iteration, and aligns
with the existing kinematic + thermodynamic + dynamic pathways which all
already use `mdotRaw_`.

Alternative (faithful HW): use `Q_pc = mdot*h_lv` (redistributed), keep
`Qcorr`, but then *also* use `mdot` (not `mdotRaw_`) in `phiStefan` and
`PCV`. This restores HW-consistency at the cost of reintroducing the
~25 % interface-position underprediction documented in the existing
in-source comments (`HardtWondra.C:1001-1014`).

---

## 4. qn / one-sided gradient reconstruction audit

### 4.1 Dead code path -- TSense

`HardtWondra.C:485-525` constructs a smoothed sensing field `TSense` and
computes `gradT = fvc::grad(TSense)`. **Neither `TSense` nor `gradT` is used
downstream.** `qn_` is assembled at `:786` from `qnSigned = kEff*dTdn_hot`
where `dTdn_hot` comes from the manual `cellCells()` neighbor search at
`:702-765`. The TSense diffusion loop (one fvc::laplacian per iteration) and
the `fvc::grad(TSense)` are pure runtime cost without effect.

**Action**: delete `TSense`, `D_Tsense`, the TSense smoothing loop, and the
`fvc::grad(TSense)` block.

### 4.2 The neighbor-search reconstruction (`:702-765`)

Logic:
- For every interface cell (`Ai > SMALL`):
  - scan face-neighbors via `cellCells()`,
  - reject neighbors with `normalProj = |dVec*nHat| < SMALL` (tangential),
  - reject neighbors with `T_nbr - Tsat < 1e-6` (cold side),
  - among survivors, **select the closest** (`if (normalProj < bestDistance)`),
  - compute `grad = (T_nbr - Tsat) / max(normalProj, floor*h_ref)`.

#### 4.2.1 Strengths

- Genuinely one-sided: excludes cold-side neighbors, so the cold-phase Tsat
  pin no longer corrupts the gradient. This is the correct response to
  diffuse-interface "both phases contaminate gradT" issue.
- Direct Dirichlet anchor at `Tsat`: gradient is `(T_hot - Tsat) / d_hot`,
  matching the sharp Stefan condition at the interface plane.
- Floor on the denominator prevents divide-by-tiny-distance singularities.

#### 4.2.2 Weaknesses (publication-critical)

| Weakness | Impact | Severity |
|---|---|---|
| **W1.** `cellCells()` returns face-neighbors only (no corner). On hexahedral meshes the search has 6 candidates; on poorly-aligned interfaces the optimal hot neighbor in the `nHat` direction may not be a face neighbor. | Direction-dependent qn bias, anisotropic Stefan speed | High |
| **W2.** Selection criterion is **closest hot neighbor**, not best-aligned. A near-tangential face-neighbor with small `normalProj` can beat a well-aligned but farther neighbor. | Spurious selection switching -> noise in qn | High |
| **W3.** Single-neighbor stencil: the gradient is a 2-point one-sided difference. Truncation error O(h). No least-squares fit, no multi-cell extrapolation. | First-order accuracy in qn -> first-order accuracy in interface speed | Critical |
| **W4.** `kinFloorCells_*h_ref` floor (default 2.0) **artificially weakens** the gradient when the chosen neighbor is closer than 2 cells. On well-resolved meshes this systematically *under*-predicts qn at exactly the cells doing the work. | Documented 25 % interface-position underprediction (per in-source comments) | Critical |
| **W5.** Boundary cells: a hot neighbor across a wall patch is excluded (cellCells does not cross patches). If hot phase abuts a wall, no qn is reconstructed there. | Wall-attached interfaces have artificially zero Stefan velocity | High |
| **W6.** Uses `T_` directly, not a smoothed field. Combined with VOF staircase and the saturation pin from the anti-trench floor, gradient can flip sign cell-to-cell. | Noise -> mdotRaw oscillation -> Helmholtz solve noise -> Qcorr noise | Medium |
| **W7.** The "10x floor" comment is stale; default `kinFloorCells_ = 2.0` is hardcoded but the comments describe the original 10x behavior. | Documentation drift; defensible parameter choice obscured. | Low (but publication-blocking) |
| **W8.** `Q_pc_thermal_` uses *the same* neighbor selection but with a different (smaller) floor. Two latent sources differing only in a magic denominator floor have no physical interpretation. | Cannot justify in publication. | Critical |

### 4.3 Recommended qn reconstruction

A defensible diffuse-interface one-sided gradient with O(h^2)-leaning behavior
that remains pure FV:

**Algorithm (per interface cell c):**

1. Build a hot-phase weight `w_n = alpha_hot[n] * pos(T[n] - Tsat - eps)` for each
   face neighbor *n*. `alpha_hot` is `alpha` or `1-alpha` depending on which side is hot.
2. Solve a **least-squares hot-side gradient** using only weighted hot
   neighbors:
   ```
   minimize Sigma_n w_n * ( (T_n - T_c) - gradT_LS * (x_n - x_c) )^2
   ```
3. Project: `dTdn = gradT_LS * nHat_c`, then clip `dTdn = max(dTdn, 0)` (cold-side
   bleed-through removed by the weights, but a safety floor remains).
4. Compute qn as `qn = kEff * dTdn`. The kEff at the interface cell is
   replaced by the *hot-side* k (k_liq for melting, k_vap for evaporation),
   not the harmonic mean -- physically the conductive flux is being carried by
   the hot phase up to the interface plane.
5. Replace the magic-number floor on `normalProj` with a tiny dimensionful
   regularizer `sqrt(h_ref^2 + (nHat*gradalpha/|gradalpha|)-^2)` that goes to zero on aligned
   interfaces -- eliminates W4 entirely.

OpenFOAM-compatible pseudocode (single-pass, drop-in for `:702-765`):

```cpp
const scalar TsatVal = T_sat_.value();
forAll(mesh_.cells(), cellI)
{
    if (interfaceArea_[cellI] < SMALL) continue;
    const vector& nh = nHat[cellI];
    const labelList& nbr = mesh_.cellCells()[cellI];
    // Hot-side weighted LS gradient
    symmTensor M = Zero;  vector b = Zero;
    forAll(nbr, k) {
        label j = nbr[k];
        scalar w = max(T_[j] - TsatVal, 0.0);     // hot-only weight
        if (w < SMALL) continue;
        vector d = mesh_.C()[j] - mesh_.C()[cellI];
        M += w * sqr(d);                          // symmetric outer product
        b += w * (T_[j] - TsatVal) * d;           // anchor: T_iface = Tsat
    }
    if (det(M) > SMALL) {
        vector gT = inv(M) & b;
        dTdn_hot[cellI] = max(gT & nh, 0.0);
    }
}
```

**Benefits**:
- Second-order via LS on whatever neighbor set is available, including more
  than one hot neighbor.
- Single floor (`det(M) > SMALL`) -- no magic mesh-units constant.
- Single gradient field, no thermal/kinematic split; eliminates the entire
  `Q_pc_thermal_` apparatus.
- Wall-adjacent cells degrade gracefully: weight goes to zero if all hot
  neighbors are wall-blocked, qn -> 0 with no spurious activation.

### 4.4 kEff inconsistency

`HardtWondra.C:436-459` (used for qn) uses a *biased* harmonic mean with
`harmonicBias_=0.3` over `alphaThermal = 1/2(1+tanh 4(alpha-1/2))`.
`kappaf()` (used for TEqn Laplacian, `:1224-1262`) uses an *unbiased*
harmonic mean over raw `alpha1_`. Two different effective conductivities for
the same physics in the same solve. **Recommendation**: collapse to a single
kEff field (raw harmonic over `alphaThermal`) and reuse it for both qn and
TEqn. If sharper thermal interface is required, sharpen `alphaThermal`, not
introduce a second field. Document `alphaThermal` sharpening with a numerical
truncation analysis in the paper.

---

## 5. TEqn audit (`TEqn.H`)

### 5.1 Current structure

```
rhoCp*dT/dt + grad*(rhoCp phi T) - grad*(k_f gradT) + (A_main + A_floor + A_hot) T
  = (A_main + A_floor + A_hot) * Tsat - Q_cond + Qcorr
```

with three implicit penalties:

| Coeff | Form | Active where | Purpose |
|---|---|---|---|
| `Acoeff` | `1/2*max(Q_pc_thermal/(T-Tsat+2.5), 0)` | T > Tsat in interface band | Latent sink (Lee-style hidden inside HW) |
| `Acoeff_floor` | `5*rhoCp/Deltat * 4alpha(1-alpha) * pos(Tsat-T)` | Diffuse band, T<Tsat | Anti-trench (cold-side floor) |
| `Acoeff_hot_floor` | `5*rhoCp/Deltat * 4alpha(1-alpha) * pos(T-Tsat) * coldSideGate` | Diffuse band, T>Tsat, cold-phase side | Anti-superheat (hot-side floor) |

### 5.2 Findings

**F1 (Critical, mathematical).** `Acoeff = Q_pc_thermal/(T-Tsat+eps_T)` re-derives
a Lee model from HW data: in cells where Acoeff is active, the implicit
penalty drags T -> Tsat at rate `Acoeff/rhoCp`. Because `Q_pc_thermal` itself
depends on `T_nbr - Tsat` (via `dTdn_hot_thermal`), the construction is
self-consistent *only when solved tightly*. PIMPLE outer-loop convergence is
not enforced; with a single outer iteration the system is effectively
explicit. **Recommendation**: either (i) make `Q_pc_thermal` explicit and
treat the latent term purely as RHS source (no Acoeff), or (ii) iterate the
phase-change update inside the PIMPLE outer loop (call `correct()` every
outer iter, not only `firstIter()` -- see Sec.6.2).

**F2 (Critical, conservation).** `Q_cond = min(Q_pc_now, 0)` peels off the
condensation branch as explicit, while `Acoeff*Tsat` provides the
evaporation branch implicitly. The split is asymmetric: condensation gets
no implicit pull-to-Tsat. For pure evaporation cases (Sheikh Stefan) Q_cond
is identically zero and this branch is dead. Document or remove for
single-direction studies.

**F3 (Critical, defensibility).** `Acoeff_floor` and `Acoeff_hot_floor` are
**numerical bandaids** for residual inconsistency in qn / mdotRaw. They use
`rhoCp/Deltat` scaling -- i.e. pure time-step-rate penalty -- which is the fingerprint
of an unstable underlying source being suppressed. Both floors should
**not appear** in a clean HW solver. Their presence indicates the qn
reconstruction (Sec.4) is releasing latent heat in cells where the temperature
state then violates Tsat bounds.

The physically defensible fix: improve qn reconstruction (LS hot-side
gradient, Sec.4.3) and let the Acoeff term alone enforce Tsat in the diffuse
band. The floors should be **off by default** and only enabled as an
emergency diagnostic, with a runtime warning when any cell triggers them.

**F4 (Medium, code-cleanup).** The dead `// REPLACE: const volScalarField
coldSideGate(...)` block at TEqn.H:51-58 should be deleted. The
`coldSideGate` constructor uses an awkward immediate `IOobject` +
post-assignment pattern; collapse to a single ternary `volScalarField`
constructor.

**F5 (Medium, conservation).** `rhoCpPhi = phi * rhoCp_f` (no phiStefan).
Comment at `alphaEqn.H:274-278` explicitly justifies this: sensible enthalpy
should not ride the Stefan velocity. **Correct decision**, but creates a
hidden conservation drift: `d(rhoCp)/dt` from `mixture.correct()` after
alphaEqn includes interface motion driven by `phiTotal=phi+phiStefan`, while
`grad*(rhoCpPhi*T) = grad*(phi*rhoCp_f*T)` uses only `phi`. The implicit assumption
is that latent heat from interface motion is captured exactly by `Q_pc` and
`Qcorr`. Add a consistency check function-object that monitors
`intrhoCp*T dV - int(...) flux - intQ_pc + intQ_cond` per timestep.

**F6 (Low, performance).** `Q_pc_now`, `Q_pc_th_now`, `Acoeff`, `Q_cond`,
`indicator_TEqn`, `Acoeff_floor`, `coldSideGate`, `Acoeff_hot_floor` are all
constructed as named `volScalarField`s every TEqn call. With AUTO_WRITE off
this is OK but it's still 6+ field allocations per outer iteration. Move all
to internal-field arithmetic with `volScalarField::Internal` to halve memory
traffic.

### 5.3 Proposed clean TEqn

```cpp
const surfaceScalarField kappaf = phaseChangePtr->kappaf();
const dimensionedScalar Tsat = phaseChangePtr->T_sat();

// Single, explicit latent source. No A_floor, no A_hot, no Q_cond split.
// Latent sink enforced via Acoeff alone, derived from HW Q_pc.
const dimensionedScalar epsDT("epsDT", dimTemperature, scalar(1.0));
const volScalarField::Internal Acoeff
(
    max(phaseChangePtr->Q_pc()()/max(T - Tsat, epsDT),
        dimensionedScalar(dimensionSet(1,-1,-3,-1,0,0,0), 0))
);

fvScalarMatrix TEqn
(
    rhoCp*fvm::ddt(T)
  + fvm::div(rhoCpPhi, T)
  - fvm::laplacian(kappaf, T)
  + fvm::Sp(Acoeff, T)
 ==
    Acoeff*Tsat
  + phaseChangePtr->Qcorr()   // becomes ==0 once Sec.3 cleanup applied
);
TEqn.relax();
TEqn.solve();
```

Anti-trench and anti-superheat floors should be **case-level diagnostics**
(printed when triggered) rather than baked-in stiffness sources. If they are
load-bearing for stability in any production case, that case fails to
demonstrate the HW formulation itself converges -- a publication-level
red flag.

---

## 6. Solver-level (`interTempFoam.C`) audit

### 6.1 `phaseChangePtr->correct()` called only on `pimple.firstIter()`

`interTempFoam.C:152-156`. Implications:
- The phase-change source is frozen across pimple correctors. With a single
  outer iteration (nOuterCorrectors=1) the energy equation uses sources
  computed from `T^n`, not the converging `T^(n+1,k)`. The implicit
  `Acoeff*(T - Tsat)` partially compensates, but does not close the
  T-feedback into qn.
- For research-grade time integration, call `correct()` every outer iter.
  Add a dictionary flag `correctPhaseChangeEveryOuter` (default false for
  backward compat) and recommend `true` in publication runs.

### 6.2 `phiTotal` updated only in outer loop, not in inner

`interTempFoam.C:158`. `phiTotal = phi + phiStefan` is computed before the
alphaEqn solve. After pEqn updates `phi`, `phiTotal` is stale until the next
outer iteration. For interface advection this is fine (alphaEqn already
ran), but it makes the rhoPhi/rhoCpPhi calculations inside alphaEqn use
pre-pEqn `phi`. This is the standard interFoam pattern, but it means the
HW Stefan velocity advecting alpha was paired with an unconverged Darcy velocity.

### 6.3 `rho`, `rhoCp` recomputed manually after pEqn

`interTempFoam.C:180-186` recomputes rho and rhoCp after pEqn solves. Then
`alphaEqnSubCycle.H:44-47` *also* recomputes them after alphaEqn. The post-
pEqn rebuild is the live one for TEqn; the post-alpha rebuild is overwritten.
**Cleanup**: remove the post-alpha rebuild from alphaEqnSubCycle.H, or remove
the post-pEqn rebuild in interTempFoam.C -- keep one. The post-pEqn version
is the right one for TEqn consistency (alpha is finalized).

### 6.4 `pcDictPtr` outer-scope hack

`createFields.H:65-94, 251-257`. The `IOdictionary` is kept alive in
`createFields.H` outer scope only so `TEqn.H` can read `coldPhaseIsHighAlpha1`.
Two clean alternatives:
1. Move `coldPhaseIsHighAlpha1` into `thermalPhaseChangeModel` base class as
   a virtual `bool isColdPhaseHighAlpha1() const` accessor. Eliminates the
   dictionary leak.
2. Read the flag once in `createFields.H` and pass it as a `const bool` --
   already done at `:252`. Then `pcDictPtr` can be a local scope variable.

Option 1 is the correct architectural choice.

---

## 7. alphaEqn / alphaSuSp audit

### 7.1 Su=0 even though `alpha1Gen` is computed

`alphaSuSp.H:8-32` sets `Su` to identically zero, but
`thermalPhaseChangeModel::alpha1Gen()` returns `-Q_pc/(rho_1*h_lv)` -- a
non-zero phase-change alpha-source. The current model relies entirely on
`phiStefan` advection for interface motion; `alpha1Gen` is computed but
*never injected into MULES*. This means:

- `alpha1Gen()` is **dead code** in the current solver. The comment in
  `alphaSuSp.H:1` ("Su = alpha1Gen: positive for condensation...") is false.
- The MULES Su/Sp boundedness machinery is bypassed for the volumetric
  phase-change source. Boundedness is enforced only through `phiStefan`
  geometry + the post-MULES clip `alpha1.max(0)/min(1)`.

This is a defensible design choice (kinematic-only interface motion is more
mesh-convergent than dual-injection Su+phiStefan), but it must be:
- documented in the paper's "model description" section,
- enforced by either (i) deleting `alpha1Gen()` from the base class or
  (ii) renaming to `alpha1GenUnused()` to prevent future re-injection.

### 7.2 phiStefan in flux compression

`alphaEqn.H:91-109`. The compression flux `phiCN = phi + phiStefan` (or
its CN blend). The compressive face flux `phir = phic*nHatf` does **not**
include phiStefan. So:
- Bulk transport: phiTotal = phi + phiStefan OK
- MULES upwind: phiCN = phiTotal OK
- Compression: phir on phi-derived phic, but applied to alpha1 (not alpha2) OK standard interFoam

Consistent.

### 7.3 phiTotal vs phiCN inside MULESCorr

`alphaEqn.H:111-171` (MULESCorr branch) uses `phiCN` (which includes
`phiStefan`) for the upwind matrix. Later, `:180-192` (corrector loop) uses
`phiTotal` for `fvc::flux(phiTotal, alpha1, alphaScheme)`. **Both paths
include phiStefan** but the variable names differ (phiCN vs phiTotal).
Cleanup: unify naming -- pass `phiTotal` everywhere, or always use `phiCN`.

### 7.4 Commented-out divergence subtraction

`alphaEqn.H:128-129`:
```cpp
// - fvm::Sp(fvc::ddt(dimensionedScalar("1", dimless, 1), mesh)
//           + fvc::div(phiCN), alpha1)
```
This is the standard interFoam non-conservative correction for compressible
flux divergence. With `divU = PCV` injected via the `==` source (`:131`),
this term is redundant *if* divU correctly cancels `fvc::div(phiCN) - fvc::div(phi)`.
**Verify by derivation**:

`grad*phiTotal = grad*phi + grad*phiStefan`. From mass continuity with phase change,
`grad*U = (1/rho_v - 1/rho_l) * mdot = PCV`. So
`grad*phi = PCV*V_cell` (in discrete form, when phi is volume flux). With
`mdot = mdotRaw`: `grad*phiStefan ~= |gradalpha|*(Ustef*nHat)*Sf ~= mdotRaw/rho_int`,
which is the **interface-localized** mass-conservation contribution.

The commented block, if uncommented, would double-correct. The current
`divU = PCV` injection is the right choice. Delete the comment.

---

## 8. pEqn / continuity audit

`pEqn.H:43-48`:
```cpp
fvScalarMatrix p_rghEqn
(
    fvm::laplacian(rAUf, p_rgh)
 ==
    fvc::div(phiHbyA) - phaseChangePtr->PCV()
);
```

### 8.1 PCV correctness

`PCV = (Q_pc/h_lv)*(1/rho_v - 1/rho_l) = mdotRaw*(1/rho_v - 1/rho_l) > 0` for
evaporation. This is the correct volumetric source from mass-flux density
difference. Sign convention matches `grad*U = PCV`. OK

### 8.2 phiStefan not in phiHbyA

`phiHbyA` does **not** include `phiStefan`. The Stefan flux only enters
alphaEqn. The dilatation enters here via `PCV`. This is the standard
diffuse-interface decomposition (see Welch & Wilson 2000, Sato & Niceno
2013) and is defensible. The pressure field accommodates the volumetric
source; the interface motion is carried by the kinematic phiStefan in
alphaEqn. Document this explicitly in Sec.3.2 of the paper.

### 8.3 Missing Rhie-Chow consistency check

If any future change adds a body force to UEqn (e.g. EHD coupling for the
user's broader research goal), that body force must also be added to
`phiHbyA`. Add a TODO comment at the phig assembly to remind future
developers.

---

## 9. Empirical regularization audit (full inventory)

| Location | Mechanism | Defensibility | Recommendation |
|---|---|---|---|
| `HardtWondra.C:283` alpha smoothing (nSmoothIter passes) | Pseudo-diffusion of alpha for gradient denoising | Defensible (Hardt & Wondra 2008); document Fo_per_iter cap | Keep, document Fo cap |
| `HardtWondra.C:307` alphaGeom extra smoothing pass | Additional smoothing for normals only | Redundant with alphaSmooth | Merge into single field |
| `HardtWondra.C:359-367` `AiMax = 1.5/h_ref` | Hard clip on |gradalpha| | Magic factor 1.5 | Derive from cell-area to volume ratio |
| `HardtWondra.C:425-434` `alphaThermal = 1/2(1+tanh 4(alpha-1/2))` | tanh sharpening for thermal kEff | Defensible (sharpens thermal interface) | Document; eliminate `harmonicBias_` |
| `HardtWondra.C:188` `harmonicBias_ = 0.3` | Bias toward liquid in qn kEff | **Indefensible** -- duplicate kEff definition (Sec.4.4) | **Remove** |
| `HardtWondra.C:486` `TSense = max(T, Tsat)` | One-sided T clip | Dead code (TSense unused for qn) | **Remove** |
| `HardtWondra.C:506` TSense diffusion (D_Tsense) | Smooth TSense | Dead code | **Remove** |
| `HardtWondra.C:561-578` `evapSwitch = 1/2(1+tanh((T-Tsat)/0.25))` | Sign-of-(T-Tsat) softening | Only used to gate Q_pc_thermal; magic DeltaT=0.25K | Defensible only if Q_pc_thermal kept; eliminate with Sec.4.3 |
| `HardtWondra.C:581-592` `condSwitch` | Cond branch sign | **Dead** -- never used | **Remove** |
| `HardtWondra.C:748` `kinFloorCells_*h_ref` denominator floor | Distance regularizer in qn search | **Indefensible** -- biases qn down on fine meshes | **Replace with LS gradient (Sec.4.3)** |
| `HardtWondra.C:782` `qnSigned *= interfaceMask` | `pos(alpha-0.01)*pos(0.99-alpha)` band mask | Defensible (excludes bulk) | Document threshold sensitivity |
| `HardtWondra.C:919` `mdot_ = clip(+/-mdotMax)` | Hard clip on Helmholtz mdot | Backstop only | Add warning when triggered |
| `HardtWondra.C:925` `mdot_ *= RelaxFac_` | Under-relaxation | Defensible | Keep |
| `HardtWondra.C:982-990` `nHatSmooth` (2 Laplacian passes) | Normal denoising | Defensible | Keep, document |
| `HardtWondra.C:1016-1021` `AiFloor_kin = 5/m` (absolute) | Division-safety in mdotRaw/Ai | Absolute floor, no mesh dependence -- good | Keep |
| `HardtWondra.C:1066-1074` `interfaceMaskF` at faces | Suppress phiStefan in bulk faces | Defensible (zero-suppression in bulk) | Keep |
| `HardtWondra.C:1133-1139` `for(i=0;i<0;...)` Q_pc smoothing | Loop count 0 -- **dead** | Dead code | **Remove** |
| `HardtWondra.C:205,1287` `betaThermal_` parameter | Read from dict, **never used** | Dead parameter | **Remove from .H/.C/dict** |
| `TEqn.H:36-43` `Acoeff_floor` (anti-trench) | Stiff penalty pulling T->Tsat (cold side) | Numerical bandaid (Sec.5.2 F3) | Disable by default; diagnostic only |
| `TEqn.H:76-84` `Acoeff_hot_floor` (anti-superheat) | Stiff penalty pulling T->Tsat (hot side, cold-phase gate) | Numerical bandaid (Sec.5.2 F3) | Disable by default |
| `TEqn.H:51-58` Commented-out `coldSideGate` block | Dead comment | **Remove** |

**Total dead code to remove**: TSense + D_Tsense + TSense smoothing loop + `fvc::grad(TSense)` + condSwitch + Q_pc smoothing loop (count=0) + betaThermal + commented-out coldSideGate block + `harmonicBias_` (and dictionary read) + dead "STEP 10/11" duplicate comment blocks. Conservative estimate: ~120 lines of dead or stale code in HardtWondra.C alone.

---

## 10. Numerical concerns (prioritized)

### 10.1 Critical (publication-blocking)

- **C1**. Helmholtz inconsistency: Q_pc uses mdotRaw, Qcorr uses mdot
  (Sec.3). Pick one consistent variant.
- **C2**. qn reconstruction is first-order, single-neighbor, with a
  mesh-dependent floor that systematically under-predicts gradient (Sec.4.2.2
  W3/W4). Replace with LS hot-side gradient (Sec.4.3).
- **C3**. Two latent fields `Q_pc_` and `Q_pc_thermal_` differing only in a
  magic denominator floor lack physical interpretation (Sec.4.2.2 W8). Collapse
  to one.
- **C4**. `Acoeff_floor` and `Acoeff_hot_floor` are stiffness penalties
  hiding upstream inconsistency. Their existence is evidence of an unstable
  qn (Sec.5.2 F3). Remove after fixing Sec.4.

### 10.2 High impact

- **H1**. `phaseChangePtr->correct()` runs only on `firstIter`. Energy
  equation cannot converge T-source feedback within an outer step (Sec.6.1).
- **H2**. Duplicate kEff definitions (Sec.4.4). Unify.
- **H3**. `alpha1Gen()` is dead infrastructure (Sec.7.1). Remove or document.
- **H4**. Dictionary-lifetime hack (`pcDictPtr` at outer scope) for one
  flag (Sec.6.4). Move flag into model base class.

### 10.3 Medium impact

- **M1**. Inconsistent naming `phiTotal` vs `phiCN` in alphaEqn (Sec.7.3).
- **M2**. Redundant rho/rhoCp rebuild between alphaEqn and post-pEqn (Sec.6.3).
- **M3**. F5 hidden conservation drift via rhoCpPhi (Sec.5.2 F5). Add
  consistency monitor.
- **M4**. Anisotropic `h_ref = min(V)/max(Sf)` correctly handles thin cells
  but is global; for unstructured meshes, local `mesh.V()/mesh.surfaceArea()`
  per-cell would be more appropriate. Defer to future work.

### 10.4 Low impact / cleanup

- **L1-L10**. Dead code removal (Sec.9 inventory).

---

## 11. Recommended implementation roadmap

**Phase 1 -- Cleanup (no physics change, ~1 day).**
- Delete dead code per Sec.9 (TSense block, condSwitch, betaThermal,
  count-0 loops, commented gates, stale STEP duplicates).
- Unify `phiTotal`/`phiCN` naming in alphaEqn.
- Remove duplicate rho/rhoCp rebuild.
- Move `coldPhaseIsHighAlpha1` into `thermalPhaseChangeModel` base class.
- Add `Q_pc_thermal_` removal / unification with `Q_pc_`.

**Phase 2 -- Helmholtz decision (1 day).**
- Decide localized vs redistributed HW (Sec.3).
- If localized: remove `mdot_` Helmholtz solve, set `Qcorr == 0`,
  rename diagnostics.
- If redistributed: switch phiStefan and PCV to use `mdot`. Document
  that this restores 25 % interface-position error; case parameters
  must be re-validated.

**Phase 3 -- qn reconstruction (2-3 days).**
- Implement LS hot-side gradient (Sec.4.3).
- Unify kEff field (Sec.4.4).
- Remove `kinFloorCells_`, `AiFloor_kin` may stay (absolute, defensible).
- Re-run Sheikh 1D Stefan validation.

**Phase 4 -- TEqn simplification (1 day).**
- Move `Acoeff_floor` and `Acoeff_hot_floor` to dictionary-controlled
  diagnostics (off by default).
- Remove `Q_cond` split if pure evaporation cases (or keep for generality
  with clear documentation).
- Replace Acoeff form with single explicit RHS option for ablation study.

**Phase 5 -- Convergence study (2 days).**
- Per-outer-iter `correct()` call (Sec.6.1).
- Conservation monitor function-object (Sec.5.2 F5).
- Mesh convergence with new qn reconstruction.

**Phase 6 -- Documentation (continuous).**
- Inline derivation of all coefficients.
- Document any remaining magic numbers with citations.
- Rebuild user manual to match cleaned solver.

---

## 12. Summary of recommended state after audit

Field inventory after cleanup:
- **HardtWondra.C**: ~3-pathway architecture preserved exactly; ~150 lines
  removed; one unified kEff, one unified qn (LS), one mdotRaw, one
  Q_pc = mdotRaw*h_lv, one phiStefan, optional Helmholtz mdot retained only
  if Qcorr remains.
- **TEqn.H**: ~30 lines (no Acoeff floors by default; single Acoeff from
  Q_pc); diagnostic flag for emergency floors.
- **alphaSuSp.H**: trivial divU=PCV-only injection; Su=0 explicitly
  documented.
- **createFields.H**: `pcDictPtr` removed; flag accessed via model base
  class.
- **interTempFoam.C**: optional per-outer-iter `correct()` flag; consistent
  rho/rhoCp rebuild point.

Result: smaller, more defensible, faster, and publishable.

---

## 13. Defence of smoothing -- what is and isn't justifiable

The audit flags many smoothing operations as suspect. Before recommending
cleanup, each must be examined on its merits. The criterion is single and
sharp:

> A regularizer is **publication-defensible** if and only if its
> characteristic length scale is set by mesh geometry alone, vanishes as
> `h -> 0`, and discretizes a continuum operator that has an honest place
> in the diffuse-interface Stefan formulation.

This is the *geometric-dependence* criterion. Anything tied to Deltat, to the
solver iteration count, or to a tuned constant with no derivation, fails.

### 13.1 Why a diffuse-interface Stefan solver needs *some* smoothing

The Hardt-Wondra formulation reinterprets the sharp-interface Stefan
condition
```
rho_l v_n h_lv = k dT/dn|_Gamma            (sharp)
```
as a volumetric source over the band where `|gradalpha| != 0`:
```
mdot(x) = |gradalpha_s|(x) * q_n(x) / h_lv      (diffuse)
```
Three field reconstructions are *required* in this formulation, and each
demands a band-limiting filter:

| Reconstructed field | Why band-limiting is required |
|---|---|
| `\|gradalpha_s\|` (interface area density) | Raw MULES `alpha` has a Heaviside-like staircase across the interface band. `gradalpha` then contains O(1/h) spikes that do not converge under mesh refinement and violate the asymptotic identity int\|gradalpha\|dV ~= A_Gamma. |
| `nHat = gradalpha/\|gradalpha\|` (interface normal) | The normal is undefined in bulk and rotates wildly on the band boundaries. Without filtering, `phiStefan = U_s * (nHat*S_f)` injects rotational components that distort the interface. |
| `q_n` and `mdot_raw` (Stefan flux & raw mass source) | Conductive flux reconstructed across the diffuse band has sub-cell variation that the energy equation cannot resolve. Without a controlled spreading of the latent sink, T develops one-cell spikes that propagate through Acoeff back into qn -- a divergent feedback. |

So three filtering operations are physically necessary. The question is
*how many distinct smoothing fields and operators are needed*, not whether
to smooth at all.

### 13.2 The diffusion-of-field-variable identity

Every smoothing operation in the current solver implements **forward-Euler
heat diffusion with a fictitious time tau**:
```
dalpha_s/dtau = D grad^2alpha_s,        alpha_s(tau=0) = alpha
```
Iterated over N explicit steps of length Deltatau. The continuum Green's function
gives
```
alpha_s(x, tau) = int G(x-x', tau) alpha(x') dx',    G = (4piDtau)^{-d/2} exp(-|x|^2/4Dtau)
```
A Gaussian filter of width
```
sigma(tau) = sqrt(2 D tau)        (1-D)
```

In the code:
```cpp
D_smooth   = Fo_per_iter * h_ref^2 / Deltat        // [m^2/s]
alpha_s <- alpha_s + D_smooth * grad^2alpha_s * Deltat             // explicit step
```
Net width after N iterations:
```
sigma^2 = 2 * D_smooth * (N*Deltat)
   = 2 * (Fo_per_iter * h_ref^2/Deltat) * (N*Deltat)
   = 2 * Fo_per_iter * N * h_ref^2
```
**The Deltat cancels exactly.** sigma depends only on `h_ref` and N. With the code's
`Fo_per_iter = alphaSmoothWidth^2/N` and the cap `Fo_per_iter <= 0.25`:
```
sigma = h_ref * sqrt(2 * alphaSmoothWidth^2 / N * N)  =  alphaSmoothWidth * h_ref * sqrt2
```
So the filter width is **purely geometric**: sigma = O(h_ref). This is the
asymptotic property that lets the diffuse-interface formulation converge
to the sharp Stefan limit as the mesh is refined.

This identity is the foundation of the cleanup proposal: **every smoothing
operation in HW *should* be geometric-only. Any that isn't is a bug.**

### 13.3 Per-operation defence

| Smoothing | Operator | Length scale | Deltat-dep? | Defensible? |
|---|---|---|---|---|
| `alpha -> alpha_s` (`nSmoothIter` passes) | Explicit `dalpha/dtau = Dgrad^2alpha` | sigma ~= `alphaSmoothWidth * h_ref` | **No** (cancels) | **Yes** -- required for `\|gradalpha\|` denoising |
| `alpha_s -> alphaGeom` (one more pass) | Same operator | sigma adds `~h_ref/sqrtN` | No | **Redundant** -- collapse into alpha_s |
| `T -> T_sense` clip + Lap | Same operator | sigma ~= 0.55*h_ref | No | **Dead** -- T_sense unused for qn |
| `nHat -> nHat_smooth` (2 passes, no Deltat) | Direct `n + lambda^2grad^2n` | sigma ~= sqrt2*lambda = sqrt2*`lambdaSmearCells*h_ref` | No | **Yes** -- required for phiStefan stability, geometric |
| `mdot_raw -> mdot` Helmholtz | Implicit `(1-lambda^2grad^2)mdot = mdot_raw` | sigma_out^2 = sigma_in^2 + lambda^2 | No | **Conditionally yes** -- see Sec.3; load-bearing only if `Q_pc=mdot*h_lv` |
| `Q_pc` extra Laplacian (count = 0) | Dead | n/a | n/a | **Dead** |
| TEqn `Acoeff_floor` | Implicit penalty | `rhoCp/Deltat` magnitude | **YES** | **NO** -- only Deltat-dependent regularizer in the entire code |
| TEqn `Acoeff_hot_floor` | Implicit penalty | `rhoCp/Deltat` magnitude | **YES** | **NO** -- same |
| qn search distance floor `kinFloorCells*h_ref` | Hard floor in denominator | `kinFloorCells * h_ref` | No | **Geometric but biases qn downward** -- replace via LS gradient (Sec.4.3) |

**Verdict**: of nine "smoothing/limiting" mechanisms inspected:
- 2 are physically required and properly mesh-set (alpha-smoothing, nHat-smoothing).
- 1 is conditionally required (Helmholtz, depends on Sec.3 choice).
- 4 are dead or redundant (alphaGeom 2nd pass, T_sense block, Q_pc smoothing,
  qn denominator floor).
- **2 are Deltat-dependent and indefensible** -- `Acoeff_floor`, `Acoeff_hot_floor`.

The "bloat" perception is correct in count but not in spirit: the
*architecture* is sound (every operation has a defensible continuum
analogue); the *implementation* duplicates, applies dead operations, and
adds two Deltat-dependent penalties that the rest of the model is structured
to avoid.

---

## 14. Why the current Deltat-dependent terms exist (and how to eliminate them)

`Acoeff_floor = C_floor * rhoCp/Deltat * 4alpha(1-alpha) * pos(Tsat - T)` provides an
implicit pull `T -> Tsat` with characteristic time `Deltat/C_floor`. As `Deltat -> 0`
the penalty becomes infinitely stiff -- pinning T to Tsat in the diffuse
band exactly. This is the *only* way to absorb spurious overshoots arising
from upstream qn noise without changing other model components.

The root cause of those overshoots is the qn reconstruction (Sec.4) -- a
single-neighbor first-order one-sided difference that switches stencil
between time steps as the interface moves, producing sub-cell `Deltaq_n`
discontinuities that release latent heat at non-physical rates. The
`rhoCp/Deltat` floor is then *not* a physics term -- it is an implicit limiter on
the *temperature jump* per timestep. This is the fingerprint of a stiff
underlying source.

**Eliminating Deltat dependence requires fixing the upstream qn**, not adding
more downstream regularizers. With the LS hot-side gradient (Sec.4.3):
- qn is continuous in time as the interface moves (LS weights vary smoothly).
- `dq_n/dt` is bounded by mesh geometry alone.
- T overshoots in the diffuse band reduce to truncation error of the
  Laplacian discretization, controllable by mesh refinement.
- The Acoeff term alone, with `Q_pc/(T-Tsat)` stiffness coefficient,
  is sufficient to enforce Tsat in the diffuse band.
- `Acoeff_floor` and `Acoeff_hot_floor` can be removed.

This converts the solver from "stabilization-by-stiffness" to
"stabilization-by-physical-consistency". The former requires
parameter tuning per case; the latter does not.

---

## 15. Clean Stefan solver architecture (proposed)

### 15.1 Design principles

1. **One filtered alpha field.** All geometric reconstructions (\|gradalpha\|, nHat,
   interfaceArea, masks) come from a single smoothed alpha with band width
   `delta = N_band * h_ref` (default `N_band = 1.5`).
2. **One filter operator.** A single implicit Helmholtz smoother
   `(I - delta^2/2 * grad^2) alpha_s = alpha` replaces multi-pass explicit smoothing.
   Equivalent in steady state to N->inf explicit passes (variance sigma^2 = delta^2),
   no stability restriction, no tau, no Deltat anywhere.
3. **One length scale delta.** The same delta controls: alpha band, nHat filter width,
   Helmholtz redistribution length (if retained), qn-search regularizer.
   Mesh-refinement convergence is then governed by a single parameter.
4. **Zero Deltat-dependent regularizers.** Every implicit penalty, every
   stabilizer, every floor, scales with `h_ref` or its derivatives only.
5. **Single qn pathway.** No thermal/kinematic split. One LS-reconstructed
   `dT/dn|_hot`, one `q_n`, one `mdot_raw`. Q_pc, phiStefan, PCV derived from
   this single chain.

### 15.2 Mathematical structure

Define the band length:
```
delta == N_band * h_ref,            N_band in [1, 2],  default 1.5
```
Filtered phase field:
```
alpha_s = (I - (delta^2/2) grad^2)^{-1} alpha          (one implicit elliptic solve)
```
Equivalently `alpha_s` is the Gaussian filter of alpha with variance delta^2, by
identity between Helmholtz smoothing and Gaussian convolution in the
continuum limit (operator splitting in radial Fourier space gives
`alpha_s = alpha/(1 + delta^2k^2/2) ~= exp(-delta^2k^2/2)` for `deltak << 1`).

Geometric primitives (all from alpha_s):
```
A_i = |gradalpha_s|                          [1/m]
nHat   = gradalpha_s / (|gradalpha_s| + eps_alpha)           eps_alpha = SMALL * 1/m
nHat_f = (filtered) face interpolation   no extra smoothing needed; alpha_s already filtered
```

Hot-side conductive gradient (LS, Sec.4.3):
```
M_c = Sigma_n w_n * (x_n - x_c)(x_n - x_c)T
b_c = Sigma_n w_n * (T_n - T_sat) * (x_n - x_c)
w_n = max(T_n - T_sat, 0)              one-sided weight
(gradT_hot)_c = M_c-^1 * b_c              (Moore-Penrose if rank-deficient)
dT/dn|_hot = max((gradT_hot * nHat), 0)
```

Effective interface conductivity (single field, hot-side biased):
```
k_eff = k_hot * alpha_hot + k_cold * (1 - alpha_hot)          arithmetic
```
where `alpha_hot = alpha_s` if cold phase has alpha=0, `1-alpha_s` otherwise. Use
arithmetic, not harmonic -- the Stefan flux is being *delivered* by the hot
phase, not transmitted through a series resistance.

Single Stefan flux chain:
```
q_n     = k_eff * dT/dn|_hot                          [W/m^2]
mdot_raw   = A_i * q_n / h_lv                            [kg/m^3/s]
mdot       = (I - delta^2 grad^2)-^1 mdot_raw                         (optional; only if Qcorr retained)
Q_pc    = mdot_raw * h_lv          (localized, default)
PCV     = mdot_raw * (1/rho_v - 1/rho_l)
U_s     = mdot_raw / (A_i + eps_A) / rho_int
eps_A     = SMALL * max(A_i)      pure zero-guard, mesh-independent in magnitude
phiStefan = (U_s nHat_f) * S_f
```

Energy equation (no floors, no condSwitch, no Q_cond split):
```
d(rhoCp T)/dt + grad*(rhoCp phi T) - grad*(k_eff,f gradT)
   + Sp(A_c, T)  =  A_c * T_sat
A_c = max( Q_pc / max(T - T_sat, eps_T), 0 )            eps_T = O(0.1) K
```
The Acoeff term alone is sufficient: in cells with Q_pc != 0, T is pulled
toward T_sat with characteristic rate `A_c / rhoCp`. There is no `rhoCp/Deltat`
penalty anywhere. Stability is achieved by physical consistency, not by
stiffness.

### 15.3 Why no Deltat appears anywhere

Trace every length scale:

| Field | Length scale | Source |
|---|---|---|
| alpha_s band sigma | delta = `N_band * h_ref` | mesh |
| nHat filter | delta (same) | mesh |
| mdot Helmholtz lambda | delta (same) | mesh |
| qn LS stencil | cellCells extent ~= h_ref | mesh |
| A_c stiffness | Q_pc/(T-T_sat) -- depends on flow state, not Deltat | physics |
| Anti-trench floor | **REMOVED** | -- |
| Anti-superheat floor | **REMOVED** | -- |
| AiFloor in U_s | absolute, eps*max(A_i) | mesh-relative |

Compare to the current code, which has `C_floor * rhoCp/Deltat` in two places.
After cleanup: zero Deltat in any source/sink coefficient. The only Deltat
remaining is in `rhoCp * dT/dt`, which is correct.

### 15.4 Asymptotic limit `h -> 0`

- delta -> 0 => alpha_s -> alpha (no smoothing), \|gradalpha\| -> delta_Gamma (surface delta).
- A_c -> inf on a measure-zero set => T = T_sat on the interface (Dirichlet).
- mdot_raw -> A_i * q_n / h_lv => `intmdot dV -> int_Gamma mdot dA` (sharp Stefan condition).
- phiStefan -> sharp jump in velocity at the interface.

This is the formal sharp-limit recovery that publication-grade Stefan
formulations must demonstrate. The current code with Deltat-dependent floors
**does not** have this limit: as h -> 0 with Deltat fixed, the floors remain at
fixed strength while the conductive gradient becomes sharper, so the floor
becomes increasingly irrelevant -- but as `h -> 0` *and* `Deltat -> 0` (CFL-locked),
the floors become infinitely stiff. The behavior depends on which limit is
taken first, which is a publication red flag.

### 15.5 Implementation footprint

The clean architecture replaces ~1300 lines of `HardtWondra.C` with ~400.
Pseudocode skeleton:

```cpp
void HardtWondra::correct()
{
    const scalar h = computeHref();
    const dimensionedScalar deltaSqr_half("d2h", dimArea, 0.5*sqr(N_band_*h));
    const dimensionedScalar deltaSqr     ("d2",  dimArea,     sqr(N_band_*h));

    // 1. ONE filter on alpha (implicit, no Deltat)
    volScalarField alphaS = alpha1_;
    fvScalarMatrix sEqn
    (
        fvm::Sp(dimensionedScalar("one", dimless, 1.0), alphaS)
      - fvm::laplacian(deltaSqr_half, alphaS)
     == alpha1_
    );
    sEqn.solve();
    alphaS.clip(0, 1);

    // 2. Geometric primitives -- all from alphaS
    volVectorField gAlpha = fvc::grad(alphaS);
    interfaceArea_ = mag(gAlpha);
    volVectorField nHat = gAlpha / (interfaceArea_ + epsA_);

    // 3. LS hot-side gradient -> qn  (Sec.4.3)
    volScalarField dTdn = leastSquaresHotGradient(T_, T_sat_, nHat);
    volScalarField kEff = k_hot_*alphaHot_ + k_cold_*(1-alphaHot_);   // arithmetic
    qn_       = kEff * dTdn;
    mdotRaw_  = interfaceArea_ * qn_ / h_lv_;

    // 4. (Optional) Helmholtz redistribution -- only if useEnthalpyCorrection
    if (useEnthalpyCorrection_)
    {
        fvScalarMatrix mEqn
        (
            fvm::Sp(dimensionedScalar("one", dimless, 1.0), mdot_)
          - fvm::laplacian(deltaSqr, mdot_) == mdotRaw_
        );
        mEqn.solve();
        Qcorr_ = -mdot_*(cp1_-cp2_)*(T_-T_sat_);
        Q_pc_  =  mdot_*h_lv_;            // consistent HW
    }
    else
    {
        Qcorr_ = dimensionedScalar(Qcorr_.dimensions(), 0);
        Q_pc_  = mdotRaw_*h_lv_;           // localized
    }

    // 5. Stefan flux  (kinematic; mesh-set epsA)
    volScalarField Ustef = mdotRaw_/(interfaceArea_ + epsA_)
                          / rhoInt(alphaS);
    phiStefan_ = fvc::interpolate(Ustef * nHat) & mesh_.Sf();
}
```

One smoothing solve. One LS gradient. One qn. One mdot chain. No dead
code, no tuned constants beyond `N_band` (geometric).

---

## 16. Generalization to arbitrary phase change

The cleaned formulation is symmetric in melting/evaporation/solidification/
condensation provided three accessors are added to the model:

```cpp
class thermalPhaseChangeModel
{
public:
    virtual bool   coldPhaseIsHighAlpha1() const = 0;     // Sec.6.4 fix
    virtual scalar h_lv_signed()           const = 0;     // h_lv for evap, +Lf for melt
    virtual scalar densityRatioFactor()    const = 0;     // (1/rho_cold - 1/rho_hot)
    ...
};
```

Then:
- **Evaporation** (water boiling): cold=liquid (alpha_1=1), `densityRatioFactor =
  1/rho_v - 1/rho_l > 0`, vapor created -> PCV > 0.
- **Condensation**: same signs reversed via `dT/dn|_hot < 0` impossibility
  (no hot neighbor with T < Tsat) -- model returns mdot=0; condensation case
  uses cold-side LS gradient instead. Add `dTdn_cold` accessor for
  bidirectional models.
- **Melting** (Sheikh Stefan): cold=solid (alpha_1=1), `densityRatioFactor =
  1/rho_l - 1/rho_s` (often small), liquid produced from solid.
- **Solidification**: hot-side flag inverts.

The single LS stencil + sign convention via `densityRatioFactor` handles
all four regimes with no branching in the inner loop. The `coldSideGate`,
`coldPhaseIsHighAlpha1`, and `Q_cond` split currently in the solver can
all be eliminated by routing through these three accessors.

---

## 17. Summary -- the defensible minimum

A publication-grade Stefan VOF solver retains, with proper derivation:

1. **One implicit alpha-filter**: `(I - delta^2/2 grad^2)alpha_s = alpha`, delta = `N_band * h_ref`.
2. **One geometric reconstruction**: `A_i = \|gradalpha_s\|`, `nHat = gradalpha_s/(\|gradalpha_s\|+eps)`.
3. **One LS hot-side gradient** for qn, no magic floors.
4. **One arithmetic hot-biased kEff** for qn (sharp Stefan delivery side).
5. **One harmonic kappaf** for TEqn Laplacian (series resistance through
   the diffuse band).
6. **One Helmholtz redistribution** for mdot -- only if `Q_pc = mdot * h_lv` with
   Qcorr enabled. Otherwise drop entirely.
7. **One Acoeff** in TEqn: `Q_pc/max(T - Tsat, eps_T)`. No floors.
8. **One length scale delta** governing all filters.
9. **Zero Deltat-dependent coefficients** anywhere outside `d/dt` operators.
10. **One phase-change correction call per outer iteration** (not just
    firstIter -- Sec.6.1).

Everything else is bloat or bandaid. The current solver contains the
correct physics buried inside this scaffolding; cleanup is removal, not
rewrite.

---

## 18. Directionality analysis of the phiStefan vector

Construction chain (HardtWondra.C):

```
gradAlpha  = grad(alphaGeom)                  // alphaGeom = doubly-smoothed alpha1
nHat       = gradAlpha / |gradAlpha|          // unit cell-centered normal
nHatSmooth = 2 Laplacian passes, then renormalize
Ustef      = mdotInterfacial / rhoInt         // scalar, >= 0 by construction
phiStefan  = Ustef_f * (nHatSmooth_f . Sf) * interfaceMaskF
```

### 18.1 Direction of nHat

`nHat = grad(alpha1)/|grad(alpha1)|` points from the alpha=0 phase toward
the alpha=1 phase. Whether this is "hot->cold" or "cold->hot" depends on
which phase alpha1 represents:

| Case | alpha1=1 | alpha1=0 | nHat points |
|---|---|---|---|
| Melting (`coldPhaseIsHighAlpha1=true`) | solid (cold) | liquid (hot) | hot -> cold |
| Evap, liquid superheated | liquid (hot) | vapor (cold) | cold -> hot |
| Evap, vapor superheated | liquid (cold) | vapor (hot) | hot -> cold |

The `coldPhaseIsHighAlpha1` flag (TEqn.H) is consulted only by the
anti-superheat gate -- the phiStefan construction never reads it. The
direction is set entirely by the local `grad(alpha)`, which is correct:
direction is a property of the alpha field, not a runtime choice.

### 18.2 Sign of Ustef

```
mdotRaw  = Ai * qn / h_lv,    qn = kEff * dTdn_hot
Ustef    = mdotRaw / [(Ai + AiFloor_kin) * rho_int]
```

The neighbor search at `HardtWondra.C:740` enforces
`superheat = T_nbr - Tsat > 1e-6`, so `dTdn_hot >= 0` always. Therefore
`qn >= 0`, `mdotRaw >= 0`, `Ustef >= 0`.

**Implication.** The solver is **uni-directional**. It models evaporation
and melting (forward phase change, hot phase consumed) but **cannot
represent condensation or solidification** -- Ustef cannot flip sign.
`Q_cond = min(Q_pc, 0)` and `condSwitch` are dead branches in any case
using this construction.

### 18.3 Direction of interface motion

Stefan contribution to `dalpha1/dt + div(alpha1 (phi + phiStefan)) = 0`:

```
div(alpha1 * Ustef * nHat) ~ Ustef * (nHat . grad(alpha1))
                           = Ustef * |grad(alpha1)|     (nHat aligns with grad)
                           >= 0
```

So `dalpha1/dt <= 0` in the band: alpha=1 phase always **shrinks**,
alpha=0 phase grows, interface translates **in the +nHat direction**.

- Melting: nHat hot->cold => interface into solid; solid alpha1=1 shrinks. **OK.**
- Evap (liquid hot): nHat cold->hot => interface into liquid; liquid alpha1=1 shrinks. **OK.**

Sign-consistency across alpha conventions is automatic because `Ustef>=0`
and `nHat` is built from `grad(alpha1)`. Direction is not user-chosen --
it is forced by the sign of the hot-side conductive flux.

### 18.4 Face-orientation convention

```
phiStefan_f = Ustef_f * (nHat_f . Sf),   Sf points owner -> neighbor
```

- positive => mass-tracking flux owner->neighbor (interface advancing
  owner-to-neighbor in alpha-transport)
- negative => opposite

Standard OpenFOAM face-flux convention; no MULES handling needed.

### 18.5 Directional weaknesses

**D1. Triple smoothing of normal direction.** `nHat` is built from
`alphaGeom` (already smoothed twice from `alpha1`) and then smoothed by
two additional Laplacian passes and renormalized. Each pass averages
normal directions across neighbors. On curved interfaces this rotates
`nHat_f` away from the true geometric normal by O(curvature *
delta_band).

**D2. Face-interpolation of a unit vector underestimates magnitude.**
`fvc::interpolate(nHatSmooth)` produces a face vector with

```
|nHat_f| = |0.5 (nHat_O + nHat_N)| = cos(theta/2)
```

where theta is the angle between cell-centered normals. phiStefan
magnitude is multiplied by cos(theta/2), under-predicting interface speed
at curvature. For theta=60deg this is 13 % error. For 1D planar fronts
(Sheikh validation) theta=0 and the error vanishes -- which is why this
defect has not surfaced.

**Mitigation.** Decouple magnitude from direction at faces:

```cpp
const surfaceVectorField nHatf_raw = fvc::interpolate(nHatSmooth);
const surfaceVectorField nHatf
(
    nHatf_raw
  / max(mag(nHatf_raw),
        dimensionedScalar("eps", dimless, SMALL))
);
phiStefan_ = fvc::interpolate(Ustef) * (nHatf & mesh_.Sf()) * interfaceMaskF;
```

One extra `mag` + `divide` per face, no extra solve.

**D3. interfaceMaskF fixed threshold.** `interfaceMaskF = min(min(alpha_f,
1-alpha_f)/0.05, 1)` uses a hard cutoff at alpha=0.05. At thin diffuse
interfaces, only one face per interface satisfies `alpha_f in (0.05,
0.95)`, producing direction-dependent step changes in phiStefan
magnitude along the interface tangent. Effect: small along-interface
waviness uncorrelated with physics. **Recommendation**: replace with
`1 - tanh^2(2*(2*alpha_f - 1))`, a smooth band centered at alpha=0.5
matched to `alphaSmoothWidth`.

**D4. Tangential leakage from smoothed nHat.** Mathematically n is purely
normal to the iso-alpha surface. After two Laplacian passes,
`nHatSmooth` is no longer exactly aligned with `grad(alpha)` and
inherits orientation from neighbors with possibly different gradients.
`nHat_f . Sf` then includes tangential-motion contributions, producing
spurious tangential interface drift. **Recommendation**: at the end of
nHat smoothing, reconstruct `nHat` from the smoothed `alphaGeom`
directly:

```cpp
nHatSmooth = grad(alphaGeom)
           / max(mag(grad(alphaGeom)),
                 dimensionedScalar("eps", dimless/dimLength, SMALL));
```

Achieves smoothness via the alpha field; guarantees purely-normal output.

**D5. Uni-directionality is hard-coded.** The `superheat > 0` filter
excludes cold neighbors entirely. Bidirectional capability requires
either (i) a parallel `dTdn_cold` reconstruction with sign-flip when
`dTdn_cold > dTdn_hot`, or (ii) the LS gradient of Sec.4.3 with signed
weights `w_n = T_n - Tsat` (which carries sign automatically). Option
(ii) is the publication-grade path: one solve handles both regimes,
Ustef sign follows naturally, `Q_cond` becomes live code.

### 18.6 Verdict

| Property | Status |
|---|---|
| Sign of Ustef | Correct for forward phase change; structurally cannot reverse |
| Direction (alignment with grad(alpha)) | Correct, automatic across alpha conventions |
| Magnitude at curvature | Under-predicted by cos(theta/2) at curved interfaces |
| Tangential purity | Degraded by nHat smoothing; restorable via grad(alphaGeom) reconstruction |
| Mask continuity | Degraded by hard alpha cutoff; restorable via smooth tanh band |

For 1-D planar validation (Sheikh) all D1-D4 vanish; the construction is
sound. For 2-D/3-D curved-interface targets, D2 and D4 should be fixed
(5-line changes each) before mesh-convergence studies. D5 is a deliberate
restriction; condensation/solidification scope requires the Sec.4.3 LS
reformulation.

---

## 19. UEqn and pEqn implementation review

### 19.1 UEqn (`UEqn.H`)

```cpp
fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevRhoReff(rho, U)
 == fvOptions(rho, U)
);
UEqn.relax();
fvOptions.constrain(UEqn);
if (pimple.momentumPredictor())
{
    solve(UEqn == fvc::reconstruct((
          mixture.surfaceTensionForce()
        - ghf*fvc::snGrad(rho)
        - fvc::snGrad(p_rgh)
    )*mesh.magSf()));
    fvOptions.correct(U);
}
```

This is **vanilla interFoam UEqn, unmodified**. Observations:

**O1. No Stefan body force.** Correct decision. The Stefan velocity is
NOT a real momentum carrier -- it is an interface-tracking auxiliary
field. Injecting a `phiStefan`-derived body force into the mixture
momentum equation would double-count the dilatation already handled by
pEqn's PCV source.

**O2. No latent-heat body force.** Correct. Volumetric expansion forces
arising from phase change enter via `PCV` in the pressure equation
(`div(U) = PCV` constraint), which produces the correct
pressure-driven outflow from the interface.

**O3. `rho` is the volume-averaged mixture density.** Built from
`alpha1*rho1 + (1-alpha1)*rho2` and updated *after* pEqn
(`interTempFoam.C:180-186`). The momentum diffusion operator
`turbulence->divDevRhoReff` uses the same alpha-weighted blending for
viscosity. This is the standard immiscible-two-phase mixture model;
no phase-change-specific modification needed.

**O4. `rhoPhi` is updated in `alphaEqn.H:273` after MULES.** This is
critical: after phase change, alpha has changed and rho has changed, so
`rhoPhi = alphaPhi10*(rho1f - rho2f) + phi*rho2f` carries the
post-phase-change momentum flux. **Correct ordering.**

**Verdict on UEqn**: Clean, no modifications required for HW phase
change. Future EHD/Marangoni extensions inject as additional RHS
`fvc::reconstruct(...)` terms; the surface-tension term placement is the
canonical Rhie-Chow-compatible insertion site.

### 19.2 pEqn (`pEqn.H`)

```cpp
rAU = 1.0/UEqn.A();
HbyA = constrainHbyA(rAU*UEqn.H(), U, p_rgh);
phiHbyA = fvc::flux(HbyA) + ddtCorr;
phig = (surfaceTensionForce() - ghf*snGrad(rho)) * rAUf * magSf;
phiHbyA += phig;

while (correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn
    (
        fvm::laplacian(rAUf, p_rgh)
     == fvc::div(phiHbyA) - phaseChangePtr->PCV()      // <-- ONLY HW modification
    );
    p_rghEqn.solve(...);
    phi = phiHbyA - p_rghEqn.flux();
    U   = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
}
```

The pEqn is **standard interFoam with a single modification**: the
`- PCV()` source on the RHS.

**Mathematical role of PCV in pEqn.**

The pressure equation is `div(U) = source`. For incompressible interFoam
this source is identically zero. With phase change:

```
div(U) = mdot * (1/rho_v - 1/rho_l) = PCV
```

This comes from summing the two phase-mass equations after dividing each
by its constant density. PCV > 0 in evaporation: vapor is being created
at lower density than the liquid it replaces, so the local fluid expands
at rate PCV. The pressure correction step computes a `p_rgh` whose
gradient produces a velocity correction such that the new face fluxes
satisfy `div(phi) = PCV * V_cell`. **Correct and necessary.**

**Consistency observations.**

**P1. `phiHbyA` does NOT contain `phiStefan`.** Stefan flux is an
interface-tracking velocity, not a momentum-carrying velocity. Including
it in `phiHbyA` would imply momentum flux at the interface from a
non-physical velocity -- a Rhie-Chow violation. **Correct exclusion.**

**P2. `phig` does NOT contain a phase-change body-force component.**
This matches O2: no Stefan body force in UEqn, no Stefan flux in
phiHbyA. The pair is Rhie-Chow-consistent.

**P3. PCV uses `Q_pc = mdotRaw*h_lv` (localized).** The PCV source is
therefore concentrated in the diffuse interface band (1-2 cells wide).
The pressure Laplacian then has a sharply-localized RHS, producing
O(1/h) pressure spikes at the interface. This is the physically-correct
behavior (mass is being created at the interface, not redistributed)
but it strains the GAMG solver on fine meshes: tighter `tolerance`
(1e-8 instead of 1e-7) and 2+ non-orthogonal correctors may be needed.

**P4. Alternative formulation: Helmholtz-smoothed PCV.** If
`Q_pc = mdot*h_lv` (Sec.3 redistributed variant), PCV is spread over
~lambda. Pressure field is smoother, easier to solve, but introduces
phantom volume sources in cells outside the physical phase-change zone
-- the same artifact Qcorr is designed to remove from the energy
equation. **Consistency rule**: if Q_pc uses mdotRaw (localized), PCV
also uses mdotRaw via the same `Q_pc` accessor. Current code does this
correctly. If Q_pc switches to `mdot`, PCV follows automatically (same
accessor) and Qcorr must be enabled.

**P5. `p_rghEqn.flux()` is used to update `phi`.** This is the standard
Rhie-Chow-consistent flux update. Phase change does not alter this
step -- the PCV source enters the pressure Poisson, the resulting
flux correction propagates to `phi`, and `phi` is automatically
divergence-balanced with `PCV` at the discrete level.

**Verdict on pEqn**: Mathematically correct and minimally modified from
vanilla interFoam. The single modification (`- PCV()`) is the canonical
diffuse-interface phase-change pressure coupling. **No issues found.**
Solver-convergence diagnostics should monitor `p_rghFinal` residuals on
fine meshes -- if residuals stall, raise `nNonOrthogonalCorrectors` or
tighten `tolerance`, not introduce empirical relaxation.

### 19.3 Joint UEqn-pEqn consistency

| Item | UEqn | pEqn | Consistent? |
|---|---|---|---|
| rho field | mixture rho | mixture rho | OK |
| Surface tension | `surfaceTensionForce()` in RHS reconstruct | `surfaceTensionForce()` in `phig` | OK (Rhie-Chow) |
| Gravity | `ghf*snGrad(rho)` in RHS reconstruct | `ghf*snGrad(rho)` in `phig` | OK |
| Stefan flux | absent (correct) | absent in `phiHbyA` (correct) | OK |
| Phase-change source | absent (correct) | `- PCV()` in Laplacian RHS | OK (physically distinct roles) |
| Update order | predicted | corrected after solve | OK |

Both equations are clean, mutually consistent, and standard. The
publication-grade audit finds no issues in this pair beyond P3's
fine-mesh pressure-solver tightening recommendation.

---

## 20. alphaEqn source structure -- Su, Sp, divU vs phiStefan and PCV

### 20.1 The user's suspected redundancy

`alphaSuSp.H` constructs three source fields:

```cpp
tmp<volScalarField> tDivU(phaseChangePtr->PCV());   // = PCV
volScalarField::Internal Su (... , 0);              // ZERO (despite comment)
volScalarField::Internal Sp (... , 0);              // ZERO
volScalarField::Internal divU (... , tDivU().primitiveField());   // = PCV
```

Injected into `alphaEqn.H:130-131`:

```cpp
fvScalarMatrix alpha1Eqn
(
    fvm::ddt(alpha1)
  + fv::gaussConvectionScheme(mesh, phiCN, upwind<scalar>(mesh, phiCN))
        .fvmDiv(phiCN, alpha1)
 == Su + fvm::Sp(Sp + divU, alpha1)
);
```

with `phiCN = phi + phiStefan` (Crank-Nicolson blend optional).

The user's concern: **Su, Sp, divU, phiStefan, and pEqn's PCV all carry
some form of "phase-change information." Is something being
double-counted?**

The answer requires deriving the conservative form rigorously.

### 20.2 Mathematical derivation

Per-phase mass conservation, constant phase densities:

```
d(alpha rho_1)/dt + div(alpha rho_1 U) = -mdot      (mass leaves phase 1)
d((1-alpha) rho_2)/dt + div((1-alpha) rho_2 U) = +mdot
```

Divide each by its density:

```
(A)  dalpha/dt + div(alpha U) = -mdot/rho_1
(B)  d(1-alpha)/dt + div((1-alpha) U) = +mdot/rho_2
```

Add:

```
0 + div(U) = mdot (1/rho_2 - 1/rho_1) = PCV          (mixture continuity)
```

Expand (A) non-conservatively:

```
dalpha/dt + U . grad(alpha) + alpha div(U) = -mdot/rho_1
dalpha/dt + U . grad(alpha) + alpha PCV   = -mdot/rho_1
```

Rearrange back to conservative form keeping the dilatation explicit:

```
dalpha/dt + div(alpha U) = -mdot/rho_1 + alpha PCV     ... (*)
```

Equation (*) is the alpha equation that must be solved. Two source
contributions appear:

- **Source A**: `-mdot/rho_1` = `alpha1Gen()` (loss of phase 1 to phase
  change).
- **Source B**: `+alpha*PCV` (volumetric dilatation of phase 1 driven by
  mixture expansion).

These two are **mathematically distinct** and **not redundant**: A is
the direct mass-loss source, B is the dilatation correction reflecting
that `div(U) != 0`.

### 20.3 Two equivalent implementations

There are two ways to discretize (*):

**Formulation B (source-based, vanilla style):**

```
phiCN = phi                            (no Stefan flux)
RHS   = Su + Sp(divU, alpha)
       = alpha1Gen() + alpha * PCV
       = -Q_pc/(rho_1 h_lv) + alpha * PCV
```

**Formulation A (kinematic, current code):**

```
phiCN = phi + phiStefan                (Stefan flux added to convection)
RHS   = Sp(divU, alpha) = alpha * PCV
Su    = 0                              (alpha1Gen absorbed into phiStefan)
```

In Formulation A, the source `-mdot/rho_1` is moved into the convection
operator as the divergence of an auxiliary flux:

```
div(alpha phiStefan) ~ alpha div(phiStefan) + phiStefan . grad(alpha)
```

In the interface band, `div(phiStefan)_cell ~ mdotRaw/rho_int * V_cell`
(by the Sec.18 cell-volume derivation), so the discrete operator
contributes approximately `-mdot/rho_int` to the alpha equation -- which
replaces the `-mdot/rho_1` source (modulo the small `rho_int` vs `rho_1`
distinction in the diffuse band).

### 20.4 Where redundancy WOULD live, and why it doesn't

**Redundancy check 1: Su and phiStefan.** If Su were `alpha1Gen()`
(non-zero) AND phiCN included phiStefan, the source `-mdot/rho_1` would
be applied twice: once explicitly via Su, once kinematically via
div(alpha phiStefan). The current code **zeroes Su** (alphaSuSp.H:11-16)
precisely to avoid this. **Not redundant by design.**

The comment in alphaSuSp.H:1-3 is **stale and misleading**:

```
// Phase-change sources for MULES (volScalarField::Internal, [1/s])
// Su  = alpha1Gen: positive for condensation, negative for evaporation
```

This describes Formulation B; the code implements Formulation A. The
comment should read:

```
// Phase-change sources for MULES.
// Su = 0: alpha1Gen is absorbed into phiStefan in alphaEqn.phiCN.
// Sp = 0: no implicit phase-change linearisation needed.
// divU = PCV: volumetric dilatation source for the alpha*div(U) term.
```

**Redundancy check 2: PCV in pEqn vs PCV in alphaEqn.** Distinct roles:

- pEqn-PCV: enforces `div(U) = PCV` on the velocity field (mixture
  continuity).
- alphaEqn-PCV (via divU): provides the `+alpha*PCV` dilatation term in
  equation (*).

The first constrains a *velocity field divergence*, the second is a
*source coefficient on alpha*. Mathematically independent operators
that happen to share the same coefficient `PCV(x)` because both
ultimately trace to the mass-conservation constraint
`div(U) = mdot * (1/rho_v - 1/rho_l)`. **Not redundant.**

**Redundancy check 3: alpha1Gen base-class implementation vs current
usage.** `thermalPhaseChangeModel::alpha1Gen()` is defined in the base
class (`thermalPhaseChangeModel.C:117-138`) and returns
`-Q_pc/(rho_1 h_lv)`. It is **never called** by the current solver
(`alphaSuSp.H` hardcodes Su=0). This is dead infrastructure -- it would
matter only if Formulation B were used, which the solver does not. See
Sec.7.1 for the architectural recommendation: either delete `alpha1Gen`
from the base class or document it as unused.

### 20.5 Quantitative consistency check

Integrating equation (*) over the interface band:

```
d/dt int_band(alpha) dV
   = -int_band div(alpha (U+U_s)) dV
     + int_band alpha PCV dV
```

The first term, by Gauss, equals the net inward flux through the band
boundary. In a closed evaporation problem, the band is interior and the
fluxes through its outer surface balance the kinematic motion of the
interface. The second term equals approximately
`alpha * mdot * (1/rho_v - 1/rho_l)`, the volumetric expansion of the
liquid (alpha=1 side) into newly created vapor space.

Mass conservation in liquid mass `M_liq = int alpha rho_1 dV`:

```
dM_liq/dt = rho_1 * int dalpha/dt dV
         = rho_1 * [-flux + int alpha PCV dV]
         = -rho_1/rho_1 * int_Gamma mdot dA            (by Sec.20.3)
           + rho_1 * (mdot * (1/rho_v - 1/rho_l))_int_band
         = -int mdot dA  + corrections from dilatation term
```

A precise audit requires the actual `phiStefan` divergence integral, but
the structural result is: the kinematic phiStefan replaces the
`-mdot/rho_1` source one-for-one, and the `alpha*PCV` term is the
volume-correction completing mass conservation.

### 20.6 What IS suspect in alphaSuSp.H / alphaEqn

Not redundancy -- but adjacent issues:

**S1. `Sp = 0` always, but the MULES API expects Sp to be the
implicit-source diagonal coefficient.** Hardcoding to zero forces MULES
to bound alpha based solely on explicit Su + alpha*divU. This is
acceptable but discards a degree of freedom that could improve
boundedness during strong phase change. **No action required**, but
document the choice.

**S2. `divU = PCV` injection via `tmp<volScalarField>` and immediate
`.primitiveField()` strip is awkward.** A direct construction:

```cpp
volScalarField::Internal divU
(
    IOobject(...),
    mesh, dimless/dimTime,
    phaseChangePtr->PCV()().primitiveField()
);
```

without the named tmp is cleaner and equivalent. **Minor cleanup.**

**S3. The commented-out term `- fvm::Sp(fvc::ddt(1) + fvc::div(phiCN),
alpha1)` (alphaEqn.H:128-129).** This would discretely add
`-(d/dt + div(phiCN))*alpha1` to the LHS, converting the conservative
`div(alpha1 phiCN)` form into a non-conservative
`alpha1 * (d/dt + div(phiCN))` correction. With `divU = PCV` on the RHS
this would **double-count the dilatation**, since the commented term
embeds `div(phiCN)` which already includes PCV. Correctly left
commented. **Delete the comment to avoid confusion.**

**S4. `tSu` allocates a `volScalarField` only to extract its
`primitiveField()` for `Su`. Wasted allocation.** Replace with direct
`Su = 0` initialization (already dimensioned correctly):

```cpp
volScalarField::Internal Su
(
    IOobject(...), mesh, dimensionedScalar("Su", dimless/dimTime, Zero)
);
```

### 20.7 Verdict

| Suspected redundancy | Verdict |
|---|---|
| Su + phiStefan both carrying phase-change source | Not redundant; Su=0 by design (Formulation A) |
| PCV in pEqn AND alphaEqn | Not redundant; mathematically distinct roles |
| alpha1Gen() base-class accessor | Dead in current solver; remove or document |
| Comment claims Su=alpha1Gen | **Stale comment**; rewrite per Sec.20.4 |
| Commented `fvm::Sp(div(phiCN), alpha1)` | Correctly inactive; remove comment |
| `tSu` temp-field allocation | Minor cleanup possible |

The architecture is **mathematically clean**: Formulation A (kinematic)
is consistently applied, no double-counting exists, PCV's two roles are
physically distinct. The perception of redundancy comes from (i) stale
comments still describing Formulation B, (ii) dead `alpha1Gen()` base-
class code, and (iii) the cosmetic noise of `tSu/tDivU` temporaries.
Cleanup is documentation and dead-code removal, **not equation
modification**.

---

## 21. Acoeff in TEqn -- what it actually does

The `Acoeff` block in `TEqn.H` is the most opaque construction in the
solver. This section breaks it down from first principles.

### 21.1 What you see in the code

```cpp
const dimensionedScalar Tsat  = phaseChangePtr->T_sat();
const dimensionedScalar eps_T("eps_T", dimTemperature, 2.5);
const scalar implicitFactor = 0.5;

const volScalarField Acoeff
(
    implicitFactor
  * max(
        Q_pc_th_now / max(T - Tsat, eps_T),
        dimensionedScalar("zero", dimensionSet(1,-1,-3,-1,0,0,0), 0)
    )
);

// ... and in the matrix:
+ fvm::Sp(Acoeff, T)
 == Acoeff*Tsat
```

(Ignoring the floor terms for now -- those are added on top of this.)

This looks like a strange combination of a divide, two `max`es, a
saturation constant, and an implicit/explicit pair that share the same
coefficient. Let me unpack it.

### 21.2 The trick being used: implicit source linearization

The physical latent-heat source in the energy equation is:

```
rho Cp dT/dt + ... - div(k grad T) = -Q_pc
```

`Q_pc = mdot * h_lv > 0` for evaporation/melting (energy absorbed by
phase change). If we just slap `-Q_pc` on the RHS as a fully explicit
source, the equation is stable only when `dt` is small enough that T
cannot overshoot Tsat in a single step. On reasonable engineering meshes
with `Q_pc` of order 10^8 W/m^3 in interface cells, this is restrictive.

The standard OpenFOAM trick is to **factor the source as A*(T - Tsat)**:

```
-Q_pc = -A * (T - Tsat)     where    A = Q_pc / (T - Tsat)
```

This identity is true at every point where `T != Tsat`. The reason for
factoring it is that **`-A*(T - Tsat)` can be split implicit/explicit**:

```
-A*(T - Tsat) = -A*T + A*Tsat
                 ^^^^   ^^^^^^
                 LHS    RHS  (both implicit in T via fvm::Sp)
```

In OpenFOAM matrix form:

```
+ fvm::Sp(A, T)     // adds +A*T to LHS (matrix diagonal)
 == ... + A*Tsat;   // adds +A*Tsat to RHS
```

After rearranging, the net contribution to the RHS is
`A*Tsat - A*T = -A*(T-Tsat) = -Q_pc`. **Mathematically equivalent to
`== -Q_pc`**, but with two stability advantages:

1. `+A*T` enlarges the diagonal of the matrix => Better-conditioned
   linear system, faster GAMG/PCG convergence.
2. If `T` overshoots above `Tsat` in one iteration, `+A*T` immediately
   pulls it back -- implicit damping. Pure-explicit `-Q_pc` cannot do
   this; it just blindly subtracts a fixed value.

This is identical in structure to radiation linearization, drag-force
linearization, and many other OpenFOAM source patterns. The pattern is:

```
explicit source S(T)  -->  fvm::Sp(dS/dT, T) == dS/dT * T_0 + S(T_0)
```

with `T_0` chosen so that `dS/dT * T_0` cancels the implicit linearisation.
Here `T_0 = Tsat`.

### 21.3 Three protective `max` operations

The naive choice `A = Q_pc / (T - Tsat)` breaks at three places:

**Issue 1: T = Tsat exactly.** Division by zero. Fixed by the inner
`max`:

```
denom = max(T - Tsat, eps_T),  eps_T = 2.5 K
```

When `T` is within 2.5 K above `Tsat`, the denominator is clipped to
`eps_T`. The Acoeff becomes `0.5 * Q_pc / 2.5`, finite. Acts as a
soft saturation in the singularity neighborhood.

**Issue 2: T < Tsat.** Then `T - Tsat < 0`. Without the `max`, the
denominator would be negative, making `A < 0`. A negative diagonal
contribution destabilises the linear solver (loss of diagonal dominance)
and physically would *push* T further away from Tsat. Fixed by the
inner `max(T - Tsat, eps_T)` clipping to `eps_T` (positive).

**Issue 3: Q_pc < 0 (condensation, hypothetically).** Then `A < 0` for
the same diagonal-dominance reason. Fixed by the outer
`max(..., 0)`, which clips negative A to zero -- in condensation cells
Acoeff is just zero, and no implicit term is added. (In the current
code Q_pc is always >= 0, but the guard is there.)

Net protected formula:

```
A = max(Q_pc / max(T - Tsat, 2.5 K),  0)
```

### 21.4 The `implicitFactor = 0.5` mystery

Now the awkward part. If we apply the full identity
`-Q_pc = -A*(T - Tsat)`, we should use `A = Q_pc/(T-Tsat)` (i.e.
`implicitFactor = 1`). The code uses `0.5`. What does this mean?

Trace the effective RHS contribution:

```
+ fvm::Sp(Acoeff, T) == Acoeff * Tsat
```

is equivalent on the RHS to:

```
Acoeff * Tsat - Acoeff * T = -Acoeff * (T - Tsat)
                            = -0.5 * Q_pc / max(T-Tsat, 2.5) * (T - Tsat)
```

Two regimes:

| Regime | T - Tsat | Acoeff*(T-Tsat) effective sink | Magnitude vs Q_pc |
|---|---|---|---|
| Strongly superheated | T - Tsat > 2.5 K | 0.5 * Q_pc | **Half** of Q_pc |
| Near-saturation | 0 < T - Tsat < 2.5 K | 0.5 * Q_pc * (T-Tsat)/2.5 | < half |
| Subcooled | T < Tsat | 0.5 * Q_pc * (T-Tsat)/2.5 < 0 | *Heat source*, pushes T up to Tsat |

So in the strongly-superheated regime, **the Acoeff term applies only
HALF of the latent heat as a sink**. The other half of `Q_pc` is **not
present** anywhere else in the equation -- there is no `- 0.5*Q_pc` on
the RHS to make up the deficit.

This is not a bug in the sense of an oversight -- it is a deliberate
**stabilisation strategy**, but it deserves to be made explicit. The
interpretation is:

**Acoeff is a "soft Dirichlet" penalty enforcing T = Tsat in the diffuse
band, NOT a literal latent-heat sink.** The literal latent heat is
delivered to the equation by the *conductive flux from the hot side
into a cell pinned at Tsat*: once T at the interface is constrained to
Tsat, the Laplacian `-div(k grad T)` automatically carries the correct
amount of energy across the interface. The Stefan condition
`k dT/dn = mdot * h_lv` then holds at the constrained surface as an
emergent property, not as an imposed source.

In this interpretation:

- The role of `Acoeff` is to **constrain T**, not to **balance latent
  heat magnitude**.
- The magnitude `0.5 * Q_pc / (T - Tsat)` need only be large enough to
  produce strong enough damping that T converges to Tsat in the band
  within a few outer iterations.
- `implicitFactor = 0.5` is then a **relaxation parameter** controlling
  how stiff the soft Dirichlet is. Larger => stiffer constraint, faster
  convergence to T=Tsat, but worse conditioning. Smaller => softer
  constraint, slower pin to Tsat, but the bulk solver is happier.

### 21.5 Why the soft-Dirichlet view is physically defensible

For a sharp-interface Stefan problem with `T = Tsat` at the interface
imposed as a Dirichlet boundary condition:

```
- div(k grad T) = 0     in each phase
T = Tsat                on Gamma
k_hot * dT/dn|_Gamma = mdot * h_lv     (Stefan condition closure)
```

The latent heat does not appear as a volumetric source. It is delivered
by the *flux through the Dirichlet boundary*: the conductive flux into
the constraint surface equals `mdot * h_lv` automatically when the
temperature field, conductivity, and `mdot` are mutually consistent.

In a diffuse-interface VOF setting we cannot impose a Dirichlet BC on a
moving surface, so we approximate it with a volumetric penalty:

```
"T = Tsat"  --(approximate)-->   + fvm::Sp(A_large, T) == A_large * Tsat
```

The larger `A_large`, the closer T sticks to Tsat in the band. Once T
is pinned, the conductive flux delivers latent heat naturally, exactly
as in the sharp-interface case. The `Acoeff = Q_pc / (T - Tsat)` choice
**adaptively sets `A_large` based on the local phase-change activity**:
where Q_pc is strong, A is large, the constraint is tight; where Q_pc
is weak, A is small, the constraint is loose. This is a smart
self-scaling choice and the reason this formulation works in practice.

### 21.6 The Acoeff_floor and Acoeff_hot_floor add-ons

The two floor terms in TEqn.H complement Acoeff:

| Coefficient | Form | Role |
|---|---|---|
| `Acoeff` | `0.5 * Q_pc / max(T-Tsat, 2.5)` | Active where Q_pc > 0; soft Dirichlet at Tsat |
| `Acoeff_floor` | `5 * rho Cp / dt * 4 alpha (1-alpha) * pos(Tsat - T)` | Active where T < Tsat in diffuse band; pulls T UP to Tsat |
| `Acoeff_hot_floor` | `5 * rho Cp / dt * 4 alpha (1-alpha) * pos(T - Tsat) * coldSideGate` | Active where T > Tsat in cold-phase side of diffuse band; pulls T DOWN to Tsat |

Same implicit `Sp(A, T) == A*Tsat` structure for each, so each is a soft
Dirichlet at Tsat. They differ in *when they activate*:

- Acoeff activates where there is positive Q_pc (phase change is
  driving heat exchange).
- Acoeff_floor activates where T has fallen below Tsat in the band
  (anti-trench: prevents sub-saturation cold spots).
- Acoeff_hot_floor activates where T exceeds Tsat on the cold side
  (anti-superheat: prevents thermal leakage into the cold phase).

The two floors are the **`dt`-dependent** terms flagged in Sec.13-14
as numerical bandaids. Their `rho Cp / dt` scaling makes them
infinitely stiff as `dt -> 0`, which is the publication red flag.

### 21.7 Summary diagram

```
TEqn LHS:
    rho Cp dT/dt + div(rho Cp phi T) - div(k grad T)
      + Sp(Acoeff,           T)    <-- soft Dirichlet @ Tsat, where Q_pc > 0
      + Sp(Acoeff_floor,     T)    <-- soft Dirichlet @ Tsat, anti-trench (cold side)
      + Sp(Acoeff_hot_floor, T)    <-- soft Dirichlet @ Tsat, anti-leakage  (hot side, cold phase)

TEqn RHS:
       Acoeff           * Tsat
     + Acoeff_floor     * Tsat
     + Acoeff_hot_floor * Tsat
     - Q_cond                       <-- condensation explicit branch (zero in evap)
     + Qcorr                        <-- enthalpy correction (HW eq. 42)
```

All three implicit penalties enforce the same constraint `T -> Tsat` in
different cells. The latent heat itself is not subtracted as a literal
source -- it is delivered emergently by the Laplacian flux into the
constrained band.

### 21.8 What to change for publication

**Recommended cleanup of Acoeff**:

1. **Remove `implicitFactor = 0.5`.** Use `A = Q_pc / max(T-Tsat, eps_T)`
   directly. The `0.5` is a relaxation parameter for which no
   justification exists in the code; documenting its origin or removing
   it improves defensibility. If empirical tuning is required for
   stability, do it via `TEqn.relax()` (which has a dictionary-set
   under-relaxation factor) rather than hardcoded scaling.

2. **Reduce `eps_T = 2.5 K` to ~0.5 K**. The denominator floor is a
   regularizer for the singularity at T=Tsat. 2.5 K is large enough to
   significantly weaken Acoeff in the typical superheat range of
   evaporation problems. A smaller value gives a tighter Dirichlet at
   the cost of stiffer matrices -- worth doing because the matrix is
   solved by GAMG which handles stiffness well.

3. **Remove `Acoeff_floor` and `Acoeff_hot_floor`** once the qn
   reconstruction (Sec.4.3) eliminates the underlying noise that
   currently makes them necessary. Their `rho Cp / dt` scaling is the
   only `dt`-dependent regularisation in the solver and is
   publication-blocking (Sec.14).

4. **Rename the variable to make intent explicit**. `Acoeff` reads as
   "some coefficient A"; `softDirichletAtTsat` documents the physical
   role.

**Pedagogical comment block to add inline**:

```cpp
// Implicit "soft Dirichlet" penalty enforcing T -> Tsat in the diffuse
// interface band. Mathematical structure: fvm::Sp(A, T) == A * Tsat is
// equivalent to a RHS source of -A*(T-Tsat). Choosing A = Q_pc/(T-Tsat)
// makes the penalty strength scale with local phase-change activity.
// The latent heat itself is NOT subtracted explicitly here -- it is
// delivered by the Laplacian -div(k grad T) flowing into the Tsat-pinned
// band, which is the diffuse-interface analogue of the sharp-interface
// Stefan-condition closure.
```

### 21.9 Connection back to the rest of the solver

The Acoeff design assumes:

- `Q_pc_thermal` is the latent activity field driving the penalty.
- The harmonic `kappaf` carries the conductive flux into the
  Tsat-pinned band.
- `phiStefan` carries the kinematic interface motion separately.
- `Qcorr` removes the artificial enthalpy if the Helmholtz redistribution
  is active.

If any of these is missing or inconsistent, Acoeff cannot achieve its
soft-Dirichlet role and the floors become necessary. **The
publication-grade fix is to make all four upstream components physically
consistent (per Sec.3 and Sec.4 recommendations), after which Acoeff
alone -- with `implicitFactor = 1.0`, `eps_T ~ 0.5 K`, no floors -- is
sufficient.**

---

## 22. Q_cond -- is it required for a generalized bidirectional solver?

**Short answer: no.** Q_cond exists only because the current `Acoeff`
construction has a `max(..., 0)` outer guard that zeroes the implicit
penalty whenever `Q_pc < 0`. With a sign-aware denominator, a single
`Acoeff` handles both evaporation/melting (Q_pc > 0) and
condensation/freezing (Q_pc < 0) uniformly, and `Q_cond` becomes
unnecessary.

### 22.1 Why the current split exists

```cpp
Acoeff = max(Q_pc / max(T - Tsat, eps_T), 0);   // zero when Q_pc < 0
Q_cond = min(Q_pc, 0);                          // explicit branch for negatives
```

- `Q_pc > 0` (evap/melt): `Acoeff > 0`, implicit pull-to-Tsat active,
  latent sink delivered via `+Sp(A, T) == A*Tsat`.
- `Q_pc < 0` (cond/freeze): `Acoeff = 0`, no implicit term; `-Q_cond`
  on RHS adds the released latent heat explicitly.

The split is symptomatic, not principled -- the sign-handling of Acoeff
is simply incomplete.

### 22.2 Bidirectional Acoeff -- one expression for both signs

Build a denominator that **carries the sign of Q_pc** and stays bounded
away from zero:

```cpp
const volScalarField sQ = sign(Q_pc);

// denom has the same sign as Q_pc; magnitude is at least eps_T
const volScalarField denom
(
    sQ * max(sQ*(T - Tsat),
             dimensionedScalar("eps_T", dimTemperature, 0.5))
);

// Bidirectional Acoeff: ratio of same-signed quantities -> always >= 0
const volScalarField Acoeff = Q_pc / denom;
```

Regime table:

| Q_pc | T - Tsat | sQ | sQ*(T-Tsat) | denom | Acoeff = Q_pc/denom | Sp(A,T) == A*Tsat effective source |
|---|---|---|---|---|---|---|
| `> 0` | `>= 0` | `+` | `T - Tsat` | `T - Tsat` | `Q_pc/(T-Tsat) > 0` | `-Q_pc` (evap sink) |
| `> 0` | `< 0` (off-eq) | `+` | `< 0` -> clipped to `0.5` | `+0.5` | `Q_pc/0.5 > 0` | weak penalty pulling T up |
| `< 0` | `<= 0` | `-` | `Tsat - T` | `T - Tsat < 0` | `Q_pc/(T-Tsat) > 0` | `-Q_pc > 0` (cond source) |
| `< 0` | `> 0` (off-eq) | `-` | `< 0` -> clipped to `-0.5` | `-0.5` | `Q_pc/(-0.5) > 0` | weak penalty pulling T down |

In every regime:

- `Acoeff >= 0` (positive diagonal preserved, solver-stable).
- Effective RHS source = `-Acoeff*(T - Tsat) = -Q_pc` (correct sign for
  both evaporation and condensation).
- Off-equilibrium cells (Q_pc and (T-Tsat) of opposite sign, indicating
  transient mismatch) get a weak penalty pulling T back toward Tsat with
  strength `Q_pc / eps_T`, which is the same regularization role the
  current `eps_T` clip plays in the unidirectional case.

`Q_cond` is no longer needed -- the implicit `Acoeff` handles
condensation by itself.

### 22.3 Why this is preferable to the current split

| Issue | Current (split) | Bidirectional |
|---|---|---|
| Condensation stability | Pure explicit `-Q_cond` -- overshoots when stiff | Implicit damping, same as evaporation |
| Code paths | Two sign-dependent branches | One expression |
| Compatibility with `Q_pc()` accessor | Requires sign branching downstream | Self-consistent |
| Dead-code status (Sec.18 D5) | `Q_cond == 0` today because qn-search forces `Q_pc >= 0` | Becomes live once qn is reformulated (Sec.4.3 LS gradient) |

The current `Q_cond` is **dead code in the live solver**: the
`HardtWondra.C:740` neighbor-search filter `superheat > 1e-6` forces
`Q_pc >= 0`. `Q_cond = min(Q_pc, 0)` is therefore identically zero in
every run today. It would only become meaningful after the Sec.4.3
LS-gradient reformulation enables bidirectional `Q_pc`. At that point,
the bidirectional `Acoeff` above is the clean replacement -- `Q_cond`
need never be reintroduced.

### 22.4 Equivalent cleaner form using `pos`/`neg` functions

If you prefer to avoid `sign()` (which has discontinuous derivative at
zero), the same result follows from a `pos/neg` split factored into the
denominator:

```cpp
const volScalarField denom
(
    pos(Q_pc) * max(T - Tsat,
                    dimensionedScalar("eps_T", dimTemperature, 0.5))
  - neg(Q_pc) * max(Tsat - T,
                    dimensionedScalar("eps_T", dimTemperature, 0.5))
);
const volScalarField Acoeff = Q_pc / denom;
```

`pos(x) = 0.5*(1+sign(x))` and `neg(x) = 0.5*(1-sign(x))` are smooth
where supported by OpenFOAM, and zero on the inactive sign. Same result,
slightly more verbose.

### 22.5 What changes in the cleaned TEqn

```cpp
// REMOVE these:
// const volScalarField Q_cond
// (
//     min(Q_pc_now, dimensionedScalar("zero", ..., 0))
// );

// REPLACE Acoeff with bidirectional version (see 22.2)
const volScalarField sQ = sign(phaseChangePtr->Q_pc()());
const volScalarField denom
(
    sQ * max(sQ*(T - Tsat), eps_T)
);
const volScalarField Acoeff = phaseChangePtr->Q_pc() / denom;

// Matrix unchanged -- one less term on RHS:
fvScalarMatrix TEqn
(
    rhoCp * fvm::ddt(T)
  + fvm::div(rhoCpPhi, T)
  - fvm::laplacian(kappaf, T)
  + fvm::Sp(Acoeff, T)         // floors removed per Sec.21.8
 ==
    Acoeff * Tsat
  + phaseChangePtr->Qcorr()    // -Q_cond branch removed
);
```

Net change: ~5 lines removed (Q_cond block), ~3 lines added
(signed denominator), one less explicit RHS term.

### 22.6 Combined with the other Sec.21 cleanups

The fully cleaned TEqn for a generalized bidirectional solver:

```cpp
const dimensionedScalar Tsat  = phaseChangePtr->T_sat();
const dimensionedScalar eps_T("eps_T", dimTemperature, 0.5);   // was 2.5
const volScalarField   Qpc   = phaseChangePtr->Q_pc();
const volScalarField   sQ    = sign(Qpc);

// Signed denominator: same sign as Qpc, magnitude >= eps_T
const volScalarField denom(sQ * max(sQ*(T - Tsat), eps_T));

// Bidirectional soft-Dirichlet penalty (was: max(..., 0) split + Q_cond)
const volScalarField softDirichletAtTsat(Qpc / denom);    // implicitFactor=1 absorbed

fvScalarMatrix TEqn
(
    rhoCp * fvm::ddt(T)
  + fvm::div(rhoCpPhi, T)
  - fvm::laplacian(phaseChangePtr->kappaf(), T)
  + fvm::Sp(softDirichletAtTsat, T)
 ==
    softDirichletAtTsat * Tsat
  + phaseChangePtr->Qcorr()
);
TEqn.relax();
TEqn.solve();
```

Net result vs. the current TEqn.H:

| Component | Current | Bidirectional clean |
|---|---|---|
| Acoeff | `0.5 * max(Q_pc/max(T-Tsat, 2.5), 0)` (uni-directional, half-magnitude) | `Q_pc / (sQ*max(sQ*(T-Tsat), 0.5))` (bidirectional, full) |
| Acoeff_floor (anti-trench) | present, `~ rho Cp / dt` | **removed** |
| Acoeff_hot_floor | present, `~ rho Cp / dt` | **removed** |
| Q_cond | `min(Q_pc, 0)` (dead today) | **removed** |
| coldSideGate logic | present (TEqn.H + createFields.H bool flag) | **removed** |
| Qcorr | present (HW eq. 42) | retained (only meaningful if Helmholtz is kept) |

Roughly 60 lines of TEqn.H collapse to 15. Zero `dt`-dependent
coefficients. Single soft-Dirichlet expression handles all four
phase-change directions (evap, cond, melt, freeze).

### 22.7 Verdict

**For a generalized phase-change solver: `Q_cond` is not required.**
It is a workaround for a sign-incomplete `Acoeff`. A signed-denominator
`Acoeff` reproduces the same physics in one expression, with implicit
damping for both signs of phase change, and removes the explicit
condensation branch entirely. Combined with the Sec.21.8 cleanups
(`implicitFactor = 1`, `eps_T = 0.5 K`, no floors) and the Sec.4.3 LS
qn reformulation (signed gradient), the solver becomes fully
bidirectional with a single soft-Dirichlet term and no dead code.
