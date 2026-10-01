# interTempFoam

A non-isothermal two-phase VOF solver with phase change for **OpenFOAM v2412**
(openfoam.com), derived from `interFoam` / `interIsoFoam`. It combines geometric VOF
(**isoAdvector**, Roenby et al. 2016, through OpenFOAM's `isoAdvection`, with PLIC
reconstruction), an energy equation, and an evaporation/condensation model (`HardtWondra`) in
which the interface is moved by the **total velocity** `u_t = u + u_PC`.

> **Status: research code.** Validated against analytic solutions for the 1D Stefan problem and a
> short axisymmetric d²-law droplet run (see [Validation](#validation)). Other case generators are
> included, but their results are not summarised here. Compare any new case with a reference before
> relying on it.

## Contents

1. [Features](#features)
2. [Notation and conventions](#notation-and-conventions)
3. [Governing equations](#governing-equations)
4. [Interface advection with phase change](#interface-advection-with-phase-change)
5. [Interface geometry and curvature](#interface-geometry-and-curvature)
6. [Phase-change model `HardtWondra`](#phase-change-model-hardtwondra)
7. [Solution algorithm](#solution-algorithm)
8. [Time-step control](#time-step-control)
9. [Repository layout](#repository-layout)
10. [Build](#build)
11. [Case set-up](#case-set-up)
12. [Output and log diagnostics](#output-and-log-diagnostics)
13. [Validation](#validation)
14. [Limitations](#limitations)
15. [References](#references)

---

## Features

- Geometric VOF with the **isoAdvector** method (Roenby et al. 2016; OpenFOAM class `isoAdvection`) and
  PLIC reconstruction (`plicRDF`); without phase change alpha stays bounded and volume is conserved
  to about 1e-11 (relative).
- Phase change through the total interface velocity `u_t = u + u_PC` (Gada & Sharma 2009, Eq. 13;
  Shaikh et al. 2016, Eq. 8), with the interface normal taken from the reconstructed distance
  function (RDF), not from `grad(alpha)`.
- `HardtWondra` phase-change model: ghost-fluid Tsat condition in the temperature equation, mass
  flux from pure-cell probes (Malan et al. 2021), shifted dilatation (Hardt & Wondra 2008).
- Interface curvature from the PLIC signed distance: `gradAlpha`, `RDF`, `heightFunction`, `constant`.
- Energy equation with a flux consistent with the geometric alpha flux (`rhoCpPhi`).
- Runs in parallel; builds against the stock v2412 libraries (no extra library).

## Notation and conventions

| Symbol | Meaning |
|---|---|
| phase 1, phase 2 | first and second name in `phases (...)`; `alpha.<phase1>` is the solved field, `alpha1 = alpha` |
| `rho_k, cp_k, nu_k, k_k` | density, specific heat, kinematic viscosity, conductivity of phase k |
| `rho = alpha rho1 + (1-alpha) rho2`, `rhoCp = alpha rho1 cp1 + (1-alpha) rho2 cp2` | mixture properties |
| `phi` | bulk face volume flux, `u` bulk velocity |
| `psi` | signed distance to the PLIC interface (`psiRDF`), **positive in phase 1** |
| `n = grad(psi)/\|grad(psi)\|` | interface normal, pointing into phase 1 |
| `a_G = \|S_PLIC\|/V` | PLIC interface area density of a cell (`aInterface`) [1/m] |
| `mdot'''` | volumetric mass-transfer rate [kg/m³/s], **> 0 converts phase 1 into phase 2** |
| `mdot'' = mdot'''/a_G` | interfacial mass flux [kg/m²/s] |
| `L = h_lv` | latent heat, phase 1 → 2 [J/kg] |
| `PCV` | dilatation source [1/s], `div(u) = PCV` |

With phase 1 the liquid, `mdot > 0` is evaporation and `mdot < 0` condensation. The equations
below are those of the code; the sections name the files that implement them.

## Governing equations

Both phases are incompressible and Newtonian with constant properties, described by one velocity and
one pressure field (single-fluid formulation).

**Continuity** (`pEqn.H`). Across the interface, with mass flux `mdot''`, the normal velocity jumps by
`mdot'' (1/rho2 - 1/rho1)`. The model spreads this into a smooth volumetric source `PCV`
(see [Shifted dilatation](#shifted-dilatation)):

$$\nabla\cdot\mathbf{u} = \mathrm{PCV}$$

**Momentum** (`UEqn.H`). With `p_rgh = p - rho g·x` and the CSF surface force:

$$\frac{\partial(\rho\mathbf{u})}{\partial t} + \nabla\cdot(\rho_\phi\,\mathbf{u}\mathbf{u})
- \left(\frac{\partial\rho}{\partial t}+\nabla\cdot\rho_\phi\right)\mathbf{u}
= -\nabla p_{rgh} - (\mathbf{g}\cdot\mathbf{x})\nabla\rho + \sigma\kappa\nabla\alpha
+ \nabla\cdot\boldsymbol{\tau}_{\mathrm{eff}}$$

`rho_phi` is the mass flux from the advection step (`rhoPhi`). The second term on the left removes the
cell-wise mass imbalance that phase change creates (a sink in interface cells and a dilatation
elsewhere). `kappa` is the curvature of [Interface geometry](#interface-geometry-and-curvature),
`sigma` comes from the stock `surfaceTensionModel` (constant or temperature dependent; only the normal
capillary force is modelled, not a tangential Marangoni stress). Turbulence uses the stock
VoF-phase-incompressible models; laminar is the tested case.

**Pressure equation** (`pEqn.H`), with `rAU = 1/A(UEqn)`, face value `rAUf`:

$$\nabla\cdot\left(r_{AU,f}\nabla p_{rgh}\right) = \nabla\cdot\phi_{HbyA} - \mathrm{PCV},\qquad
\phi_{HbyA} = \phi[\mathbf{H}/A] + \left(\sigma\kappa\,\mathrm{snGrad}\,\alpha
- (\mathbf{g}\cdot\mathbf{x})\,\mathrm{snGrad}\,\rho\right) r_{AU,f}\,|\mathbf{S}|$$

followed by the usual flux and velocity corrections. The total dilatation does not vanish when the
densities differ, so the domain needs an open boundary through which the generated volume leaves.

**Energy** (`TEqn.H`). For the temperature, in conservative form with a continuity correction:

$$\frac{\partial(\rho c_p T)}{\partial t} + \nabla\cdot(\phi_{\rho c_p}\,T)
- \left(\frac{\partial(\rho c_p)}{\partial t} + \nabla\cdot\phi_{\rho c_p}\right)T
- \nabla\cdot(k_f\nabla T) + \mathcal{S}_{\Gamma}(T) = 0$$

- `phi_rhoCp` (`rhoCpPhi`) is the enthalpy-capacity flux built from the same geometric alpha flux as
  the mass flux: `rhoCpPhi = (rho1 cp1 - rho2 cp2) alphaPhi_T + rho2 cp2 phi` with the phase-change part
  of the carrier flux removed (next section).
- `k_f` is the model's face conductivity: for `HardtWondra` the phase conductivity `k1` or `k2`, **zero
  on faces whose two cells lie on different sides of the interface**; for `noPhaseChange` an arithmetic
  mixture of `k_k = rho_k nu_k cp_k / Pr_k`.
- `S_Gamma(T)` is the interface term of the model: the implicit ghost-fluid Dirichlet condition
  `a (T - T_sat)` at the interface ([below](#temperature-ghost-fluid-condition)). It removes the latent
  heat `L * mdot` from the interface cells by construction. No other latent-heat source exists.

## Interface advection with phase change

Implemented in `alphaEqn.H`, `alphaEqnSubCycle.H`, `alphaPCVelocity.H`, `createFields.H`.

The interface is advected with the total velocity (Gada & Sharma 2009, Eq. 13), `u_t = u + u_PC`:

$$\frac{\partial\alpha}{\partial t} + \mathbf{u}_t\cdot\nabla\alpha = 0,\qquad
\mathbf{u}_{PC} = \frac{\dot m''}{\rho_1}\,\mathbf{n},\quad
\mathbf{n} = \frac{\nabla\psi}{|\nabla\psi|}$$

isoAdvector (`isoAdvection`) solves the flux form, so with the face fluxes `phiPC = u_PC·S` and
`phiT = phi + phiPC` the equation actually solved is

$$\frac{\partial\alpha}{\partial t} + \nabla\cdot(\alpha\,\phi_T) = \alpha\,\nabla\cdot\phi_{PC}$$

- The right-hand side is required. Without it the flux form only redistributes volume
  (`div(alpha phiPC)` integrates to a boundary term). It is passed to `advect()` as the explicit
  source `Su = clip(alpha, 0, 1) div(phiPC)` (`Sp = 0`).
- In a pure cell the face value equals the cell value, the two terms cancel exactly, and the volume
  change takes place only where `alpha_f` differs from the cell value, i.e. at the interface. Overshoots
  of interface cells are handled by the `isoAdvection` bounding (conservative limiting of the fluxes,
  then clipping and snapping).
- The bulk term `alpha div(u)` is **not** included; the dilatation enters through `PCV` in the
  pressure equation, where the model places it in non-interface cells of the lighter phase.

**Coefficient `1/rho1`.** `HardtWondra` shifts the dilatation off the interface, so the bulk velocity at
the interface is the liquid velocity and the interface moves at `mdot''/rho1` relative to it. The
Gada & Sharma coefficient `0.5 (1/rho1 + 1/rho2)` assumes the dilatation sits at the interface;
it was tried in the 1D Stefan problem and gave +26.8 % error in the vapour-layer thickness against
+4.6 % for `1/rho1`, so only `1/rho1` is kept.

**Face and cell fields** (`alphaPCVelocity.H`):

1. `mdot''' = -rho1 alpha1Gen` (`alpha1Gen = -mdot'''/rho1`, see the model) and
   `mdot'' = mdot'''/a_G` in cells that carry the source (cells without a PLIC facet use `a_G = 1/h`).
2. `n` from `fvc::grad(psiRDF)`, normalised, interpolated to the faces and normalised again.
3. Face flux `phiPC = (mdot''_f/rho1)(n_f·S_f)`, where `mdot''_f` is the mean of the neighbouring cells
   that carry a source (weights of the linear interpolation, ratio of the interpolated source to the
   interpolated mask). Non-coupled boundary faces get zero.
4. Cell velocity `U_PC = (mdot''/rho1) n`, used by `isoAdvection` to time the isoface crossing of faces.

`isoAdvection` is built on `phiT` and `UT = U + U_PC` (`createFields.H`); both are refreshed before
every `advect()` call. On a moving mesh `U` is made relative first, as in `interIsoFoam`.

**Mass and enthalpy fluxes.** `getRhoPhi` uses `phiT`, so the carrier part is removed again:

```
rhoPhi   = getRhoPhi(rho1, rho2)       - rho2   * phiPC
rhoCpPhi = getRhoPhi(rho1 cp1, rho2 cp2) - rho2 cp2 * phiPC
```

Globally `d(rho)/dt + div(rhoPhi) = (rho1-rho2) alpha div(phiPC) + rho2 PCV` sums to zero:
`-(rho1-rho2) mdot/rho1 + rho2 mdot (1/rho2 - 1/rho1) = 0`. The balance is exact over the domain but
not cell by cell, because `PCV` sits away from the interface; this is what the non-conservative term
in the momentum equation compensates.

**Source mode.** `phaseChangeAdvection source;` replaces all of the above with the previous scheme
(`alphaSuSp.H`): the explicit source `alpha1Gen` is distributed over neighbouring cells within the
phase-1 volume budget of a time step (up to 5 passes) and `phiPC = 0`. It is kept for comparison.

**Sub-cycling and outer iterations.** With `nAlphaSubCycles > 1` the mass and enthalpy fluxes are
time-averaged over the sub-cycles. With `nOuterCorrectors > 1` the second and later outer iterations
re-advect alpha from the start-of-step field with `0.5 (U + U_prev)` and `0.5 (phi + phi_prev)`.

## Interface geometry and curvature

Implemented in `thermalPhaseChangeModels/interfaceCurvatureITF.[CH]`. It is evaluated once per outer
iteration from the PLIC reconstruction of the **advected** alpha, and provides `psiRDF`, `aInterface`,
`kappaI` and the surface force.

**Interface cells.** A cell is a genuine interface cell if it carries a PLIC facet and has a face
neighbour on the other side of `alpha = 0.5`. Cells with a facet but no such neighbour ("wisps") are
counted and ignored. `a_G = |S_PLIC|/V` on genuine cells.

**Signed distance `psi`** on the interface cells and three point-neighbour rings around them:
- ring 0: plane distance of the cell centre to its own facet;
- ring `L` = 1..3: from the facets carried by rings `< L`. Within the closest facet distance plus one
  cell size, the distance is the weighted plane distance
  `psi = sum(w (x - c_j)·n_j) / sum(w)`, `w = ((x - c_j)/|x - c_j| · n_j)^2`
  (the weighting of `reconstructedDistanceFunction`); the closest facet is carried on to the next ring;
- the sign is set by a vote over pure band cells so that `psi > 0` in phase 1.

**Curvature models** (`curvatureModel` in `transportProperties`):

| model | curvature |
|---|---|
| `gradAlpha` (default) | stock `interfaceProperties`: `kappa = -div(n_hat_f)` with the stock `surfaceTensionForce` |
| `RDF` | `kappa_j = -div(grad(psi)/\|grad(psi)\|)` on rings 0–1 (least-squares gradient). Interface cell: distance-weighted mean over its stencil of `kappa_j/(1 + kappa_j psi_j/n_c)`, weights `1/(\|psi_j\| + 0.1 h)`, `n_c` = number of curvature directions (1 for planar 2D, 2 for axisymmetric and 3D); stencil values with denominator < 0.5 are skipped. Ring-1 cells take the value of the nearest interface cell. |
| `heightFunction` | 3 × 7 column height function with the non-uniform three-point second derivative, for planar 2D and axisymmetric wedge meshes (including the hoop term and the axis mirror); cells without a valid stencil keep the `RDF` value. In 3D it is the `RDF` value everywhere. |
| `constant` | prescribed `constantCurvature` on rings 0–1 (verification) |

**Surface force.** For all models but `gradAlpha`,
`F_f = sigma_f kappa_f snGrad(alpha)` with `kappa_f = (w_P kappa_P + w_N kappa_N)/(w_P + w_N)`,
`w = alpha(1-alpha) + 1e-6`.

**Capillary time-step limit** (Brackbill–Kothe–Zemach):
`dt_c = sqrt(rho_avg h_min^3 / (2 pi sigma_max))`, `h_min` the smallest in-plane cell size
(`V` over the largest face area, wedge and empty faces excluded).

## Phase-change model `HardtWondra`

Implemented in `thermalPhaseChangeModels/HardtWondra.[CH]` on top of the base class
`thermalPhaseChangeModel` (runtime selection through `phaseChangeProperties`). In this model
conduction is cut at the interface, the interface temperature is `T_sat`, and the mass flux follows
from the heat flux into the interface.

### Which side a cell is on

A cell is in phase 1 if `psi > 1e-6 h`; if `|psi|` is smaller (interface through the cell centre) the
side is decided by `alpha >= 0.5 - 1e-8`.

### Temperature: ghost-fluid condition

On every face `f` joining cells `P` and `N` on different sides, the position of the interface along the
cell-centre line is

```
theta_P = psi_P / (psi_P - psi_N)        (from alpha if psi has no sign change)
theta_N = 1 - theta_P,   both limited below by thetaMin
```

The conduction through `f` is removed (`k_f = 0`) and each cell gets an implicit Dirichlet term,
added to the left side of the temperature equation (Gibou et al. 2002):

$$a_P (T_P - T_{sat}),\qquad a_P = \frac{k_P\,|S_f|}{\theta_P\,d_{PN}}$$

with `k_P` the conductivity of the phase of `P` and `d_PN` the centre distance; likewise for `N`.
Coupled (processor) faces are treated the same way with the neighbour values. A heat term
`heatReturn` from the previous step (below) is added explicitly to the source.

### Mass flux (after the temperature is converged)

`correctAfterT()` runs at the last outer iteration and gives the `mdot'''` used in the **next** step
(one-step lag).

*Ghost-fluid estimate.* The heat leaving the two sides of each crossing face is
`Q_f = a_P (T_P - T_sat) + a_N (T_N - T_sat)` [W], `mdot_f = Q_f/L`. It is assigned to the interface
cells of the pair (split by `theta` if both are interface cells, otherwise to the one that is; with
neither, to the phase-1 cell). This is exactly energy consistent but uses mixed-cell temperatures
divided by `theta`.

*Probe estimate (default, `fluxMethod probe`, Malan et al. 2021).* For every interface cell the heat flux
into the interface is taken from **pure** cells (`alpha < 1e-3` or `> 1 - 1e-3`, not interface cells):

- probe value `G = (T - T_sat)/d` with `d = max(|psi|, 0.05 h)`; if a farther pure cell of the same phase
  exists (`d_2 - d_1 > 0.3 h`) the quadratic through `(0, T_sat)`, `(d_1, T_1)`, `(d_2, T_2)` is used
  for a second-order gradient;
- for each phase the probes in the interface cell's stencil are averaged with weights
  `w = |n·dx| / |dx|^3` (`n` the PLIC normal, `dx` the cell-to-cell vector);
- `q = k1 G_1 + k2 G_2` [W/m²], `mdot'' = q/L`, `mdot''' = mdot'' a_G`;
- an interface cell with probes on only one side keeps its ghost-fluid value.

The heat drained by the Dirichlet term but not turned into mass,
`L (mdot_GF V - mdot_probe V)`, is returned to the probe cells in the next step through `heatReturn`
(weights `w k_phase`), so the temperature equation loses exactly `L mdot`.

### Sources derived from `mdot'''`

- Alpha source (used in source mode; in velocity mode it defines `mdot''` via `alphaPCVelocity.H`):
  `alpha1Gen = -mdot'''/rho1` [1/s].
- `Q_pc = L mdot'''` [W/m³] (diagnostic).

### Shifted dilatation

`calcPCV()` follows Hardt & Wondra (2008): `mdot'''` is smoothed with a Helmholtz filter

$$s - \lambda^2\nabla^2 s = \dot m''',\qquad \lambda = \texttt{lambdaCells}\cdot\bar h$$

(`bar h` the mean cell size of the source cells; solved with `fvSolution` entry `mdotSmoothHW` if present,
else PCG/DIC to 1e-14). The dilatation is placed only in **non-interface cells of the lighter phase**:

$$\mathrm{PCV}_i = N\,s_i\left(\frac{1}{\rho_2}-\frac{1}{\rho_1}\right),\qquad
N = \frac{\sum \dot m'''V}{\sum_{\text{light, non-interface}} s\,V}$$

so that `sum(PCV V) = sum(mdot''' V)(1/rho2 - 1/rho1)`. `PCV = 0` for equal densities or no mass
transfer. The normalisation assumes one connected interface; several interface regions give a warning.

### Maximum interface speed

`maxInterfaceSpeed() = max |mdot'''|/(rho1 a_G)` (with `a_G = 1/h` where no facet exists) feeds the
phase-change time-step limit.

### `noPhaseChange`

No mass or energy transfer: `Q_pc = 0`, `PCV = 0`, `alpha1Gen = 0`, conduction with the `Pr`-based
face conductivity. Used when `constant/phaseChangeProperties` is absent.

## Solution algorithm

Per time step (`interTempFoam.C`):

```
setDeltaT                      Co, alpha Co, capillary, phase-change limits
PIMPLE loop
  mesh.update()                (only for dynamic meshes)
  alphaEqnSubCycle
      u_PC, phiPC, UT, phiT    (alphaPCVelocity.H; uses mdot''' and psi of the previous step)
      isoAdvection::advect     PLIC at start-of-step alpha, flux, bounding
      rhoPhi, rhoCpPhi, rho, rhoCp
  surf.reconstruct()           PLIC of the advected alpha
  mixture.correct()
  interface.correct()          genuine cells, psi (3 rings), a_G, curvature
  RDF.constructRDF()           RDF band for the solver's own distance field
  UEqn.H                       momentum predictor (optional)
  pEqn.H  (nCorrectors)        pressure with PCV, flux/velocity correction
  update rho, rhoCp
  TEqn.H                       energy, including ghost-fluid Tsat terms
  correctAfterT()   [final outer iteration]   mdot''', PCV for the next step
  turbulence correct
```

`Phase-change coupling` is explicit in time: alpha, momentum and energy of step `n+1` use `mdot'''` and
`PCV` evaluated from the converged temperature of step `n`.

## Time-step control

With `adjustTimeStep yes`, the new step is limited by (`setDeltaT.H`):

- `maxCo` (flow Courant number) and `maxAlphaCo` (interface Courant number, computed from
  `phi + phiPC`) with the usual damped increase (at most 1.2 per step);
- the capillary limit `maxCapillaryNum * dt_c` (default `maxCapillaryNum = 1`); a warning is issued if a
  fixed step exceeds it;
- the phase-change limit `maxPhaseChangeCo * h_min / max(|mdot''|/rho1)` (default `maxPhaseChangeCo =
  0.25`);
- `maxDeltaT`.

Local time stepping (`localEuler`) is not supported (isoAdvection and the phase-change source need a
uniform step); the solver stops with an error if it is selected.

## Repository layout

```
applications/solvers/multiphase/interTempFoam/interTempFoam/
├── interTempFoam.C           main loop
├── createFields.H            fields, advector (phiT, UT, phiPC), phase-change model, keywords
├── alphaControls.H, alphaEqnSubCycle.H, alphaEqn.H, alphaCourantNo.H
├── alphaPCVelocity.H         u_PC, phiPC, Su (default mode)
├── alphaSuSp.H               explicit alpha source (phaseChangeAdvection source)
├── UEqn.H, pEqn.H, TEqn.H, correctPhi.H, initCorrectPhi.H, rhofs.H, setDeltaT.H
├── thermalPhaseChangeModels/ compiled into the solver
│   ├── thermalPhaseChangeModel.[CH], newThermalPhaseChangeModel.C   base class + selector
│   ├── HardtWondra.[CH]                                             phase-change model
│   ├── noPhaseChange.[CH]                                           disabled model
│   └── interfaceCurvatureITF.[CH]                                   psiRDF, a_G, curvature, surface force
├── validation/               case generators and run scripts (see Validation)
└── report.md, report.pdf
```

## Build

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc      # or your v2412 environment
cd applications/solvers/multiphase/interTempFoam/interTempFoam
wmake
```

The executable is installed to `$FOAM_USER_APPBIN`; `Make/options` uses `$(LIB_SRC)` only. The warnings
that appear come from OpenFOAM's own headers (`-Woverloaded-virtual`).

## Case set-up

Case generators in `validation/` write complete working cases (`gen.py`, with `Allrun` scripts);
copy a generator to a run directory before use, they write into their own directory.

### `constant/transportProperties`

```
phases (water vapour);                  // phase 1 first: the liquid for evaporation
water  { transportModel Newtonian; nu 2.94e-7; rho 958; cp 4216; Pr 1.75; }
vapour { transportModel Newtonian; nu 2.05e-5; rho 0.6; cp 2080; Pr 1.02; }
sigma           0.059;
curvatureModel  heightFunction;         // gradAlpha | RDF | heightFunction | constant
// constantCurvature 1/R;               // only with curvatureModel constant
```

`cp` and `Pr` are required per phase. `Pr` only sets the conductivity of `noPhaseChange`;
`HardtWondra` uses `k1`, `k2`.

### `constant/phaseChangeProperties` (optional)

Without the file the model is `noPhaseChange`.

| keyword | default | meaning |
|---|---|---|
| `model` | – | `HardtWondra` or `noPhaseChange` |
| `T_sat` | – | saturation temperature [K] |
| `h_lv` | – | latent heat, phase 1 → 2 [J/kg] |
| `k1`, `k2` | – | conductivity of phase 1, 2 [W/m/K] (`HardtWondra`) |
| `thetaMin` | 0.05 | lower limit of the interface position `theta` on a crossing face |
| `lambdaCells` | 1 | dilatation smoothing length in cell sizes |
| `fluxMethod` | `probe` | `probe` or `ghostFluid` |
| `DilatationSource` | true | include `PCV` in the pressure equation |
| `PhaseFractionSource` | true | include `alpha1Gen` (source mode) / `mdot'''` in the alpha advection |

### `system/fvSolution`

```
solvers
{
    "alpha.water.*"
    {
        reconstructionScheme plicRDF;
        plicRDFCoeffs { tol 1e-6; relTol 0.1; iterations 5; interpolateNormal true; }
        nAlphaBounds   3;   snapTol 1e-12;   clip true;      // stock isoAdvection controls
        nAlphaCorr     1;   nAlphaSubCycles 1;
        cAlpha         1;                    // read by interfaceProperties, unused by isoAdvection
        phaseChangeAdvection velocity;       // velocity (default) | source
    }
    p_rgh   { solver PCG; preconditioner DIC; tolerance 1e-10; relTol 0.01; }
    p_rghFinal { $p_rgh; relTol 0; }
    "(U|T).*" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0; }
    // mdotSmoothHW  { ... }                // optional, for the dilatation smoothing
}
PIMPLE { momentumPredictor no; nOuterCorrectors 1; nCorrectors 3; nNonOrthogonalCorrectors 0;
         pRefCell 0; pRefValue 0; }
```

### `system/fvSchemes` (minimum)

```
ddtSchemes      { default Euler; }
gradSchemes     { default Gauss linear; }
divSchemes
{
    div(rhoPhi,U)   Gauss linear;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
    div(rhoCpPhi,T) Gauss limitedLinear 1;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes   { default corrected; }
```

### `system/controlDict` (solver-specific keys)

```
application      interTempFoam;
adjustTimeStep   yes;
maxCo            0.2;
maxAlphaCo       0.2;       // required
maxCapillaryNum  0.5;       // default 1
maxPhaseChangeCo 0.25;      // default 0.25
maxDeltaT        1e-3;
```

### Initial and boundary fields

`0/alpha.<phase1>`, `0/U`, `0/p_rgh`, `0/T` (all required). For a closed or nearly closed domain make
sure there is an open boundary (e.g. `inletOutlet`/`totalPressure`) where the generated volume can leave.
`constant/g`, `constant/turbulenceProperties` (e.g. `simulationType laminar`) as for `interFoam`.

### Running

```bash
blockMesh && interTempFoam                         # serial
decomposePar && mpirun -np 4 interTempFoam -parallel && reconstructPar
```

## Output and log diagnostics

Written fields (besides `alpha.*`, `U`, `p_rgh`, `T`): `psiRDF` (signed distance), `kappaI`
(interface curvature), `mdot` (`mdot'''`) and `PCV`, for `HardtWondra`.

Log lines that are useful for checking a run:

| line | meaning |
|---|---|
| `alphaPCVelocity: max\|u_PC\| ..., requested ... m3/s` | largest phase-change velocity; requested volume rate `sum(alpha1Gen V)` |
| `alpha source [m3/s]: requested ..., implicit form ..., applied ...` | requested rate, rate handed to `advect()` and rate actually applied (after bounding) |
| `isoAdvection: After conservative bounding: min/max(alpha)` | boundedness before clipping |
| `interfaceCurvatureITF: interface cells = ..., wisps ignored = ...` | genuine interface cells and isolated facets |
| `interfaceCurvatureITF: kappa ... min/mean/max` | curvature over interface cells |
| `HardtWondra: ghost-fluid faces = ..., clipped ..., min raw theta` | crossing faces and theta clipping |
| `HardtWondra: probe flux: mdotTotal = ... (ghost-fluid estimate ...)` | probe versus ghost-fluid mass rate |
| `HardtWondra: sum(PCV V) = ..., target ...` | dilatation integral versus its target |
| `Total mass sum(rho V)` | total mass |

## Validation

Case generators for these cases are in `validation/` (run scripts `run*.sh`, one per case):
static droplet (spurious currents and Laplace jump), 1D Stefan melting, 1D Stefan vaporisation and sucking
interface, axisymmetric Scriven bubble growth, axisymmetric d²-law droplet, bubble rise in a superheated
liquid, 2D film boiling (Welch–Wilson type) and the Sun et al. (2014) film-boiling set-up.

Results verified for the current code:

| Case | Set-up | Result |
|---|---|---|
| 1D Stefan vaporisation | 50 cells, 10 s, water/steam | vapour-layer thickness +4.63 % against the analytic solution at t = 10 s |
| d²-law droplet | axisymmetric wedge, R0/h = 8, 2 ms | evaporation constant K −0.08 % against the middle of the analytic range; mass rate −0.54 % |

The d²-law run is short (D² changes by about 6e-4 in it), so it checks the instantaneous evaporation rate
rather than a developed d²-law. The explicit-source path (`phaseChangeAdvection source`) gives +4.61 % and
−0.50 % on the same two cases. The static-droplet generator (no phase change) is the check for spurious
currents and the Laplace jump of a chosen curvature model at your resolution.

## Limitations

- `HardtWondra` is the only phase-change model besides `noPhaseChange`; condensation (`mdot < 0`) follows
  from the sign of the heat flux but has not been validated.
- The mass transfer uses the converged temperature of the previous step (one-step lag).
- The velocity form needs a model whose `mdot'''` sits in interface cells.
- The dilatation normalisation assumes one connected interface (warning otherwise); with different phase
  densities the domain needs an outlet.
- The coefficient of `u_PC` is fixed to `1/rho1` (see above); a model that places the dilatation at the
  interface would need `0.5 (1/rho1 + 1/rho2)`.
- Spurious currents of the CSF force are not removed; check them with the static droplet at your resolution.
- Local time stepping is not supported. Mesh motion and turbulence are inherited from `interFoam` and are
  not validated here. Parallel runs were tested on up to 4 processors.

## References

- Gada & Sharma (2009), Numer. Heat Transfer B 56(4), 307–322.
- Shaikh, Sharma & Bhardwaj (2016), Int. J. Heat Mass Transfer 96, 458–473.
- Hardt & Wondra (2008), J. Comput. Phys. 227, 5871–5895.
- Gibou, Chen, Fedkiw (2002), J. Comput. Phys. 176, 205–227.
- Malan et al. (2021), J. Comput. Phys. 426, 109920.
- Roenby, Bredmose, Jasak (2016), isoAdvector, Royal Society Open Science 3, 160405.
