# interTempFoam

A non-isothermal two-phase VOF solver for OpenFOAM v2412, extending `interFoam` with a coupled energy equation for heat transfer across immiscible fluid interfaces.

---

## Overview

`interTempFoam` solves for two incompressible, immiscible fluids using a Volume of Fluid (VOF) phase-fraction based interface capturing approach, with an additional energy equation to account for convective and conductive heat transfer. The solver is suited for problems such as:

- Boiling and condensation (without phase change)
- Thermal stratification in two-phase systems
- Droplet heating/cooling
- Conjugate heat transfer across fluid interfaces
- Thermocapillary (Marangoni) flow (with temperature-dependent surface tension)

---

## Governing Equations

The solver resolves the following coupled system:

### Phase fraction (VOF)

```
∂α/∂t + ∇·(αU) + ∇·[α(1−α)Uᵣ] = 0
```

where `α` is the phase fraction of fluid 1, and `Uᵣ` is the compression velocity at the interface. Solved using MULES (Multidimensional Universal Limiter with Explicit Solution) with optional sub-cycling.

### Continuity

```
∇·U = 0
```

### Momentum

```
∂(ρU)/∂t + ∇·(ρUU) = −∇p* − g·x∇ρ + ∇·(μ_eff ∇U) + f_σ
```

where `p*` is the modified pressure (excluding hydrostatic), `f_σ` is the surface tension body force via the CSF model, and `μ_eff` includes turbulent contributions.

### Energy

```
∂(ρCpT)/∂t + ∇·(ρCpφ T) − ∇·(κ ∇T) = 0
```

where:
- `ρCp = α·ρ₁·Cp₁ + (1−α)·ρ₂·Cp₂` — volume-weighted mixture heat capacity
- `κ = α·κ₁ + (1−α)·κ₂` — mixture thermal conductivity (derived from Prandtl numbers)
- `ρCpφ` — enthalpy flux at cell faces, computed consistently with the phase-fraction flux

The enthalpy flux `rhoCpPhi` is updated alongside `rhoPhi` after each alpha equation solve, including within sub-cycles, ensuring thermodynamic consistency with the VOF advection.

---

## Repository Structure

```
interTempFoam/
├── applications/
│   └── solvers/
│       └── multiphase/
│           └── interTempFoam/
│               └── interTempFoam/
│                   ├── interTempFoam.C       # Main solver
│                   ├── createFields.H        # Field declarations (U, p_rgh, T, rhoCp, rhoCpPhi)
│                   ├── TEqn.H                # Energy equation
│                   ├── UEqn.H                # Momentum equation
│                   ├── pEqn.H                # Pressure equation
│                   ├── alphaEqn.H            # VOF equation (MULES, rhoCpPhi update)
│                   ├── alphaEqnSubCycle.H    # Sub-cycling with rhoCpPhi time-averaging
│                   ├── alphaCourantNo.H      # Alpha-based CFL check
│                   ├── createAlphaFluxes.H   # Alpha flux field setup
│                   ├── correctPhi.H          # Flux correction for moving meshes
│                   ├── initCorrectPhi.H
│                   ├── rhofs.H               # Face-interpolated densities
│                   ├── setDeltaT.H           # Adaptive time-stepping
│                   ├── setRDeltaT.H          # LTS reciprocal delta-T
│                   ├── alphaSuSp.H
│                   ├── Make/
│                   │   ├── files
│                   │   └── options
│                   ├── interMixingFoam/      # Three-phase variant (not primary)
│                   └── overInterDyMFoam/     # Overset mesh variant (not primary)
│
└── src/
    └── transportModels/
        ├── Allwmake                          # Build script for all transport libs
        ├── incompressible/                   # Extended two-phase mixture model
        │   └── incompressibleTwoPhaseMixture/
        │       ├── incompressibleTwoPhaseMixture.H  # Adds cp, Pr, kappa accessors
        │       └── incompressibleTwoPhaseMixture.C
        ├── immiscibleIncompressibleTwoPhaseMixture/  # VOF mixture + interface
        ├── interfaceProperties/              # Surface tension, contact angle models
        │   └── surfaceTensionModels/
        │       └── temperatureDependent/     # σ(T) model for Marangoni flows
        ├── twoPhaseMixture/
        ├── twoPhaseProperties/               # Contact angle boundary conditions
        ├── compressible/
        └── geometricVoF/                     # isoAdvection schemes
```

---

## Prerequisites

- **OpenFOAM v2412** (openfoam.com release). Other versions may require adaptation.
- A working OpenFOAM environment with `WM_PROJECT_DIR` and `FOAM_USER_APPBIN` set.
- GCC or compatible C++14 compiler.

Source your OpenFOAM environment before building:

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
# or wherever your installation lives, e.g.:
# source $HOME/OpenFOAM/OpenFOAM-v2412/etc/bashrc
```

---

## Compilation

Build in two steps: first the custom transport model libraries, then the solver binary. Both must be compiled from within the repository root.

```bash
cd /path/to/interTempFoam
```

### Step 1 — Build transport model libraries

```bash
cd src/transportModels
./Allwmake
cd ../..
```

This compiles (in order):
1. `twoPhaseMixture`
2. `interfaceProperties` (includes temperature-dependent surface tension)
3. `twoPhaseProperties`
4. `incompressible` (extended with `cp`, `Pr`, `kappa` accessors)
5. `compressible`
6. `immiscibleIncompressibleTwoPhaseMixture`
7. `geometricVoF`

Libraries are installed to `$FOAM_USER_LIBBIN`.

### Step 2 — Build the solver

```bash
cd applications/solvers/multiphase/interTempFoam/interTempFoam
wmake
```

The `interTempFoam` executable is installed to `$FOAM_USER_APPBIN`.

### Clean rebuild (if needed)

```bash
# Clean transport models
cd src/transportModels
./Allwmake -clean
./Allwmake
cd ../..

# Clean and rebuild solver
cd applications/solvers/multiphase/interTempFoam/interTempFoam
wclean
wmake
```

### Verify the build

```bash
which interTempFoam
interTempFoam --help
```

---

## Running a Case

### Required initial fields (`0/`)

| Field   | Type            | Description                          |
|---------|-----------------|--------------------------------------|
| `U`     | `volVectorField`| Velocity [m/s]                       |
| `p_rgh` | `volScalarField`| Modified pressure [Pa] (p − ρgh)     |
| `alpha.water` | `volScalarField` | Phase fraction [0,1]          |
| `T`     | `volScalarField`| Temperature [K]                      |

### `constant/transportProperties`

In addition to the standard `interFoam` entries, `interTempFoam` requires thermal properties for each phase:

```cpp
phases (water air);

water
{
    transportModel  Newtonian;
    nu              1e-06;       // kinematic viscosity [m²/s]
    rho             1000;        // density [kg/m³]
    cp              4182;        // specific heat capacity [J/(kg·K)]
    Pr              7.0;         // Prandtl number [-]
}

air
{
    transportModel  Newtonian;
    nu              1.48e-05;
    rho             1.0;
    cp              1005;
    Pr              0.71;
}

sigma           0.07;            // surface tension [N/m]
```

Thermal conductivity is computed internally as `κ = μ·Cp/Pr`.

### `constant/turbulenceProperties`

```cpp
simulationType  laminar;
```

Or specify a turbulence model as in standard `interFoam`.

### `system/fvSchemes` (energy-relevant entries)

```cpp
ddtSchemes
{
    default         Euler;
}

divSchemes
{
    div(rhoCpPhi,T) Gauss linearUpwind grad(T);
    // ...standard interFoam div schemes...
}
```

### Running

```bash
interTempFoam
# or in parallel:
decomposePar
mpirun -np 4 interTempFoam -parallel
reconstructPar
```

---

## Key Features

- **Consistent enthalpy flux**: `rhoCpPhi` is derived from the same alpha flux used for `rhoPhi`, guaranteeing energy conservation consistency with the VOF advection.
- **Sub-cycle support**: When `nAlphaSubCycles > 1`, `rhoCpPhi` is time-averaged over sub-cycles identically to `rhoPhi`.
- **Temperature-dependent surface tension**: The `temperatureDependent` surface tension model (`interfaceProperties/surfaceTensionModels/temperatureDependent/`) allows Marangoni-driven flows.
- **PIMPLE loop**: Supports outer correctors, pressure-velocity coupling, and turbulence correction within the standard PIMPLE framework.
- **Dynamic mesh**: Inherits `interFoam`'s dynamic mesh / adaptive re-meshing capability.
- **LTS (Local Time Stepping)**: Supported for steady-state acceleration.

---

## Differences from `interFoam`

| Feature | `interFoam` | `interTempFoam` |
|---------|-------------|-----------------|
| Energy equation | No | Yes — `fvm::ddt(rhoCp,T) + fvm::div(rhoCpPhi,T) - fvm::laplacian(kappaf,T)` |
| `rhoCpPhi` field | No | Yes — enthalpy flux, updated with `rhoPhi` |
| `rhoCp` field | No | Yes — volume-weighted mixture heat capacity |
| `transportProperties` | `nu`, `rho` | Adds `cp`, `Pr` per phase |
| Temperature-dependent σ | No (base OF model) | Available via `temperatureDependent` model |

---

## References

- Natali, A. et al. (2014). *Non-isothermal two-phase flow simulation using VOF method.* — basis for the enthalpy flux formulation.
- Liu, J. et al. (2011). *Heat transfer in two-phase VOF simulations.* — reference for `rhoCpPhi` consistency.
- OpenFOAM `interFoam` solver — base implementation this solver extends.
- Rusche, H. (2002). *Computational Fluid Dynamics of Dispersed Two-Phase Flows at High Phase Fractions.* PhD thesis, Imperial College London.

---

## License

This project is distributed under the GNU General Public License v3 — the same license as OpenFOAM. See the individual source file headers for copyright notices.
