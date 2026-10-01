#!/usr/bin/env python3
"""
Build a self-contained OpenFOAM case (no Python needed to run it) for the
2D film boiling problem of

  D. Sun, J. Xu, Q. Chen (2014), Modeling of the evaporation and condensation
  phase-change problems with FLUENT, Numer. Heat Transfer B 66:326-342,
  Sec. 3.2, Figs. 6-10.

Setup (Sun et al. 2014):
  vapour rho 5, mu 0.005, k 1, cp 200; liquid rho 200, mu 0.1;
  h_fg 1e4 J/kg, sigma 0.1 N/m, g 9.81 m/s2
  lambda0 = 2 pi sqrt(3 sigma/((rho_l - rho_v) g))          (eq. 31)
  domain lambda0/2 x 3 lambda0, grid 100 x 600 (their grid-independent one)
  bottom: no-slip wall, T_w = T_sat + 5 K; sides: symmetry;
  top: pressure outlet; U = 0, liquid at T_sat, vapour T linear from the
  wall to the interface y = lambda0/128 (4 + cos(2 pi x/lambda0)) (eq. 32).
  Liquid k: equal to the vapour one (their model assumption 1, Sec. 4);
  liquid cp (irrelevant, liquid stays saturated) = 200; T_sat not given in
  the paper (only differences enter): 500 K.
Reference result: time-averaged wall Nusselt number (eq. 33)
  Nu = 2/lambda0 int_0^{lambda0/2} l0/(T_w - T_sat) dT/dy|_0 dx,
  l0 = sqrt(sigma/((rho_l - rho_v) g))                      (eq. 34)
  Sun et al.: 1.83; Klimenko (1981) correlation: 1.91 (computed here).

Usage: make_case.py <outdir> [nx=100] [endTime=3.0] [nWrites=50] [nProcs=8]
"""
import math
import os
import sys

out = sys.argv[1]
nx = int(sys.argv[2]) if len(sys.argv) > 2 else 100
tEnd = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0
nWrites = int(sys.argv[4]) if len(sys.argv) > 4 else 50
nProcs = int(sys.argv[5]) if len(sys.argv) > 5 else 8

rhoL, muL, kL, cpL = 200.0, 0.1, 1.0, 200.0
rhoV, muV, kV, cpV = 5.0, 0.005, 1.0, 200.0
SIGMA, L, G = 0.1, 1e4, 9.81
TSAT, DT = 500.0, 5.0
LAM = 2*math.pi*math.sqrt(3*SIGMA/((rhoL - rhoV)*G))
L0 = math.sqrt(SIGMA/((rhoL - rhoV)*G))
ny = 6*nx
X, Y = LAM/2, 3*LAM
h = X/nx

# Klimenko (1981), Int. J. Heat Mass Transfer 24:69: Nu = 0.19 (Gr Pr)^(1/3) f1,
# f1 = 0.89 Ja^(1/3) for Ja = h_fg/(cp_v dT) > 0.71
Gr = rhoV*(rhoL - rhoV)*G*L0**3/muV**2
Pr = muV*cpV/kV
Ja = L/(cpV*DT)
f1 = 0.89*Ja**(1/3) if Ja > 0.71 else 1.0
NU_KLIMENKO = 0.19*(Gr*Pr)**(1/3)*f1

HDR = r"""/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Version:  v2106
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
"""


def foam(cls, obj, loc=None):
    s = HDR + 'FoamFile\n{\n    version     2.0;\n    format      ascii;\n'
    s += f'    class       {cls};\n'
    if loc:
        s += f'    location    "{loc}";\n'
    s += f'    object      {obj};\n}}\n' + '// ' + '* '*37 + '//\n\n'
    return s


END = '\n// ' + '*'*73 + ' //\n'


def w(rel, text):
    p = os.path.join(out, rel)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, 'w').write(text)


def y0(x):
    return LAM/128*(4 + math.cos(2*math.pi*x/LAM))


def nonuniform(vals):
    return 'nonuniform List<scalar>\n%d\n(\n%s\n)' % (
        len(vals), '\n'.join('%.12g' % v for v in vals))


# ---- initial fields, blockMesh cell order: i (x) fastest, then j (y)
NS = 64
alpha, T = [], []
for j in range(ny):
    for i in range(nx):
        xc, yc = (i + 0.5)*h, (j + 0.5)*h
        if (j + 1)*h < LAM/128*3 - 1e-12:
            vap = 1.0
        elif j*h > LAM/128*5 + 1e-12:
            vap = 0.0
        else:
            vap = 0.0
            for k in range(NS):
                xs = i*h + (k + 0.5)*h/NS
                vap += min(max((y0(xs) - j*h)/h, 0.0), 1.0)
            vap /= NS
        alpha.append(1 - vap)
        yi = y0(xc)
        T.append(TSAT + DT*(1 - yc/yi) if yc < yi else TSAT)


def field(cls, name, dims, internal, wall, outlet):
    return foam(cls, name, '0') + f"""dimensions      {dims};

internalField   {internal};

boundaryField
{{
    wall
    {{
        {wall}
    }}
    outlet
    {{
        {outlet}
    }}
    left
    {{
        type            symmetryPlane;
    }}
    right
    {{
        type            symmetryPlane;
    }}
    frontAndBack
    {{
        type            empty;
    }}
}}
""" + END


Tw = TSAT + DT
w('0/U', field('volVectorField', 'U', '[0 1 -1 0 0 0 0]', 'uniform (0 0 0)',
               'type            noSlip;',
               'type            pressureInletOutletVelocity;\n        value           uniform (0 0 0);'))
w('0/p_rgh', field('volScalarField', 'p_rgh', '[1 -1 -2 0 0 0 0]', 'uniform 0',
                   'type            fixedFluxPressure;\n        value           uniform 0;',
                   'type            prghTotalPressure;\n        p0              uniform 0;\n        value           uniform 0;'))
w('0/T', field('volScalarField', 'T', '[0 0 0 1 0 0 0]', nonuniform(T),
               f'type            fixedValue;\n        value           uniform {Tw};',
               f'type            inletOutlet;\n        inletValue      uniform {TSAT};\n        value           uniform {TSAT};'))
w('0/alpha.liquid', field('volScalarField', 'alpha.liquid', '[0 0 0 0 0 0 0]', nonuniform(alpha),
                          'type            zeroGradient;',
                          'type            inletOutlet;\n        inletValue      uniform 1;\n        value           uniform 1;'))

# ---- mesh
w('system/blockMeshDict', foam('dictionary', 'blockMeshDict', 'system') + f"""// lambda0 = {LAM*1e3:.6f} mm; domain lambda0/2 x 3 lambda0 (Sun et al. 2014, Fig. 6b)
scale   1;

vertices
(
    (0 0 0) ({X:.10g} 0 0) ({X:.10g} {Y:.10g} 0) (0 {Y:.10g} 0)
    (0 0 {h:.10g}) ({X:.10g} 0 {h:.10g}) ({X:.10g} {Y:.10g} {h:.10g}) (0 {Y:.10g} {h:.10g})
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({nx} {ny} 1) simpleGrading (1 1 1)
);

boundary
(
    wall
    {{
        type wall;
        faces ((0 1 5 4));
    }}
    outlet
    {{
        type patch;
        faces ((3 7 6 2));
    }}
    left
    {{
        type symmetryPlane;
        faces ((0 4 7 3));
    }}
    right
    {{
        type symmetryPlane;
        faces ((1 2 6 5));
    }}
    frontAndBack
    {{
        type empty;
        faces ((0 3 2 1) (4 5 6 7));
    }}
);
""" + END)

# ---- properties
w('constant/g', foam('uniformDimensionedVectorField', 'g', 'constant')
  + f'dimensions      [0 1 -2 0 0 0 0];\nvalue           (0 -{G} 0);\n' + END)
w('constant/turbulenceProperties', foam('dictionary', 'turbulenceProperties', 'constant')
  + 'simulationType  laminar;\n' + END)
w('constant/transportProperties', foam('dictionary', 'transportProperties', 'constant') + f"""// Sun, Xu & Chen (2014) Numer. Heat Transfer B 66:326, Sec. 3.2
phases (liquid vapour);

liquid
{{
    transportModel  Newtonian;
    nu              {muL/rhoL:.8g};      // mu 0.1 / rho 200
    rho             {rhoL};
    cp              {cpL};              // not given; liquid stays saturated
    Pr              {muL*cpL/kL:.8g};    // k_l = k_v = 1 (model assumption 1)
}}

vapour
{{
    transportModel  Newtonian;
    nu              {muV/rhoV:.8g};      // mu 0.005 / rho 5
    rho             {rhoV};
    cp              {cpV};
    Pr              {muV*cpV/kV:.8g};    // k 1
}}

sigma           {SIGMA};

// Height-function curvature (2D planar), balanced-force CSF
curvatureModel  heightFunction;
""" + END)
w('constant/phaseChangeProperties', foam('dictionary', 'phaseChangeProperties', 'constant') + f"""// Sharp ghost-fluid phase change; phase 1 = liquid, phase 2 = vapour
model           HardtWondra;

T_sat           {TSAT};     // not given in the paper; only differences enter
h_lv            {L:g};      // [J/kg]
k1              {kL};       // liquid [W/m/K]
k2              {kV};       // vapour [W/m/K]
""" + END)

# ---- numerics
w('system/fvSchemes', foam('dictionary', 'fvSchemes', 'system') + """ddtSchemes
{
    default         Euler;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    div(rhoPhi,U)   Gauss linearUpwind grad(U);
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
    div(rhoCpPhi,T) Gauss limitedLinear 1;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}
""" + END)
w('system/fvSolution', foam('dictionary', 'fvSolution', 'system') + """solvers
{
    "alpha.liquid.*"
    {
        reconstructionScheme plicRDF;
        plicRDFCoeffs
        {
            tol         1e-6;
            relTol      0.1;
            iterations  5;
            interpolateNormal true;
        }
        nAlphaBounds    3;
        snapTol         1e-12;
        clip            true;
        nAlphaCorr      1;
        nAlphaSubCycles 1;
        cAlpha          1;     // read by interfaceProperties, unused by isoAdvector
    }

    "pcorr.*"
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-10;
        relTol          0;
    }

    p_rgh
    {
        solver          GAMG;
        smoother        DIC;
        tolerance       1e-9;
        relTol          0.01;
    }

    p_rghFinal
    {
        $p_rgh;
        relTol          0;
    }

    "(U|T).*"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-10;
        relTol          0;
    }

    "mdotSmoothHW.*"
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-14;
        relTol          0;
    }
}

PIMPLE
{
    momentumPredictor   no;
    nOuterCorrectors    1;
    nCorrectors         3;
    nNonOrthogonalCorrectors 0;
}
""" + END)

cap = math.sqrt(0.5*(rhoL + rhoV)*h**3/(2*math.pi*SIGMA))
w('system/controlDict', foam('dictionary', 'controlDict', 'system') + f"""application     interTempFoam;

startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {tEnd:g};

// Capillary (Brackbill-Kothe-Zemach) limit sqrt(rho_avg h^3/(2 pi sigma))
// = {cap:.4e} s for h = {h*1e3:.4f} mm; with maxCapillaryNum 0.5 the solver
// runs at ~{0.5*cap:.2e} s (about {tEnd/(0.5*cap):.0f} steps to endTime).
deltaT          {0.5*cap:.4e};
adjustTimeStep  yes;
maxCo           0.2;
maxAlphaCo      0.2;
maxCapillaryNum 0.5;
maxPhaseChangeCo 0.25;
maxDeltaT       1e-3;

writeControl    adjustable;
writeInterval   {tEnd/nWrites:g};     // {nWrites} writes
purgeWrite      0;
writeFormat     binary;
writePrecision  10;
writeCompression off;
timeFormat      general;
timePrecision   8;
runTimeModifiable yes;

functions
{{
    // vapour volume V(t) = domain volume - volIntegrate(alpha.liquid)
    liquidVolume
    {{
        type            volFieldValue;
        libs            (fieldFunctionObjects);
        regionType      all;
        operation       volIntegrate;
        fields          (alpha.liquid);
        writeFields     false;
        writeControl    timeStep;
        writeInterval   20;
    }}

    // wall temperature gradient for the Nusselt number (eq. 33):
    // Nu = -l0/(T_w - T_sat) * areaAverage(grad(T)_y) on the wall
    gradT
    {{
        type            grad;
        libs            (fieldFunctionObjects);
        field           T;
        result          grad(T);
        executeControl  timeStep;
        executeInterval 20;
        writeControl    writeTime;
    }}
    wallGradT
    {{
        type            surfaceFieldValue;
        libs            (fieldFunctionObjects);
        regionType      patch;
        name            wall;
        operation       areaAverage;
        fields          (grad(T));
        writeFields     false;
        writeControl    timeStep;
        writeInterval   20;
    }}
}}
""" + END)
w('system/decomposeParDict', foam('dictionary', 'decomposeParDict', 'system') + f"""numberOfSubdomains {nProcs};

// Horizontal slabs (built-in method, needs no scotch/metis)
method          simple;

coeffs
{{
    n           (1 {nProcs} 1);
}}
""" + END)

w('Allrun', f"""#!/bin/sh
cd "${{0%/*}}" || exit                                # Run from this directory
. ${{WM_PROJECT_DIR:?}}/bin/tools/RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------
# 2D film boiling, Sun, Xu & Chen (2014), Sec. 3.2. Initial fields in 0/ are
# precomputed for the {nx} x {ny} mesh (blockMesh cell order). Keep 0/ intact.

runApplication blockMesh
runApplication checkMesh
runApplication decomposePar
runParallel $(getApplication)
runApplication reconstructPar

#------------------------------------------------------------------------------
""")
w('Allclean', """#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory
. ${WM_PROJECT_DIR:?}/bin/tools/CleanFunctions      # Tutorial clean functions
#------------------------------------------------------------------------------
# keeps the precomputed 0/ fields
rm -rf processor* postProcessing log.* constant/polyMesh
foamListTimes -rm 2>/dev/null || true

#------------------------------------------------------------------------------
""")
os.chmod(os.path.join(out, 'Allrun'), 0o755)
os.chmod(os.path.join(out, 'Allclean'), 0o755)

w('post_Nu.py', f'''#!/usr/bin/env python3
"""Nusselt number (Sun et al. 2014 eq. 33) and V/V0 from postProcessing.
Nu = -l0/dT * areaAverage(dT/dy) on the wall; reference: Sun 1.83,
Klimenko 1.91."""
import glob, sys
L0, DT, AREA, VDOM = {L0!r}, {DT!r}, {X*h!r}, {X*Y*h!r}
rows = []
for f in sorted(glob.glob('postProcessing/wallGradT/*/surfaceFieldValue.dat')):
    for l in open(f):
        if l.startswith('#') or not l.strip():
            continue
        t = l.split()
        v = [float(x) for x in l.replace('(', ' ').replace(')', ' ').split()]
        rows.append((v[0], -L0/DT*v[2]))
rows.sort()
vol = []
for f in sorted(glob.glob('postProcessing/liquidVolume/*/volFieldValue.dat')):
    for l in open(f):
        if l.startswith('#') or not l.strip():
            continue
        t, a = map(float, l.split()[:2])
        vol.append((t, VDOM - a))
vol.sort()
tmin = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
late = [n for t, n in rows if t >= tmin]
with open('Nu_history.csv', 'w') as f:
    f.write('t,Nu\\n')
    for t, n in rows:
        f.write(f'{{t}},{{n}}\\n')
print(f'samples {{len(rows)}}, t_end {{rows[-1][0] if rows else 0}}')
if late:
    print(f'time-averaged Nu (t >= {{tmin}} s) = {{sum(late)/len(late):.4f}}'
          f'   (Sun et al. 2014: 1.83, Klimenko: {NU_KLIMENKO:.3f})')
if vol:
    print(f'V/V0 at t_end = {{vol[-1][1]/vol[0][1]:.4f}}')
''')

open(os.path.join(out, 'README'), 'w').write(f"""2D film boiling - Sun, Xu & Chen (2014), Numer. Heat Transfer B 66:326, Sec. 3.2

Solver: interTempFoam (isoAdvector + plicRDF, height-function curvature,
sharp ghost-fluid phase change, model HardtWondra).

Mesh {nx} x {ny}, h = {h*1e3:.4f} mm, lambda0 = {LAM*1e3:.4f} mm,
domain lambda0/2 x 3 lambda0. endTime {tEnd:g} s, {nWrites} writes
(writeInterval {tEnd/nWrites:g} s). Capillary time-step limit
{cap:.3e} s -> dt ~ {0.5*cap:.2e} s with maxCapillaryNum 0.5.

Run:       ./Allrun            (decomposition: {nProcs} slabs in y, system/decomposeParDict)
Restart:   startFrom latestTime is set; rerun the solver step only.
Result:    python3 post_Nu.py [tmin=1.0]   -> time-averaged wall Nu
           reference: Sun et al. 1.83; Klimenko (1981) correlation {NU_KLIMENKO:.3f}
""")
print(f'case {out}: {nx}x{ny}, h={h*1e3:.4f} mm, lambda0={LAM*1e3:.4f} mm, '
      f'dt_cap={cap:.3e}, Nu_Klimenko={NU_KLIMENKO:.3f}, writeInterval={tEnd/nWrites:g}')
