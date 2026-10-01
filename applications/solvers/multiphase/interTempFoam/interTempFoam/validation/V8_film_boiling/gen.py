#!/usr/bin/env python3
"""
V8: two-dimensional film boiling (Rayleigh-Taylor), planar.

Setup and properties ("fluid B") from arXiv:2407.09969 (momentum balance
correction, one-fluid VOF boiling), Sec. 6 and Table 3, following Guo et al.
(2011) Numer. Heat Transfer A 59:857, Sun et al. (2014) Numer. Heat
Transfer B 66:326, Boyd & Ling (2023) Comput. Fluids 254:105807:
  rho_L 200, rho_G 5, mu_L 0.1, mu_G 0.005, k_L = k_G = 1, cp_L = cp_G = 200,
  sigma 0.1, h_LV 1e4, T_sat = 1 K, T_wall = T_sat + 5 K, g = 9.81 (-y)
  lambda_d = 2 pi sqrt(3 sigma/((rho_L - rho_G) g)) = 78.6844 mm
  domain [0, lambda_d/2] x [0, lambda_d]; symmetry at x = 0, lambda_d/2;
  no-slip wall at y = 0; outflow at y = lambda_d
  film y0(x) = lambda_d/128 (4 + cos(2 pi x/lambda_d)), linear T in film,
  liquid at T_sat.
Metrics: V/V0 and wall Nusselt number
  Nu = 2/lambda_d int_0^{lambda_d/2} l0/(T_w - T_sat) dT/dy|_{y=0} dx,
  l0 = sqrt(sigma/((rho_L - rho_G) g)),
vs time; time average vs Berenson (1961), J. Heat Transfer 83:351:
  Nu_B = 0.425 [rho_G (rho_L - rho_G) g h_LV l0^3/(mu_G k_G dT)]^(1/4).

Usage: gen.py <nx> [init]     (ny = 2 nx; M3 of the paper: nx = 200)
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
rhoL, rhoG, muL, muG = 200.0, 5.0, 0.1, 0.005
kL, kG, cpL, cpG = 1.0, 1.0, 200.0, 200.0
SIGMA, L, TSAT, DT, G = 0.1, 1e4, 1.0, 5.0, 9.81
LAM = 2*math.pi*math.sqrt(3*SIGMA/((rhoL - rhoG)*G))
L0 = math.sqrt(SIGMA/((rhoL - rhoG)*G))
TEND = 0.5
NU_BERENSON = 0.425*(rhoG*(rhoL - rhoG)*G*L*L0**3/(muG*kG*DT))**0.25


def y0(x):
    return LAM/128*(4 + math.cos(2*math.pi*x/LAM))


def foam(cls, o):
    hdr = open(os.path.join(HERE, 'system', 'fvSchemes')).read()
    hdr = hdr[:hdr.index('FoamFile')]
    return hdr + ('FoamFile\n{\n    version     2.0;\n    format      '
                  'ascii;\n    class       %s;\n    object      %s;\n}\n\n'
                  % (cls, o))


def field(cls, name, dims, internal, b):
    return foam(cls, name) + f"""dimensions      {dims};
internalField   {internal};
boundaryField
{{
    wall            {{ {b[0]} }}
    outlet          {{ {b[1]} }}
    left            {{ type symmetryPlane; }}
    right           {{ type symmetryPlane; }}
    frontAndBack    {{ type empty; }}
}}
"""


def mesh(nx):
    ny = 2*nx
    X, Y = LAM/2, LAM
    dz = X/nx
    w = lambda rel, s: open(os.path.join(HERE, rel), 'w').write(s)
    w('system/blockMeshDict', foam('dictionary', 'blockMeshDict') + f"""convertToMeters 1;
vertices
(
    (0 0 0) ({X} 0 0) ({X} {Y} 0) (0 {Y} 0)
    (0 0 {dz}) ({X} 0 {dz}) ({X} {Y} {dz}) (0 {Y} {dz})
);
blocks ( hex (0 1 2 3 4 5 6 7) ({nx} {ny} 1) simpleGrading (1 1 1) );
boundary
(
    wall   {{ type wall; faces ((0 1 5 4)); }}
    outlet {{ type patch; faces ((3 7 6 2)); }}
    left   {{ type symmetryPlane; faces ((0 4 7 3)); }}
    right  {{ type symmetryPlane; faces ((1 2 6 5)); }}
    frontAndBack {{ type empty; faces ((0 3 2 1) (4 5 6 7)); }}
);
""")
    os.makedirs(os.path.join(HERE, '0'), exist_ok=True)
    Tw = TSAT + DT
    w('0/U', field('volVectorField', 'U', '[0 1 -1 0 0 0 0]', 'uniform (0 0 0)',
      ('type noSlip;', 'type pressureInletOutletVelocity; value uniform (0 0 0);')))
    w('0/p_rgh', field('volScalarField', 'p_rgh', '[1 -1 -2 0 0 0 0]', 'uniform 0',
      ('type fixedFluxPressure; value uniform 0;',
       'type prghTotalPressure; p0 uniform 0; value uniform 0;')))
    w('0/T', field('volScalarField', 'T', '[0 0 0 1 0 0 0]', f'uniform {TSAT}',
      (f'type fixedValue; value uniform {Tw};',
       f'type inletOutlet; inletValue uniform {TSAT}; value uniform {TSAT};')))
    w('0/alpha.liquid', field('volScalarField', 'alpha.liquid', '[0 0 0 0 0 0 0]',
      'uniform 1', ('type zeroGradient;',
                    'type inletOutlet; inletValue uniform 1; value uniform 1;')))
    w('constant/g', foam('uniformDimensionedVectorField', 'g')
      + f'dimensions      [0 1 -2 0 0 0 0];\nvalue           (0 -{G} 0);\n')
    w('constant/transportProperties', foam('dictionary', 'transportProperties') + f"""// fluid B, arXiv:2407.09969 Table 3
phases (liquid vapour);
liquid
{{
    transportModel  Newtonian;
    nu              {muL/rhoL:.8g};
    rho             {rhoL};
    cp              {cpL};
    Pr              {muL*cpL/kL:.8g};
}}
vapour
{{
    transportModel  Newtonian;
    nu              {muG/rhoG:.8g};
    rho             {rhoG};
    cp              {cpG};
    Pr              {muG*cpG/kG:.8g};
}}
sigma           {SIGMA};
curvatureModel  heightFunction;
""")
    w('constant/phaseChangeProperties', foam('dictionary', 'phaseChangeProperties') + f"""model           HardtWondra;
T_sat           {TSAT};
h_lv            {L};
k1              {kL};
k2              {kG};
""")
    w('system/controlDict', foam('dictionary', 'controlDict') + f"""application     interTempFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {TEND};
deltaT          1e-6;
writeControl    adjustable;
writeInterval   0.025;
writeFormat     ascii;
writePrecision  10;
timeFormat      general;
timePrecision   8;
runTimeModifiable no;
adjustTimeStep  yes;
maxCo           0.2;
maxAlphaCo      0.2;
maxCapillaryNum 0.5;
maxPhaseChangeCo 0.25;
maxDeltaT       1e-3;

functions
{{
    liquidVolume
    {{
        type            volFieldValue;
        libs            (fieldFunctionObjects);
        regionType      all;
        operation       volIntegrate;
        fields          (alpha.liquid);
        writeFields     false;
        writeControl    timeStep;
        writeInterval   10;
    }}
}}
""")
    with open(os.path.join(HERE, 'case.info'), 'w') as f:
        f.write(f'nx {nx}\nlambda {LAM}\nl0 {L0}\nNuB {NU_BERENSON}\n')
    print(f'V8 nx={nx} lambda_d={LAM*1e3:.4f} mm (paper 78.6844) '
          f'h={X/nx*1e3:.4f} mm Nu_Berenson={NU_BERENSON:.3f}')


def init(nx):
    ny = 2*nx
    X = LAM/2
    h = X/nx
    NS = 50
    al, T = [], []
    for j in range(ny):
        for i in range(nx):
            xc, yc = (i + 0.5)*h, (j + 0.5)*h
            vap = 0.0
            for k in range(NS):
                x = i*h + (k + 0.5)*h/NS
                vap += min(max((y0(x) - j*h)/h, 0.0), 1.0)
            vap /= NS
            al.append(1 - vap)
            yi = y0(xc)
            T.append(TSAT + DT*(1 - yc/yi) if yc < yi else TSAT)
    for name, vals in (('alpha.liquid', al), ('T', T)):
        p = os.path.join(HERE, '0', name)
        s = open(p).read()
        body = 'nonuniform List<scalar> %d\n(\n%s\n)\n;' % (
            len(vals), '\n'.join('%.12g' % v for v in vals))
        s = re.sub(r'internalField\s+[^;]*;', 'internalField   ' + body, s,
                   count=1, flags=re.S)
        open(p, 'w').write(s)
    print(f'init: {len(vals)} cells')


if __name__ == '__main__':
    n = int(sys.argv[1])
    (init if len(sys.argv) > 2 and sys.argv[2] == 'init' else mesh)(n)
