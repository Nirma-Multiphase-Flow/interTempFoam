#!/usr/bin/env python3
"""
V7: vapour bubble rising and growing in superheated ethanol, axisymmetric
(3D 5-deg wedge about the y axis; y = up).

Experiment of Florschuetz, Henry & Rashid Khan (1969), Int. J. Heat Mass
Transfer 12:1465, as simulated by Bures & Sato and by Long, Pan & Zaleski,
EBIT, arXiv:2402.13677, Sec. 3.5 (setup, properties, Table 1 from there):
  ethanol, 1 atm: rho_l 757, mu_l 4.29e-4, k_l 0.154, cp_l 3000,
                  rho_v 1.435, mu_v 1.04e-5, k_v 0.02, cp_v 1830,
                  T_sat 351.45 K, T_inf 354.55 K, h_lg 9.63e5, sigma 0.018
  domain r in [0, 4 mm], y in [0, 20 mm]; bottom and side no-slip walls,
  outflow at the top; g = 9.81 m/s2 in -y.
  initial bubble D0 = 420 um centred at y = 1 mm, the Scriven solution at
  t = 5.6 ms; T from the Scriven profile; run to t = 0.0856 s.
Metric: normalised radius R/(2 beta sqrt(alpha_l)) [sqrt(s)] vs t (t
counted from the Scriven time of the initial bubble, as in EBIT), compared
with the experiment (EBIT Table 1).

Mesh: uniform h in the core r < 2 mm, y < 12 mm (the bubble rises ~7 mm and
reaches R ~ 1.2 mm), geometrically graded outside. Walls: T zeroGradient.

Usage: gen.py <h_um> [init]
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

rhoL, muL, kL, cpL = 757.0, 4.29e-4, 0.154, 3000.0
rhoV, muV, kV, cpV = 1.435, 1.04e-5, 0.02, 1830.0
Tsat, Tinf, L, SIGMA = 351.45, 354.55, 9.63e5, 0.018
dT = Tinf - Tsat
alphaL = kL/(rhoL*cpL)
eps = 1 - rhoV/rhoL
R0, Y0 = 210e-6, 1e-3
RDOM, YDOM = 4e-3, 20e-3
RCORE, YCORE = 2e-3, 12e-3
TEND = 0.0856

EXPERIMENT = [  # t [s], R/(2 beta sqrt(alpha_l)) [sqrt(s)]  (EBIT Table 1)
    (0.0147, 0.122), (0.0239, 0.157), (0.0371, 0.197), (0.0474, 0.256),
    (0.0586, 0.306), (0.0667, 0.364), (0.0755, 0.425),
]


def integral(beta, a, n=4000):
    if a >= 1:
        return 0.0
    f = lambda x: 0.0 if x >= 1 else math.exp(
        -beta*beta*((1 - x)**-2 - 2*eps*x - 1))
    hh = (1 - a)/n
    s = f(a) + f(1)
    for i in range(1, n):
        s += (4 if i % 2 else 2)*f(a + i*hh)
    return s*hh/3


def beta_root():
    lhs = rhoL*cpL*dT/(rhoV*(L + (cpL - cpV)*dT))
    g = lambda b: 2*b*b*integral(b, 0.0) - lhs
    lo, hi = 0.5, 40.0
    for _ in range(100):
        mid = 0.5*(lo + hi)
        if g(lo)*g(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5*(lo + hi)


BETA = beta_root()
T0 = (R0/(2*BETA))**2/alphaL      # Scriven time of the initial bubble


def T_scriven(r, R):
    if r <= R:
        return Tsat
    return Tinf - 2*BETA**2*(rhoV*(L + (cpL - cpV)*dT)/(rhoL*cpL)) \
        * integral(BETA, 1 - R/r, 400)


def foam(o):
    hdr = open(os.path.join(HERE, 'system', 'fvSchemes')).read()
    hdr = hdr[:hdr.index('FoamFile')]
    return hdr + ('FoamFile\n{\n    version     2.0;\n    format      '
                  'ascii;\n    class       dictionary;\n    object      '
                  '%s;\n}\n\n' % o)


def outer_cells(h, Lout):
    """geometric stretching from h: number of cells and total expansion"""
    n = max(4, int(math.ceil(Lout/(4*h))))
    n = min(n, 40)
    lo, hi = 1.0000001, 3.0
    for _ in range(200):
        q = 0.5*(lo + hi)
        if h*(q**n - 1)/(q - 1) > Lout:
            hi = q
        else:
            lo = q
    return n, q**(n - 1)


def mesh(hum):
    h = hum*1e-6
    nr = int(round(RCORE/h))
    ny = int(round(YCORE/h))
    nro, gr = outer_cells(h, RDOM - RCORE)
    nyo, gy = outer_cells(h, YDOM - YCORE)
    th = math.radians(2.5)
    c, s = RDOM*math.cos(th), RDOM*math.sin(th)
    fr, fy = RCORE/RDOM, YCORE/YDOM
    body = f"""convertToMeters 1;
mergeType points;
vertices
(
    (0 0 0) ({c} 0 {-s}) ({c} {YDOM} {-s}) (0 {YDOM} 0)
    ({c} 0 {s}) ({c} {YDOM} {s})
);
blocks
(
    hex (0 1 2 3 0 4 5 3) ({nr + nro} {ny + nyo} 1)
    simpleGrading
    (
        (({fr:.8f} {nr} 1) ({1 - fr:.8f} {nro} {gr:.6f}))
        (({fy:.8f} {ny} 1) ({1 - fy:.8f} {nyo} {gy:.6f}))
        1
    )
);
boundary
(
    axis   {{ type empty; faces ((0 3 3 0)); }}
    bottom {{ type wall; faces ((0 1 4 0)); }}
    side   {{ type wall; faces ((1 2 5 4)); }}
    outlet {{ type patch; faces ((3 2 5 3)); }}
    wedge1 {{ type wedge; faces ((0 1 2 3)); }}
    wedge2 {{ type wedge; faces ((0 4 5 3)); }}
);
"""
    open(os.path.join(HERE, 'system', 'blockMeshDict'), 'w').write(
        foam('blockMeshDict') + body)
    open(os.path.join(HERE, 'system', 'controlDict'), 'w').write(
        foam('controlDict') + f"""
application     interTempFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {TEND - T0:.8g};
deltaT          1e-7;
writeControl    adjustable;
writeInterval   0.005;
purgeWrite      0;
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
        fields          (alpha.ethanol);
        writeFields     false;
        writeControl    timeStep;
        writeInterval   10;
    }}
    Umax
    {{
        type            fieldMinMax;
        libs            (fieldFunctionObjects);
        mode            magnitude;
        fields          (U);
        writeControl    timeStep;
        writeInterval   10;
    }}
}}
""")
    with open(os.path.join(HERE, 'case.info'), 'w') as f:
        f.write(f'h {h}\nbeta {BETA:.8f}\nt0 {T0:.8g}\nalphaL {alphaL:.8g}\n'
                f'nr {nr}\nny {ny}\n')
    print(f'V7 h={hum} um core {nr}x{ny} outer {nro}x{nyo} grading {gr:.2f}/{gy:.2f} '
          f'beta={BETA:.5f} t0_Scriven={T0*1e3:.3f} ms (EBIT: 5.6 ms)')


def init():
    from post import read_field
    C = read_field(os.path.join(HERE, '0', 'C'), vector=True)
    info = dict(l.split() for l in open(os.path.join(HERE, 'case.info')))
    h = float(info['h'])
    NS = 100
    al, T = [], []
    for c in C:
        r = math.hypot(c[0], c[2])
        y = c[1]
        rs = math.hypot(r, y - Y0)
        if abs(rs - R0) < 2*h:
            r0 = math.floor(r/h + 1e-9)*h
            y0 = math.floor(y/h + 1e-9)*h
            num = den = 0.0
            for k in range(NS):
                rr = r0 + (k + 0.5)*h/NS
                yc = math.sqrt(R0*R0 - rr*rr) if rr < R0 else 0.0
                lo, hi = max(y0, Y0 - yc), min(y0 + h, Y0 + yc)
                if hi > lo:
                    num += rr*(hi - lo)
                den += rr*h
            vap = num/den
        else:
            vap = 1.0 if rs < R0 else 0.0
        al.append(1 - vap)
        T.append(T_scriven(rs, R0) if rs < 10*R0 else Tinf)
    for name, vals in (('alpha.ethanol', al), ('T', T)):
        p = os.path.join(HERE, '0', name)
        s = open(p).read()
        body = 'nonuniform List<scalar> %d\n(\n%s\n)\n;' % (
            len(vals), '\n'.join('%.12g' % v for v in vals))
        s = re.sub(r'internalField\s+[^;]*;', 'internalField   ' + body, s,
                   count=1, flags=re.S)
        open(p, 'w').write(s)
    print(f'init: {len(C)} cells, T range {min(T):.4f}..{max(T):.4f}')


if __name__ == '__main__':
    if len(sys.argv) > 2 and sys.argv[2] == 'init':
        init()
    else:
        mesh(float(sys.argv[1]))
