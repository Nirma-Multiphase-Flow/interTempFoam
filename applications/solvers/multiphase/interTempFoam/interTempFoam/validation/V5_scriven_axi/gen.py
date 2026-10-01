#!/usr/bin/env python3
"""
V5: Scriven (1959) vapour bubble growth in superheated water, axisymmetric
(3D 5-deg wedge about the y axis, symmetry plane y = 0: half bubble).

Setup of Tanguy et al. (2014) JCP 264:1-22, values as in
basilisk.fr/sandbox/ecipriano/run/scrivenproblem.c (see
constant/transportProperties): Ja = 3, R0 = 1 mm, grown to Rf = 2 mm.
Domain 6 mm (Basilisk uses 12 mm; the thermal layer is < 0.3 mm, the
bubble ends at 2 mm).

Usage: gen.py <R0/h> [mode]
  mode 'mesh' : blockMeshDict + controlDict (before blockMesh)
  mode 'init' : alpha.water and T in 0/ (needs 0/C)
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

# Properties (Tanguy et al. 2014)
rhoL, rhoV = 958.0, 0.59
cpL, cpV = 4216.0, 2034.0
kL = 0.6
L = 2.257e6
Tsat = 373.0
Ja = 3.0
dT = L*rhoV*Ja/(rhoL*cpL)
Tbulk = Tsat + dT
alphaL = kL/(rhoL*cpL)
eps = 1 - rhoV/rhoL
R0, Rf = 1e-3, 2e-3
Ldom = 6e-3


def integral(beta, a, n=4000):
    # int_a^1 exp(-beta^2((1-x)^-2 - 2 eps x - 1)) dx, Simpson
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
    lo, hi = 0.5, 20.0
    for _ in range(100):
        mid = 0.5*(lo + hi)
        if g(lo)*g(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5*(lo + hi)


BETA = beta_root()
TSHIFT = (R0/(2*BETA))**2/alphaL
TEND = Rf*Rf/(4*BETA*BETA*alphaL) - TSHIFT


def R_exact(t):
    return 2*BETA*math.sqrt(alphaL*(t + TSHIFT))


def T_exact(r, R):
    if r <= R:
        return Tsat
    return Tbulk - 2*BETA**2*(rhoV*(L + (cpL - cpV)*dT)/(rhoL*cpL)) \
        * integral(BETA, 1 - R/r, 400)


def mesh(n):
    h = R0/n
    N = int(round(Ldom/h))
    th = math.radians(2.5)
    c, s = Ldom*math.cos(th), Ldom*math.sin(th)
    body = f"""convertToMeters 1;
mergeType points;
vertices
(
    (0 0 0) ({c} 0 {-s}) ({c} {Ldom} {-s}) (0 {Ldom} 0)
    ({c} 0 {s}) ({c} {Ldom} {s})
);
blocks ( hex (0 1 2 3 0 4 5 3) ({N} {N} 1) simpleGrading (1 1 1) );
boundary
(
    axis   {{ type empty; faces ((0 3 3 0)); }}
    bottom {{ type symmetryPlane; faces ((0 1 4 0)); }}
    outlet {{ type patch; faces ((1 2 5 4) (3 2 5 3)); }}
    wedge1 {{ type wedge; faces ((0 1 2 3)); }}
    wedge2 {{ type wedge; faces ((0 4 5 3)); }}
);
"""
    hdr = open(os.path.join(HERE, 'system', 'fvSchemes')).read()
    hdr = hdr[:hdr.index('FoamFile')]
    foam = lambda o: hdr + ('FoamFile\n{\n    version     2.0;\n    format      '
        'ascii;\n    class       dictionary;\n    object      %s;\n}\n\n' % o)
    open(os.path.join(HERE, 'system', 'blockMeshDict'), 'w').write(
        foam('blockMeshDict') + body)
    nOut = 10
    open(os.path.join(HERE, 'system', 'controlDict'), 'w').write(foam('controlDict') + f"""
application     interTempFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {TEND:.8g};
deltaT          1e-6;
writeControl    adjustable;
writeInterval   {TEND/nOut:.8g};
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
        fields          (alpha.water);
        writeFields     false;
        writeControl    timeStep;
        writeInterval   1;
    }}
    Umax
    {{
        type            fieldMinMax;
        libs            (fieldFunctionObjects);
        mode            magnitude;
        fields          (U);
        writeControl    timeStep;
        writeInterval   1;
    }}
}}
""")
    with open(os.path.join(HERE, 'case.info'), 'w') as f:
        f.write(f'Rh {n}\nh {h}\nbeta {BETA:.8f}\ntshift {TSHIFT:.8g}\n'
                f'tend {TEND:.8g}\nTbulk {Tbulk:.8f}\n')
    print(f'V5 R0/h={n} h={h:.3e} N={N} beta={BETA:.6f} tshift={TSHIFT:.5g} '
          f'tend={TEND:.5g} Tbulk={Tbulk:.6f}')


def init(n):
    from post import read_field
    h = R0/n
    C = read_field(os.path.join(HERE, '0', 'C'), vector=True)
    NS = 200
    al, T = [], []
    for c in C:
        r = math.hypot(c[0], c[2])
        y = c[1]
        rs = math.hypot(r, y)
        if abs(rs - R0) < 2*h:
            r0 = math.floor(r/h + 1e-9)*h
            y0 = math.floor(y/h + 1e-9)*h
            num = den = 0.0
            for k in range(NS):
                rr = r0 + (k + 0.5)*h/NS
                yc = math.sqrt(R0*R0 - rr*rr) if rr < R0 else 0.0
                lo, hi = max(y0, -yc), min(y0 + h, yc)
                if hi > lo:
                    num += rr*(hi - lo)
                den += rr*h
            vap = num/den
        else:
            vap = 1.0 if rs < R0 else 0.0
        al.append(1 - vap)
        T.append(T_exact(rs, R0))
    for name, vals in (('alpha.water', al), ('T', T)):
        p = os.path.join(HERE, '0', name)
        s = open(p).read()
        body = 'nonuniform List<scalar> %d\n(\n%s\n)\n;' % (
            len(vals), '\n'.join('%.12g' % v for v in vals))
        s = re.sub(r'internalField\s+[^;]*;', 'internalField   ' + body, s,
                   count=1, flags=re.S)
        s = s.replace('373.98934', '%.8f' % Tbulk)
        open(p, 'w').write(s)
    print(f'init: {len(C)} cells, T range {min(T):.5f}..{max(T):.5f}')


if __name__ == '__main__':
    n = int(sys.argv[1])
    mode = sys.argv[2] if len(sys.argv) > 2 else 'mesh'
    (mesh if mode == 'mesh' else init)(n)
