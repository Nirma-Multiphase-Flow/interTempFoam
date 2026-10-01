#!/usr/bin/env python3
"""
V6: water droplet at T_sat evaporating in superheated vapour (d2 law),
axisymmetric (3D 5-deg wedge about the y axis, symmetry plane y = 0).

Quasi-steady solution with Stefan flow (Godsave 1953; Spalding 1953; e.g.
Law, Prog. Energy Combust. Sci. 8 (1982) 171), vapour held at T_inf on a
sphere of radius Rinf:
    B      = cp_v (T_inf - T_sat)/L
    mdot   = 4 pi (k_v/cp_v) ln(1 + B)/(1/R - 1/Rinf)          [kg/s]
    d(D^2)/dt = -8 k_v ln(1 + B)/(rho_l cp_v)/(1 - R/Rinf)      (d2 law)
    T(r)   = T_sat + (L/cp_v)(exp(zeta (1/R - 1/r)) - 1),
             zeta = mdot cp_v/(4 pi k_v)                        (r >= R)
The box has no spherical outer boundary; the correction uses the radius
of the sphere of equal distance to the nearest outer face (Rinf = Lbox)
and, as a bound, Rinf = sqrt(2) Lbox (corners). The difference is logged.

Properties: Tanguy et al. (2014) water/vapour (constant/), T_sat 373 K,
T_inf 393 K (superheat 20 K, as run/bubble_evap2_*), R0 = 50 um.

Usage: gen.py <R0/h> [mode]   mode 'mesh' | 'init'
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

rhoL, rhoV = 958.0, 0.59
cpV = 2034.0
kV = 0.026
L = 2.257e6
Tsat, Tinf = 373.0, 393.0
R0 = 50e-6
LBOX = 16*R0
B = cpV*(Tinf - Tsat)/L
SIGMA = 0.05891


def mdot(R, Rinf):
    return 4*math.pi*(kV/cpV)*math.log(1 + B)/(1/R - 1/Rinf)


def K(R, Rinf):
    return 8*kV*math.log(1 + B)/(rhoL*cpV)/(1 - R/Rinf)


def T_qs(r, R, Rinf):
    if r <= R:
        return Tsat
    zeta = mdot(R, Rinf)*cpV/(4*math.pi*kV)
    return Tsat + (L/cpV)*(math.exp(zeta*(1/R - 1/r)) - 1)


def foam(o):
    hdr = open(os.path.join(HERE, 'system', 'fvSchemes')).read()
    hdr = hdr[:hdr.index('FoamFile')]
    return hdr + ('FoamFile\n{\n    version     2.0;\n    format      '
                  'ascii;\n    class       dictionary;\n    object      '
                  '%s;\n}\n\n' % o)


def mesh(n, tEnd):
    h = R0/n
    # uniform core [0, 1.5 R0] with h, then geometric stretching to LBOX
    nCore = int(round(1.5*R0/h))
    nOut = int(round(1.2*n))
    # expansion ratio last/first cell of the outer block, first cell ~ h
    Lout = LBOX - 1.5*R0
    # solve for ratio q per cell: h (q^nOut - 1)/(q - 1) = Lout
    lo, hi = 1.0000001, 2.0
    for _ in range(200):
        q = 0.5*(lo + hi)
        if h*(q**nOut - 1)/(q - 1) > Lout:
            hi = q
        else:
            lo = q
    grad = q**(nOut - 1)
    th = math.radians(2.5)
    c, s = LBOX*math.cos(th), LBOX*math.sin(th)
    fc = 1.5*R0/LBOX
    g = f'(({fc:.8f} {nCore} 1) ({1 - fc:.8f} {nOut} {grad:.6f}))'
    body = f"""convertToMeters 1;
mergeType points;
vertices
(
    (0 0 0) ({c} 0 {-s}) ({c} {LBOX} {-s}) (0 {LBOX} 0)
    ({c} 0 {s}) ({c} {LBOX} {s})
);
blocks
(
    hex (0 1 2 3 0 4 5 3) ({nCore + nOut} {nCore + nOut} 1)
    simpleGrading ({g} {g} 1)
);
boundary
(
    axis   {{ type empty; faces ((0 3 3 0)); }}
    bottom {{ type symmetryPlane; faces ((0 1 4 0)); }}
    outlet {{ type patch; faces ((1 2 5 4) (3 2 5 3)); }}
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
endTime         {tEnd:.8g};
deltaT          1e-9;
writeControl    adjustable;
writeInterval   {tEnd/5:.8g};
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
        f.write(f'Rh {n}\nh {h}\ntend {tEnd}\nB {B}\n'
                f'K_Rinf_L {K(R0, LBOX)}\nK_Rinf_sqrt2L {K(R0, math.sqrt(2)*LBOX)}\n'
                f'mdot_Rinf_L {mdot(R0, LBOX)}\n'
                f'mdot_Rinf_sqrt2L {mdot(R0, math.sqrt(2)*LBOX)}\n')
    print(f'V6 R0/h={n} h={h:.3e} cells/dir={nCore + nOut} grading={grad:.3f} '
          f'B={B:.5f} K={K(R0, LBOX):.4e}..{K(R0, math.sqrt(2)*LBOX):.4e} '
          f'tend={tEnd:g}')


def init(n):
    from post import read_field
    h = R0/n
    C = read_field(os.path.join(HERE, '0', 'C'), vector=True)
    NS = 200
    al, T = [], []
    Rinf = LBOX
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
            liq = num/den
        else:
            liq = 1.0 if rs < R0 else 0.0
        al.append(liq)
        T.append(T_qs(rs, R0, Rinf))
    for name, vals in (('alpha.water', al), ('T', T)):
        p = os.path.join(HERE, '0', name)
        s = open(p).read()
        body = 'nonuniform List<scalar> %d\n(\n%s\n)\n;' % (
            len(vals), '\n'.join('%.12g' % v for v in vals))
        s = re.sub(r'internalField\s+[^;]*;', 'internalField   ' + body, s,
                   count=1, flags=re.S)
        open(p, 'w').write(s)
    print(f'init: {len(C)} cells, T range {min(T):.4f}..{max(T):.4f}')


if __name__ == '__main__':
    n = int(sys.argv[1])
    mode = sys.argv[2] if len(sys.argv) > 2 else 'mesh'
    if mode == 'mesh':
        mesh(n, float(sys.argv[3]) if len(sys.argv) > 3 else 2e-3)
    else:
        init(n)
