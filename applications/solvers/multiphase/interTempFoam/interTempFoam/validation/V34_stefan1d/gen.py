#!/usr/bin/env python3
"""
V3 (Stefan, vapour-side flux) and V4 (sucking interface), 1D along x, one
row of cells (sides symmetryPlane, frontAndBack empty), sigma on.

Setups and properties from the Basilisk sandbox cases of E. Cipriano
(basilisk.fr/sandbox/ecipriano/run/stefanproblem.c, suckingproblem.c),
after Malan et al., JCP 426 (2021) 109920 / Welch & Wilson, JCP 160 (2000)
662 / Tanguy et al., JCP 264 (2014) 1:
  stefan : L0 = 10 mm, liquid (phase 1) at Tsat = 373 K for x < x_i,
           vapour layer 0.3225 mm at the hot wall x = L0 (T = 383 K);
           liquid leaves through x = 0 (p = 0, T = Tsat).
           delta(t) = 2 lam sqrt(a_v t), lam e^lam^2 erf lam = St/sqrt(pi)
  sucking: L0 = 10 mm, vapour layer 0.476 mm at the wall x = 0 (T = Tsat
           = 373.15 K), superheated liquid, T = 378.15 K and p = 0 at x = L0.
           delta(t) = 2 beta sqrt(a_v t), beta from Welch & Wilson (2000).
The analytical solution is time-shifted to the initial layer thickness;
the initial temperatures are the exact profiles at that time.

Usage: gen.py <stefan|sucking> <N>
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
HDR = open(os.path.join(HERE, 'system', 'fvSchemes')).read()
HDR = HDR[:HDR.index('FoamFile')]
L0 = 10e-3

P = {
    'stefan': dict(rhoL=958.0, rhoV=0.6, muL=2.82e-4, muV=1.23e-5,
                   kL=0.68, kV=0.025, cpL=4216.0, cpV=2080.0, L=2.256e6,
                   sigma=0.059, Tsat=373.0, Tw=383.0, delta0=10e-3 - 9.6775e-3,
                   tEnd=10.0),
    'sucking': dict(rhoL=958.4, rhoV=0.597, muL=2.80e-4, muV=1.26e-5,
                    kL=0.679, kV=0.025, cpL=4216.0, cpV=2030.0, L=2.26e6,
                    sigma=0.0059, Tsat=373.15, Tb=378.15, delta0=0.476e-3,
                    tEnd=0.2),
}


def bisect(f, lo, hi):
    for _ in range(200):
        mid = 0.5*(lo + hi)
        if f(lo)*f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5*(lo + hi)


def growth(case):
    p = P[case]
    aV = p['kV']/(p['rhoV']*p['cpV'])
    aL = p['kL']/(p['rhoL']*p['cpL'])
    if case == 'stefan':
        St = p['cpV']*(p['Tw'] - p['Tsat'])/p['L']
        g = lambda l: l*math.exp(l*l)*math.erf(l) - St/math.sqrt(math.pi)
        return bisect(g, 1e-6, 2.0), aV, aL
    rr = p['rhoV']/p['rhoL']
    def g(b):
        return math.exp(b*b)*math.erf(b)*(b - (p['Tb'] - p['Tsat'])*p['cpV']*p['kL']
            *math.sqrt(aV)*math.exp(-b*b*rr*rr*aV/aL)
            /(p['L']*p['kV']*math.sqrt(math.pi*aL)
              *math.erfc(b*rr*math.sqrt(aV/aL))))
    return bisect(g, 1e-3, 5.0), aV, aL


def delta(case, t):
    b, aV, _ = growth(case)
    return 2*b*math.sqrt(aV*t)


def tshift(case):
    b, aV, _ = growth(case)
    return (P[case]['delta0']/(2*b))**2/aV


def T_exact(case, x, t):
    p = P[case]
    b, aV, aL = growth(case)
    d = delta(case, t)
    if case == 'stefan':
        xi = L0 - x                               # distance from hot wall
        if xi >= d:
            return p['Tsat']
        return p['Tw'] + (p['Tsat'] - p['Tw'])/math.erf(b)*math.erf(
            xi/(2*math.sqrt(aV*t)))
    if x <= d:
        return p['Tsat']
    rr = p['rhoV']/p['rhoL']
    return p['Tb'] - (p['Tb'] - p['Tsat'])/math.erfc(b*rr*math.sqrt(aV/aL)) \
        * math.erfc(x/(2*math.sqrt(aL*t)) + b*(rr - 1)*math.sqrt(aV/aL))


def foam(cls, obj):
    return HDR + ('FoamFile\n{\n    version     2.0;\n    format      ascii;\n'
                  '    class       %s;\n    object      %s;\n}\n\n' % (cls, obj))


def field(cls, name, dims, internal, left, right):
    return foam(cls, name) + f"""dimensions      {dims};
internalField   {internal};
boundaryField
{{
    left            {{ {left} }}
    right           {{ {right} }}
    bottom          {{ type symmetryPlane; }}
    top             {{ type symmetryPlane; }}
    frontAndBack    {{ type empty; }}
}}
"""


def main():
    case, N = sys.argv[1], int(sys.argv[2])
    p = P[case]
    h = L0/N
    t0 = tshift(case)
    os.makedirs(os.path.join(HERE, '0'), exist_ok=True)
    # mesh
    open(os.path.join(HERE, 'system', 'blockMeshDict'), 'w').write(
        foam('dictionary', 'blockMeshDict') + f"""convertToMeters 1;
vertices
(
    (0 0 0) ({L0} 0 0) ({L0} {h} 0) (0 {h} 0)
    (0 0 {h}) ({L0} 0 {h}) ({L0} {h} {h}) (0 {h} {h})
);
blocks ( hex (0 1 2 3 4 5 6 7) ({N} 1 1) simpleGrading (1 1 1) );
boundary
(
    left  {{ type patch; faces ((0 4 7 3)); }}
    right {{ type patch; faces ((1 2 6 5)); }}
    bottom {{ type symmetryPlane; faces ((0 1 5 4)); }}
    top {{ type symmetryPlane; faces ((3 7 6 2)); }}
    frontAndBack {{ type empty; faces ((0 3 2 1) (4 5 6 7)); }}
);
""")
    # cell-centred initial fields
    xs = [(i + 0.5)*h for i in range(N)]
    if case == 'stefan':
        xi = L0 - p['delta0']
        alpha = [min(max((xi - (x - h/2))/h, 0.0), 1.0) for x in xs]
        bc = dict(
            U=('type pressureInletOutletVelocity; value uniform (0 0 0);',
               'type noSlip;'),
            p=('type fixedValue; value uniform 0;',
               'type fixedFluxPressure; value uniform 0;'),
            T=(f'type fixedValue; value uniform {p["Tsat"]};',
               f'type fixedValue; value uniform {p["Tw"]};'),
            a=('type inletOutlet; inletValue uniform 1; value uniform 1;',
               'type zeroGradient;'))
    else:
        xi = p['delta0']
        alpha = [min(max(((x + h/2) - xi)/h, 0.0), 1.0) for x in xs]
        bc = dict(
            U=('type noSlip;',
               'type pressureInletOutletVelocity; value uniform (0 0 0);'),
            p=('type fixedFluxPressure; value uniform 0;',
               'type fixedValue; value uniform 0;'),
            T=(f'type fixedValue; value uniform {p["Tsat"]};',
               f'type fixedValue; value uniform {p["Tb"]};'),
            a=('type zeroGradient;',
               'type inletOutlet; inletValue uniform 1; value uniform 1;'))
    # Exact cell averages of the initial temperature: the sucking problem
    # converts liquid enthalpy into vapour volume with the amplification
    # rho_l cp_l dT/(rho_v L) ~ 15, and at t0 the liquid thermal layer is
    # thinner than a cell; centre sampling put ~0.75 mm of spurious vapour
    # into the solution (V4 +21..+40 %).
    NS = 400
    T = [sum(T_exact(case, x - h/2 + (k + 0.5)*h/NS, t0) for k in range(NS))/NS
         for x in xs]
    lst = lambda v: 'nonuniform List<scalar> %d\n(\n%s\n)\n' % (
        len(v), '\n'.join('%.12g' % x for x in v))
    w = lambda n, s: open(os.path.join(HERE, '0', n), 'w').write(s)
    w('U', field('volVectorField', 'U', '[0 1 -1 0 0 0 0]',
                 'uniform (0 0 0)', *bc['U']))
    w('p_rgh', field('volScalarField', 'p_rgh', '[1 -1 -2 0 0 0 0]',
                     'uniform 0', *bc['p']))
    w('T', field('volScalarField', 'T', '[0 0 0 1 0 0 0]', lst(T), *bc['T']))
    w('alpha.water', field('volScalarField', 'alpha.water',
                           '[0 0 0 0 0 0 0]', lst(alpha), *bc['a']))
    # properties
    open(os.path.join(HERE, 'constant', 'transportProperties'), 'w').write(
        foam('dictionary', 'transportProperties') + f"""// {case}: see gen.py for sources
phases (water vapour);
water
{{
    transportModel  Newtonian;
    nu              {p['muL']/p['rhoL']:.8g};
    rho             {p['rhoL']};
    cp              {p['cpL']};
    Pr              {p['muL']*p['cpL']/p['kL']:.8g};
}}
vapour
{{
    transportModel  Newtonian;
    nu              {p['muV']/p['rhoV']:.8g};
    rho             {p['rhoV']};
    cp              {p['cpV']};
    Pr              {p['muV']*p['cpV']/p['kV']:.8g};
}}
sigma           {p['sigma']};
curvatureModel  heightFunction;
""")
    open(os.path.join(HERE, 'constant', 'phaseChangeProperties'), 'w').write(
        foam('dictionary', 'phaseChangeProperties') + f"""model           HardtWondra;
T_sat           {p['Tsat']};
h_lv            {p['L']};
k1              {p['kL']};
k2              {p['kV']};
""")
    tEnd = p['tEnd']
    open(os.path.join(HERE, 'system', 'controlDict'), 'w').write(
        foam('dictionary', 'controlDict') + f"""application     interTempFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {tEnd};
deltaT          1e-7;
writeControl    adjustable;
writeInterval   {tEnd/10};
writeFormat     ascii;
writePrecision  12;
timeFormat      general;
timePrecision   8;
runTimeModifiable no;
adjustTimeStep  yes;
maxCo           0.2;
maxAlphaCo      0.2;
maxCapillaryNum 0.5;
maxPhaseChangeCo 0.25;
maxDeltaT       1e-2;

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
}}
""")
    b, aV, aL = growth(case)
    with open(os.path.join(HERE, 'case.info'), 'w') as f:
        f.write(f'case {case}\nN {N}\nh {h}\ngrowth {b}\ntshift {t0}\n'
                f'tEnd {tEnd}\nA {h*h}\n')
    print(f'{case}: N={N} h={h:.4e} growth={b:.8f} tshift={t0:.5g} tEnd={tEnd}')


if __name__ == '__main__':
    main()
