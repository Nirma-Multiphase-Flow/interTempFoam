#!/usr/bin/env python3
"""
Generate the geometry-/resolution-dependent files of V1 (static droplet).

Usage: gen.py <planar|wedge> <R/h> [nSteps]

Droplet D = 0.8 (R = 0.4) centred at the origin (Popinet 2009, Sec. 5.1
parameters; the box size 4R x 4R is our choice, walls are 2R from the
interface).  planar: 2D Cartesian, 4R x 4R, cylinder droplet.  wedge:
axisymmetric 5 deg wedge about the y axis, r in [0, 2R], y in [-2R, 2R],
spherical droplet.

Time step: fixed, half the Brackbill-Kothe-Zemach capillary limit
dt_c = sqrt(rho_avg h^3 / (2 pi sigma)) (Brackbill et al., JCP 100, 1992),
so the baseline and the new code run the same step sequence.
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
HDR = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Version:  v2412
    \\\\  /    A nd           | Website:  www.openfoam.com
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      {0};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""

R, sigma, rho = 0.4, 1.0, 1.0


def write(rel, obj, body):
    with open(os.path.join(HERE, rel), 'w') as f:
        f.write(HDR.format(obj) + body)


def main():
    geom = sys.argv[1]
    n = int(sys.argv[2])
    nSteps = int(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] else 1000
    h = R/n
    L = 2*R
    dt = 0.5*math.sqrt(rho*h**3/(2*math.pi*sigma))

    if geom == 'planar':
        N = int(round(2*L/h))
        dz = h
        body = f"""convertToMeters 1;
vertices
(
    ({-L} {-L} 0) ({L} {-L} 0) ({L} {L} 0) ({-L} {L} 0)
    ({-L} {-L} {dz}) ({L} {-L} {dz}) ({L} {L} {dz}) ({-L} {L} {dz})
);
blocks ( hex (0 1 2 3 4 5 6 7) ({N} {N} 1) simpleGrading (1 1 1) );
boundary
(
    walls {{ type wall; faces ((0 4 7 3) (1 2 6 5) (0 1 5 4) (3 7 6 2)); }}
    frontAndBack {{ type empty; faces ((0 3 2 1) (4 5 6 7)); }}
);
"""
        alpha = f"""field           alpha.liquid;
type            cylinder;
radius          {R};
direction       (0 0 1);
origin          (0 0 {dz/2});
"""
    elif geom == 'wedge':
        Nr = int(round(L/h))
        Ny = int(round(2*L/h))
        th = math.radians(2.5)
        c, s = L*math.cos(th), L*math.sin(th)
        body = f"""convertToMeters 1;
mergeType points;
vertices
(
    (0 {-L} 0) ({c} {-L} {-s}) ({c} {L} {-s}) (0 {L} 0)
    ({c} {-L} {s}) ({c} {L} {s})
);
blocks ( hex (0 1 2 3 0 4 5 3) ({Nr} {Ny} 1) simpleGrading (1 1 1) );
boundary
(
    axis   {{ type empty; faces ((0 3 3 0)); }}
    walls  {{ type wall;  faces ((1 2 5 4) (3 2 5 3) (0 1 4 0)); }}
    wedge1 {{ type wedge; faces ((0 1 2 3)); }}
    wedge2 {{ type wedge; faces ((0 4 5 3)); }}
);
"""
        alpha = f"""field           alpha.liquid;
type            sphere;
radius          {R};
origin          (0 0 0);
"""
    elif geom == 'translate':
        # Pure-advection check (Gate 2): same droplet on the wedge mesh,
        # frozen uniform axial velocity U = (0 1 0), sigma = 0, moved 2.5 R
        # in t = 1.  dt from Co = 0.2.
        Nr = int(round(L/h))
        Ny = int(round((2*L + 2.5*R)/h))
        th = math.radians(2.5)
        c, s_ = L*math.cos(th), L*math.sin(th)
        y1 = L + 2.5*R
        body = f"""convertToMeters 1;
mergeType points;
vertices
(
    (0 {-L} 0) ({c} {-L} {-s_}) ({c} {y1} {-s_}) (0 {y1} 0)
    ({c} {-L} {s_}) ({c} {y1} {s_})
);
blocks ( hex (0 1 2 3 0 4 5 3) ({Nr} {Ny} 1) simpleGrading (1 1 1) );
boundary
(
    axis   {{ type empty; faces ((0 3 3 0)); }}
    walls  {{ type patch; faces ((1 2 5 4) (3 2 5 3) (0 1 4 0)); }}
    wedge1 {{ type wedge; faces ((0 1 2 3)); }}
    wedge2 {{ type wedge; faces ((0 4 5 3)); }}
);
"""
        alpha = f"""field           alpha.liquid;
type            sphere;
radius          {R};
origin          (0 0 0);
"""
        dt = 0.2*h
        nSteps = int(round(1.0/dt))
    else:
        raise SystemExit('geom must be planar, wedge or translate')

    write('system/blockMeshDict', 'blockMeshDict', body)
    write('system/setAlphaFieldDict', 'setAlphaFieldDict', alpha)
    write('system/controlDict', 'controlDict', f"""application     interTempFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {nSteps*dt:.10g};
deltaT          {dt:.10g};
writeControl    timeStep;
writeInterval   {nSteps};
purgeWrite      0;
writeFormat     ascii;
writePrecision  12;
timeFormat      general;
timePrecision   10;
runTimeModifiable no;
adjustTimeStep  no;
maxCo           1;
maxAlphaCo      1;
maxDeltaT       1;

functions
{{
    Umax
    {{
        type            fieldMinMax;
        libs            (fieldFunctionObjects);
        mode            magnitude;
        fields          (U);
        writeControl    timeStep;
        writeInterval   1;
    }}
    liquidVolume
    {{
        type            volFieldValue;
        libs            (fieldFunctionObjects);
        regionType      all;
        operation       volIntegrate;
        fields          (alpha.liquid);
        writeFields     false;
        writeControl    timeStep;
        writeInterval   1;
    }}
}}
""")
    with open(os.path.join(HERE, 'case.info'), 'w') as f:
        f.write(f'geom {geom}\nRh {n}\nh {h}\ndt {dt}\nnSteps {nSteps}\n')
    print(f'{geom} R/h={n} h={h} dt={dt:.4e} endTime={nSteps*dt:.4e}')


if __name__ == '__main__':
    main()
