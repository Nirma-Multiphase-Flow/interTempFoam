#!/usr/bin/env python3
"""
Exact initial volume fractions for V1 (replaces setAlphaField, whose
approximate cut-cell integration limits the curvature accuracy at the level
of ~2%, cf. Scheufler & Roenby 2023, Sec. 4.2 "imperfect initialization").

planar: circle R centred at 0, area fraction of each h x h cell, integrated
        column-wise in x (n sub-columns) with the exact chord in y.
wedge : sphere R centred at 0, volume fraction of each (r, y) cell of the
        3D wedge, weight r (volume element r dr dtheta dy), integrated in r
        with the exact chord in y.
translate: as wedge.

Needs 0/C (postProcess -func writeCellCentres -time 0) and case.info.
"""
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from post import read_field, info, HERE  # noqa: E402

R = 0.4
NSUB = 400


def chord(x):
    return math.sqrt(R*R - x*x) if abs(x) < R else 0.0


def frac_planar(x0, y0, h):
    s = 0.0
    for k in range(NSUB):
        x = x0 + (k + 0.5)*h/NSUB
        yc = chord(x)
        lo, hi = max(y0, -yc), min(y0 + h, yc)
        if hi > lo:
            s += hi - lo
    return s/(NSUB*h)


def frac_wedge(r0, y0, h):
    num = den = 0.0
    for k in range(NSUB):
        r = r0 + (k + 0.5)*h/NSUB
        yc = chord(r)
        lo, hi = max(y0, -yc), min(y0 + h, yc)
        if hi > lo:
            num += r*(hi - lo)
        den += r*h
    return num/den


def main():
    ci = info()
    geom = ci['geom']
    h = float(ci['h'])
    C = read_field(os.path.join(HERE, '0', 'C'), vector=True)
    vals = []
    for c in C:
        if geom == 'planar':
            x0 = math.floor(c[0]/h + 1e-9)*h
            y0 = math.floor(c[1]/h + 1e-9)*h
            v = frac_planar(x0, y0, h) if abs(math.hypot(c[0], c[1]) - R) < 2*h \
                else (1.0 if math.hypot(c[0], c[1]) < R else 0.0)
        else:
            r = math.hypot(c[0], c[2])
            r0 = math.floor(r/h + 1e-9)*h
            y0 = math.floor(c[1]/h + 1e-9)*h
            v = frac_wedge(r0, y0, h) if abs(math.hypot(r, c[1]) - R) < 2*h \
                else (1.0 if math.hypot(r, c[1]) < R else 0.0)
        vals.append(v)
    p = os.path.join(HERE, '0', 'alpha.liquid')
    s = open(p).read()
    body = 'nonuniform List<scalar> %d\n(\n%s\n)\n;' % (
        len(vals), '\n'.join('%.15g' % v for v in vals))
    s = re.sub(r'internalField\s+[^;]*;', 'internalField   ' + body, s,
               count=1, flags=re.S)
    open(p, 'w').write(s)
    print('init_alpha: %s, %d cells, liquid volume fraction sum %.10g'
          % (geom, len(vals), sum(vals)))


if __name__ == '__main__':
    main()
