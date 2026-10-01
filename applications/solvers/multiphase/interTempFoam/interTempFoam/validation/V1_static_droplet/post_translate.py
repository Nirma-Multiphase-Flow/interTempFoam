#!/usr/bin/env python3
"""
Gate 2 pure-advection metrics: sphere R = 0.4 translated by (0, 1, 0) t on
the wedge mesh, frozen flow.

  vol_relErr   = integral(alpha) / integral(alpha_0) - 1
  shape_L1     = sum |alpha - alpha_exact| V / sum alpha_exact V, alpha_exact
                 from 8x8 sub-sampling of each (r, y) cell with weight r
                 (exact for the wedge volume element)
  bounds       = min(alpha), max(alpha) - 1
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from post import read_field, last_time, info, HERE  # noqa: E402

R = 0.4


def main():
    ci = info()
    h = float(ci['h'])
    t = last_time()
    tdir = os.path.join(HERE, t)
    C = read_field(os.path.join(tdir, 'C'), vector=True)
    V = read_field(os.path.join(tdir, 'V'))
    a = read_field(os.path.join(tdir, 'alpha.liquid'))
    a0 = read_field(os.path.join(HERE, '0', 'alpha.liquid'))
    yc = float(t)
    n = 8
    err = 0.0
    ex_tot = 0.0
    for c, v, x in zip(C, V, a):
        r0 = math.hypot(c[0], c[2]) - 0.5*h
        y0 = c[1] - 0.5*h
        num = den = 0.0
        for i in range(n):
            r = r0 + (i + 0.5)*h/n
            for j in range(n):
                y = y0 + (j + 0.5)*h/n
                den += r
                if r*r + (y - yc)**2 < R*R:
                    num += r
        ex = num/den if den > 0 else 0.0
        err += abs(x - ex)*v
        ex_tot += ex*v
    vol0 = sum(x*v for x, v in zip(a0, V))
    vol1 = sum(x*v for x, v in zip(a, V))
    lines = [
        f'Rh = {ci["Rh"]}', f't_end = {t}',
        f'vol_relErr = {vol1/vol0 - 1}',
        f'shape_L1 = {err/ex_tot}',
        f'min_alpha = {min(a)}', f'max_alpha_minus_1 = {max(a) - 1}',
    ]
    open(os.path.join(HERE, 'translate_metrics.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
