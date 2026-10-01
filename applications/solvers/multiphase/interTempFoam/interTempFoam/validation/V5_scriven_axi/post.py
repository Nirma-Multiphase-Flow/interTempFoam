#!/usr/bin/env python3
"""
V5 metrics: bubble radius R(t) vs Scriven, R = 2 beta sqrt(alpha_l (t + t0)).

The case is a 5-deg wedge slice (3D) of a half bubble (symmetry plane y=0):
  V_bubble = 2 * (2 pi/theta) * (V_domain - integral(alpha.water))
  R_eff = (3 V_bubble/(4 pi))^(1/3)
Writes V5_metrics.txt and V5_radius.csv.
"""
import glob
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def read_field(path, vector=False):
    txt = open(path).read()
    m = re.search(r'internalField\s+(uniform|nonuniform)', txt)
    if m.group(1) == 'uniform':
        return None
    m2 = re.compile(r'List<\w+>\s*(\d+)\s*\(').search(txt, m.end())
    n = int(m2.group(1))
    body = txt[m2.end():]
    if vector:
        vals = re.findall(r'\(([^()]*)\)', body)[:n]
        return [tuple(float(x) for x in v.split()) for v in vals]
    toks = body[:body.index(')')].split()
    return [float(t) for t in toks[:n]]


def main():
    import gen
    info = dict(l.split() for l in open(os.path.join(HERE, 'case.info')))
    Vd = sum(read_field(os.path.join(HERE, 'constant', 'V')) or [])
    theta = math.radians(5.0)
    rows = []
    for f in glob.glob(os.path.join(HERE, 'postProcessing', 'liquidVolume',
                                    '*', 'volFieldValue.dat')):
        for line in open(f):
            if line.startswith('#') or not line.strip():
                continue
            t, v = map(float, line.split()[:2])
            Vb = 2*(2*math.pi/theta)*(Vd - v)
            rows.append((t, (3*Vb/(4*math.pi))**(1/3)))
    rows.sort()
    tend = float(info['tend'])
    out = open(os.path.join(HERE, 'V5_radius.csv'), 'w')
    out.write('t,R_num,R_exact,relErr\n')
    errs = []
    for t, R in rows:
        Rx = gen.R_exact(t)
        out.write(f'{t},{R},{Rx},{(R - Rx)/Rx}\n')
        if t >= 0.2*tend:
            errs.append(abs(R - Rx)/Rx)
    t, R = rows[-1]
    lines = [
        f'R0/h = {info["Rh"]}', f'beta = {info["beta"]}',
        f't_end = {t}', f'R_num = {R}', f'R_exact = {gen.R_exact(t)}',
        f'relErr_end = {(R - gen.R_exact(t))/gen.R_exact(t)}',
        f'maxRelErr_after_0.2tend = {max(errs) if errs else float("nan")}',
    ]
    open(os.path.join(HERE, 'V5_metrics.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
