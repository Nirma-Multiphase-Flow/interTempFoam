#!/usr/bin/env python3
"""
V7 metric: normalised bubble radius R/(2 beta sqrt(alpha_l)) vs time since
nucleation t_abs = t + t0 (t0 = Scriven time of the initial bubble),
compared with Florschuetz et al. (1969) (EBIT Table 1). Full bubble volume
from the 5-deg wedge slice: V = (2 pi/theta)(V_domain - integral alpha).
Writes V7_metrics.txt, V7_radius.csv.
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
    t0 = float(info['t0'])
    norm = 2*gen.BETA*math.sqrt(gen.alphaL)
    Vd = sum(read_field(os.path.join(HERE, 'constant', 'V')))
    fac = 2*math.pi/math.radians(5.0)
    rows = []
    for f in glob.glob(os.path.join(HERE, 'postProcessing', 'liquidVolume',
                                    '*', 'volFieldValue.dat')):
        for line in open(f):
            if line.startswith('#') or not line.strip():
                continue
            t, v = map(float, line.split()[:2])
            R = (3*fac*(Vd - v)/(4*math.pi))**(1/3)
            rows.append((t + t0, R/norm))
    rows.sort()
    with open(os.path.join(HERE, 'V7_radius.csv'), 'w') as f:
        f.write('t_abs,R_norm\n')
        for t, r in rows:
            f.write(f'{t},{r}\n')
    lines = [f'h = {info["h"]}', f'beta = {gen.BETA}', f't0 = {t0}',
             't_abs  R_norm_num  R_norm_exp  relDiff']
    errs = []
    for te, re_ in gen.EXPERIMENT:
        if rows[-1][0] < te:
            break
        rn = min(rows, key=lambda r: abs(r[0] - te))[1]
        errs.append(abs(rn/re_ - 1))
        lines.append(f'{te}  {rn:.4f}  {re_:.3f}  {rn/re_ - 1:+.4f}')
    if errs:
        lines.append(f'mean |relDiff| = {sum(errs)/len(errs):.4f}')
    open(os.path.join(HERE, 'V7_metrics.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
