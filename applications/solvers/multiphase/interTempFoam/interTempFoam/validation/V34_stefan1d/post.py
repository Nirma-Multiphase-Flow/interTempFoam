#!/usr/bin/env python3
"""V3/V4 metric: vapour-layer thickness delta = L0 - V_liquid/A vs exact."""
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gen  # noqa: E402

info = dict(l.split() for l in open(os.path.join(HERE, 'case.info')))
case, A, t0 = info['case'], float(info['A']), float(info['tshift'])
rows = []
for f in glob.glob(os.path.join(HERE, 'postProcessing', 'liquidVolume', '*',
                                'volFieldValue.dat')):
    for line in open(f):
        if line.startswith('#') or not line.strip():
            continue
        t, v = map(float, line.split()[:2])
        rows.append((t, gen.L0 - v/A))
rows.sort()
tEnd = rows[-1][0]
err = [abs(d/gen.delta(case, t + t0) - 1) for t, d in rows if t >= 0.2*tEnd]
t, d = rows[-1]
lines = [f'case = {case}', f'N = {info["N"]}', f't_end = {t}',
         f'delta_num = {d}', f'delta_exact = {gen.delta(case, t + t0)}',
         f'relErr_end = {d/gen.delta(case, t + t0) - 1}',
         f'maxRelErr_after_0.2tend = {max(err)}']
open(os.path.join(HERE, 'V34_metrics.txt'), 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
