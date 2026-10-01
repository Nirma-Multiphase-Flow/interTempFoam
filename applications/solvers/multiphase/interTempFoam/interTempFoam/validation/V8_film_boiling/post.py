#!/usr/bin/env python3
"""
V8 metrics: V/V0 (vapour volume) and wall Nusselt number vs time,
Nu = 2/lambda int l0/(Tw - Tsat) (Tw - T_P)/(h/2) dx over the bottom cells
(arXiv:2407.09969 eq. 47), time average over t >= 0.1 s vs Berenson.
Writes V8_metrics.txt, V8_history.csv.
"""
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gen  # noqa: E402

info = dict(l.split() for l in open(os.path.join(HERE, 'case.info')))
nx = int(info['nx'])
h = gen.LAM/2/nx
area = (gen.LAM/2)*h          # domain x-width times dz = h


def readT(t):
    s = open(os.path.join(HERE, t, 'T')).read()
    m = re.compile(r'List<scalar>\s*(\d+)\s*\(').search(s)
    return list(map(float, s[m.end():].split(')')[0].split()[:nx]))


vol = []
for f in glob.glob(os.path.join(HERE, 'postProcessing', 'liquidVolume', '*',
                                'volFieldValue.dat')):
    for line in open(f):
        if line.startswith('#') or not line.strip():
            continue
        t, v = map(float, line.split()[:2])
        vol.append((t, (gen.LAM/2)*gen.LAM*h - v))
vol.sort()
V0 = vol[0][1]
Tw = gen.TSAT + gen.DT
nus = []
for d in os.listdir(HERE):
    try:
        t = float(d)
    except ValueError:
        continue
    if t <= 0 or not os.path.exists(os.path.join(HERE, d, 'T')):
        continue
    T = readT(d)
    nu = sum(gen.L0/gen.DT*(Tw - Tp)/(h/2) for Tp in T)/nx
    nus.append((t, nu))
nus.sort()
with open(os.path.join(HERE, 'V8_history.csv'), 'w') as f:
    f.write('t,Nu\n')
    for t, n in nus:
        f.write(f'{t},{n}\n')
    f.write('t,V_over_V0\n')
    for t, v in vol[::50]:
        f.write(f'{t},{v/V0}\n')
late = [n for t, n in nus if t >= 0.1]
lines = [f'nx = {nx}', f't_end = {vol[-1][0]}', f'V/V0(t_end) = {vol[-1][1]/V0}',
         f'Nu_mean(t>=0.1) = {sum(late)/len(late) if late else float("nan")}',
         f'Nu_Berenson = {gen.NU_BERENSON}']
lines += [f'Nu({t}) = {n:.4f}' for t, n in nus]
open(os.path.join(HERE, 'V8_metrics.txt'), 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
