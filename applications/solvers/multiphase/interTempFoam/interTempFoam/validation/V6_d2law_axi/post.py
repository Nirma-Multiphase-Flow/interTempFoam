#!/usr/bin/env python3
"""
V6 metrics: evaporation rate and d2-law slope vs the quasi-steady law
(gen.py). Wedge (5 deg) slice of a half droplet:
  full-droplet volume V = 2 (2 pi/theta) integral(alpha.water),
  D = (6 V/pi)^(1/3); slope K = -d(D^2)/dt by least squares over t >= t_end/2
  mdot_full = 2 (2 pi/theta) mdotTotal (solver log)
Writes V6_metrics.txt, V6_history.csv.
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
    fac = 2*(2*math.pi/math.radians(5.0))
    rows = []
    for f in glob.glob(os.path.join(HERE, 'postProcessing', 'liquidVolume',
                                    '*', 'volFieldValue.dat')):
        for line in open(f):
            if line.startswith('#') or not line.strip():
                continue
            t, v = map(float, line.split()[:2])
            rows.append((t, (6*fac*v/math.pi)**(2/3)))
    rows.sort()
    md = []
    for lg in glob.glob(os.path.join(HERE, 'log.interTempFoam')):
        t = 0.0
        for line in open(lg):
            if line.startswith('Time = '):
                t = float(line.split()[2])
            m = re.search(r'mdotTotal = ([-0-9.eE+]+)', line)
            if m:
                md.append((t, fac*float(m.group(1))))
    tend = rows[-1][0]
    sel = [(t, d2) for t, d2 in rows if t >= 0.5*tend]
    n = len(sel)
    mt = sum(t for t, _ in sel)/n
    md2 = sum(d for _, d in sel)/n
    slope = sum((t - mt)*(d - md2) for t, d in sel)/sum((t - mt)**2 for t, _ in sel)
    Klo, Khi = float(info['K_Rinf_sqrt2L']), float(info['K_Rinf_L'])
    mlo, mhi = float(info['mdot_Rinf_sqrt2L']), float(info['mdot_Rinf_L'])
    mlate = [m for t, m in md if t >= 0.5*tend]
    mavg = sum(mlate)/len(mlate) if mlate else float('nan')
    with open(os.path.join(HERE, 'V6_history.csv'), 'w') as f:
        f.write('t,D2\n')
        for t, d2 in rows:
            f.write(f'{t},{d2}\n')
    Kmid = 0.5*(Klo + Khi)
    mmid = 0.5*(mlo + mhi)
    lines = [
        f'R0/h = {info["Rh"]}', f't_end = {tend}',
        f'K_num = {-slope}',
        f'K_law = {Klo} .. {Khi} (Rinf = sqrt2 L .. L)',
        f'K_relErr_vs_mid = {(-slope - Kmid)/Kmid}',
        f'mdot_num_mean(t>=tend/2) = {mavg}',
        f'mdot_law = {mlo} .. {mhi}',
        f'mdot_relErr_vs_mid = {(mavg - mmid)/mmid}',
        f'D2_relChange = {rows[-1][1]/rows[0][1] - 1}',
    ]
    open(os.path.join(HERE, 'V6_metrics.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
