#!/usr/bin/env python3
"""
V1 metrics (static droplet, no gravity, no phase change).

  Ca_s      = mu max|U| / sigma   (at the last step, and max over the run)
  dp error  = (p_in - p_out)/dp_exact - 1, dp_exact = sigma/R (planar),
              2 sigma/R (axisymmetric); p_in = volume average of p_rgh over
              r < 0.5 R, p_out over r > 1.5 R (Popinet 2009, Sec. 5.1).
  kappa err = mean and max of |kappa - kappa_exact|/kappa_exact over
              interface cells (0.01 < alpha < 0.99), if the solver wrote a
              curvature field (interfaceProperties:K or kappaI).
  volume    = relative change of integral(alpha dV).

Writes V1_metrics.txt (key = value lines) next to this script.
"""
import glob
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
R, sigma, mu = 0.4, 1.0, 8.16497e-3


def info():
    d = {}
    for line in open(os.path.join(HERE, 'case.info')):
        k, v = line.split()
        d[k] = v
    return d


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


def last_time():
    ts = []
    for d in os.listdir(HERE):
        try:
            t = float(d)
        except ValueError:
            continue
        if t > 0 and os.path.isdir(os.path.join(HERE, d)):
            ts.append((t, d))
    return max(ts)[1]


def umax_series():
    fs = glob.glob(os.path.join(HERE, 'postProcessing', 'Umax', '*',
                                'fieldMinMax.dat'))
    out = []
    for f in fs:
        hdr = None
        for line in open(f):
            if line.startswith('#'):
                hdr = line[1:].split()
                continue
            tok = line.split()
            if hdr and 'max' in hdr:
                # location columns are "(x y z)" -> 3 tokens each; rebuild
                vals = re.sub(r'\([^)]*\)', 'LOC', line).split()
                out.append((float(vals[0]), float(vals[hdr.index('max')])))
    out.sort()
    return out


def main():
    ci = info()
    geom = ci['geom']
    t = last_time()
    tdir = os.path.join(HERE, t)
    C = read_field(os.path.join(tdir, 'C'), vector=True)
    V = read_field(os.path.join(tdir, 'V'))
    p = read_field(os.path.join(tdir, 'p_rgh'))
    a = read_field(os.path.join(tdir, 'alpha.liquid'))
    a0 = read_field(os.path.join(HERE, '0', 'alpha.liquid'))
    Uf = read_field(os.path.join(tdir, 'U'), vector=True)

    def rad(c):
        return math.hypot(c[0], c[1])

    pin = sum(pi*v for pi, v, c in zip(p, V, C) if rad(c) < 0.5*R)
    vin = sum(v for v, c in zip(V, C) if rad(c) < 0.5*R)
    pout = sum(pi*v for pi, v, c in zip(p, V, C) if rad(c) > 1.5*R)
    vout = sum(v for v, c in zip(V, C) if rad(c) > 1.5*R)
    dp = pin/vin - pout/vout
    dpx = sigma/R if geom == 'planar' else 2*sigma/R
    kx = 1/R if geom == 'planar' else 2/R

    umax_end = max(math.sqrt(u[0]**2 + u[1]**2 + u[2]**2) for u in Uf)
    ser = umax_series()
    umax_run = max(u for _, u in ser) if ser else float('nan')

    vol0 = sum(x*v for x, v in zip(a0, V))
    vol1 = sum(x*v for x, v in zip(a, V))

    res = {
        'geom': geom, 'Rh': ci['Rh'], 'nSteps': ci['nSteps'], 't_end': t,
        'Ca_s_end': mu*umax_end/sigma,
        'Ca_s_maxRun': mu*umax_run/sigma,
        'dp_num': dp, 'dp_exact': dpx, 'dp_relErr': dp/dpx - 1,
        'vol_relChange': vol1/vol0 - 1,
    }
    for kname in ('kappaI', 'interfaceProperties:K'):
        kp = os.path.join(tdir, kname)
        if os.path.exists(kp):
            K = read_field(kp)
            if K is None:
                break
            e = [abs(k - kx)/kx for k, x in zip(K, a) if 0.01 < x < 0.99]
            if e:
                res['kappa_field'] = kname
                res['kappa_meanRelErr'] = sum(e)/len(e)
                res['kappa_maxRelErr'] = max(e)
                res['kappa_nCells'] = len(e)
            break
    lines = [f'{k} = {v}' for k, v in res.items()]
    open(os.path.join(HERE, 'V1_metrics.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    return 0


if __name__ == '__main__':
    sys.exit(main())
