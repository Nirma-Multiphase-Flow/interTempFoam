#!/usr/bin/env python3
"""
V2 post-processing: 1D one-phase Stefan (Neumann) melting.

Analytical solution (Carslaw & Jaeger; Alexiades & Solomon, "Mathematical
Modeling of Melting and Freezing Processes", 1993, Sec. 2.1):
    St = cp (Tw - Tm) / L
    lambda exp(lambda^2) erf(lambda) = St / sqrt(pi)
    x_s(t) = 2 lambda sqrt(a t),     a = k / (rho cp)
    T(x,t) = Tw - (Tw - Tm) erf(x / (2 sqrt(a t))) / erf(lambda),  x < x_s
Properties: Shaikh et al. (2016) paraffin benchmark, as in constant/.

The case starts with the front at x0 = 2.8e-3 m (setFieldsDict) and the
analytical linearised profile; the similarity solution reaches x0 at
t0 = x0^2 / (4 lambda^2 a).  Simulation time t is compared with the
analytical solution at t + t0.

Writes: V2_metrics.txt, V2_interface.csv, V2_interface.png (if matplotlib).
"""
import glob
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

Tw, Tm = 350.0, 313.0
k, rho, cp, L = 0.21, 750.0, 2400.0, 175000.0
x0 = 2.8e-3
a = k/(rho*cp)
St = cp*(Tw - Tm)/L


def root():
    f = lambda l: l*math.exp(l*l)*math.erf(l) - St/math.sqrt(math.pi)
    lo, hi = 1e-6, 3.0
    for _ in range(200):
        mid = 0.5*(lo + hi)
        if f(lo)*f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5*(lo + hi)


lam = root()
t0 = x0*x0/(4*lam*lam*a)


def xs(t):
    return 2*lam*math.sqrt(a*(t + t0))


def Tan(x, t):
    return Tw - (Tw - Tm)*math.erf(x/(2*math.sqrt(a*(t + t0))))/math.erf(lam)


def read_set(path):
    rows = []
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        rows.append([float(v) for v in line.split()])
    return rows


def main():
    files = glob.glob(os.path.join(HERE, 'postProcessing', '*', '*',
                                   'xLine_*.xy'))
    if not files:
        print('no sampling output found')
        return 1
    header = None
    data = []
    for f in files:
        t = float(os.path.basename(os.path.dirname(f)))
        if t <= 0:
            continue
        name = os.path.basename(f)[len('xLine_'):-3]
        cols = name.split('_')
        # field names may contain '.', never '_' except as separator
        rows = read_set(f)
        ia = cols.index('alpha.solid') + 1
        iT = cols.index('T') + 1
        # interface: first alpha.solid = 0.5 crossing from the hot wall
        xi = None
        for r0, r1 in zip(rows[:-1], rows[1:]):
            a0, a1 = r0[ia], r1[ia]
            if (a0 - 0.5)*(a1 - 0.5) <= 0 and a0 != a1:
                xi = r0[0] + (0.5 - a0)*(r1[0] - r0[0])/(a1 - a0)
                break
        # liquid T error, normalised by (Tw - Tm), for x < 0.9 x_s(analytic)
        xa = xs(t)
        errT = 0.0
        for r in rows:
            if 0 < r[0] < 0.9*xa:
                errT = max(errT, abs(r[iT] - Tan(r[0], t))/(Tw - Tm))
        data.append((t, xi, xa, errT))
    data.sort()
    with open(os.path.join(HERE, 'V2_interface.csv'), 'w') as fh:
        fh.write('t_s,x_num_m,x_an_m,relErr_x,maxErrT_norm\n')
        for t, xi, xa, eT in data:
            e = (xi - xa)/xa if xi is not None else float('nan')
            fh.write(f'{t},{xi},{xa},{e},{eT}\n')
    tend, xi, xa, eT = data[-1]
    errs = [abs((d[1] - d[2])/d[2]) for d in data
            if d[1] is not None and d[0] >= 0.1*tend]
    Terrs = [d[3] for d in data if d[0] >= 0.1*tend]
    lines = [
        f'lambda = {lam:.6f}  St = {St:.5f}  t0 = {t0:.2f} s',
        f't_end = {tend}  x_num = {xi}  x_an = {xa:.6e}',
        f'relErr_x(t_end) = {(xi - xa)/xa:.4%}',
        f'max |relErr_x| over t >= 0.1 t_end = {max(errs):.4%}',
        f'max normalised T error (liquid, t >= 0.1 t_end) = {max(Terrs):.4%}',
        f'n_samples = {len(data)}',
    ]
    open(os.path.join(HERE, 'V2_metrics.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        ts = [d[0]/3600 for d in data]
        plt.plot(ts, [d[1]*1e3 for d in data], 'o', ms=2, label='interTempFoam')
        plt.plot(ts, [d[2]*1e3 for d in data], '-', label='Neumann')
        plt.xlabel('t [h]')
        plt.ylabel('x_s [mm]')
        plt.legend()
        plt.savefig(os.path.join(HERE, 'V2_interface.png'), dpi=120)
    except Exception as exc:
        print('plot skipped:', exc)
    return 0


if __name__ == '__main__':
    sys.exit(main())
