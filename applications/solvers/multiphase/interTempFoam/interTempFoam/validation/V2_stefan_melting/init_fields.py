#!/usr/bin/env python3
"""
Row-uniform initial fields for V2 (replaces setFields + setExprFields).

The front x0 = 2.8e-3 m lies exactly on the cell centres of column 2
(h = 0.28/250 = 1.12e-3, x0 = 2.5 h), so the box test x < x0 of setFields /
setExprFields was decided by floating-point round-off per row and the
initial front was a random staircase across the 125 rows. Here:
  alpha.solid = exact area fraction of x > x0 in each cell (0.5 in column 2)
  T = 350 - 37 x/x0 for x < x0 (setExprFieldsDict formula), 313 otherwise
Needs 0/C.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
x0, h = 2.8e-3, 0.28/250


def read_C(path):
    txt = open(path).read()
    m = re.compile(r'List<vector>\s*(\d+)\s*\(').search(txt)
    n = int(m.group(1))
    vals = re.findall(r'\(([^()]*)\)', txt[m.end():])[:n]
    return [tuple(float(v) for v in s.split()) for s in vals]


def write(name, vals):
    p = os.path.join(HERE, '0', name)
    s = open(p).read()
    body = 'nonuniform List<scalar> %d\n(\n%s\n)\n;' % (
        len(vals), '\n'.join('%.15g' % v for v in vals))
    s = re.sub(r'internalField\s+[^;]*;', 'internalField   ' + body, s,
               count=1, flags=re.S)
    open(p, 'w').write(s)


C = read_C(os.path.join(HERE, '0', 'C'))
alpha, T = [], []
for c in C:
    lo, hi = c[0] - h/2, c[0] + h/2
    alpha.append(min(max((hi - x0)/h, 0.0), 1.0))
    T.append(350 - 37*c[0]/x0 if c[0] < x0 - 1e-12 else 313.0)
write('alpha.solid', alpha)
write('T', T)
print('init_fields: %d cells, partial cells %d' %
      (len(C), sum(1 for a in alpha if 0 < a < 1)))
