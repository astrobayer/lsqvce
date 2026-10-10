"""Does the block path refuse a free (rank-deficient) network, or solve it with an
arbitrary datum? Mirror of BlockAccumulator(qr) with nl = 0 (kind B):
  loc = QRSystem(nc); loc.add(A,l); blocks -> N = R1'R1; merge(chol(N) | eig); solve().
Also the NEQ backend: chol(N) of the accumulated normals, refusal only on chol failure.
Many random networks/seeds: count how often the refusal fires.
"""
import numpy as np
from m1_qr_roundtrip import QRSystem

def run(A, l):
    nc = A.shape[1]
    loc = QRSystem(nc); loc.add(A, l)
    ic = np.arange(nc); il = np.arange(nc, nc)
    _, _, Ncc, _, ncr, ltl, m = loc.blocks(ic, il)
    g = QRSystem(nc); path = g.merge(Ncc, ncr, ltl, m)
    x, st = g.solve()
    # NEQ backend
    N = A.T @ A
    try:
        np.linalg.cholesky(N); neq = 'solved'
    except np.linalg.LinAlgError:
        neq = 'refused'
    return path, st, x, neq

rs = np.random.RandomState(0)
stats = {}
examples = []
for trial in range(2000):
    u = rs.randint(4, 12)
    # random connected levelling network: a ring plus random chords
    pairs = [(i, (i + 1) % u) for i in range(u)]
    for _ in range(rs.randint(1, u)):
        i, j = rs.choice(u, 2, replace=False); pairs.append((i, j))
    A = np.zeros((len(pairs), u))
    for k, (i, j) in enumerate(pairs):
        A[k, i] = -1; A[k, j] = 1
    A *= rs.uniform(0.5, 3.0)            # scale
    h = rs.randn(u) * 10
    l = A @ h + 0.01 * rs.randn(len(pairs))
    path, st, x, neq = run(A, l)
    key = (path.split('(')[0], st, neq)
    stats[key] = stats.get(key, 0) + 1
    if st == 'ok' and len(examples) < 3:
        examples.append((u, x.mean(), np.abs(x).max()))
for k, v in sorted(stats.items(), key=lambda kv: -kv[1]):
    print('merge path %-5s  qr solve %-10s  neq %-8s : %4d / 2000' % (k + (v,)))
print('examples of a FREE network "solved" by the qr path (u, mean(x), max|x|):')
for e in examples: print('   u=%2d  mean(x)=%+.3e  max|x|=%.3e' % e)
