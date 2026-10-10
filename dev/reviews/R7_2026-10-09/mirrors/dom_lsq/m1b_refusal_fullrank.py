"""Full-rank but ill-conditioned commons through the accumulator's QR route:
how often is a problem that QR-in-isolation solves to 1e-8 REFUSED (singularR via the
eig range-space fold) or returned with no correct digits?"""
import numpy as np
from m1_qr_roundtrip import accumulator_qr, qr_isolated
m = 800
for k in (7.5, 8, 8.5, 9):
    out = {}; errs = []; errq = []
    for seed in range(200):
        rs = np.random.RandomState(1000 + seed)
        U, _ = np.linalg.qr(rs.randn(m, 8)); V, _ = np.linalg.qr(rs.randn(8, 8))
        A = U @ np.diag(np.logspace(0, -k, 8)) @ V.T; B = np.zeros((m, 0))
        xt = rs.randn(8); l = A @ xt + 1e-10 * rs.randn(m)
        xr = np.linalg.lstsq(A, l, rcond=None)[0]
        xa, sa, pa = accumulator_qr(A, B, l); xq, sq = qr_isolated(A, B, l)
        key = pa.split('(')[0] + '/' + sa; out[key] = out.get(key, 0) + 1
        if xa is not None: errs.append(np.linalg.norm(xa - xr) / np.linalg.norm(xr))
        errq.append(np.linalg.norm(xq - xr) / np.linalg.norm(xr))
    print('cond(A)=1e%.1f: accumulator outcomes %s; median rel err accumulated %.1e, QR isolated %.1e'
          % (k, out, np.median(errs) if errs else np.nan, np.median(errq)))
