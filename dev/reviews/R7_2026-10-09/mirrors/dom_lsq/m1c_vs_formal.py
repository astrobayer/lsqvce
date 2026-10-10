"""Numerical error of the accumulator-QR route measured in units of the FORMAL error
(the brief's yardstick): err / sigma_x along the worst parameter."""
import numpy as np
from m1_qr_roundtrip import accumulator_qr, qr_isolated
m = 800
print('%6s %8s | %-28s | %-12s' % ('cond', 'noise', 'max|dx|/sigma_x accumulator', 'QR isolated'))
for k in (5, 6, 7, 8):
    for noise in (1e-4, 1e-7, 1e-10):
        ra = []; rq = []
        for seed in range(20):
            rs = np.random.RandomState(50 + seed)
            U, _ = np.linalg.qr(rs.randn(m, 8)); V, _ = np.linalg.qr(rs.randn(8, 8))
            A = U @ np.diag(np.logspace(0, -k, 8)) @ V.T
            l = A @ rs.randn(8) + noise * rs.randn(m)
            xr = np.linalg.lstsq(A, l, rcond=None)[0]
            sx = noise * np.sqrt(np.diag(np.linalg.pinv(A.T @ A)))
            xa, sa, _ = accumulator_qr(A, np.zeros((m, 0)), l); xq, _ = qr_isolated(A, np.zeros((m, 0)), l)
            if xa is not None: ra.append(np.max(np.abs(xa - xr) / sx))
            rq.append(np.max(np.abs(xq - xr) / sx))
        print('%6.0e %8.0e | median %9.2e  max %9.2e    | %9.2e' % (10.**k, noise, np.median(ra), np.max(ra), np.max(rq)))
