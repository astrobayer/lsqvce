"""
v9verify.py -- Python reference for the functions moved into +lsqvce this
session (selectArAcv, residualAcv, collectResiduals, whitenessCorr,
normalizeSatList, resolveSatSpecs, private satKeyOf). Mirror the numerics here,
then assert the same properties in MATLAB (tests/accum/test_v9.m).

Covered numerically:
    selectArAcv         -- reflection-ladder sweep recovers a known AR order
    whitenessCorr          -- white below band, correlated above
    residualAcv -- D^-1 v of coloured-by-D white noise -> white ACV
All green => the MATLAB tolerances in test_v9.m are achievable.
"""
import numpy as np
from numpy.linalg import norm
from scipy.signal import lfilter

rng = np.random.default_rng(9)
OK = []
def check(name, cond, extra=""):
    OK.append(bool(cond)); print(f"[{'PASS' if cond else 'FAIL'}] {name}   {extra}")

# ---- toolbox-consistent primitives -----------------------------------------
def levinson(r, p):
    r = np.asarray(r, float); a = np.zeros(p+1); a[0] = 1; E = r[0]; refl = np.zeros(p)
    for m in range(1, p+1):
        k = -(r[m] + sum(a[i]*r[m-i] for i in range(1, m))) / E
        an = a.copy()
        for i in range(1, m):
            an[i] = a[i] + k*a[m-i]
        an[m] = k; a = an; refl[m-1] = k; E *= (1 - k*k)
    return a, E, refl

def armaACV(a, b, n):
    a = np.asarray(a, float); b = np.asarray(b, float); a = a/a[0]; b = b/b[0]
    nh = max(4*n, 512); imp = np.zeros(nh); imp[0] = 1; h = lfilter(b, a, imp)
    nfft = 1 << int(np.ceil(np.log2(2*nh))); H = np.fft.rfft(h, nfft)
    R = np.fft.irfft(np.abs(H)**2, nfft); return R[:n]

# ---- 1. selectArAcv core: order recovery -----------------------------------
def selectArAcv(acv, nData, pMax=120, tolA=1e-3, plat=0.02):
    acv = np.asarray(acv, float); nLag = len(acv); r0 = acv[0]
    pMax = min(nLag-1, pMax)
    if r0 <= 0 or pMax < 1:
        return np.array([1.0]), max(r0, 1e-300), 0
    _, _, refl = levinson(acv, pMax)
    pMax = len(refl)
    normAcv = max(norm(acv), 1e-300)
    err = np.full(pMax+1, np.inf)
    err[0] = norm(np.concatenate([[r0], np.zeros(nLag-1)]) - acv) / normAcv  # p=0
    prev = err[0]; knee = 0; pStar = -1
    if err[0] <= tolA:
        pStar = 0
    for pk in range(1, pMax+1):
        if pStar >= 0:
            break
        ak, s2k, _ = levinson(acv, pk)
        rFit = s2k * armaACV(ak, [1.0], nLag)
        err[pk] = norm(rFit - acv) / normAcv
        if knee == 0 and (prev - err[pk]) <= plat*prev:
            knee = pk-1
        prev = err[pk]
        if err[pk] <= tolA:
            pStar = pk
    if pStar < 0:
        pStar = knee if knee >= 1 else int(np.argmin(err))
    if pStar == 0:
        return np.array([1.0]), r0, 0
    a, s2, _ = levinson(acv, pStar)
    return a, s2, pStar

for (a0, s0, want) in [([1, -1.6, 0.9, -0.25], 0.7, 3),
                       ([1, -0.6], 1.0, 1),
                       ([1, -1.2, 0.72, -0.11, 0.05], 0.5, 4)]:
    r = s0**2 * armaACV(a0, [1.0], 200)
    a, s2, p = selectArAcv(r, 5000)
    check(f"selectArAcv recovers AR({want})", p == want, f"got p={p}")
    if p == want:
        check(f"  AR({want}) coeffs+sigma", norm(a - np.array(a0)) < 1e-6 and abs(np.sqrt(s2)-s0) < 1e-6)

# white input -> p = 0
rw = np.zeros(50); rw[0] = 1.0
_, _, pw = selectArAcv(rw, 5000)
check("selectArAcv white -> p=0", pw == 0, f"got p={pw}")

# ---- 2. whitenessCorr ----------------------------------------------------------
def whitenessCorr(x, maxLag):
    x = np.asarray(x, float) - np.mean(x); n = len(x); lag0 = x @ x
    if n < 2 or lag0 == 0:
        return 0.0
    return max(abs(x[:n-l] @ x[l:]) / lag0 for l in range(1, min(maxLag, n-1)+1))

w = rng.standard_normal(4000)
check("whitenessCorr white below band", whitenessCorr(w, 5) < 3/np.sqrt(4000),
      f"{whitenessCorr(w,5):.4f} < {3/np.sqrt(4000):.4f}")
ar1 = lfilter([1.0], [1, -0.8], rng.standard_normal(4000))   # strong AR(1)
check("whitenessCorr correlated above band", whitenessCorr(ar1, 5) > 3/np.sqrt(4000),
      f"{whitenessCorr(ar1,5):.4f}")
check("whitenessCorr in [0,1]", 0 <= whitenessCorr(ar1, 5) <= 1)

# ---- 3. residualAcv: D^-1 (D w) -> white ACV -----------------------
# pos-like: 3 axes, per-epoch 3x3 D. Build residuals v = chol(D) w, then the
# de-scaled series chol(D)^-1 v = w must have a white ACV.
d, nE, maxLag = 3, 600, 40
# a fixed but non-trivial per-epoch covariance (same each epoch here)
L = np.array([[1.0, 0, 0], [0.3, 0.8, 0], [-0.2, 0.1, 0.6]])
Dcov = L @ L.T
acc = np.zeros(maxLag+1); nSamp = 0
for _ in range(4):                       # a few "arcs"
    w = rng.standard_normal((d, nE))
    v = L @ w                            # colour by D (instantaneous only)
    desc = np.linalg.solve(L, v)         # chol(D)^-1 v  == w
    for ax in range(d):
        s = desc[ax] - desc[ax].mean(); m = len(s)
        for lag in range(0, min(maxLag, m-1)+1):
            acc[lag] += s[:m-lag] @ s[lag:]
        nSamp += m
acv = acc / nSamp
r0, r1 = acv[0], abs(acv[1])/acv[0]
check("descaled ACV ~ white (r0~1)", abs(r0 - 1.0) < 0.1, f"r0={r0:.3f}")
check("descaled ACV ~ white (|r1/r0| small)", r1 < 0.1, f"|r1/r0|={r1:.3f}")
check("descaled nSamp = arcs*axes*nE", nSamp == 4*d*nE, f"nSamp={nSamp}")

print("\n=== v9verify:", f"{sum(OK)}/{len(OK)} PASS ===")
import sys; sys.exit(0 if all(OK) else 1)
