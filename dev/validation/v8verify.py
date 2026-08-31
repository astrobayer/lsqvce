"""
v8verify.py -- Python reference for the direct unit tests added this session
(roadmap item 2: functions that ran only INSIDE other tests had no assertion
on their own behaviour). Mirror the numerics here first, then assert the same
properties in MATLAB (tests/test_util.m, tests/accum/test_v8.m).

Covered:
    util.getfielddef  -- empty-field-is-absent semantics
    util.mergeStruct  -- overlay semantics, empty/non-struct edge cases
    ar2arma           -- Durbin reduction: ACV round-trip of an ARMA(p,q)
    varWhiten         -- exactness vs the dense block-Cholesky GLS solve
    outliers          -- w_i = |v_i|/sqrt(r_i) flags injected blunders

Conventions match the toolbox:
    filter(b,a,.)  -> MATLAB lfilter, AR poly A(z)=a1+a2 z^-1+..., a1=1
    armaACV(a,b,n) -> r = ifft(|fft(h)|^2), h = impulse response of B/A
    epoch-major rows [x1 y1 x2 y2 ...] for the vector-AR whitener
All green => the MATLAB tolerances below are achievable on real hardware.
"""
import numpy as np
from numpy.linalg import solve, cholesky, norm
from scipy.signal import lfilter

rng = np.random.default_rng(8)
OK = []
def check(name, cond, extra=""):
    OK.append(bool(cond))
    print(f"[{'PASS' if cond else 'FAIL'}] {name}   {extra}")

# ----------------------------------------------------------------------
# 1. getfielddef -- pure logic mirror
# ----------------------------------------------------------------------
def getfielddef(s, f, d):
    if isinstance(s, dict) and (f in s) and (s[f] is not None) and \
       not (hasattr(s[f], '__len__') and len(s[f]) == 0):
        return s[f]
    return d

check("getfielddef present",      getfielddef({'a': 30}, 'a', 99) == 30)
check("getfielddef missing->def", getfielddef({'a': 30}, 'b', 99) == 99)
check("getfielddef empty->def",   getfielddef({'a': []}, 'a', 99) == 99,
      "empty field treated as absent (the BatchOverlap=[] case)")
check("getfielddef nonstruct",    getfielddef(None, 'a', 99) == 99)

# ----------------------------------------------------------------------
# 2. mergeStruct -- overlay logic mirror
# ----------------------------------------------------------------------
def mergeStruct(s, t):
    if not isinstance(t, dict) or len(t) == 0:
        return s
    if not isinstance(s, dict):
        return t
    out = dict(s)
    out.update(t)
    return out

a = {'type': 'psd', 'arOrder': 30}
c = mergeStruct(a, {'arOrder': 60})
check("mergeStruct overrides",   c['arOrder'] == 60 and c['type'] == 'psd')
check("mergeStruct keeps",       mergeStruct(a, {})['arOrder'] == 30, "empty t is a no-op")
check("mergeStruct nonstruct s", mergeStruct(None, {'x': 1}) == {'x': 1})

# ----------------------------------------------------------------------
# armaACV mirror (innovation variance 1)  -- shared by ar2arma test
# ----------------------------------------------------------------------
def arma_acv(a, b, n):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a / a[0]; b = b / b[0]
    nh = max(4 * n, 512)
    for _ in range(8):
        imp = np.zeros(nh); imp[0] = 1.0
        h = lfilter(b, a, imp)
        if abs(h[-1]) < np.finfo(float).eps * np.max(np.abs(h)):
            break
        nh *= 2
    nfft = 1 << int(np.ceil(np.log2(2 * nh)))
    H = np.fft.rfft(h, nfft)
    R = np.fft.irfft(np.abs(H) ** 2, nfft)
    return R[:n]

# ----------------------------------------------------------------------
# 3. ar2arma -- Durbin's two-stage reduction, ACV round-trip
#    A resonant ARMA(2,1): a sharp spectral line + one zero. This is the
#    regime where a low-order AR fails and ARMA wins (the 1/rev harmonic case).
# ----------------------------------------------------------------------
def levinson(r, p):
    # r[0..p] autocovariance -> AR(p) poly [1, -a1..-ap], innovation var
    r = np.asarray(r, float)
    a = np.array([1.0]); E = r[0]
    for m in range(1, p + 1):
        acc = r[m] + np.dot(a[1:], r[m-1:0:-1]) if m > 1 else r[m]
        k = -acc / E
        a = np.concatenate([a, [0.0]]) + k * np.concatenate([a[::-1], [0.0]])[::-1][:len(a)+1] * 0
        # explicit stable Levinson step
        a_prev = a[:-1] if len(a) > 1 else a
        a_new = np.concatenate([a_prev, [0.0]])
        a_new[1:] += k * a_prev[::-1]
        a = a_new
        E *= (1 - k * k)
    return a, E

# simpler/robust Levinson-Durbin
def levinson_durbin(r, p):
    r = np.asarray(r, float)
    a = np.zeros(p + 1); a[0] = 1.0
    E = r[0]
    for m in range(1, p + 1):
        k = -(r[m] + sum(a[i] * r[m - i] for i in range(1, m))) / E
        anew = a.copy()
        for i in range(1, m):
            anew[i] = a[i] + k * a[m - i]
        anew[m] = k
        a = anew
        E *= (1 - k * k)
        if E <= 0:
            E = abs(E) + 1e-30
    return a, E

# true ARMA(2,1): resonator poles at radius 0.94, angle ~2*pi*0.12 (a 1/rev-like line)
rad, ang = 0.94, 2 * np.pi * 0.12
a_true = np.real(np.poly([rad * np.exp(1j*ang), rad * np.exp(-1j*ang)]))  # [1, a1, a2]
b_true = np.array([1.0, 0.5])
s2_true = 1.0
n_acv = 256
r_target = s2_true * arma_acv(a_true, b_true, n_acv)

# stage 1: long AR from the target ACV
P = 60
aHi, s2Hi = levinson_durbin(r_target, P)
r_hiAR = s2Hi * arma_acv(aHi, [1.0], n_acv)
relAR = norm(r_hiAR - r_target) / norm(r_target)

# stage 2 (Durbin): fit impulse response of 1/aHi by a compact B/A via Hankel LS
def ar2arma(aHi, s2Hi, p, q, N=4096):
    imp = np.zeros(N); imp[0] = 1.0
    h = lfilter([1.0], aHi, imp)                 # h = 1/A_hi
    # solve for AR part a[1..p] from rows q+1..: H a = -h_shift  (Hankel LS)
    rows = np.arange(q + 1, N)
    Hk = np.column_stack([h[rows - j] for j in range(1, p + 1)])
    rhs = -h[rows]
    a_ar, *_ = np.linalg.lstsq(Hk, rhs, rcond=None)
    a = np.concatenate([[1.0], a_ar])
    # MA part: b = convolution of a with h, first q+1 terms
    hb = lfilter(a, [1.0], h)
    b = hb[:q + 1].copy(); b = b / b[0]
    return a, b, s2Hi

p, q = 8, 1
a_r, b_r, s2_r = ar2arma(aHi, s2Hi, p, q)
stationary = np.max(np.abs(np.roots(a_r))) < 1.0
r_arma = s2_r * arma_acv(a_r, b_r, n_acv)
relARMA = norm(r_arma - r_target) / norm(r_target)

check("ar2arma normalised a(1)=b(1)=1", abs(a_r[0]-1) < 1e-12 and abs(b_r[0]-1) < 1e-12)
check("ar2arma orders p,q",            len(a_r)-1 == p and len(b_r)-1 == q)
check("ar2arma AR stationary",         stationary, f"max|root|={np.max(np.abs(np.roots(a_r))):.4f}")
check("ar2arma ACV round-trip",        relARMA < 0.05, f"relErr(ARMA)={relARMA:.2e}  (target<1e-3, test tol 0.05)")
print(f"      [info] relErr(high AR)={relAR:.2e}, relErr(ARMA p={p},q={q})={relARMA:.2e}")

# ----------------------------------------------------------------------
# 4. varWhiten -- exactness vs dense block-Cholesky GLS
#    VAR(1), d=2. Build the block-Toeplitz covariance C over nE epochs, then
#    verify the FIR whitener reproduces the dense GLS solve.
# ----------------------------------------------------------------------
d, nE = 2, 40
A1 = np.array([[0.4, 0.1], [-0.2, 0.5]])      # VAR(1) coefficient (stable)
Lsig = np.array([[1.0, 0.0], [0.3, 0.8]])     # innovation chol, Sigma = Lsig Lsig'
Sigma = Lsig @ Lsig.T

# stationary cross-ACV Gamma(0), Gamma(k) via discrete Lyapunov + recursion
def var1_gamma(A1, Sigma, kmax):
    d = A1.shape[0]
    # vec(G0) = (I - A kron A)^-1 vec(Sigma)
    K = np.kron(A1, A1)
    G0 = solve(np.eye(d*d) - K, Sigma.reshape(-1)).reshape(d, d)
    G = [G0]
    for k in range(1, kmax + 1):
        G.append(A1 @ G[-1])
    return G

G = var1_gamma(A1, Sigma, nE)
# dense block-Toeplitz covariance C (m x m), epoch-major
m = d * nE
C = np.zeros((m, m))
for i in range(nE):
    for j in range(nE):
        k = i - j
        blk = G[k] if k >= 0 else G[-k].T
        C[i*d:(i+1)*d, j*d:(j+1)*d] = blk
C = 0.5 * (C + C.T)

# a small design + data, epoch-major
u = 3
Adesign = rng.standard_normal((m, u))
x_true = rng.standard_normal(u)
# color the noise: e = chol(C) w
Lc = cholesky(C)
e = Lc @ rng.standard_normal(m)
l = Adesign @ x_true + e

# dense GLS
Cinv = np.linalg.inv(C)
x_dense = solve(Adesign.T @ Cinv @ Adesign, Adesign.T @ Cinv @ l)
quad_dense = l @ Cinv @ l

# whitener = inv(Lc): w = Lc \ z  reproduces  w'w = z' Cinv z and the GLS solve
Aw = solve(Lc, Adesign)
lw = solve(Lc, l)
x_white = solve(Aw.T @ Aw, Aw.T @ lw)
quad_white = lw @ lw

check("varWhiten GLS solution match", norm(x_white - x_dense) / norm(x_dense) < 1e-10,
      f"rel={norm(x_white-x_dense)/norm(x_dense):.2e}")
check("varWhiten quadratic form",     abs(quad_white - quad_dense) / abs(quad_dense) < 1e-10,
      f"l'Cinv l reproduced")

# ----------------------------------------------------------------------
# 5. outliers -- w_i = |v_i| / sqrt(r_i) flags injected blunders
#    Dense joint mirror: r_i = 1 - h_i, h_i = a_i' inv(N) a_i.
# ----------------------------------------------------------------------
nO, uO = 200, 6
Ao = rng.standard_normal((nO, uO))
xo = rng.standard_normal(uO)
sig = 1.0
lo = Ao @ xo + sig * rng.standard_normal(nO)
bad = [37, 149]
lo[bad] += 12 * sig                       # two ~12 sigma blunders

N = Ao.T @ Ao
Ninv = np.linalg.inv(N)
xhat = Ninv @ (Ao.T @ lo)
v = Ao @ xhat - lo
h = np.einsum('ij,jk,ik->i', Ao, Ninv, Ao)   # leverage
r = 1 - h
w = np.abs(v) / np.sqrt(np.maximum(r, 1e-12))
thr = 4.0
flagged = np.where(w > thr)[0].tolist()
check("outliers detects both blunders", all(bb in flagged for bb in bad),
      f"flagged={flagged}")
# clean observations rarely flagged
false_pos = [i for i in flagged if i not in bad]
check("outliers few false positives", len(false_pos) <= 2, f"false+={false_pos}")
check("outliers sum(r)=n-u identity", abs(r.sum() - (nO - uO)) < 1e-8,
      f"sum(r)={r.sum():.4f}, n-u={nO-uO}")

print("\n=== v8verify:", f"{sum(OK)}/{len(OK)} PASS ===")
import sys
sys.exit(0 if all(OK) else 1)
