#!/usr/bin/env python3
"""WP-126 mirror: lsqvce.kernel.ar2arma (Durbin's method) -- the orientation of
the Hankel regression, and the cases the MATLAB tests rest on.

Durbin: fit a long AR to the spectrum, take its impulse response h = 1/A_hi(z),
and fit B(z)/A(z) to it. For n > q the ARMA impulse response obeys

    h(n) + a_1 h(n-1) + ... + a_p h(n-p) = 0,

so with the Hankel rows H(i,:) = h(q+i .. q+i+p) the LAST column is the target
and a_k multiplies column p+1-k. The code before WP-126 regressed the FIRST
column on the later ones -- the recursion run backwards -- which returns the
reciprocal roots (a stationary 0.94 becomes 1.064), so the stationarity guard
fired on almost every input and the reduction fell back to the long AR.

Cases (printed and asserted):
  1. the round-trip of test_signalCore>testAr2armaRoundTrip -- a resonant
     ARMA(2,1) (poles 0.94 at 0.12 cycles, b = [1 0.5]) through AR(60),
     armaP 2, q 1: the old orientation gives |root| 1.064, the fixed one the
     true a and b to ~1e-14;
  2. the fallback case of testAr2armaNonInvertibleMaFallsBackToTheAr -- an
     ARMA(1,2) (a = [1 -0.9], b = [1 1.2 0.81]) reduced to ARMA(1,1): a
     stationary AR part (0.940) and a NON-INVERTIBLE MA part (root 1.160).

Run:  python3 dev/validation/wp126_ar2arma_mirror.py
Created 2026-10-12 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 5.5 -- WP-126.
"""
import numpy as np
from scipy.linalg import hankel, solve_toeplitz
from scipy.signal import lfilter


def arma_acv(a, b, n, nimp=20000):
    """Autocovariance (innovation variance 1) of the ARMA b(z)/a(z), lags 0..n-1."""
    h = lfilter(b, a, np.r_[1.0, np.zeros(nimp - 1)])
    return np.array([h[:nimp - k] @ h[k:] for k in range(n)])


def levinson(r, p):
    """AR(p) polynomial [1 a_1 .. a_p] and innovation variance from the ACV r."""
    a = np.r_[1.0, -solve_toeplitz(r[:p], r[1:p + 1])]
    return a, r[:p + 1] @ a


def ar2arma(a_hi, q, p, n, fixed=True):
    """The MATLAB function's arithmetic, MATLAB indices shifted by one."""
    h = lfilter([1.0], a_hi, np.r_[1.0, np.zeros(n - 1)])
    H = hankel(h[q:n - p], h[n - p - 1:])          # h(q+1:end-p), h(end-p:end)
    if fixed:
        aa = -np.linalg.lstsq(H[:, :-1][:, ::-1], H[:, -1], rcond=None)[0]
    else:
        aa = -np.linalg.lstsq(H[:, 1:], H[:, 0], rcond=None)[0]
    a = np.r_[1.0, aa]
    hb = lfilter(a, [1.0], h)
    b = np.r_[1.0, hb[1:q + 1]] / hb[0]
    return a, b


def maxroot(c):
    return max(abs(np.roots(c))) if len(c) > 1 and np.any(c[1:]) else 0.0


def main():
    rad, ang = 0.94, 2 * np.pi * 0.12
    a_true = np.real(np.poly([rad * np.exp(1j * ang), rad * np.exp(-1j * ang)]))
    b_true = np.array([1.0, 0.5])
    r = arma_acv(a_true, b_true, 256)
    a_hi, s2 = levinson(r, 60)
    a_old, _ = ar2arma(a_hi, 1, 2, 4096, fixed=False)
    a, b = ar2arma(a_hi, 1, 2, 4096)
    rel = np.linalg.norm(s2 * arma_acv(a, b, 256) - r) / np.linalg.norm(r)
    print(f"1. old orientation: max|root a| = {maxroot(a_old):.4f} (true 0.94)")
    print(f"   fixed: |a - a_true| = {abs(a - a_true).max():.1e}, "
          f"|b - b_true| = {abs(b - b_true).max():.1e}, ACV relErr = {rel:.1e}")
    assert maxroot(a_old) > 1
    assert abs(a - a_true).max() < 1e-12 and abs(b - b_true).max() < 1e-12 and rel < 1e-12

    r2 = arma_acv(np.array([1.0, -0.9]), np.array([1.0, 1.2, 0.81]), 256)
    a_hi2, _ = levinson(r2, 60)
    a2, b2 = ar2arma(a_hi2, 1, 1, 1024)
    print(f"2. ARMA(1,2) as ARMA(1,1): a = {np.round(a2, 4)}, b = {np.round(b2, 4)}; "
          f"|root a| = {maxroot(a2):.3f}, |root b| = {maxroot(b2):.3f}")
    assert maxroot(a2) < 1 and maxroot(b2) > 1.1


if __name__ == "__main__":
    main()
