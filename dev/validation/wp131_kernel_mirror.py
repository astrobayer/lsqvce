#!/usr/bin/env python3
"""WP-131 mirror: three kernel findings of the second whole-code review.

(6) redundancyLags>sigmaInvSubdiagSums took Sigma for Toeplitz when columns 1,
    2 and ceil(N/2) agreed with a Toeplitz matrix -- a per-row sigma after
    N/2 passed, and the sums came from inv(toeplitz(c1)): trace(Sigma^-1) 40
    instead of 31.1 (N = 40, sigma = 3 from row 31). A deterministic probe
    y, Sigma y against toeplitz(c1) y, now refuses it; a true Toeplitz Sigma
    still passes.
(7) armaACV stopped the impulse response when ONE sample h(end) was below
    eps max|h|; a sparse AR polynomial is exactly 0 there long before it has
    decayed: r(0) of a = [1 0 0.998] was 160.5 instead of 250.25. The tail
    energy over max(p,q)+1 samples stops at the right length.
(13) crossCorr tested the covariance R(tau) against the band z/sqrt(nE) of a
    correlation: white noise at sigma 2 failed every time, a leftover AR(1)
    at sigma 0.5 was never seen. The correlation Rn = R / sqrt(R0_ii R0_jj)
    makes the verdict scale-free.

Run:  python3 dev/validation/wp131_kernel_mirror.py
Created 2026-10-13 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 5.5 -- WP-131.
"""
import numpy as np
import scipy.linalg as sl
import scipy.signal as ss

EPS = np.finfo(float).eps


def is_toep(W, probe):
    N = W.shape[0]
    c1, c2 = W.T @ W[:, 0], W.T @ W[:, 1]
    jm = max(3, int(np.ceil(N / 2)))
    cm = W.T @ W[:, jm - 1]
    tol = 1e-10 * abs(c1).max()
    ok = (max(abs(c2[1:] - c1[:-1])) <= tol and abs(c2[0] - c1[1]) <= tol
          and max(abs(cm[jm - 1:] - c1[:N - jm + 1])) <= tol and max(abs(cm[:jm] - c1[jm - 1::-1])) <= tol)
    if ok and probe:
        y = np.sin(np.sqrt(2) * np.arange(1, N + 1) ** 1.5)
        ok = np.linalg.norm(W.T @ (W @ y) - sl.toeplitz(c1) @ y) <= \
            1e-10 * np.linalg.norm(abs(W.T) @ (abs(W) @ abs(y)))
    return ok


def arma_r0(a, b, tail):
    nh = 512
    for _ in range(8):
        imp = np.zeros(nh)
        imp[0] = 1
        h = ss.lfilter(b, a, imp)
        K = min(nh, max(len(a), len(b)))
        stop = (np.sum(h[-K:] ** 2) < EPS ** 2 * np.sum(h ** 2)) if tail else (abs(h[-1]) < EPS * abs(h).max())
        if stop or nh >= 2 ** 20:
            break
        nh *= 2
    return np.sum(h ** 2)


def verdict(w, d, corr):
    nE = len(w) // d
    W = w[:nE * d].reshape(nE, d)
    mx = 10
    R = np.array([W[:nE - t].T @ W[t:] / (nE - t) for t in range(mx + 1)])
    if corr:
        sd = np.sqrt(np.diag(R[0]))
        R = R / np.outer(sd, sd)
    band = 1.959964 / np.sqrt(nE)
    cross0 = np.abs(R[0] - np.diag(np.diag(R[0]))).max()
    diagT = max(np.abs(np.diag(R[t])).max() for t in range(1, mx + 1))
    crossT = max(np.abs(R[t] - np.diag(np.diag(R[t]))).max() for t in range(1, mx + 1))
    return 'ok' if max(cross0, diagT, crossT) <= 3 * band else ('ar-order-low' if diagT > 3 * band else 'other')


def main():
    N = 40
    s = np.ones(N)
    s[30:] = 3
    W = np.diag(s)
    print(f'(6) per-row sigma: old Toeplitz {is_toep(W, False)}, new {is_toep(W, True)} '
          f'(trace inv(toeplitz(c1)) {np.trace(np.linalg.inv(sl.toeplitz(W.T @ W[:, 0]))):.1f}, '
          f'true {np.trace(np.linalg.inv(W.T @ W)):.1f})')
    assert is_toep(W, False) and not is_toep(W, True)
    Wt = sl.cholesky(sl.toeplitz(0.7 ** np.arange(N)))
    assert is_toep(Wt, False) and is_toep(Wt, True)
    print('    a true Toeplitz Sigma passes both')
    for a in ([1, 0, 0.998], [1, 0, 0, 0, -0.995]):
        old, new = arma_r0(a, [1], False), arma_r0(a, [1], True)
        print(f'(7) a = {a}: r0 old {old:.2f}, new {new:.2f}')
    assert abs(arma_r0([1, 0, 0.998], [1], True) - 1 / (1 - 0.998 ** 2)) < 1e-6 * 250
    rng = np.random.default_rng(7)
    res = {}
    for sig in (1.0, 2.0):
        fails_old = sum(verdict(sig * rng.standard_normal(12000), 3, False) != 'ok' for _ in range(50))
        fails_new = sum(verdict(sig * rng.standard_normal(12000), 3, True) != 'ok' for _ in range(50))
        res[sig] = (fails_old, fails_new)
        print(f'(13) white, sigma {sig}: non-ok old {fails_old}/50, new {fails_new}/50')
    assert res[2.0][0] == 50 and res[2.0][1] <= 3
    e = ss.lfilter([1], [1, -0.1], 0.5 * rng.standard_normal((4000, 3)), axis=0).ravel()   # per axis
    print(f'     AR(1) 0.1 leftover at sigma 0.5: old {verdict(e, 3, False)}, new {verdict(e, 3, True)}')
    assert verdict(e, 3, False) == 'ok' and verdict(e, 3, True) == 'ar-order-low'


if __name__ == '__main__':
    main()
