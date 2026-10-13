#!/usr/bin/env python3
"""WP-134 mirror: the kernel findings of the third whole-code review (§12.68).

(9) robustWeights tested a near-tie MAD against 1e3*eps of the LARGEST
    deviation -- the worst blunder. One residual beyond ~4.5e12 MAD flipped
    the scale to the std that the same blunder inflates: Huber's estimate
    moved by 2e8 and 8-sigma rows kept factor 1. The reference is now the
    75th percentile of the deviations (the bulk's spread; a blunder in the
    upper quarter cannot move it). The R7-21 fixtures (51 % exact ties, 51 %
    near-ties at 1e-15) still fall back.
(low) armaACV started the impulse response at max(4n, 512) and doubled at
    most 7 times: for n <= 128 the cap was 65 536 samples, not 2^20, and a
    near-unit root was truncated n-dependently (|root| 0.99999: r(0) 27 %
    low at n = 2). The length now follows the decay alone, from 512 to 2^20
    for any n, with one FFT of nh + n - 1 points (no wrap); the cap's
    truncation rho^(2 nh) is reported. The speed point armaACV(4n)/armaACV(n)
    (n = 2000) now times FFTs of 16384 / 4096 points instead of 65536 / 16384.

Run:  python3 dev/validation/wp134_kernel_mirror.py
Created 2026-10-13 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 5.5 -- WP-134.
"""
import numpy as np
import scipy.signal as ss

EPS = np.finfo(float).eps


def nextpow2(x):
    return int(np.ceil(np.log2(x))) if x > 1 else 0


def arma_acv_old(a, b, n):
    nh = max(4 * n, 512)
    for _ in range(8):
        imp = np.zeros(nh)
        imp[0] = 1
        h = ss.lfilter(b, a, imp)
        K = min(nh, max(len(a), len(b)))
        if np.sum(h[-K:] ** 2) < EPS ** 2 * np.sum(h ** 2) or nh >= 2 ** 20:
            break
        nh *= 2
    nfft = 2 ** nextpow2(2 * nh)
    R = np.fft.irfft(np.abs(np.fft.rfft(h, nfft)) ** 2, nfft)
    return R[:n], nh, nfft


def arma_acv_new(a, b, n):
    nh = 512
    K = max(len(a), len(b))
    while True:
        imp = np.zeros(nh)
        imp[0] = 1
        h = ss.lfilter(b, a, imp)
        if np.sum(h[-min(nh, K):] ** 2) < EPS ** 2 * np.sum(h ** 2) or nh >= 2 ** 20:
            break
        nh *= 2
    rts = np.roots(a)
    rho = max(np.abs(rts)) if len(rts) else 0.0
    warned = rho ** (2 * nh) > EPS
    nfft = 2 ** nextpow2(nh + n - 1)
    R = np.fft.irfft(np.abs(np.fft.rfft(h, nfft)) ** 2, nfft)
    return R[:n], nh, nfft, warned


def robust_scale(vw, ref):
    med = np.median(vw)
    dev = np.abs(vw - med)
    s = 1.4826 * np.median(dev)
    if ref == 'max':
        r = dev.max()
    else:
        d = np.sort(dev)
        r = d[int(np.ceil(0.75 * len(d))) - 1]
    if not np.isfinite(s) or s <= 1e3 * EPS * r:
        return np.std(vw, ddof=1), 'std'
    return s, 'mad'


def huber_loc(y, ref, iters=50):
    x = np.median(y)
    for _ in range(iters):
        v = y - x
        s, path = robust_scale(v, ref)
        w = np.abs(v) / s
        fac = np.where(w <= 1.345, 1.0, 1.345 / np.maximum(w, 1e-300))
        x = np.sum(fac * y) / np.sum(fac)
    v = y - x
    s, path = robust_scale(v, ref)
    w = np.abs(v) / s
    fac = np.where(w <= 1.345, 1.0, 1.345 / np.maximum(w, 1e-300))
    return x, s, path, fac


def main():
    ok = True
    print('== armaACV: AR(1) r(0) against 1/(1 - phi^2) ==')
    for phi in (0.9, 0.9999, 0.99995, 0.99999):
        r0 = 1 / (1 - phi ** 2)
        for n in (2, 100, 1000):
            ro, nho, _ = arma_acv_old([1, -phi], [1], n)
            rn, nhn, _, wn = arma_acv_new([1, -phi], [1], n)
            eo, en = ro[0] / r0 - 1, rn[0] / r0 - 1
            print(f'  phi {phi:<8} n {n:>5}: old {eo:+.2e} (nh {nho:>7})   new {en:+.2e} '
                  f'(nh {nhn:>7}{", warned" if wn else ""})')
            # tol: the dominant pole's tail rho^(2 nh) bounds the truncation;
            # 1e-12 relative above it for the FFT roundoff on r(0) ~ 1/(1-phi^2)
            bound = phi ** (2 * nhn) + 1e-12
            ok &= abs(en) <= bound
            ok &= wn == (phi ** (2 * nhn) > EPS)
    print('== the same law for any n (r(1), phi 0.99999) ==')
    r1 = [arma_acv_new([1, -0.99999], [1], n)[0][1] for n in (2, 100, 1000)]
    spread = (max(r1) - min(r1)) / r1[0]
    print(f'  r(1) over n = 2, 100, 1000: relative spread {spread:.1e}')
    ok &= spread < 1e-12                       # tol: FFT roundoff of three lengths
    print('== ARMA(1,1), sparse AR, narrow AR(2) line ==')
    r = arma_acv_new([1, -0.9], [1, 0.5], 3)[0]
    r0 = (1 + 2 * 0.9 * 0.5 + 0.25) / (1 - 0.81)
    print(f'  ARMA(1,1) r0 {r[0]:.12f} vs {r0:.12f}, r1 {r[1]:.12f} vs {0.9 * r0 + 0.5:.12f}')
    ok &= abs(r[0] / r0 - 1) < 1e-13 and abs(r[1] / (0.9 * r0 + 0.5) - 1) < 1e-13
    rs = arma_acv_new([1, 0, 0.998], [1], 2)[0]
    print(f'  [1 0 0.998] r0 {rs[0]:.6f} vs {1 / (1 - 0.998 ** 2):.6f}')
    ok &= abs(rs[0] * (1 - 0.998 ** 2) - 1) < 1e-12
    rho, th = 0.99999, 0.3
    a2 = [1, -2 * rho * np.cos(th), rho ** 2]
    ro2 = arma_acv_old(a2, [1], 12)[0][0]
    rn2, nh2, _, w2 = arma_acv_new(a2, [1], 12)
    rref = arma_acv_new(a2, [1], 5000)[0][0]
    print(f'  AR(2) line rho {rho}: old r0 {ro2:.6e}, new {rn2[0]:.6e} (nh {nh2}, warned {w2}), '
          f'n = 5000: {rref:.6e}')
    ok &= abs(rn2[0] / rref - 1) < 1e-12
    print('== speed point: the polynomial of test_speedBudget, n = 2000 / 8000 ==')
    a = [1, -1.6, 0.9, -0.25]
    for n in (2000, 8000):
        _, nho, nfo = arma_acv_old(a, [1], n)
        _, nhn, nfn, _ = arma_acv_new(a, [1], n)
        print(f'  n {n}: old nh {nho} nfft {nfo};  new nh {nhn} nfft {nfn}')
    ok &= np.allclose(arma_acv_old(a, [1], 2000)[0], arma_acv_new(a, [1], 2000)[0],
                      rtol=0, atol=1e-13 * arma_acv_new(a, [1], 1)[0][0])
    print('== robustWeights: one extreme blunder ==')
    rng = np.random.default_rng(134)
    y = np.concatenate([rng.standard_normal(994), 8 + 0.1 * rng.standard_normal(5), [0.0]])
    for blunder in (4e12, 5e12, 1e13):
        y[-1] = blunder
        for ref in ('max', 'q75'):
            x, s, path, fac = huber_loc(y, ref)
            print(f'  blunder {blunder:.0e} ref {ref}: path {path}, s {s:.3g}, x {x:+.3g}, '
                  f'8-sigma factor {np.median(fac[994:999]):.3f}, blunder factor {fac[-1]:.1e}')
            if ref == 'q75':
                ok &= path == 'mad' and abs(s - 1) < 0.15 and fac[-1] < 1e-11
    print('== the R7-21 fixtures still fall back ==')
    vw = np.concatenate([np.zeros(51), np.linspace(-1, 1, 44), [20, -20, 20, -20, 20]])
    near = vw.copy()
    near[:51] = 1e-15 * np.arange(1, 52)
    for name, v in (('ties', vw), ('near-ties', near)):
        _, p = robust_scale(v, 'q75')
        print(f'  {name}: path {p}')
        ok &= p == 'std'
    _, p = robust_scale(np.linspace(-3, 3, 101), 'q75')
    ok &= p == 'mad'
    print('ALL CHECKS PASS' if ok else 'A CHECK FAILED')
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
