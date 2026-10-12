"""WP-120 mirror (§12.48 residue 2): the size and power of the pooled whiteness test.

The joint estimator tests a stream's law-whitened residuals with the POOLED Ljung-Box
statistic (R7-38): Q = N(N+2) sum_k rho_k^2 / N_k over lags 1..K, K = 20 fixed. The law
is a spectrum on Noise.Grid bins (default 60), so it resolves structure down to
~ Nyquist/Grid in frequency -- correlation out to ~ 2 Grid lags in time. A misfit the
law leaves at long lags (a narrow-band line, a low-frequency excess) is invisible at
K = 20 on 2160-epoch arcs.

Measured (user decision 2026-10-12: K = 20 STAYS, the blind spot recorded): pooled LB at
K = 20 / 60 / 120 / 200 / 540 and the two-scale test (K = 20 and K = 2 Grid = 120,
Bonferroni) on weak misfits, 8 arcs x 2160 epochs (and 9 x 256 at K = 20 / 60):
  white     e ~ N(0,1)                                                      (size)
  ar1       short-lag misfit: AR(1), phi = 0.025
  line      narrow band: a random-phase sinusoid, period 200 epochs, amplitude 0.04
  lowfreq   low-frequency excess: white + 0.05 * AR(1) phi = 0.995 (unit variance)
  echo150   correlation at lag 150 ONLY: (z_i + 0.08 z_{i-150}) / sqrt(1 + 0.08^2)
A sinusoid or a low-frequency excess shows at the short lags too, so a larger K does not
help it; a larger K dilutes the short-lag power; only an isolated long-lag correlation
needs K beyond its lag.
"""
import numpy as np
from scipy.special import gammaincc

def lag_sums(e, K):
    n = e.size
    nfft = 1 << int(np.ceil(np.log2(2 * n - 1)))
    F = np.fft.rfft(e, nfft)
    r = np.fft.irfft(F * np.conj(F), nfft)[:K + 1]
    s = np.zeros(K + 1); s[:min(K + 1, n)] = r[:min(K + 1, n)]
    return s

def pooled_p(arcs, K):
    s = sum(lag_sums(e, K) for e in arcs)
    Nk = sum(np.maximum(e.size - np.arange(K + 1), 0) for e in arcs).astype(float)
    m, r0 = Nk[0], s[0]
    ok = Nk[1:] > 0
    rho = s[1:] / r0
    Q = m * (m + 2) * np.sum(rho[ok] ** 2 / Nk[1:][ok])
    return gammaincc(ok.sum() / 2, Q / 2)

def ar1(rs, n, phi):
    e = rs.randn(n); x = np.empty(n); x[0] = e[0] / np.sqrt(1 - phi ** 2)
    for i in range(1, n): x[i] = phi * x[i - 1] + e[i]
    return x * np.sqrt(1 - phi ** 2)

def scenario(rs, kind, n):
    w = rs.randn(n)
    if kind == 'white':   return w
    if kind == 'ar1':     return ar1(rs, n, 0.025)
    if kind == 'line':    return w + 0.04 * np.sqrt(2) * np.cos(2 * np.pi * np.arange(n) / 200 + rs.uniform(0, 2 * np.pi))
    if kind == 'lowfreq': return w + 0.05 * ar1(rs, n, 0.995)
    if kind == 'echo150':
        z = rs.randn(n + 150)
        return (z[150:] + 0.08 * z[:-150]) / np.sqrt(1 + 0.08 ** 2)

if __name__ == '__main__':
    rs = np.random.RandomState(7)
    print('rejection rate at alpha 0.05 (white = size, others = power); two = K 20 + K 120, Bonferroni')
    for na, n, Ks, T in [(8, 2160, (20, 60, 120, 200, 540), 300), (9, 256, (20, 60), 400)]:
        print('\n%d arcs x %d epochs, %d trials' % (na, n, T))
        print('  %-8s ' % '' + ' '.join('%7s' % ('K%d' % K) for K in Ks) + ('     two' if 120 in Ks else ''))
        for kind in ('white', 'ar1', 'line', 'lowfreq', 'echo150'):
            rej = np.zeros(len(Ks)); two = 0
            for _ in range(T):
                arcs = [scenario(rs, kind, n) for _ in range(na)]
                ps = [pooled_p(arcs, K) for K in Ks]
                rej += np.array(ps) < 0.05
                if 120 in Ks: two += min(2 * min(ps[0], ps[Ks.index(120)]), 1) < 0.05
            print('  %-8s ' % kind + ' '.join('%7.3f' % r for r in rej / T) + ('  %6.3f' % (two / T) if 120 in Ks else ''))
