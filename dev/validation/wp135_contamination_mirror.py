#!/usr/bin/env python3
"""WP-135 mirror: the contamination ladder's n-dependent calibration
(third whole-code review, finding 19; user decision: a measurement).

The ladder read fhat -- the share of robustly standardised residuals beyond
2.5 -- against a null band calibrated at n ~ 300 (tol = P + 4 * 1.30 *
sqrt(P(1-P)/n), P = 2 Phi(-2.5)) and called fhat <= 0.005 "collapsed".
Measured here on clean data (the diagnose port: LS residuals, raw scores
v / sqrt(qv), the MAD scale):

  * the mean of fhat is P (1 + 27/n), not P: the MAD of a small sample is
    biased low, more residuals pass 2.5 robust sigmas;
  * around that mean the count K = n fhat is over-dispersed by a CONSTANT
    1.30^2 for n = 20..2000 (the 1.30 the code carried, measured at n ~ 300
    against the asymptotic mean, was the large-n end of the same thing);
  * so K is negative binomial: mean mu = n P (1 + 27/n), variance 1.69 mu.

Calibration (both tails from that law, n = the scored residuals):
  fhatTol(n)  = kHi / n, kHi the smallest k with P(K > k) <= 0.25 %  (the
                false-alarm rate the code always claimed for k = 4);
  collapsed   = P(K <= n fhat) <= 1 % -- fhat is implausibly low for clean
                data. Below n ~ 450 even fhat = 0 is not (P(K = 0) > 1 %): a
                MAD captured by the blunders cannot be told from a sigma0
                declared too small, and the ladder says so instead of
                "breakdown".
  The order: breakdown (backstop big AND collapsed), then clean (fhat <=
  fhatTol, whatever the backstop -- a uniformly inflated residual is the
  global test's finding), then refuse / severe / moderate.

Before / after (clean data, 4000+ draws per cell):
  correct sigma, any alarm   n = 30: 3.2 % -> <= 0.3 %
  sigma 3.5x off, breakdown  n = 30: ~50 % -> 0 %;  n = 500: 8.6 % -> <= 1 %

WP-137 (the code review of WP-135): ladder_new is now the MATLAB ladder in
full -- 'unknown' without a robust scale, and 'indeterminate' (backstop big,
fhat in its clean band, a collapse undetectable at this n) before 'clean'.
The WP-135 port lacked both, so its "after" rates were measured on a ladder
the code does not run, and its 55 % column read 0 % breakdown below n ~ 450
without saying where those runs went: to 'indeterminate', which snoop then
passed in silence. The table now counts both levels, and the check asserts
that a 55 % contamination is never 'clean' (breakdown -> snoop refuses,
indeterminate -> snoop warns, WP-137).

Run:  python3 dev/validation/wp135_contamination_mirror.py
Created 2026-10-14 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 5.5 -- WP-135.
Updated 2026-10-15 (CEST, UTC+02:00), WP-137: ladder_new is the MATLAB ladder in
full ('unknown', 'indeterminate'); the 55 % cell counts both. With assistance
from Claude (Anthropic), model Claude Opus 5.5.
"""
import math
import numpy as np

P = 2 * 0.5 * (1 + math.erf(-2.5 / math.sqrt(2)))      # 0.012419...
PHI = 1.30 ** 2                                          # measured over-dispersion
BIAS = 27.0                                              # measured: mean = P (1 + 27/n)
A_HI, A_LO = 0.0025, 0.01
Z_HI = 2.807033768343811                                 # normInv(1 - 0.0025)


def nb_law(n):
    mu = n * P * (1 + BIAS / n)
    return mu, mu / (PHI - 1)                            # mean, NB size r


def nb_cdf_upto(kmax, n):
    """P(K <= k) for k = 0..kmax (log-space recursion; the MATLAB port)."""
    mu, r = nb_law(n)
    q = mu / (r + mu)
    lp = r * math.log1p(-q)                              # log pmf(0)
    cdf = np.empty(kmax + 1)
    c = 0.0
    for k in range(kmax + 1):
        if k > 0:
            lp += math.log((r + k - 1) / k * q)
        c += math.exp(lp)
        cdf[k] = c
    return cdf


def fhat_tol(n):
    mu, _ = nb_law(n)
    if mu > 200:                                         # the normal limit (exact enough)
        return (mu + Z_HI * math.sqrt(PHI * mu)) / n
    cdf = nb_cdf_upto(int(math.ceil(mu + 12 * math.sqrt(PHI * mu) + 20)), n)
    return int(np.argmax(cdf >= 1 - A_HI)) / n


def collapsed(fhat, n):
    mu, _ = nb_law(n)
    k = int(round(fhat * n))
    if mu > 200:
        z = (k + 0.5 - mu) / math.sqrt(PHI * mu)
        return 0.5 * math.erfc(-z / math.sqrt(2)) <= A_LO
    return nb_cdf_upto(k, n)[k] <= A_LO


def diag_port(n, u, rng, c=1.0, contam=0.0):
    A = rng.standard_normal((n, u))
    e = c * rng.standard_normal(n)
    if contam > 0:
        k = int(round(contam * n))
        idx = rng.choice(n, k, replace=False)
        e[idx] += 20 * np.sign(rng.standard_normal(k)) * (1 + rng.random(k))
    Q, _ = np.linalg.qr(A)
    v = e - Q @ (Q.T @ e)
    qv = 1 - np.sum(Q * Q, axis=1)
    raw = v / np.sqrt(qv)
    rs = 1.4826 * np.median(np.abs(raw - np.median(raw)))
    fhat = np.mean(np.abs(raw / rs) > 2.5)
    r = n - u
    backstop = math.sqrt(v @ v / r)                       # sigma0 = 1 declared
    big = backstop > max(3, 2 * (1 + 3 / math.sqrt(2 * r)))
    return fhat, big


def ladder_old(fhat, big, n):
    tol = P + 4 * 1.30 * math.sqrt(P * (1 - P) / n)
    if fhat <= 0.005 and big:
        return 'breakdown'
    if fhat <= tol and not big:
        return 'clean'
    if fhat >= 0.25:
        return 'refuse'
    if fhat >= 0.10:
        return 'severe'
    return 'moderate' if fhat > tol else 'clean'


def collapse_detectable(n):
    """P(K = 0) <= 1 %: fhat = 0 is a collapse only where clean data rarely show it."""
    mu, r = nb_law(n)
    return math.exp(r * math.log1p(-mu / (r + mu))) <= A_LO


def ladder_new(fhat, big, n, mad_available=True):
    """verdict>contaminationOf, every level (WP-137: 'unknown', 'indeterminate')."""
    if not mad_available:
        return 'unknown'
    if big and collapsed(fhat, n):
        return 'breakdown'
    if fhat <= fhat_tol(n):
        return 'indeterminate' if big and not collapse_detectable(n) else 'clean'
    if fhat >= 0.25:
        return 'refuse'
    return 'severe' if fhat >= 0.10 else 'moderate'


def main():
    ok = True
    rng = np.random.default_rng(135)
    print('== the law against the measured null (clean, sigma right) ==')
    for n in (20, 30, 50, 100, 200, 500, 1000, 2000):
        u = min(8, max(3, n // 10))
        T = 20000 if n <= 200 else 6000
        f = np.array([diag_port(n, u, rng)[0] for _ in range(T)])
        mu, _ = nb_law(n)
        p0 = nb_cdf_upto(0, n)[0] if mu <= 200 else 0.0
        fa = np.mean(f > fhat_tol(n))
        print(f'  n {n:5d}: mean {f.mean():.4f} (law {mu / n:.4f}); P(fhat=0) {np.mean(f == 0):.4f} '
              f'(law {p0:.4f}); tol {fhat_tol(n):.4f} -> false alarms {100 * fa:.2f} %')
        # tol: the 0.25 % target over T draws has a binomial sd of
        # sqrt(0.0025/T) <= 0.07 %; 0.45 % is the target plus ~3 sd
        ok &= fa <= 0.0045
        # tol: the mean law against the sample mean, sd(fhat)/sqrt(T) ~ 2e-4
        ok &= abs(f.mean() - mu / n) <= 0.1 * mu / n
    print('== before / after on the review\'s cells (4000 draws) ==')
    for n, u in ((30, 4), (50, 5), (100, 5), (200, 8), (500, 8)):
        res = {}
        for name, c, ct in (('clean', 1.0, 0.0), ('sigma 3.5x', 3.5, 0.0),
                            ('5 % blunders', 1.0, 0.05), ('55 % blunders', 1.0, 0.55)):
            cnt_o, cnt_n = {}, {}
            for _ in range(4000):
                fh, big = diag_port(n, u, rng, c, ct)
                lo, ln = ladder_old(fh, big, n), ladder_new(fh, big, n)
                cnt_o[lo] = cnt_o.get(lo, 0) + 1
                cnt_n[ln] = cnt_n.get(ln, 0) + 1
            res[name] = (cnt_o, cnt_n)
        alarm = lambda c: sum(v for k, v in c.items() if k != 'clean') / 4000
        bd = lambda c: c.get('breakdown', 0) / 4000
        ind = lambda c: c.get('indeterminate', 0) / 4000
        det = lambda c: sum(v for k, v in c.items() if k in ('moderate', 'severe', 'refuse')) / 4000
        print(f'  n {n:4d}: clean any alarm {100 * alarm(res["clean"][0]):5.2f} % -> '
              f'{100 * alarm(res["clean"][1]):5.2f} %;  sigma 3.5x breakdown '
              f'{100 * bd(res["sigma 3.5x"][0]):5.1f} % -> {100 * bd(res["sigma 3.5x"][1]):4.1f} %;  '
              f'5 % blunders detected {100 * det(res["5 % blunders"][0]):5.1f} % -> '
              f'{100 * det(res["5 % blunders"][1]):5.1f} %;  55 % breakdown '
              f'{100 * bd(res["55 % blunders"][0]):5.1f} % -> {100 * bd(res["55 % blunders"][1]):5.1f} % '
              f'(+ indeterminate {100 * ind(res["55 % blunders"][1]):5.1f} %);  sigma 3.5x '
              f'indeterminate {100 * ind(res["sigma 3.5x"][1]):5.1f} %')
        # tol: binomial sd of a rate over 4000 draws <= 0.25 % at these levels
        ok &= alarm(res['clean'][1]) <= 0.006 and bd(res['sigma 3.5x'][1]) <= 0.015
        # WP-137: a 55 % contamination is never 'clean' -- breakdown (snoop
        # refuses) or indeterminate (snoop warns); tol: as above
        ok &= bd(res['55 % blunders'][1]) + ind(res['55 % blunders'][1]) >= 0.99
    print('ALL CHECKS PASS' if ok else 'A CHECK FAILED')
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
