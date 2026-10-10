"""WP-118 part 2 (R7-38): size of the joint loop's pooled Ljung-Box test under a white null.
Per-arc lag sums s_k = sum_t e_t e_{t+k} are pooled over a stream's arcs (as sweepArcs does).
  shipped : Q = m(m+2) sum rho_k^2 / (m-k),  m = TOTAL count (each arc loses k pairs, not m)
  pooled  : Q = N(N+2) sum rho_k^2 / N_k,    N_k = sum_a max(m_a - k, 0) (the pairs actually summed)
K = 20 lags (lags with N_k = 0 drop out, and the dof with them). alpha 0.05, 2000 trials."""
import numpy as np
from scipy.stats import chi2
def lagsums(e, K):
    n = len(e); s = np.zeros(K + 1)
    for k in range(min(K + 1, n)): s[k] = e[:n-k] @ e[k:]
    return s
def run(nArc, m, K=20, trials=2000, seed=0):
    rng = np.random.default_rng(seed); rej = {'shipped': 0, 'pooled': 0}; Qs = {'shipped': [], 'pooled': []}
    for _ in range(trials):
        S = np.zeros(K + 1); N = 0; Nk = np.zeros(K + 1)
        for a in range(nArc):
            e = rng.normal(size=m); S += lagsums(e, K); N += m
            Nk += np.maximum(m - np.arange(K + 1), 0)
        rho = S[1:] / S[0]; k = np.arange(1, K + 1)
        Q = N * (N + 2) * np.sum(rho**2 / (N - k)); Qs['shipped'].append(Q)
        rej['shipped'] += chi2.sf(Q, K) < 0.05
        ok = Nk[1:] > 0; Kp = ok.sum()
        Q = N * (N + 2) * np.sum(rho[ok]**2 / Nk[1:][ok]); Qs['pooled'].append(Q)
        rej['pooled'] += chi2.sf(Q, Kp) < 0.05
    return {k: (rej[k] / trials, np.mean(Qs[k])) for k in rej}
for nArc, m in [(9, 256), (8, 2160), (64, 30), (200, 12), (1, 500)]:
    r = run(nArc, m, trials=1000 if m > 1000 else 2000)
    print(f'{nArc:4d} x {m:5d}: shipped size {r["shipped"][0]:.3f} E[Q] {r["shipped"][1]:6.2f} | pooled size {r["pooled"][0]:.3f} E[Q] {r["pooled"][1]:6.2f}  (K_eff {min(20, m-1)})')
