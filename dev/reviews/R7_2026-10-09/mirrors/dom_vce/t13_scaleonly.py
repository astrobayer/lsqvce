"""Control-flow mirror of adjustVCEjoint lines 471-489 + 599 for a run in which NO
stream has estimate "law" (demo_lsqvce case 5: two white groups with per-row sigma
vectors -> streamPolicy makes them "scale"; group 2 truly 3x noisier than declared).
lawCols all false -> groups = {} -> pvMax = 0 -> pv.maxDev = 0 <= Tol -> converged
on pass 1.  The returned sol is pass 1's (accumulated BEFORE the level update)."""
import numpy as np
rs = np.random.default_rng(2026)
n, u = 2000, 8
A = rs.normal(size=(n, u)); xt = rs.normal(size=u); sg = 0.5 + rs.random(n)
g = np.r_[np.ones(n//2), 2*np.ones(n//2)]
sgT = sg.copy(); sgT[g == 2] *= 3
l = A@xt + sgT*rs.normal(size=n)
def solve(level):
    s = sg*np.sqrt(np.where(g == 1, level[0], level[1]))
    Ab = A/s[:, None]; lb = l/s
    N = Ab.T@Ab; x = np.linalg.solve(N, Ab.T@lb); v = lb - Ab@x
    H = Ab@np.linalg.solve(N, Ab.T); r = 1 - np.diag(H)
    s0sq = v@v/(n-u)
    lev = [ (v[g == k]@v[g == k])/r[g == k].sum() * level[k-1] for k in (1, 2)]   # Helmert (exact, = alpha2*level)
    return x, np.sqrt(s0sq), np.sqrt(np.diag(np.linalg.inv(N))*s0sq), lev
lev0 = [1.0, 1.0]
x1, s01, se1, lev1 = solve(lev0)                    # pass 1: accumulated under the DECLARED weights
pvMax = 0.0                                         # no law group
converged = pvMax <= 1e-2
print('pass 1: converged =', converged, ' (maxDev =', pvMax, ')')
print('  returned sol: sigma0 = %.3f   (computed under level [1 1])' % s01)
print('  returned noise: level = [%.3f %.3f]  -> R.history.sigma ratio = %.3f (truth 3)' % (lev1[0], lev1[1], np.sqrt(lev1[1]/lev1[0])))
print('  Manual reconstruction declared*sqrt(level)*R.sigma0 / declared: group1 %.3f (truth 1), group2 %.3f (truth 3)'
      % (np.sqrt(lev1[0])*s01, np.sqrt(lev1[1])*s01))
x2, s02, se2, _ = solve(lev1)
print('what one more solve under the estimated levels gives: sigma0 = %.3f;  |x - xt|/se : returned %s  re-solved %s'
      % (s02, np.round(np.abs(x1-xt)/se1, 2), np.round(np.abs(x2-xt)/se2, 2)))
print('formal se returned / re-solved (mean): %.3f' % np.mean(se1/se2))
