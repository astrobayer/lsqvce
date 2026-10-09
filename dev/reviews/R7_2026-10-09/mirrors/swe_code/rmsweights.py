# Mirror of adjustVCEjoint>epochResiduals (d-axis epoch -> RMS over its d rows)
# feeding lsqvce.kernel.robustWeights (s = 1.4826*MAD about the median, w = vw/s NOT centred).
import numpy as np
def robustWeights(vw, infl, tune=None):
    med = np.median(vw); s = 1.4826*np.median(np.abs(vw-med))
    if not (s > 0 and np.isfinite(s)):
        s = np.std(vw, ddof=1)
        if not (s > 0 and np.isfinite(s)): s = 1.0
    w = vw/s; aw = np.abs(w)
    if infl == 'huber':
        c = 1.345 if tune is None else tune; fac = np.ones_like(aw); m = aw > c; fac[m] = c/aw[m]
    elif infl == 'igg3':
        k0, k1 = (1.5, 3.0) if tune is None else tune
        fac = np.ones_like(aw); mid = (aw > k0) & (aw <= k1)
        fac[mid] = (k0/aw[mid])*((k1-aw[mid])/(k1-k0))**2; fac[aw > k1] = 0
    else:
        c = 4.685 if tune is None else tune; fac = np.zeros_like(aw); m = aw < c; fac[m] = (1-(aw[m]/c)**2)**2
    return fac, s
rng = np.random.default_rng(1)
n = 200000
for d in (1, 2, 3, 6, 9, 12):
    v = rng.standard_normal((n, d))       # CLEAN, whitened, unit level
    u = v[:, 0] if d == 1 else np.sqrt(np.mean(v**2, axis=1))
    for infl in ('huber', 'igg3', 'bisquare'):
        fac, s = robustWeights(u, infl)
        print('d=%d %-8s s=%.3f  mean factor=%.3f  zeroed=%.1f%%  factor<0.5=%.1f%%' %
              (d, infl, s, fac.mean(), 100*np.mean(fac == 0), 100*np.mean((fac > 0) & (fac < 0.5))))
