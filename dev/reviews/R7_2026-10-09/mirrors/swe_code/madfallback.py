# robustWeights' silent MAD -> std -> 1 fallback: when does it change the answer?
import numpy as np, importlib.util, sys
spec = importlib.util.spec_from_file_location('m', '/tmp/claude-0/-home-claude/de66562f-20a2-5cf1-9c5b-d39655ca4eeb/scratchpad/mirrors/swe_code/rmsweights_lib.py')
def robustWeights(vw, infl):
    med = np.median(vw); s = 1.4826*np.median(np.abs(vw-med)); path='MAD'
    if not (s > 0 and np.isfinite(s)):
        s = np.std(vw, ddof=1); path='std'
        if not (s > 0 and np.isfinite(s)): s = 1.0; path='1'
    aw = np.abs(vw/s)
    if infl == 'igg3':
        fac = np.ones_like(aw); mid=(aw>1.5)&(aw<=3); fac[mid]=(1.5/aw[mid])*((3-aw[mid])/1.5)**2; fac[aw>3]=0
    else:
        fac = np.ones_like(aw); m = aw>1.345; fac[m]=1.345/aw[m]
    return fac, s, path
rng = np.random.default_rng(7)
# quantised residuals (resolution 1 unit, true sigma 0.4 unit): >50% round to 0 -> MAD = 0
import sys
BL = float(sys.argv[1]); n = 10000; v = np.round(0.4*rng.standard_normal(n))
bl = rng.choice(n, 500, replace=False); v[bl] += BL
for infl in ('huber','igg3'):
    fac, s, path = robustWeights(v, infl)
    print('%-5s scale path=%s s=%.3f  blunder factor median=%.3f  blunders zeroed=%d/500  clean mean factor=%.3f'
          % (infl, path, s, np.median(fac[bl]), np.sum(fac[bl]==0), np.mean(np.delete(fac, bl))))
print('fraction exactly at median:', np.mean(v == np.median(v)))
