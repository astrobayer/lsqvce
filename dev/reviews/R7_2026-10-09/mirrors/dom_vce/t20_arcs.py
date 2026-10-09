import numpy as np
from t3_robust_level import igg3
rng = np.random.default_rng(1)
# (a) robustWeights arithmetic on one stream whose arc 2 of 3 is 3x noisier (declared per-arc factor, NOT a blunder):
m = 256*9; arc = np.repeat(np.arange(9) % 3, 256)
v = rng.normal(size=m)*np.where(arc == 1, 3.0, 1.0)
s = 1.4826*np.median(np.abs(v - np.median(v))); w = igg3(v/s)
for a in range(3):
    print(f'(a) arc type {a} (sigma {3 if a == 1 else 1}): mean IGG3 factor {w[arc == a].mean():.3f}, zeroed {np.mean(w[arc == a] == 0)*100:.1f} %')
# (b) what the per-arc factor does to x when it never reaches the whitener: mean of 9 arcs, one arc sigma k
for k in (3.0, 30.0):
    sig2 = np.ones(9); sig2[4] = k**2
    eq = sig2.sum()/81; opt = 1/np.sum(1/sig2)
    print(f'(b) one of 9 arcs {k:.0f}x worse: variance of the equal-weight estimate / Ellmer-weighted = {eq/opt:.1f}')
