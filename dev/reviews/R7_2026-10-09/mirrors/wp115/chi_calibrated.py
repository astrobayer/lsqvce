"""WP-115 (R7-08): mirror of adjustVCEjoint>chiCalibrated -- the chi-calibrated
multi-axis epoch statistic -- against the RMS it replaces, through the same
uncentred IGG3 standardisation as kernel.robustWeights. Tier B evidence of the
mechanism, not a MATLAB result."""
import numpy as np
from scipy.special import gammaincc, erfcinv
from scipy.stats import chi2

rng = np.random.default_rng(1)


def igg3(vw, k0=1.5, k1=3.0):
    med = np.median(vw); s = 1.4826 * np.median(abs(vw - med)); aw = abs(vw / s)
    f = np.ones_like(aw); mid = (aw > k0) & (aw <= k1)
    f[mid] = (k0 / aw[mid]) * ((k1 - aw[mid]) / (k1 - k0)) ** 2; f[aw > k1] = 0
    return f


def chi_cal(V):
    d = V.shape[1]; q = (V ** 2).sum(1); sg = np.sign(V.sum(1)); sg[sg == 0] = 1
    s2 = np.median(q / chi2.median(d))
    return sg * np.minimum(np.sqrt(2) * erfcinv(gammaincc(d / 2, q / (2 * s2))), 40)


for d in (1, 3, 6, 12):
    V = rng.standard_normal((200000, d)) * 2.7          # an arbitrary level
    z = V[:, 0] if d == 1 else chi_cal(V)
    line = 'd=%2d  chi-cal: sd(z)/level=%.3f zeroed=%.4f mean factor=%.3f' % (
        d, z.std() / (2.7 if d == 1 else 1), np.mean(igg3(z) == 0), igg3(z).mean())
    if d > 1:
        line += '   | old RMS: zeroed=%.3f' % np.mean(igg3(np.sqrt((V ** 2).mean(1))) == 0)
    print(line)
V = rng.standard_normal((20000, 3)); V[:2000, 0] += 8   # 10 % blunders, 8 sigma on one axis
f = igg3(chi_cal(V))
print('3-axis, 10%% blunders of 8 sigma: blunders zeroed %.3f, clean zeroed %.4f'
      % (np.mean(f[:2000] == 0), np.mean(f[2000:] == 0)))
