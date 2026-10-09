import numpy as np
from scipy.linalg import toeplitz
from kern import dctI, lawSpecs
L = 60
# (a) every DCT-I ortho basis column (= the covariance function a single positive amplitude implies)
neg = 0
for j in range(L):
    S = np.zeros(L); S[j] = 1.0
    c = dctI(S)                                   # psdVCE: Cnew = dctI(Shat, Inverse) = X*Shat (ortho)
    ev = np.linalg.eigvalsh(toeplitz(c)).min()
    neg += ev < -1e-12
print(f'(a) {neg}/{L} single-bin POSITIVE spectra give an INDEFINITE Toeplitz(C) on {L} lags')
# (b) the zero-lag weight: ortho DCT-I vs the symmetric cosine series of a Toeplitz ACF
C = 0.9**np.arange(L)
fN, Sl = lawSpecs(C); Sd = dctI(C)
k = np.ones(L); k[0] = k[-1] = 1/np.sqrt(2)
print('(b) DCT-ortho S_j = sqrt(2/(L-1)) k_j [C0/sqrt2 + sum_{1..L-2} Cn cos + C_{L-1}(-1)^j/sqrt2];  lawSpecs/4 = C0/2 + sum_{1..L-1} Cn cos')
# (c) random positive Shat -> C = X Shat: how often is Toeplitz(C) indefinite / |C1| > C0
rng = np.random.default_rng(0); bad = 0; badc1 = 0
for t in range(5000):
    Sh = np.exp(rng.normal(size=L))
    c = dctI(Sh)
    bad += np.linalg.eigvalsh(toeplitz(c)).min() < 0
    badc1 += abs(c[1]) > c[0]
print(f'(c) random positive (log-normal) Shat: Toeplitz(C) indefinite in {bad}/5000, |C1|>C0 in {badc1}/5000')
# (d) smooth near-truth spectra: AR(1) phi=0.9 DCT amplitudes perturbed by +-30 % per bin (one-pass VCE noise level)
bad = 0; worst = 0
Sd = dctI(0.9**np.arange(L))
print('    DCT of AR(1) 0.9 truncated at 60 lags, min/max = %.3e' % (Sd.min()/Sd.max()))
for t in range(5000):
    Sh = np.maximum(Sd*np.exp(0.3*rng.normal(size=L)), 1e-12)
    c = dctI(Sh); ev = np.linalg.eigvalsh(toeplitz(c)).min()
    bad += ev < 0; worst = min(worst, ev/c[0])
print(f'(d) AR(1) 0.9 spectrum with 30 % per-bin noise: Toeplitz(C) indefinite in {bad}/5000 (worst min-eig/C0 {worst:.3f})')
