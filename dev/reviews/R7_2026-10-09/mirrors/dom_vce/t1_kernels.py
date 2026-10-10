import numpy as np, sys
from scipy.linalg import toeplitz, cholesky
from kern import *
rng = np.random.default_rng(1)
# ---- (1) redundancyLags Toeplitz path + FFT lag sums vs dense trace(R V_n), exact probes (UiZ = U^-1 full basis)
N, P, nl = 200, 7, 2
acv = 0.6**np.arange(N); acv[0] += 0.3
Sig = toeplitz(acv); W = cholesky(Sig)            # upper, Sig = W'W
A = rng.normal(size=(N, P)); B = rng.normal(size=(N, nl)); l = rng.normal(size=N)
Abar = np.linalg.solve(W.T, A); Bbar = np.linalg.solve(W.T, B); lbar = np.linalg.solve(W.T, l)
Nll = Bbar.T@Bbar; Nlc = Bbar.T@Abar; G = np.linalg.solve(Nll, Nlc)
Nred = Abar.T@Abar - Nlc.T@G
U = cholesky(Nred)
UiZ_exact = np.linalg.solve(U, np.eye(P))          # Zbar = I  -> exact trace
At = Abar - Bbar@G
Ll = cholesky(Nll, lower=True)
x = np.linalg.solve(Nred, At.T@lbar - 0) ; # reduced rhs: At' (lbar - Bbar Nll^-1 Bbar' lbar)
Pb = np.eye(N) - Bbar@np.linalg.solve(Nll, Bbar.T)
x = np.linalg.solve(Nred, Abar.T@Pb@lbar)
y = np.linalg.solve(Nll, Bbar.T@(lbar - Abar@x))
ebar = lbar - Abar@x - Bbar@y
s, om = redundancyLags(At, W, UiZ_exact, ebar, 30, Bbar, Ll)
Si = np.linalg.inv(Sig); AB = np.c_[A, B]
Rd = Si - Si@AB@np.linalg.solve(AB.T@Si@AB, AB.T)@Si
eh = W.T@ebar
def Vn(n):
    V = np.eye(N, k=n) + np.eye(N, k=-n) if n else np.eye(N); return V
sd = np.array([np.trace(Rd@Vn(n)) for n in range(30)])
od = np.array([eh@Si@Vn(n)@Si@eh for n in range(30)])
print('s_n kernel vs dense trace(R V_n): max rel err %.2e' % (np.max(np.abs(s-sd))/np.max(np.abs(sd))))
print('Omega_n kernel vs dense         : max rel err %.2e' % (np.max(np.abs(om-od))/np.max(np.abs(od))))
H = Si@AB@np.linalg.solve(AB.T@Si@AB, AB.T)@Si
print('trace of WHICH matrix: s_0 = %.4f, tr(R)=%.4f, tr(Sigma^-1 H-part)=%.4f' % (s[0], np.trace(Rd), np.trace(H)))
# unbiasedness of E[Omega_n] = s_n under the model (Monte Carlo)
M = 4000; acc = np.zeros(30)
for k in range(M):
    e = W.T@rng.normal(size=N)
    lb = np.linalg.solve(W.T, e)
    xx = np.linalg.solve(Nred, Abar.T@Pb@lb); yy = np.linalg.solve(Nll, Bbar.T@(lb - Abar@xx))
    _, o = redundancyLags(At, W, UiZ_exact, lb - Abar@xx - Bbar@yy, 30, Bbar, Ll)
    acc += o
acc /= M
print('E[Omega_n] (MC, %d) vs s_n: lags 0..5' % M, np.round(acc[:6], 3), np.round(s[:6], 3))
# ---- (2) Hutchinson with sharedProbe normalisation: unbiased, variance vs closed form
Z = 100; reps = 3000; est = []
for r in range(reps):
    raw = np.sign(rng.normal(size=(P, Z))); raw[raw == 0] = 1
    UiZ = np.linalg.solve(U, raw/np.sqrt(Z))
    s2, _ = redundancyLags(At, W, UiZ, None, 4, Bbar, Ll)
    est.append(s2)
est = np.array(est)
T = np.linalg.solve(U.T, (np.linalg.inv(W)@At).T@(np.linalg.inv(W)@At))@np.linalg.inv(U)   # u x u matrix traced for s_0
T = (T+T.T)/2
varCF = 2*(np.sum(T**2) - np.sum(np.diag(T)**2))/Z
print('Hutchinson s_0: mean %.5f exact %.5f; sd MC %.4f closed-form %.4f' % (est[:,0].mean(), s[0], est[:,0].std(), np.sqrt(varCF)))
