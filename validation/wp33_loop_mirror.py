"""wp33_loop_mirror.py — faithful Python mirror of the adjustVCEjoint
closed loop (declared covariance fed back, EXACT s_n instead of probes).
Written 2026-08-19 (WP-33 r2) as the F-15 evidence; EXTENDED 2026-08-19
(WP-40 r2) with the post-WP-35 floor update, the damping exponent, the
Robbins-Monro schedule ('rm2': gamma_t = 2/(2+it)) and Polyak tail
averaging. MEASURED GRIDS (8 seeds each, corr over bins 2:end):
  undamped floor:  min -0.21 (8x96) / 0.08 (16x96) / -0.66 (16x128) / 0.97 (16x256)
  rm2 decay:       min  0.87 / 0.90 / 0.89 / 0.91
  rm2 + tail4:     min  0.89 / 0.92 / 0.91 / 0.93
  tail4 alone:     does NOT rescue the undamped collapse
Also measured: the WP-35 floor does NOT cure F-15 (some seeds worsen);
the ZeroFreqRule and nLag-resolution suspects are exonerated (shifting,
not stabilising). Conclusion: stochastic-approximation instability of
the undamped multiplicative step -> the decay schedule is the cure,
shipped opt-in as adjustVCEjoint StepSchedule/AvgTail + psdVCE Gamma.
Run: python3 wp33_loop_mirror.py (numpy + scipy). Model Claude Fable 5."""
# Faithful Python mirror of the adjustVCEjoint closed loop on the coloured
# fixture (exact s_n via dense Rbar, no probes): find the honest config with
# corr > 0.95, and quantify the DC-rule effect.
import numpy as np
from scipy.linalg import toeplitz, cholesky, solve_toeplitz

def dctI_ortho(x):
    x = np.atleast_2d(x.T).T if x.ndim == 1 else x
    n = x.shape[0]
    k = np.ones(n); k[0] = k[-1] = 1/np.sqrt(2)
    j = np.arange(n)
    X = np.cos(np.pi*np.outer(j, j)/(n-1)) * np.outer(k, k) * np.sqrt(2/(n-1))
    return (X @ x)

def arma_acv(ar, nlag, nsim=20000, rng=None):
    rng = rng or np.random.default_rng(0)
    e = rng.normal(size=nsim+500)
    x = np.zeros(nsim+500)
    for i in range(len(ar)-1, nsim+500):
        x[i] = e[i] - sum(ar[k]*x[i-k] for k in range(1, len(ar)))
    x = x[500:]
    return np.array([x[:nsim-l] @ x[l:] / (nsim-l) for l in range(nlag)])

def theoretical_acv(ar, nlag):
    # exact ACV via long AR impulse response
    h = np.zeros(4000); h[0] = 1
    for i in range(1, 4000):
        h[i] = -sum(ar[k]*h[i-k] for k in range(1, min(len(ar), i+1)))
    return np.array([h[:4000-l] @ h[l:] for l in range(nlag)])

def run(nArc=8, m=96, nLag=10, iters=8, dc_rule=True, seed=1,
        ar=(1,-0.75,0.4), nc=3, nl=2, update='abs', gamma=1.0, schedule=None, avgTail=0):
    rng = np.random.default_rng(seed)
    acv_t = theoretical_acv(np.array(ar), m)
    SigT = toeplitz(acv_t)                      # true noise covariance
    Wt = cholesky(SigT)                         # upper
    A = [rng.normal(size=(m, nc)) for _ in range(nArc)]
    B = [rng.normal(size=(m, nl)) for _ in range(nArc)]
    xT, yT = rng.normal(size=nc), rng.normal(size=nl)
    L = [a@xT + b@yT + Wt.T@rng.normal(size=m) for a, b in zip(A, B)]

    C = np.r_[1.0, 0.3*np.ones(nLag-1)]        # the driver's default C0
    tail = []
    for it in range(iters):
        # declared covariance from C: upsampled psd route (as jointAccMapped r2c)
        Sh = np.maximum(dctI_ortho(C).ravel(), 0)
        fC = np.arange(nLag)/(nLag-1)
        fD = np.linspace(0, 1, 64)
        SD = np.maximum(np.interp(fD, fC, Sh), 1e-10)
        # psd -> acv on m lags (inverse DCT-I of dense spectrum, ortho pair)
        acv_d = dctI_ortho(SD).ravel()
        # map the 64-lag acv to m lags: dense-spectrum acv has 64 lags; extend 0
        acv_full = np.zeros(m); acv_full[:min(64, m)] = acv_d[:min(64, m)]
        Sig = toeplitz(acv_full)
        w, V = np.linalg.eigh((Sig+Sig.T)/2)
        Sig = V@np.diag(np.maximum(w, 1e-6*np.max(w)))@V.T
        Sig = (Sig+Sig.T)/2
        W = cholesky(Sig + 1e-10*np.eye(m))
        # joint solve over all arcs (whitened)
        Aw = [np.linalg.solve(W.T, a) for a in A]
        Bw = [np.linalg.solve(W.T, b) for b in B]
        lw = [np.linalg.solve(W.T, l) for l in L]
        Ncc = sum(a.T@a for a in Aw); Nll = sum(b.T@b for b in Bw)
        Nlc = sum(b.T@a for a, b in zip(Aw, Bw))
        ncv = sum(a.T@l for a, l in zip(Aw, lw))
        nlv = sum(b.T@l for b, l in zip(Bw, lw))
        G = np.linalg.solve(Nll, Nlc)
        Nred = Ncc - Nlc.T@G
        x = np.linalg.solve(Nred, ncv - Nlc.T@np.linalg.solve(Nll, nlv))
        y = np.linalg.solve(Nll, nlv - Nlc@x)
        Ll = cholesky((Nll+Nll.T)/2, lower=True)
        Sinv = np.linalg.inv(Sig)
        Xw = np.linalg.inv(W)                   # W^-1 (upper)
        om = np.zeros((nLag, nArc)); sn = np.zeros((nLag, nArc))
        for k in range(nArc):
            At = Aw[k] - Bw[k]@G
            e = lw[k] - Aw[k]@x - Bw[k]@y       # ebar (whitened residual)
            et = np.linalg.solve(W, e)          # etilde
            om[0, k] = et@et
            Rt = np.linalg.solve(W, At)         # W^-1 AtilBar (exact, no probes)
            K = np.linalg.solve(W, np.linalg.solve(Ll, Bw[k].T).T)
            Rb = Sinv - Rt@np.linalg.solve(Nred, Rt.T) - K@K.T
            sn[0, k] = np.trace(Rb)
            for n in range(1, nLag):
                om[n, k] = 2*(et[n:]@et[:-n])
                sn[n, k] = 2*np.sum(np.diag(Rb, -n))
        OmT = dctI_ortho(C[:, None]*om); sT = dctI_ortho(C[:, None]*sn)
        Omj = OmT.sum(1); sj = sT.sum(1)
        alpha2 = np.where(np.abs(sj) > 0, Omj/sj, 1.0)
        Sprior = dctI_ortho(C).ravel()
        if update == 'abs':                       # PRE-WP-35 (the F-15 grid)
            Shat = np.maximum(np.abs(alpha2*Sprior), 1e-6)
        else:                                     # POST-WP-35: FLOOR, never abs
            a2 = np.where(alpha2 > 0, alpha2, np.nan)
            a2 = np.nan_to_num(a2, nan=1.0)       # nonpositive: hold (floored)
            g = gamma
            if schedule == 'rm':   g = 1.0/(1+it)          # Robbins-Monro decay
            if schedule == 'rm2':  g = 2.0/(2+it)
            Shat = np.maximum((a2**g)*Sprior, 1e-6)
        if dc_rule: Shat[0] = Shat[1]
        C = dctI_ortho(Shat).ravel()
        if avgTail: tail.append(np.log(Shat))
    if avgTail and len(tail) >= avgTail:
        Shat = np.exp(np.mean(tail[-avgTail:], axis=0))    # Polyak tail average
    St = dctI_ortho(theoretical_acv(np.array(ar), nLag)).ravel()
    Stc = St.copy()
    if dc_rule: Stc[0] = Stc[1]
    def cc(a, b): return np.corrcoef(a, b)[0, 1]
    return cc(Shat, St), cc(Shat, Stc), cc(Shat[1:], St[1:])

if __name__ == "__main__":
    import numpy as _np
    for cfg in [dict(nArc=8, m=96, iters=12), dict(nArc=16, m=96, iters=12),
                dict(nArc=16, m=128, iters=12), dict(nArc=16, m=256, iters=12)]:
        raw = [run(seed=s, update='floor', **cfg)[2] for s in range(1, 9)]
        stb = [run(seed=s, update='floor', schedule='rm2', avgTail=4, **cfg)[2]
               for s in range(1, 9)]
        print(f"{cfg['nArc']}x{cfg['m']:3d}  undamped min {min(raw):+.3f}   "
              f"rm2+tail4 min {min(stb):+.3f}")
