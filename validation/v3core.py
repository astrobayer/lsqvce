"""
v3 feature mirror. Validate BEFORE transcribing to MATLAB.

Core insight (measured): trace(Sxx * A_g'A_g) = sum_i ||L^-1 a_i||^2, exact,
via one blocked triangular solve. Gives per-group redundancy AND per-observation
leverage h_ii in the same pass. Hutchinson was a mistake (20-50% error on
realistically ill-conditioned normals).
"""
import numpy as np
from numpy.linalg import solve, norm, qr, lstsq
from scipy.linalg import cholesky, solve_triangular, toeplitz
from scipy.signal import lfilter

# ============================================================ leverage core
def joint_leverage(A, B, Lc, Ll, G, SG, GSG):
    """
    h_ii of the JOINT system [x_c ; y_s] for one arc, WITHOUT forming inv(N).

    Block-arrow inverse:
        inv(N)_cc   = Sxx
        inv(N)_ll,s = inv(Nll) + G Sxx G'
        inv(N)_lc,s = -G Sxx
    so
        h_i = a'Sxx a  +  b'(inv(Nll) + G Sxx G')b  -  2 a'Sxx G' b

    Lc  = chol(Ncc_red, lower)      SG  = Sxx @ G'      (nc x nl)
    Ll  = chol(Nll,     lower)      GSG = G @ Sxx @ G'  (nl x nl)
    """
    Ya = solve_triangular(Lc, A.T, lower=True)       # nc x m
    t1 = np.sum(Ya**2, axis=0)                       # a' Sxx a

    Yb = solve_triangular(Ll, B.T, lower=True)       # nl x m
    t2 = np.sum(Yb**2, axis=0) + np.sum((B @ GSG) * B, axis=1)

    t3 = -2.0 * np.sum((A @ SG) * B, axis=1)         # -2 a' Sxx G' b
    return t1 + t2 + t3


def leverage_blocks(Ncc_red, Nll, Nlc):
    """Precompute the small matrices needed by joint_leverage, once per system."""
    Lc  = cholesky(Ncc_red, lower=True)
    Ll  = cholesky(Nll,     lower=True)
    G   = solve(Nll, Nlc)                            # nl x nc
    SG  = solve_triangular(Lc.T, solve_triangular(Lc, G.T, lower=True))  # Sxx G'
    GSG = G @ SG                                     # nl x nl
    return Lc, Ll, G, SG, GSG


# ============================================================ regularisation
class CommonPrior:
    """
    Basis-agnostic prior on the common parameters:  K x ~ x0,  weight 1/sig_K^2.

    K = I                       -> plain damping
    K = diag(1./sigma_l)        -> Kaula (you build sigma_l; toolbox never sees 'SH')
    K = graph Laplacian         -> mascon smoothness
    K = chol(inv(C_signal))     -> full signal-covariance prior
    K = [0 I 0]                 -> regularise a subset only
    """
    def __init__(self, K, x0=None, sigma=1.0):
        self.K = np.atleast_2d(K)
        self.nk, self.nc = self.K.shape
        self.x0 = np.zeros(self.nk) if x0 is None else np.asarray(x0).ravel()
        self.sigma = float(sigma)

    def normals(self):
        w = 1.0 / self.sigma
        Kw = self.K * w
        return Kw.T @ Kw, Kw.T @ (self.x0 * w)

    def redundancy(self, Lc):
        """
        r_K = nk - (1/sig^2) trace(Sxx K'K) = nk - sum_j ||L^-1 k_j||^2 / sig^2
        DoF absorbed by the prior = nk - r_K.
        """
        w2 = 1.0 / self.sigma**2
        Y = solve_triangular(Lc, self.K.T, lower=True)     # nc x nk
        absorbed = w2 * np.sum(Y**2)
        return self.nk - absorbed, absorbed

    def vtv(self, x):
        v = (self.K @ x - self.x0) / self.sigma
        return v @ v


# ============================================================ separable multi-axis
def separable_whiten(A, B, l, a_ar, s_ar, r_acv, sig_epoch, d, first):
    """
    Separable covariance:  C = D (R_t kron I_d) D

      1. de-scale each EPOCH by its (d x d) Cholesky factor  -> D^-1
      2. filter each AXIS independently down epochs           -> F_t kron I_d

    Rows are EPOCH-MAJOR: [x1 y1 z1 x2 y2 z2 ...],  m = nEpoch * d.

    sig_epoch : (nEpoch, d)      per-epoch, per-axis sigma      -> diagonal D_k
             or (nEpoch, d, d)   per-epoch full covariance      -> D_k = chol(C_k)
    a_ar, s_ar, r_acv : the SCALAR temporal AR model, shared by all axes
                        (or pass a list of d models for per-axis temporal behaviour)

    Cost: O(nEpoch*d^3) for the de-scale + O(m*p) for the filter. No penalty
    for multi-axis: the same as the scalar cost.
    """
    m = l.shape[0]
    nE = m // d
    assert m == nE * d, "rows must be a whole number of epochs"

    # ---- 1. de-scale per epoch --------------------------------------------
    def descale(X):
        if X.size == 0:
            return X
        Xr = X.reshape(nE, d, -1)                    # epoch, axis, col
        if sig_epoch.ndim == 2:                      # diagonal
            Xr = Xr / sig_epoch[:, :, None]
        else:                                        # full d x d per epoch
            for k in range(nE):
                Lk = cholesky(sig_epoch[k], lower=True)
                Xr[k] = solve_triangular(Lk, Xr[k], lower=True)
        return Xr.reshape(m, -1)

    Ad = descale(A.copy() if A.size else A)
    Bd = descale(B.copy() if B.size else B)
    ld = descale(l.reshape(m, 1).copy()).ravel()

    # ---- 2. filter each axis down epochs ----------------------------------
    p = len(a_ar) - 1

    def tfilt(X):
        if X.size == 0:
            return X
        nCol = X.shape[1] if X.ndim == 2 else 1
        Xr = X.reshape(nE, d * nCol)                 # epoch-major -> filter down axis 0
        Yr = lfilter(a_ar, [1.0], Xr, axis=0) / s_ar
        if p > 0 and first:
            nS = min(p, nE)
            Gm = toeplitz(r_acv[:nS])
            Wp = solve_triangular(cholesky(Gm, lower=True), np.eye(nS), lower=True)
            Yr[:nS] = Wp @ Xr[:nS]
        return Yr.reshape(nE * d, nCol) if X.ndim == 2 else Yr.ravel()

    Aw = tfilt(Ad) if A.size else A
    Bw = tfilt(Bd) if B.size else B
    lw = tfilt(ld.reshape(nE * d, 1)).ravel()
    return Aw, Bw, lw


# ============================================================ cross-correlation
def cross_corr(w, d, maxLag=20):
    """
    Empirical R(tau) of the whitened residuals, reshaped to (nEpoch, d).

    Under a CORRECT separable model, R(0) = I and R(tau>0) = 0.

      R(0) off-diagonal significant  -> cross-axis correlation at zero lag.
                                        Your per-epoch covariance is wrong (or
                                        you passed a diagonal sigma when you have
                                        the full 3x3). Fix D, not the AR model.
      R(tau>0) diagonal nonzero      -> AR order too low. Raise p.
      R(tau>0) OFF-diagonal nonzero  -> genuine cross-axis TEMPORAL correlation.
                                        Separability fails; you need a VAR.

    Significance band: +-1.96/sqrt(nEpoch) per element.
    """
    nE = len(w) // d
    W = w[:nE * d].reshape(nE, d)
    R = np.zeros((maxLag + 1, d, d))
    for tau in range(maxLag + 1):
        R[tau] = W[:nE - tau].T @ W[tau:] / (nE - tau)
    band = 1.96 / np.sqrt(nE)

    diag_ok  = np.abs(np.diag(R[0]) - 1).max() < 5 * band
    cross0   = np.abs(R[0] - np.diag(np.diag(R[0]))).max()
    diagT    = max(np.abs(np.diag(R[t])).max() for t in range(1, maxLag + 1))
    crossT   = max(np.abs(R[t] - np.diag(np.diag(R[t]))).max()
                   for t in range(1, maxLag + 1))

    verdict, advice = 'ok', 'separable model is adequate'
    if cross0 > 3 * band:
        verdict, advice = 'cross-axis-zero-lag', \
            'per-epoch covariance is wrong -- supply the full dxd, not a diagonal sigma'
    elif crossT > 3 * band:
        verdict, advice = 'needs-var', \
            'genuine cross-axis temporal correlation -- separability fails, use a VAR'
    elif diagT > 3 * band:
        verdict, advice = 'ar-order-low', 'raise the AR order'

    return dict(R=R, band=band, verdict=verdict, advice=advice,
                cross0=cross0, diagT=diagT, crossT=crossT)


# ============================================================ AR order selection
def select_ar_order(f, S, dt, grid=None, acvTol=1e-3, checkLags=512, nData=None):
    """
    Select the AR order on ACV FIDELITY, not on AIC/BIC.

    AIC/BIC target PREDICTION error; we want DECORRELATION fidelity. They are
    not the same, and the difference matters: AIC happily picks a low order that
    predicts well on average while leaving structure in the ACF tail -- which is
    exactly what corrupts the normal equations.
    """
    from math import log
    if grid is None:
        grid = [2, 5, 10, 20, 30, 40, 60, 80, 120, 160]
    if nData is None:
        nData = 10000

    r_ref = psd2acv(f, S, dt, checkLags)
    out = []
    best = None
    for p in grid:
        a, s2 = levinson(psd2acv(f, S, dt, p + 1), p)
        r_fit = s2 * arma_acv(a, [1.0], checkLags)
        e = norm(r_fit - r_ref) / norm(r_ref)
        rts = np.abs(np.roots(a)) if len(a) > 1 else np.array([0.0])
        rec = dict(p=len(a) - 1, acvRelErr=e,
                   aic=nData * log(s2) + 2 * p,
                   bic=nData * log(s2) + p * log(nData),
                   maxRoot=rts.max(), sigma=np.sqrt(s2))
        out.append(rec)
        if best is None or e < best['acvRelErr']:
            best = rec
        if e < acvTol:
            break

    # --- guards ---------------------------------------------------------
    warn = []
    if best['maxRoot'] > 0.999:
        warn.append('near-unit AR root: the PSD probably has an unmodelled trend '
                    'or a pole at f=0. Detrend before estimating it.')
    if len(out) >= 2:
        impr = out[-1]['acvRelErr'] / max(out[-2]['acvRelErr'], 1e-300)
        if impr > 0.9 and best['acvRelErr'] > acvTol:
            warn.append('improvement has plateaued: AR cannot represent this '
                        'spectrum (sharp lines?). Use ARMA or psdMode=toeplitz.')
    verdict = ('good' if best['acvRelErr'] < acvTol else
               'marginal' if best['acvRelErr'] < 10 * acvTol else 'bad')
    return dict(p=best['p'], acvRelErr=best['acvRelErr'], verdict=verdict,
                curve=out, warnings=warn, sigma=best['sigma'])


# ============================================================ primitives (v2)
def psd2acv(f, S, dt, n):
    f = np.asarray(f).ravel(); S = np.asarray(S).ravel()
    tau = np.arange(n)[:, None] * dt
    return np.trapezoid(S[None, :] * np.cos(2 * np.pi * tau * f[None, :]), f, axis=1)


def levinson(r, p):
    r = np.asarray(r).ravel(); a = np.array([1.0]); s2 = r[0]
    for k in range(1, p + 1):
        kap = -(r[k] + a[1:] @ r[k - 1:0:-1]) / s2
        a = np.concatenate([a, [0.]]) + kap * np.concatenate([[0.], a[::-1]])
        s2 *= (1 - kap**2)
        if s2 <= 0:
            a = a[:-1]; s2 /= (1 - kap**2); break
    return a, s2


def arma_acv(a, b, n):
    a = np.atleast_1d(a); b = np.atleast_1d(b)
    nh = max(4 * n, 512)
    imp = np.zeros(nh); imp[0] = 1
    h = lfilter(b, a, imp)
    R = np.fft.irfft(np.abs(np.fft.rfft(h, 2 * nh))**2, 2 * nh)
    return R[:n]
