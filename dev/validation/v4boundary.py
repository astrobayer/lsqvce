"""
Boundary-bias measurement for the DENSE covariance path (ROADMAP 2.2).

The filter path has ZERO batch-boundary error once overlap >= p. The dense
path truncates the correlation at the batch edge: per-batch chol(C_bb) treats
each batch as independent, discarding the cross-batch covariance. The toolbox
says "choose BatchSize >> correlation length" -- this script measures the
number behind that advice.

EXPERIMENT
  AR(1) noise with phi = exp(-1/L)  ->  correlation length L samples (1/e).
  Whiten in batches of size b with NO overlap, per-batch dense Cholesky of the
  TRUE covariance (so the only error is the boundary truncation).
  Reference: exact whole-arc whitening.
  Report: relative error in x, and the sigma0 distortion, versus b/L.

The parameter error and the sigma0 distortion are what a user actually pays.
"""
import numpy as np
from numpy.linalg import solve, norm
from scipy.linalg import cholesky, solve_triangular, toeplitz
from scipy.signal import lfilter

rng = np.random.default_rng(2026)

def run(L, m=4096, nc=6, nSeed=32):
    phi = np.exp(-1.0 / L)
    s_inn = np.sqrt(1 - phi**2)                      # unit process variance
    lag = np.arange(m)
    rows = []
    for b in [max(2, int(round(f * L))) for f in (0.5, 1, 2, 4, 8, 16, 32)]:
        if b > m:
            continue
        dx_form, es0, err_ref, err_bat = [], [], [], []
        for sd in range(nSeed):
            r = np.random.default_rng(1000 * sd + L)
            A = r.standard_normal((m, nc))
            xt = r.standard_normal(nc)
            e = lfilter([1.0], [1.0, -phi], s_inn * r.standard_normal(m + 200))[200:]
            l = A @ xt + e

            # exact reference: AR(1) FIR whitener over the whole arc IS chol^-1
            Aw = lfilter([1.0, -phi], [1.0], A, axis=0) / s_inn
            lw = lfilter([1.0, -phi], [1.0], l) / s_inn
            Nr = Aw.T @ Aw
            xr = solve(Nr, Aw.T @ lw)
            wr = lw - Aw @ xr
            s0r = np.sqrt(wr @ wr / (m - nc))
            sig_form = np.sqrt(np.trace(solve(Nr, np.eye(nc))))  # formal |dx|

            # batched dense, NO overlap: per-batch chol of the true covariance
            N = np.zeros((nc, nc)); n = np.zeros(nc); vtv = 0.0
            for i0 in range(0, m, b):
                i1 = min(i0 + b, m)
                k = i1 - i0
                C = toeplitz(phi ** lag[:k])
                Lc = cholesky(C, lower=True)
                Ab = solve_triangular(Lc, A[i0:i1], lower=True)
                lb = solve_triangular(Lc, l[i0:i1], lower=True)
                N += Ab.T @ Ab; n += Ab.T @ lb; vtv += lb @ lb
            xb = solve(N, n)
            s0b = np.sqrt(max(vtv - n @ xb, 0) / (m - nc))

            dx_form.append(norm(xb - xr) / sig_form)  # vs the FORMAL error
            es0.append(abs(s0b / s0r - 1))
            err_ref.append(norm(xr - xt))             # true-error, exact GLS
            err_bat.append(norm(xb - xt))             # true-error, batched
        # efficiency loss: RMS true error of the batched estimator vs exact GLS
        effloss = (np.sqrt(np.mean(np.square(err_bat))) /
                   np.sqrt(np.mean(np.square(err_ref))) - 1)
        rows.append((b, b / L, float(np.median(dx_form)), float(np.median(es0)),
                     float(effloss)))
    return rows

if __name__ == '__main__':
    print('Boundary bias, dense path, batches WITHOUT overlap')
    print('AR(1), correlation length L (1/e), m = 4096, nc = 6, 32 realisations\n')
    print('  |dx|/formal : boundary-induced parameter shift vs the FORMAL error')
    print('  eff. loss   : RMS true-error inflation of the batched estimator\n')
    for L in (16, 64):
        print(f'  L = {L} samples')
        print(f'    {"b":>6s} {"b/L":>6s} {"|dx|/formal":>12s} '
              f'{"sigma0 dist.":>13s} {"eff. loss":>10s}')
        for b, bl, exf, es, ef in run(L):
            print(f'    {b:6d} {bl:6.1f} {exf:12.4f} {es:13.3e} {ef:10.2%}')
        print()
