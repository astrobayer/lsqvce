# adjustVCEjoint>oneArc on a 'var' arc: nAx = d > 1 -> Sig = a.Cfcn(1:mE) with mE = m/d EPOCHS,
# but covModel's var Cfcn = varDense(Av, Sg, d, mrows) interprets its argument as ROWS (block-Toeplitz
# over mrows/d epochs, axes interleaved). Mirror both for VAR(1), d = 2, A = 0.5 I, Sigma = I.
import numpy as np
d, a = 2, 0.5
def Gamma(k):  # VAR(1) diagonal: Gamma(k) = a^k/(1-a^2) I
    return (a**k/(1-a*a))*np.eye(d)
def varDense(mrows):
    if mrows % d: raise ValueError('lsqvce:covModel:notWholeEpochs (mrows=%d, d=%d)' % (mrows, d))
    nE = mrows//d; C = np.zeros((mrows, mrows))
    for i in range(nE):
        for j in range(nE):
            k = i-j; B = Gamma(k) if k >= 0 else Gamma(-k).T
            C[i*d:(i+1)*d, j*d:(j+1)*d] = B
    return C
for nE in (6, 7):
    m = nE*d; mE = m//d
    want = np.array([[a**abs(i-j)/(1-a*a) for j in range(mE)] for i in range(mE)])  # per-axis temporal cov over epochs
    try:
        got = varDense(mE)
        print('nE=%d: Sig used by oneArc (%dx%d) == per-axis temporal covariance? %s; max|diff| = %.3f; lag-1 entry used %.3f vs true %.3f'
              % (nE, mE, mE, np.allclose(got, want), np.abs(got-want).max(), got[0,1], want[0,1]))
    except ValueError as e:
        print('nE=%d: oneArc -> %s (an internal id, not a refusal naming var+noise)' % (nE, e))
