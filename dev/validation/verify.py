"""Mirror of the MATLAB toolbox semantics; runs the same assertions as
tests/test_lsqAccum.m so the tolerances and logic are validated."""
import numpy as np
from numpy.linalg import norm, eig, lstsq, cond
from scipy.linalg import toeplitz, cholesky, solve_triangular, qr, hankel
from scipy.signal import lfilter

# ============================================================ primitives
def psd2acv(f, S, dt, n):
    f = np.asarray(f).ravel(); S = np.asarray(S).ravel()
    tau = np.arange(n)[:, None]*dt
    return np.trapezoid(S[None, :]*np.cos(2*np.pi*tau*f[None, :]), f, axis=1)

def levinson(r, p):
    r = np.asarray(r).ravel(); a = np.array([1.0]); s2 = r[0]
    assert s2 > 0
    for k in range(1, p+1):
        kap = -(r[k] + a[1:] @ r[k-1:0:-1]) / s2
        a = np.concatenate([a, [0.]]) + kap*np.concatenate([[0.], a[::-1]])
        s2 *= (1-kap**2)
        if s2 <= 0:
            a = a[:-1]; s2 /= (1-kap**2); break
    return a, s2

def arma_acv(a, b, n):
    a = np.atleast_1d(a); b = np.atleast_1d(b)
    nh = max(4*n, 512)
    imp = np.zeros(nh); imp[0] = 1
    h = lfilter(b, a, imp)
    R = np.fft.irfft(np.abs(np.fft.rfft(h, 2*nh))**2, 2*nh)
    return R[:n]

def decorr(C, tol=1e-10):
    from scipy.linalg import eigh
    C = (C+C.T)/2
    try:
        L = cholesky(C, lower=True)
        return solve_triangular(L, np.eye(C.shape[0]), lower=True)
    except Exception:
        # symmetric eig -> V orthonormal -> W = diag(1/sqrt(d)) @ V.T
        d, V = eigh(C)
        dmax = d.max()
        if dmax <= 0:
            raise ValueError('covariance has no positive eigenvalue')
        d = np.maximum(d, tol*dmax)
        return np.diag(1/np.sqrt(d)) @ V.T

def arma_whiten(a, b, s, r, A, B, l, first):
    a = np.atleast_1d(a); b = np.atleast_1d(b)
    p = len(a)-1
    f = lambda X: lfilter(a, b, X, axis=0)/s if X.size else X
    Aw, Bw, lw = f(A), f(B), f(l)
    if p > 0 and first:
        n = l.shape[0]
        nS = min(p, n)
        G = toeplitz(r[:nS])
        Wp = solve_triangular(cholesky(G, lower=True), np.eye(nS), lower=True)
        if A.size: Aw[:nS] = Wp @ A[:nS]
        if B.size: Bw[:nS] = Wp @ B[:nS]
        lw[:nS] = Wp @ np.asarray(l)[:nS]
    return Aw, Bw, lw

def cov_factory(cm, t, opt):
    m = len(t); dt = np.median(np.diff(t)) if m > 1 else 1
    Cf = Wf = None; fm = 0
    ty = cm['type']
    if ty == 'diag':
        s = cm['sigma']; s = np.full(m, s) if np.isscalar(s) else s
        Cf = lambda idx, o: np.diag(s[idx]**2)
        def Wf(idx, A, B, l, o):
            sc = s[idx]
            return (A/sc[:, None] if A.size else A,
                    B/sc[:, None] if B.size else B,
                    np.asarray(l)/(sc if np.asarray(l).ndim == 1 else sc[:, None]))
    elif ty == 'arma':
        a = np.atleast_1d(cm['ar']); b = np.atleast_1d(cm.get('ma', [1.]))
        s = cm['sigma']; r = s**2*arma_acv(a, b, max(len(a), 2)); fm = len(a)-1
        Wf = lambda idx, A, B, l, o: arma_whiten(a, b, s, r, A, B, l, idx[0] == 0)
        Cf = lambda idx, o: s**2*toeplitz(arma_acv(a, b, len(idx)))
    elif ty == 'psd':
        md = opt.get('psdMode', 'ar'); p = opt.get('arOrder', 30)
        f = np.asarray(cm['f']); S = np.asarray(cm['S'])
        S = np.maximum(S, 1e-10*S.max())
        if md == 'toeplitz':
            nmax = min(opt.get('BatchSize', 20000)+opt.get('BatchOverlap', 0), m)
            r = psd2acv(f, S, dt, nmax)
            Cf = lambda idx, o: toeplitz(r[:len(idx)])
        else:
            a, s2 = levinson(psd2acv(f, S, dt, p+1), p)
            s = np.sqrt(s2); rr = psd2acv(f, S, dt, max(len(a), 2)); fm = len(a)-1
            Wf = lambda idx, A, B, l, o: arma_whiten(a, [1.], s, rr, A, B, l, idx[0] == 0)
            Cf = lambda idx, o: s2*toeplitz(arma_acv(a, [1.], len(idx)))
    return Cf, Wf, fm

# ============================================================ backends
class NormalEq:
    def __init__(s, d): s.d = d; s.N = np.zeros((d, d)); s.n = np.zeros(d); s.ltl = 0.; s.m = 0
    def add(s, A, l): s.N += A.T@A; s.n += A.T@l; s.ltl += l@l; s.m += l.size
    def blocks(s, ic, il):
        return (s.N[np.ix_(il, il)], s.N[np.ix_(il, ic)], s.N[np.ix_(ic, ic)],
                s.n[il], s.n[ic], s.ltl, s.m)
    def merge(s, R): s.N += R[0]; s.n += R[1]; s.ltl += R[2]; s.m += R[3]
    def solve(s):
        R = cholesky(s.N)
        x = solve_triangular(R, solve_triangular(R.T, s.n, lower=True))
        vtv = s.ltl - s.n@x; dof = s.m - s.d
        Ri = solve_triangular(R, np.eye(s.d))
        return x, (vtv/dof)*(Ri@Ri.T), vtv/dof, dof, vtv

class QRSys:
    def __init__(s, d): s.d = d; s.R = np.zeros((0, d+1)); s.m = 0
    def add(s, A, l):
        M = np.vstack([s.R, np.hstack([A, l[:, None]])])
        s.R = np.triu(qr(M, mode='r')[0])[:min(M.shape[0], s.d+1)]
        s.m += l.size
    def blocks(s, ic, il):
        d = s.d
        R1 = np.zeros((d, d)); rr = np.zeros(d)
        k = min(s.R.shape[0], d)
        R1[:k] = s.R[:k, :d]; rr[:k] = s.R[:k, -1]
        N = R1.T@R1; nn = R1.T@rr
        ltl = rr@rr
        if s.R.shape[0] > d: ltl += s.R[d, -1]**2
        return (N[np.ix_(il, il)], N[np.ix_(il, ic)], N[np.ix_(ic, ic)],
                nn[il], nn[ic], ltl, s.m)
    def merge(s, R):
        N, n, ltl, dof = R
        try: Rc = cholesky(N)
        except Exception:
            w, V = eig((N+N.T)/2); Rc = np.diag(np.sqrt(np.maximum(np.real(w), 0)))@V.T
        rhs = solve_triangular(Rc.T, n, lower=True)
        # ltl = ||rhs||^2 + rho^2 : the residual mass rho must be carried as an
        # extra pseudo-observation row, otherwise vtv (and sigma0) collapse.
        rho2 = ltl - rhs @ rhs
        blk = np.hstack([Rc, rhs[:, None]])
        if rho2 > 0:
            res_row = np.zeros((1, s.d+1)); res_row[0, -1] = np.sqrt(rho2)
            blk = np.vstack([blk, res_row])
        M = np.vstack([s.R, blk])
        s.R = np.triu(qr(M, mode='r')[0])[:min(M.shape[0], s.d+1)]
        s.m += dof
    def solve(s):
        d = s.d; R1 = s.R[:d, :d]; rr = s.R[:d, -1]
        x = solve_triangular(R1, rr)
        vtv = s.R[d, -1]**2 if s.R.shape[0] > d else 0.
        dof = s.m - d
        Ri = solve_triangular(R1, np.eye(d))
        return x, (vtv/dof)*(Ri@Ri.T), vtv/dof, dof, vtv

# ============================================================ regularisation
def diagnostics(Nll, opt):
    w, V = eig((Nll+Nll.T)/2); w = np.sort(np.real(w))[::-1]
    d = dict(eig=w, condNll=w[0]/max(w[-1], 1e-300),
             rankDef=int(np.sum(w < opt['TruncTol']*w[0])), spd=bool(np.all(w > 0)))
    d['status'] = 'ok'
    if d['condNll'] > opt['CondWarn']: d['status'] = 'ill-conditioned'
    if d['condNll'] > opt['CondAbort'] or not d['spd']: d['status'] = 'rank-deficient'
    d['regUsed'] = 'none'; d['lambdaUsed'] = 0.
    return d

def gcv_lambda(N, n, L):
    from scipy.linalg import eigh
    w, V = eigh((N+N.T)/2, L.T@L + 1e-16*np.eye(N.shape[0]))
    b = V.T@n; dm = w.max()
    if dm <= 0: return 0.
    grid = np.logspace(np.log10(dm*1e-12), np.log10(dm), 80)
    g = [norm(gk/(w+gk)*b)**2/max(np.sum(gk/(w+gk))**2, 1e-300) for gk in grid]
    return grid[int(np.argmin(g))]

def regularise(Nll, nl, dg, opt):
    p = Nll.shape[0]; L = np.eye(p); reg = opt['Reg']
    if reg == 'auto':
        reg = {'ok': 'none', 'ill-conditioned': 'tikhonov'}.get(dg['status'], 'tsvd')
    Nr = Nll.copy()
    if reg == 'tikhonov':
        lam = opt['RegLambda'] or gcv_lambda(Nll, nl, L)
        Nr = Nll + lam*(L.T@L); dg['regUsed'] = 'tikhonov'; dg['lambdaUsed'] = lam
    elif reg == 'lm':
        lam = opt['RegLambda'] or 1e-8
        Nr = Nll + lam*np.diag(np.diag(Nll)); dg['regUsed'] = 'lm'; dg['lambdaUsed'] = lam
    elif reg == 'tsvd':
        w, V = eig((Nll+Nll.T)/2); w = np.real(w); V = np.real(V)
        fl = opt['TruncTol']*w.max(); w = np.maximum(w, fl)
        Nr = V@np.diag(w)@V.T; dg['regUsed'] = 'tsvd'; dg['lambdaUsed'] = fl
    return (Nr+Nr.T)/2, nl, dg

# ============================================================ accumulator
class Accumulator:
    def __init__(s, nc, backend='qr', BatchSize=20000, BatchOverlap=None,
                 Reg='none', RegLambda=0., TruncTol=1e-10, CondWarn=1e8,
                 CondAbort=1e14, UserOpt=None):
        s.nc = nc; s.backend = backend
        s.opt = dict(Reg=Reg, RegLambda=RegLambda, TruncTol=TruncTol,
                     CondWarn=CondWarn, CondAbort=CondAbort,
                     BatchSize=BatchSize, BatchOverlap=BatchOverlap)
        s.uo0 = dict(UserOpt or {})
        s.sys = (NormalEq if backend == 'normal' else QRSys)(nc)
        s.diag = []; s.store = {}

    def _new(s, n): return (NormalEq if s.backend == 'normal' else QRSys)(n)

    def addSystem(s, sid, nl, dataFcn, days):
        loc = s._new(s.nc+nl)
        uo = dict(s.uo0); uo.update(sysID=sid, nc=s.nc, nl=nl,
                                    BatchSize=s.opt['BatchSize'])
        if s.opt['BatchOverlap'] is not None:
            uo['BatchOverlap'] = s.opt['BatchOverlap']
        for d in days:
            meta = dataFcn('meta', d, None, uo)
            for ot in meta['obsTypes']:
                uo['obsType'] = ot
                for a in meta['arcs'][ot]:
                    uo['day'] = d; uo['arc'] = a
                    arc = dataFcn(d, ot, a, uo)
                    if arc is None or arc['l'].size == 0: continue
                    s._accum(loc, arc, uo)
        red, dg, keep = s._preelim(loc, nl, sid)
        s.sys.merge(red); dg['sysID'] = sid; s.diag.append(dg)
        keep.update(nl=nl, dataFcn=dataFcn, days=days, uo=uo)
        s.store[sid] = keep

    def _accum(s, acc, arc, uo):
        m = arc['A'].shape[0]
        bs = uo.get('BatchSize', 20000)
        nov = uo.get('BatchOverlap')
        if nov is None: nov = arc.get('filterMemory', 0)
        for i0 in range(0, m, bs):
            i1 = min(i0+bs, m); j0 = max(0, i0-nov)
            idx = np.arange(j0, i1); keep = slice(i0-j0, i1-j0)
            A, B, l = arc['A'][idx], arc['B'][idx], arc['l'][idx]
            if arc.get('Wfcn'):
                A, B, l = arc['Wfcn'](idx, A, B, l, uo)
            elif arc.get('Cfcn'):
                W = decorr(arc['Cfcn'](idx, uo)); A, B, l = W@A, W@B, W@l
            acc.add(np.hstack([A[keep], B[keep]]), l[keep])

    def _preelim(s, loc, nl, sid):
        nc = s.nc; ic = np.arange(nc); il = np.arange(nc, nc+nl)
        Nll, Nlc, Ncc, nlr, ncr, ltl, m = loc.blocks(ic, il)
        dg = diagnostics(Nll, s.opt)
        Nll_r, nl_r, dg = regularise(Nll, nlr, dg, s.opt)
        RHS = np.column_stack([Nlc, nl_r])
        try:
            # Only intercept true singularity; ill-conditioning is what the QR
            # backend exists to handle, so do not divert it to pinv.
            if not np.all(np.isfinite(Nll_r)) or np.linalg.cond(Nll_r) > 1/np.finfo(float).eps:
                raise np.linalg.LinAlgError('numerically singular')
            Z = np.linalg.solve(Nll_r, RHS)
        except np.linalg.LinAlgError:
            # Never crash: fall back to the pseudo-inverse and tell the user.
            Z = np.linalg.pinv(Nll_r, rcond=s.opt['TruncTol']) @ RHS
            dg['regUsed'] = 'pinv-fallback'
            dg['lambdaUsed'] = s.opt['TruncTol']
        Z = np.asarray(Z)
        Ncc_r = Ncc - Nlc.T@Z[:, :nc]
        ncc_r = ncr - Nlc.T@Z[:, -1]
        ltl_r = ltl - nl_r@Z[:, -1]
        Ncc_r = (Ncc_r+Ncc_r.T)/2
        return ((Ncc_r, ncc_r, ltl_r, m-nl), dg,
                dict(Nll_r=Nll_r, Nlc=Nlc, nl_r=nl_r, mObs=m))

    def solve(s):
        x, Cxx, s0sq, dof, vtv = s.sys.solve()
        return dict(x=x, Cxx=Cxx, s0sq=s0sq, dof=dof, vtv=vtv,
                    sigma0=np.sqrt(s0sq), diag=s.diag)

    def backSubstitute(s, sid, x):
        st = s.store[sid]
        return np.linalg.solve(st['Nll_r'], st['nl_r'] - st['Nlc']@x)

    def localCovariance(s, sid, Cxx, s0sq):
        st = s.store[sid]
        G = np.linalg.solve(st['Nll_r'], st['Nlc'])
        C = s0sq*np.linalg.inv(st['Nll_r']) + G@Cxx@G.T
        return (C+C.T)/2

    def residuals(s, sol):
        out = {'sys': [], 'vtv': 0., 'dof': 0}
        npar = len(sol['x']); mtot = 0; vtot = 0.
        for sid, st in s.store.items():
            y = s.backSubstitute(sid, sol['x'])
            S = dict(sysID=sid, y=y, vtv=0., m=0, w=[])
            uo = st['uo']
            for d in st['days']:
                meta = st['dataFcn']('meta', d, None, uo)
                for ot in meta['obsTypes']:
                    uo['obsType'] = ot
                    for a in meta['arcs'][ot]:
                        uo['day'] = d; uo['arc'] = a
                        arc = st['dataFcn'](d, ot, a, uo)
                        if arc is None or arc['l'].size == 0: continue
                        v = arc['l'] - arc['A']@sol['x'] - arc['B']@y
                        w = s._whiten_resid(arc, v, uo)
                        S['vtv'] += w@w; S['m'] += v.size; S['w'].append(w)
            S['w'] = np.concatenate(S['w']) if S['w'] else np.array([])
            S['rms'] = np.sqrt(S['vtv']/max(S['m'], 1))
            vtot += S['vtv']; mtot += S['m']; npar += st['nl']
            out['sys'].append(S)
        out['vtv'] = vtot; out['dof'] = mtot-npar
        out['s0sq'] = vtot/max(out['dof'], 1)
        return out

    def _whiten_resid(s, arc, v, uo):
        m = v.size
        bs = uo.get('BatchSize', 20000)
        nov = uo.get('BatchOverlap')
        if nov is None: nov = arc.get('filterMemory', 0)
        w = np.zeros(m)
        for i0 in range(0, m, bs):
            i1 = min(i0+bs, m); j0 = max(0, i0-nov)
            idx = np.arange(j0, i1); keep = slice(i0-j0, i1-j0)
            z = np.zeros((len(idx), 0))
            if arc.get('Wfcn'):
                _, _, wb = arc['Wfcn'](idx, z, z, v[idx], uo)
            elif arc.get('Cfcn'):
                wb = decorr(arc['Cfcn'](idx, uo))@v[idx]
            else:
                wb = v[idx]
            w[i0:i1] = wb[keep]
        return w

# ============================================================ synth loader
_CACHE = {}
def synthLoader(cfg, sat, day, obsType, arcID, opt):
    if day == 'meta':
        return {'obsTypes': cfg['obsTypes'],
                'arcs': {ot: list(range(1, cfg['arcsPerDay']+1)) for ot in cfg['obsTypes']}}
    key = (sat, day, obsType, cfg['tag'])
    if key not in _CACHE:
        _CACHE.clear()
        _CACHE[key] = _buildDay(cfg, sat, day, obsType)
    D = _CACHE[key]
    r = D['arcIdx'][arcID-1]
    out = dict(A=D['A'][r], B=D['B'][r], l=D['l'][r], t=D['t'][r])
    Cf, Wf, fm = cov_factory(cfg['cov'][obsType], out['t'], opt)
    out.update(Cfcn=Cf, Wfcn=Wf, filterMemory=fm)
    return out

def _buildDay(cfg, sat, day, obsType):
    M = cfg['arcsPerDay']*cfg['mArc']
    rs = np.random.default_rng(abs(hash((cfg['seed'], sat, day, obsType))) % (2**31))
    t = np.arange(M)*cfg['dt']
    A = rs.standard_normal((M, cfg['nc']))
    B = rs.standard_normal((M, cfg['nl']))
    if obsType == 'kbr': A[:, :3] *= 0.1
    if cfg.get('collinear'):
        B[:, -1] = B[:, -2] + cfg['collinearEps']*rs.standard_normal(M)
    e = _noise(rs, cfg['cov'][obsType], M, cfg['dt'])
    l = A@cfg['x_c'] + B@cfg['y'][sat] + e
    return dict(t=t, A=A, B=B, l=l,
                arcIdx=[np.arange(a*cfg['mArc'], (a+1)*cfg['mArc'])
                        for a in range(cfg['arcsPerDay'])])

def _noise(rs, cm, M, dt):
    nb = 2000
    if cm['type'] == 'diag':
        return cm['sigma']*rs.standard_normal(M)
    if cm['type'] == 'arma':
        w = cm['sigma']*rs.standard_normal(M+nb)
        return lfilter(np.atleast_1d(cm.get('ma', [1.])), np.atleast_1d(cm['ar']), w)[nb:]
    if cm['type'] == 'psd':
        a, s2 = levinson(psd2acv(cm['f'], cm['S'], dt, 201), 200)
        w = np.sqrt(s2)*rs.standard_normal(M+nb)
        return lfilter([1.], a, w)[nb:]

# ============================================================ test harness
def baseCfg(**kw):
    rs = np.random.default_rng(42)
    cfg = dict(nc=6, nl=4, sats=[1, 2, 3], obsTypes=['pos', 'kbr'],
               arcsPerDay=2, mArc=400, dt=5, seed=42, tag='base',
               collinear=False, collinearEps=1e-7)
    cfg['x_c'] = rs.standard_normal(6)
    cfg['y'] = {s: rs.standard_normal(4) for s in cfg['sats']}
    cfg['cov'] = {'pos': dict(type='arma', ar=[1, -1.6, .9, -.25], ma=[1.], sigma=.7),
                  'kbr': dict(type='diag', sigma=.2)}
    cfg.update(kw); cfg['tag'] = str(sorted(kw.keys()))+str(cfg['collinearEps'])
    return cfg

def runIt(cfg, days, **kw):
    acc = Accumulator(cfg['nc'], **kw)
    for s in cfg['sats']:
        acc.addSystem(s, cfg['nl'], lambda d, o, a, u, s=s: synthLoader(cfg, s, d, o, a, u), days)
    return acc, acc.solve()

def assembleFull(cfg, days):
    ns = len(cfg['sats']); Af = []; lf = []
    uo = dict(BatchSize=10**9)
    for k, s in enumerate(cfg['sats']):
        for d in days:
            meta = synthLoader(cfg, s, 'meta', None, None, uo)
            for ot in meta['obsTypes']:
                for a in meta['arcs'][ot]:
                    arc = synthLoader(cfg, s, d, ot, a, uo)
                    m = arc['A'].shape[0]; idx = np.arange(m)
                    if arc['Wfcn']:
                        Aw, Bw, lw = arc['Wfcn'](idx, arc['A'], arc['B'], arc['l'], uo)
                    else:
                        W = decorr(arc['Cfcn'](idx, uo))
                        Aw, Bw, lw = W@arc['A'], W@arc['B'], W@arc['l']
                    blk = np.zeros((m, cfg['nc']+ns*cfg['nl']))
                    blk[:, :cfg['nc']] = Aw
                    blk[:, cfg['nc']+k*cfg['nl']:cfg['nc']+(k+1)*cfg['nl']] = Bw
                    Af.append(blk); lf.append(lw)
    return np.vstack(Af), np.concatenate(lf)

# ---------------------------------------------------------------- tests
R = []
def T(name, fn):
    try:
        fn(); R.append((name, 'PASS', ''))
    except AssertionError as e:
        R.append((name, 'FAIL', str(e)[:70]))
    except Exception as e:
        R.append((name, 'ERROR', f'{type(e).__name__}: {e}'[:70]))

def t_psd2acv():
    a = [1, -1.6, .9, -.25]; s = .7; dt = 5
    r_ana = s**2*arma_acv(a, [1.], 64)
    f = np.linspace(1e-9, .5/dt, 8192)
    H = 1/(np.exp(-1j*2*np.pi*np.outer(f, np.arange(4))*dt)@np.array(a))
    S = s**2*dt*2*np.abs(H)**2
    e = norm(psd2acv(f, S, dt, 64)-r_ana)/norm(r_ana)
    assert e < 2e-3, f'acv err {e:.2e}'

def t_levinson():
    a = np.array([1, -1.6, .9, -.25]); s = .7
    ah, s2 = levinson(s**2*arma_acv(a, [1.], 8), 3)
    assert norm(ah-a) < 1e-8, f'a err {norm(ah-a):.2e}'
    assert abs(np.sqrt(s2)-s)/s < 1e-8

def t_filter_eq_chol():
    a = [1, -1.6, .9, -.25]; s = .7; n = 300
    r = s**2*arma_acv(a, [1.], n)
    X = np.random.default_rng(3).standard_normal((n, 5))
    Y1 = decorr(toeplitz(r))@X
    z = np.zeros((n, 0))
    _, _, Y2 = arma_whiten(a, [1.], s, r, z, z, X, True)
    e = norm(Y2.T@Y2 - Y1.T@Y1)/norm(Y1.T@Y1)
    assert e < 1e-9, f'YtY err {e:.2e}'

def t_arma_acv():
    r = arma_acv([1, -.9], [1, .5], 5)
    phi, th = .9, .5
    r0 = (1+2*phi*th+th**2)/(1-phi**2); r1 = phi*r0+th
    assert abs(r[0]-r0)/r0 < 1e-6 and abs(r[1]-r1)/r1 < 1e-6

def t_recover_white():
    cfg = baseCfg()
    cfg['cov'] = {'pos': dict(type='diag', sigma=.05), 'kbr': dict(type='diag', sigma=.05)}
    cfg['tag'] = 'white'
    _, sol = runIt(cfg, [1, 2, 3], backend='qr')
    e = np.max(np.abs(sol['x']-cfg['x_c']))
    assert e < 5e-3, f'max err {e:.2e}'

def t_sigma0():
    cfg = baseCfg()
    _, sol = runIt(cfg, [1, 2, 3, 4, 5], backend='qr')
    assert abs(sol['sigma0']-1) < .06, f"sigma0={sol['sigma0']:.4f}"

def t_backends_agree():
    cfg = baseCfg()
    _, s1 = runIt(cfg, [1, 2, 3], backend='normal')
    _, s2 = runIt(cfg, [1, 2, 3], backend='qr')
    e = norm(s2['x']-s1['x'])/norm(s1['x'])
    assert e < 1e-7, f'{e:.2e}'

def t_preelim_exact():
    cfg = baseCfg(); days = [1, 2]
    _, sol = runIt(cfg, days, backend='qr')
    Af, lf = assembleFull(cfg, days)
    xf = lstsq(Af, lf, rcond=None)[0][:cfg['nc']]
    e = norm(sol['x']-xf)/norm(xf)
    assert e < 1e-8, f'{e:.2e}'

def t_backsub():
    cfg = baseCfg(); days = [1, 2]
    acc, sol = runIt(cfg, days, backend='qr')
    Af, lf = assembleFull(cfg, days)
    xf = lstsq(Af, lf, rcond=None)[0]
    for k, s in enumerate(cfg['sats']):
        y = acc.backSubstitute(s, sol['x'])
        yr = xf[cfg['nc']+k*cfg['nl']:cfg['nc']+(k+1)*cfg['nl']]
        e = norm(y-yr)/norm(yr)
        assert e < 1e-7, f'sat{s} {e:.2e}'

def t_overlap():
    cfg = baseCfg(); days = [1, 2]
    _, ref = runIt(cfg, days, backend='qr', BatchSize=10**9)
    _, ok = runIt(cfg, days, backend='qr', BatchSize=100, BatchOverlap=3)
    _, bad = runIt(cfg, days, backend='qr', BatchSize=100, BatchOverlap=0)
    eo = norm(ok['x']-ref['x'])/norm(ref['x'])
    eb = norm(bad['x']-ref['x'])/norm(ref['x'])
    assert eo < 1e-6, f'overlap-ok err {eo:.2e}'
    assert eb > 10*eo, f'no-overlap not worse: {eb:.2e} vs {eo:.2e}'
    globals()['_ov'] = (eo, eb)

def t_auto_overlap():
    cfg = baseCfg()
    _, ref = runIt(cfg, [1, 2], backend='qr', BatchSize=10**9)
    _, aut = runIt(cfg, [1, 2], backend='qr', BatchSize=100)
    e = norm(aut['x']-ref['x'])/norm(ref['x'])
    assert e < 1e-6, f'{e:.2e}'

def t_batchsize_invariance():
    cfg = baseCfg(); x = None
    for bs in [50, 137, 400, 10**9]:
        _, s = runIt(cfg, [1, 2], backend='qr', BatchSize=bs, BatchOverlap=3)
        if x is None: x = s['x']
        else:
            e = norm(s['x']-x)/norm(x)
            assert e < 1e-6, f'bs={bs} {e:.2e}'

def t_psd_vs_toeplitz():
    dt = 5; a = [1, -1.6, .9, -.25]; s = .7
    f = np.linspace(1e-9, .5/dt, 4096)
    H = 1/(np.exp(-1j*2*np.pi*np.outer(f, np.arange(4))*dt)@np.array(a))
    cfg = baseCfg()
    cfg['cov'] = {'pos': dict(type='psd', f=f, S=s**2*dt*2*np.abs(H)**2),
                  'kbr': dict(type='diag', sigma=.2)}
    cfg['tag'] = 'psd'
    _, sA = runIt(cfg, [1, 2], backend='qr', BatchSize=400, BatchOverlap=20,
                  UserOpt=dict(psdMode='ar', arOrder=20))
    _, sT = runIt(cfg, [1, 2], backend='qr', BatchSize=400, BatchOverlap=20,
                  UserOpt=dict(psdMode='toeplitz'))
    e = norm(sA['x']-sT['x'])/norm(sT['x'])
    assert e < 1e-3, f'AR vs Toeplitz {e:.2e}'
    globals()['_psd'] = e

def t_qr_beats_normal():
    # isolated backend comparison: pipeline test is masked by regularisation
    from scipy.linalg import qr as sqr, solve_triangular, cholesky
    rng = np.random.default_rng(5); m, n = 800, 8
    U, _ = np.linalg.qr(rng.standard_normal((m, n)))
    V, _ = np.linalg.qr(rng.standard_normal((n, n)))
    A = U @ np.diag(np.logspace(0, -6, n)) @ V.T
    l = A @ rng.standard_normal(n) + 1e-10*rng.standard_normal(m)
    xr = lstsq(A, l, rcond=None)[0]
    aN = NormalEq(n); aN.add(A, l)
    aQ = QRSys(n);    aQ.add(A, l)
    eN = norm(aN.solve()[0]-xr)/norm(xr); eQ = norm(aQ.solve()[0]-xr)/norm(xr)
    assert eQ < eN, f'qr {eQ:.2e} not better than normal {eN:.2e}'
    assert eQ < 1e-10, f'qr {eQ:.2e}'
    globals()['_cond'] = (cond(A), eN, eQ)

def t_diag_flags():
    cfg = baseCfg(collinear=True, collinearEps=1e-9)
    acc, _ = runIt(cfg, [1, 2], backend='qr')
    assert acc.diag[0]['condNll'] > 1e8, f"cond {acc.diag[0]['condNll']:.2e}"
    assert acc.diag[0]['status'] != 'ok'

def t_tikhonov():
    cfg = baseCfg(collinear=True, collinearEps=1e-10)
    acc, sol = runIt(cfg, [1, 2], backend='qr', Reg='tikhonov', RegLambda=1e-3)
    assert np.all(np.isfinite(sol['x']))
    assert acc.diag[0]['regUsed'] == 'tikhonov'

def t_auto_reg():
    cfg = baseCfg()
    acc, _ = runIt(cfg, [1, 2], Reg='auto')
    assert acc.diag[0]['regUsed'] == 'none', acc.diag[0]['regUsed']
    cfg2 = baseCfg(collinear=True, collinearEps=1e-11)
    acc2, _ = runIt(cfg2, [1, 2], Reg='auto')
    assert acc2.diag[0]['regUsed'] != 'none'

def t_tsvd_exact_deficiency():
    cfg = baseCfg(collinear=True, collinearEps=0.0)
    acc, sol = runIt(cfg, [1, 2], backend='qr', Reg='tsvd')
    assert np.all(np.isfinite(sol['x']))
    assert acc.diag[0]['rankDef'] >= 1, f"rankDef {acc.diag[0]['rankDef']}"

def t_resid_sigma0():
    cfg = baseCfg()
    acc, sol = runIt(cfg, [1, 2, 3], backend='qr')
    res = acc.residuals(sol)
    e = abs(np.sqrt(res['s0sq'])-sol['sigma0'])/sol['sigma0']
    assert e < 1e-8, f'sigma0 mismatch {e:.2e}'

def t_resid_white():
    cfg = baseCfg()
    acc, sol = runIt(cfg, [1, 2, 3, 4, 5], backend='qr')
    res = acc.residuals(sol)
    w = res['sys'][0]['w']
    assert abs(np.std(w)-1) < .12, f'std {np.std(w):.4f}'
    rho = (w[:-1]@w[1:])/(w@w)
    assert abs(rho) < .1, f'rho1 {rho:.4f}'
    globals()['_res'] = (np.std(w), rho)

def t_local_cov_pd():
    cfg = baseCfg()
    acc, sol = runIt(cfg, [1, 2], backend='qr')
    C = acc.localCovariance(cfg['sats'][0], sol['Cxx'], sol['s0sq'])
    assert norm(C-C.T) < 1e-12
    assert np.min(np.real(eig(C)[0])) > 0

def t_useropt_forwarded():
    cfg = baseCfg(); seen = []
    def probe(d, ot, a, uo):
        assert uo.get('myFlag') == 77, 'opt not forwarded'
        seen.append(1)
        return synthLoader(cfg, 1, d, ot, a, uo)
    acc = Accumulator(cfg['nc'], UserOpt=dict(myFlag=77))
    acc.addSystem(1, cfg['nl'], probe, [1])
    assert seen

def t_empty_arc_skipped():
    cfg = baseCfg()
    def gappy(d, ot, a, uo):
        out = synthLoader(cfg, 1, d, ot, a, uo)
        if d != 'meta' and a == 2:
            out['A'] = np.zeros((0, cfg['nc'])); out['B'] = np.zeros((0, cfg['nl']))
            out['l'] = np.zeros(0)
        return out
    acc = Accumulator(cfg['nc'])
    acc.addSystem(1, cfg['nl'], gappy, [1, 2])
    sol = acc.solve()
    assert np.all(np.isfinite(sol['x']))

for nm, fn in [
    ('psd2acv matches analytic AR',        t_psd2acv),
    ('levinson recovers AR coeffs',        t_levinson),
    ('FIR filter == Cholesky whitening',   t_filter_eq_chol),
    ('armaACV vs analytic ARMA(1,1)',      t_arma_acv),
    ('recovers truth (white noise)',       t_recover_white),
    ('sigma0 unbiased (~1)',               t_sigma0),
    ('backends agree (well-cond)',         t_backends_agree),
    ('preelimination == joint solve',      t_preelim_exact),
    ('back-substitution recovers locals',  t_backsub),
    ('batch overlap exactness',            t_overlap),
    ('auto overlap from filterMemory',     t_auto_overlap),
    ('batch-size invariance',              t_batchsize_invariance),
    ('PSD->AR == dense Toeplitz',          t_psd_vs_toeplitz),
    ('QR beats normal (ill-cond)',         t_qr_beats_normal),
    ('diagnostics flag ill-conditioning',  t_diag_flags),
    ('tikhonov regularisation runs',       t_tikhonov),
    ('auto reg selects by status',         t_auto_reg),
    ('tsvd handles exact rank deficiency', t_tsvd_exact_deficiency),
    ('residuals reproduce sigma0',         t_resid_sigma0),
    ('whitened residuals are white',       t_resid_white),
    ('local covariance is PD',             t_local_cov_pd),
    ('UserOpt reaches loader',             t_useropt_forwarded),
    ('empty arc skipped',                  t_empty_arc_skipped),
]:
    T(nm, fn)

print('='*78)
npass = sum(1 for r in R if r[1] == 'PASS')
for nm, st, msg in R:
    mark = 'ok  ' if st == 'PASS' else 'FAIL'
    print(f'  [{mark}] {nm:<38s} {msg}')
print('='*78)
print(f'  {npass}/{len(R)} passed')
print('='*78)
if '_ov' in globals():
    print(f"batch overlap:  ov=p -> {_ov[0]:.2e}   ov=0 -> {_ov[1]:.2e}  ({_ov[1]/_ov[0]:.0f}x worse)")
if '_cond' in globals():
    print(f"ill-cond cond(A)={_cond[0]:.2e}:  normal {_cond[1]:.2e}   qr {_cond[2]:.2e}  ({_cond[1]/_cond[2]:.1f}x better)")
if '_psd' in globals():
    print(f"PSD->AR(20) vs dense Toeplitz normals: {_psd:.2e}")
if '_res' in globals():
    print(f"whitened residuals: std={_res[0]:.4f}  rho(1)={_res[1]:+.4f}")
