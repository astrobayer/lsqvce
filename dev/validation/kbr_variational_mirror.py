"""Mirror of template/kbr_variational.m (WP-122): the physics of a variational-equation
range-rate builder, checked against oracles that were not built here.

  1. normalized solid harmonics  Vb_nm + i Wb_nm = (R/r)^(n+1) Pbar_nm(sin phi) e^(i m lam)
     by the Cartesian (pole-safe) recursion  -- vs scipy.special.lpmv + factorials
  2. the derivative operators Dx, Dy, Dz (sparse, degree n -> n+1) -- vs central
     finite differences of the harmonics
  3. acceleration of a C00 + C20 field -- vs the closed-form J2 acceleration;
     gravity gradient -- vs finite differences of the acceleration
  4. ABM8 PECE (exact rational coefficients) at 5 s -- vs the Kepler orbit, one day
  5. variational partials (state transition, SH, ACC bias/scale) -- vs re-integration
  6. range-rate partials -- vs finite differences of the computed range-rate
  7. closed loop: simulate a perturbed field, recover it by least squares
Run: python3 dev/validation/kbr_variational_mirror.py
"""
import numpy as np
from fractions import Fraction
from scipy.special import lpmv
from math import factorial, sqrt

GM, RE = 3.986004415e14, 6378136.3
OMEGA_E = 7.2921150e-5

# ------------------------------------------------------------------ index layout
def coeff_index(N):
    """(n, m, cs) for degree 0..N: per degree, C_n0..C_nn then S_n1..S_nn."""
    out = []
    for n in range(N + 1):
        out += [(n, m, 0) for m in range(n + 1)] + [(n, m, 1) for m in range(1, n + 1)]
    return out

def pos_of(N):
    return {k: i for i, k in enumerate(coeff_index(N))}

# ------------------------------------------------------------------ 1. harmonics
def harmonics(r, L):
    """Vb, Wb (L+1 x L+1) at the Cartesian point r (TRF)."""
    x, y, z = r
    rr2 = x * x + y * y + z * z
    q = RE / rr2
    V = np.zeros((L + 1, L + 1)); W = np.zeros((L + 1, L + 1))
    V[0, 0] = RE / np.sqrt(rr2)
    for m in range(0, L + 1):
        if m > 0:
            c = sqrt((2 * m + 1) / ((1 if m == 1 else 2) * m))
            V[m, m] = c * q * (x * V[m - 1, m - 1] - y * W[m - 1, m - 1])
            W[m, m] = c * q * (x * W[m - 1, m - 1] + y * V[m - 1, m - 1])
        for n in range(m + 1, L + 1):
            a = sqrt((2 * n - 1) * (2 * n + 1) / ((n - m) * (n + m)))
            V[n, m] = a * z * q * V[n - 1, m]
            W[n, m] = a * z * q * W[n - 1, m]
            if n - 2 >= m:
                b = sqrt((2 * n + 1) * (n + m - 1) * (n - m - 1) / ((2 * n - 3) * (n + m) * (n - m)))
                V[n, m] -= b * RE * q * V[n - 2, m]
                W[n, m] -= b * RE * q * W[n - 2, m]
    return V, W

def hvec(V, W, N):
    return np.array([V[n, m] if cs == 0 else W[n, m] for (n, m, cs) in coeff_index(N)])

def harmonics_ref(r, L):
    """Oracle: scipy's associated Legendre functions (Condon-Shortley phase removed)."""
    x, y, z = r
    rn = np.linalg.norm(r); sphi = z / rn; lam = np.arctan2(y, x)
    V = np.zeros((L + 1, L + 1)); W = np.zeros((L + 1, L + 1))
    for n in range(L + 1):
        for m in range(n + 1):
            Nnm = sqrt((1 if m == 0 else 2) * (2 * n + 1) * factorial(n - m) / factorial(n + m))
            P = (-1) ** m * lpmv(m, n, sphi) * Nnm
            V[n, m] = (RE / rn) ** (n + 1) * P * np.cos(m * lam)
            W[n, m] = (RE / rn) ** (n + 1) * P * np.sin(m * lam)
    return V, W

# ------------------------------------------------------------------ 2. operators
def ratios(n, m):
    up = sqrt((1 if m == 0 else 2) * (2 * n + 1) * (n + m + 1) * (n + m + 2) / (2 * (2 * n + 3)))
    same = sqrt((2 * n + 1) * (n + m + 1) / ((2 * n + 3) * (n - m + 1)))
    down = sqrt(2 * (2 * n + 1) / ((1 if m == 1 else 2) * (2 * n + 3) * (n - m + 2) * (n - m + 1))) if m >= 1 else 0.0
    return up, same, down

def operators(N):
    """Dx, Dy, Dz: (N+2)^2 x (N+1)^2, R * d/dx of a degree-N expansion in the degree-(N+1) basis."""
    src, dst = pos_of(N), pos_of(N + 1)
    D = [np.zeros(((N + 2) ** 2, (N + 1) ** 2)) for _ in range(3)]
    for (n, m, cs), j in src.items():
        up, same, down = ratios(n, m)
        k = (n - m + 2) * (n - m + 1)
        V = lambda nn, mm: dst[(nn, mm, 0)]
        Wi = lambda nn, mm: dst[(nn, mm, 1)] if mm > 0 else None
        def put(Di, idx, val):
            if idx is not None:
                Di[idx, j] += val
        if cs == 0:
            if m == 0:
                put(D[0], V(n + 1, 1), -up)
                put(D[1], Wi(n + 1, 1), -up)
            else:
                put(D[0], V(n + 1, m + 1), -0.5 * up); put(D[0], V(n + 1, m - 1), 0.5 * k * down)
                put(D[1], Wi(n + 1, m + 1), -0.5 * up); put(D[1], Wi(n + 1, m - 1), -0.5 * k * down)
            put(D[2], V(n + 1, m), -(n - m + 1) * same)
        else:
            put(D[0], Wi(n + 1, m + 1), -0.5 * up); put(D[0], Wi(n + 1, m - 1), 0.5 * k * down)
            put(D[1], V(n + 1, m + 1), 0.5 * up); put(D[1], V(n + 1, m - 1), 0.5 * k * down)
            put(D[2], Wi(n + 1, m), -(n - m + 1) * same)
    return D

# ------------------------------------------------------------------ 3. the field
class Field:
    """A field c (degree N, normalized, in coeff_index order): its acceleration and
    gradient operators precomputed once, evaluated per point with ONE harmonic sweep."""
    def __init__(self, c, N, est):
        self.N = N; self.c = np.asarray(c, float)
        D1, D2 = operators(N), operators(N + 1)
        self.ac = [D1[i] @ self.c for i in range(3)]                       # degree N+1
        self.gc = [[D2[i] @ self.ac[j] for j in range(3)] for i in range(3)]  # degree N+2
        self.est = np.asarray(est, int)                                     # estimated indices
        self.Dest = [D1[i][:, self.est] for i in range(3)]                  # partial operators
    def eval(self, r_trf):
        V, W = harmonics(r_trf, self.N + 2)
        h1 = hvec(V, W, self.N + 1); h2 = hvec(V, W, self.N + 2)
        s1, s2 = GM / RE ** 2, GM / RE ** 3
        a = s1 * np.array([h1 @ self.ac[i] for i in range(3)])
        G = s2 * np.array([[h2 @ self.gc[i][j] for j in range(3)] for i in range(3)])
        P = s1 * np.array([h1 @ self.Dest[i] for i in range(3)])            # 3 x np
        return a, G, P

def j2_accel(r, C20):
    """Closed form: point mass + J2 (J2 = -sqrt(5) C20bar)."""
    J2 = -sqrt(5) * C20
    x, y, z = r; rn = np.linalg.norm(r)
    f = 1.5 * J2 * GM * RE ** 2 / rn ** 5
    a0 = -GM * r / rn ** 3
    return a0 + f * np.array([x * (5 * z * z / rn ** 2 - 1), y * (5 * z * z / rn ** 2 - 1), z * (5 * z * z / rn ** 2 - 3)])

def check_core():
    rs = np.random.RandomState(1)
    pts = [np.array([6.8e6, 1.2e6, -2.0e6]), np.array([10.0, -20.0, 6.87e6]),       # near the pole
           rs.randn(3) * 4e6 + np.array([0, 0, 1e6])]
    L = 20
    err = 0
    for p in pts:
        V, W = harmonics(p, L); Vr, Wr = harmonics_ref(p, L)
        err = max(err, np.abs(V - Vr).max() / np.abs(Vr).max(), np.abs(W - Wr).max() / np.abs(Vr).max())
    print('1. harmonics vs scipy lpmv (L 20, 3 points incl. 1e-6 rad off the pole): max rel %.1e' % err)
    N = 12; D = operators(N); h = 1.0
    worst = 0
    for p in pts:
        V0, W0 = harmonics(p, N + 1); hv1 = hvec(V0, W0, N + 1)
        for i in range(3):
            e = np.zeros(3); e[i] = h
            Vp, Wp = harmonics(p + e, N); Vm, Wm = harmonics(p - e, N)
            fd = (hvec(Vp, Wp, N) - hvec(Vm, Wm, N)) / (2 * h) * RE
            an = hv1 @ D[i]
            worst = max(worst, np.abs(an - fd).max() / np.abs(an).max())
    print('2. operators R d/dx_i vs central differences (N 12, h 1 m): max rel %.1e' % worst)
    pos = pos_of(4); c = np.zeros(25); c[pos[(0, 0, 0)]] = 1.0; c[pos[(2, 0, 0)]] = -4.841694573e-4
    F = Field(c, 4, [pos[(2, 0, 0)]])
    worst = 0; wg = 0
    for p in pts:
        a, G, _ = F.eval(p)
        worst = max(worst, np.abs(a - j2_accel(p, c[pos[(2, 0, 0)]])).max() / np.linalg.norm(a))
        Gfd = np.zeros((3, 3))
        for j in range(3):
            e = np.zeros(3); e[j] = 1.0
            Gfd[:, j] = (F.eval(p + e)[0] - F.eval(p - e)[0]) / 2
        wg = max(wg, np.abs(G - Gfd).max() / np.abs(G).max())
    print('3. acceleration C00 + C20 vs closed-form J2: max rel %.1e;  gradient vs FD of a: max rel %.1e' % (worst, wg))
    return err, worst, wg

# ------------------------------------------------------------------ 4. integrator
AB8 = np.array([434241, -1152169, 2183877, -2664477, 2102243, -1041723, 295767, -36799]) / 120960.0
AM8 = np.array([36799, 139849, -121797, 123133, -88547, 41499, -11351, 1375]) / 120960.0

def rot_trf2crf(t):
    """Default Earth rotation: about z, theta = OMEGA_E * t (t in seconds of the
    reference epoch). r_crf = R r_trf; v_crf = R v_trf + Rdot r_trf."""
    th = OMEGA_E * t; c, s = np.cos(th), np.sin(th)
    R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])
    Rd = OMEGA_E * np.array([[-s, -c, 0], [c, -s, 0], [0, 0, 0.0]])
    return R, Rd

class Dyn:
    """Force model of one satellite: field + accelerometer (Q(t), acc(t)) with
    calibration b0, s0; column layout of Psi: [IC 6 | SH np | bias 3 | scale 3]."""
    def __init__(self, field, accfun, b0=np.zeros(3), s0=np.ones(3)):
        self.f = field; self.accfun = accfun; self.b0 = b0; self.s0 = s0
        self.np = len(field.est); self.ncol = 6 + self.np + 6
    def rhs(self, t, Y):
        r, v = Y[:3], Y[3:6]; Psi = Y[6:].reshape(6, self.ncol)
        R, _ = rot_trf2crf(t)
        a_t, G_t, P_t = self.f.eval(R.T @ r)
        Q, acc = self.accfun(t)
        a = R @ a_t + Q @ (self.b0 + self.s0 * acc)
        G = R @ G_t @ R.T
        Fp = np.zeros((3, self.ncol))
        Fp[:, 6:6 + self.np] = R @ P_t
        Fp[:, 6 + self.np:9 + self.np] = Q
        Fp[:, 9 + self.np:] = Q * acc[None, :]
        dPsi = np.vstack([Psi[3:6], G @ Psi[0:3] + Fp])
        return np.concatenate([v, a, dPsi.ravel()])

def integrate(dyn, t0, y0, h, nstep, sub=20):
    """ABM8 PECE on the grid t0 + k h, k = 0..nstep; RK4 (h/sub) start. Returns
    Y (nstep+1 x len) at the grid."""
    Y0 = np.concatenate([y0, np.hstack([np.eye(6), np.zeros((6, dyn.ncol - 6))]).ravel()])
    Ys = [Y0]; Fs = [dyn.rhs(t0, Y0)]
    Y = Y0.copy(); hs = h / sub
    for k in range(1, min(8, nstep + 1)):                 # start: RK4 substeps
        for j in range(sub):
            t = t0 + (k - 1) * h + j * hs
            k1 = dyn.rhs(t, Y); k2 = dyn.rhs(t + hs / 2, Y + hs / 2 * k1)
            k3 = dyn.rhs(t + hs / 2, Y + hs / 2 * k2); k4 = dyn.rhs(t + hs, Y + hs * k3)
            Y = Y + hs / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        Ys.append(Y.copy()); Fs.append(dyn.rhs(t0 + k * h, Y))
    for k in range(8, nstep + 1):
        hist = Fs[-1:-9:-1]                                # f_n, f_n-1, ..., f_n-7
        Yp = Ys[-1] + h * sum(b * f for b, f in zip(AB8, hist))
        fp = dyn.rhs(t0 + k * h, Yp)
        Yc = Ys[-1] + h * (AM8[0] * fp + sum(b * f for b, f in zip(AM8[1:], hist[:7])))
        Ys.append(Yc); Fs.append(dyn.rhs(t0 + k * h, Yc))
    return np.array(Ys)

def kepler(r0, v0, dt):
    """Oracle: two-body propagation by Kepler's equation (elliptic)."""
    rn = np.linalg.norm(r0); vv = v0 @ v0
    a = 1 / (2 / rn - vv / GM); n = sqrt(GM / a ** 3)
    ecosE = 1 - rn / a; esinE = (r0 @ v0) / sqrt(GM * a)
    E0 = np.arctan2(esinE, ecosE); e = np.hypot(ecosE, esinE)
    M = E0 - esinE + n * dt; E = M
    for _ in range(50):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    dE = E - E0
    f = 1 - a / rn * (1 - np.cos(dE)); g = dt - (dE - np.sin(dE)) / n
    return f * r0 + g * v0

def leo_state(rs, inc=89.0, h=490e3, u=0.0, node=0.0):
    a = RE + h; vc = sqrt(GM / a); i = np.radians(inc)
    Rn = np.array([[np.cos(node), -np.sin(node), 0], [np.sin(node), np.cos(node), 0], [0, 0, 1]])
    p = np.array([np.cos(u), np.sin(u) * np.cos(i), np.sin(u) * np.sin(i)])
    q = np.array([-np.sin(u), np.cos(u) * np.cos(i), np.cos(u) * np.sin(i)])
    return np.concatenate([Rn @ p * a, Rn @ q * vc])

def test_field(N):
    """Deterministic field: GM point mass, C20, Kaula-like n>=3 (seeded)."""
    pos = pos_of(N); c = np.zeros((N + 1) ** 2)
    c[pos[(0, 0, 0)]] = 1.0; c[pos[(2, 0, 0)]] = -4.841694573e-4
    rs = np.random.RandomState(7)
    for (n, m, cs), j in pos.items():
        if n >= 2 and not (n == 2 and m == 0 and cs == 0):
            c[j] = 1e-5 / n ** 2 * rs.randn() * (0.3 if n == 2 else 1.0)
    return c

def no_acc(t):
    return np.eye(3), np.zeros(3)

def check_integrator():
    pos = pos_of(2); c = np.zeros(9); c[pos[(0, 0, 0)]] = 1.0
    dyn = Dyn(Field(c, 2, [pos[(2, 0, 0)]]), no_acc)
    y0 = leo_state(None)
    out = []
    for h, nst in ((5.0, 17280),):
        Y = integrate(dyn, 0.0, y0, h, nst)
        rk = kepler(y0[:3], y0[3:6], h * nst)
        out.append(np.linalg.norm(Y[-1, :3] - rk))
        print('4. ABM8 PECE, h %.0f s, one day, point mass: |r - Kepler| = %.2e m' % (h, out[-1]))
    return out

# ------------------------------------------------------------------ 5. variational partials
W_ORB = 2 * np.pi / 5650.0

def sat_acc(phase):
    """Synthetic accelerometer + attitude: Q(t) a slowly turning SRF->CRF rotation,
    acc(t) drag-like along x with once- and twice-per-rev variations (m/s^2)."""
    def f(t):
        th = W_ORB * t + phase
        c, s = np.cos(th), np.sin(th)
        Q = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]]) @ np.array([[1, 0, 0], [0, np.cos(0.3), -np.sin(0.3)], [0, np.sin(0.3), np.cos(0.3)]])
        acc = np.array([-1.0e-7 * (1 + 0.3 * np.cos(th)), 2e-8 * np.sin(th), 5e-8 * np.cos(2 * th)])
        return Q, acc
    return f

def check_partials(N=8, nstep=1080, h=5.0):
    pos = pos_of(N); c = test_field(N)
    est = [pos[k] for k in coeff_index(N) if k[0] >= 2]
    y0 = leo_state(None, u=0.4)
    def run(cc, yy, b0=np.zeros(3), s0=np.ones(3)):
        return integrate(Dyn(Field(cc, N, est), sat_acc(0.2), b0, s0), 0.0, yy, h, nstep)
    Y = run(c, y0); npar = len(est); ncol = 6 + npar + 6
    Psi = Y[-1, 6:].reshape(6, ncol)
    res = {}
    for name, col, delta in (('IC x', 0, 1.0), ('IC vy', 4, 1e-3)):
        yy = y0.copy(); yy[col] += delta
        d = run(c, yy)[-1, :6] - Y[-1, :6]
        res[name] = np.linalg.norm(Psi[:, col] * delta - d) / np.linalg.norm(d)
    j = est.index(pos[(5, 3, 1)]); cc = c.copy(); cc[pos[(5, 3, 1)]] += 1e-9
    d = run(cc, y0)[-1, :6] - Y[-1, :6]
    res['S53'] = np.linalg.norm(Psi[:, 6 + j] * 1e-9 - d) / np.linalg.norm(d)
    d = run(c, y0, b0=np.array([0, 1e-8, 0]))[-1, :6] - Y[-1, :6]
    res['bias y'] = np.linalg.norm(Psi[:, 6 + npar + 1] * 1e-8 - d) / np.linalg.norm(d)
    d = run(c, y0, s0=np.array([1.01, 1, 1]))[-1, :6] - Y[-1, :6]
    res['scale x'] = np.linalg.norm(Psi[:, 9 + npar] * 0.01 - d) / np.linalg.norm(d)
    print('5. Psi(t_end) * delta vs re-integration, 1.5 h arc, degree %d: ' % N +
          ', '.join('%s %.1e' % kv for kv in res.items()))
    return res

# ------------------------------------------------------------------ 6./7. range-rate rows, closed loop
def rr_and_partials(YC, YD):
    """rho_dot and its partials w.r.t. the 6-states of C and D, per epoch."""
    dr = YD[:, :3] - YC[:, :3]; dv = YD[:, 3:6] - YC[:, 3:6]
    rho = np.linalg.norm(dr, axis=1); e = dr / rho[:, None]
    rd = np.sum(e * dv, axis=1)
    hvec_ = (dv - rd[:, None] * e) / rho[:, None]
    J = np.hstack([hvec_, e])                       # d rho_dot / d (r_D, v_D); C gets -J
    return rd, J

def arg_lat(Y):
    r, v = Y[:, :3], Y[:, 3:6]
    hh = np.cross(r, v); hh /= np.linalg.norm(hh, axis=1)[:, None]
    nod = np.cross(np.array([0, 0, 1.0]), hh); nod /= np.linalg.norm(nod, axis=1)[:, None]
    return np.arctan2(np.sum(np.cross(nod, r) * hh, axis=1), np.sum(nod * r, axis=1))

def design(YC, YD, npar, t):
    """Rows of one arc. Columns: [SH np | C: IC 6, bias 3, scale 3 | D: same | kbr 4]."""
    ncol = 6 + npar + 6; m = len(t)
    PC = YC[:, 6:].reshape(m, 6, ncol); PD = YD[:, 6:].reshape(m, 6, ncol)
    nu = npar + 12 + 12 + 4
    def satcols(P, off):
        out = np.zeros(P.shape[:2] + (nu,))
        out[..., :npar] = P[..., 6:6 + npar]
        out[..., off:off + 6] = P[..., :6]
        out[..., off + 6:off + 12] = P[..., 6 + npar:]
        return out
    SC, SD = satcols(PC, npar), satcols(PD, npar + 12)
    rd, J = rr_and_partials(YC, YD)
    Akbr = np.einsum('mk,mkn->mn', J, SD) - np.einsum('mk,mkn->mn', J, SC)
    tm = t - t.mean(); u = arg_lat(YC)
    Akbr[:, -4:] = np.column_stack([np.ones(m), tm / 3600.0, np.cos(u), np.sin(u)])
    Apc = SC[:, :3, :].reshape(3 * m, nu); Apd = SD[:, :3, :].reshape(3 * m, nu)
    return rd, Akbr, Apc, Apd

def closed_loop(N=8, nstep=17280, h=5.0, iterations=2):
    """Simulate a perturbed field and ACC calibration, recover them; each further
    iteration re-linearizes at the previous solution (field, ACC b/s, ICs) -- the
    a-priori the template takes from the user."""
    pos = pos_of(N); c_ref = test_field(N)
    est = [pos[k] for k in coeff_index(N) if k[0] >= 2]; npar = len(est)
    rs = np.random.RandomState(11)
    dc = np.zeros_like(c_ref); dc[est] = 2e-9 * rs.randn(npar) / np.array([coeff_index(N)[j][0] for j in est]) ** 1.0
    b_t = {'C': np.array([3e-9, -2e-9, 1e-9]), 'D': np.array([-1e-9, 2.5e-9, 0.5e-9])}
    s_t = {'C': np.array([1.02, 0.98, 1.01]), 'D': np.array([0.99, 1.03, 1.0])}
    sep = 220e3 / (RE + 490e3)
    y0 = {'C': leo_state(None, u=0.4), 'D': leo_state(None, u=0.4 - sep)}
    acc = {'C': sat_acc(0.2), 'D': sat_acc(0.2 - sep)}
    truth = {k: integrate(Dyn(Field(c_ref + dc, N, est), acc[k], b_t[k], s_t[k]), 0.0, y0[k], h, nstep) for k in 'CD'}
    t = np.arange(nstep + 1) * h
    obs_pos = {k: truth[k][:, :3] for k in 'CD'}
    rd_obs, _ = rr_and_partials(truth['C'], truth['D'])
    # a-priori: the reference field, b = 0, s = 1, ICs = the observed state at arc start
    ap = {'c': c_ref.copy(), 'b': {k: np.zeros(3) for k in 'CD'}, 's': {k: np.ones(3) for k in 'CD'},
          'y0': {k: truth[k][0, :6].copy() for k in 'CD'}}
    errs = []
    for it in range(iterations):
        ref = {k: integrate(Dyn(Field(ap['c'], N, est), acc[k], ap['b'][k], ap['s'][k]), 0.0, ap['y0'][k], h, nstep) for k in 'CD'}
        rd_ref, Akbr, Apc, Apd = design(ref['C'], ref['D'], npar, t)
        if it == 0:
            j = est.index(pos[(6, 2, 0)]); cc = ap['c'].copy(); cc[pos[(6, 2, 0)]] += 1e-9
            pert = {k: integrate(Dyn(Field(cc, N, est), acc[k]), 0.0, ap['y0'][k], h, nstep) for k in 'CD'}
            drd = rr_and_partials(pert['C'], pert['D'])[0] - rr_and_partials(
                *[integrate(Dyn(Field(ap['c'], N, est), acc[k]), 0.0, ap['y0'][k], h, nstep) for k in 'CD'])[0]
            e6 = np.abs(Akbr[:, j] * 1e-9 - drd).max() / np.abs(drd).max()
            print('6. range-rate partial (C62) vs finite difference: max rel %.1e (|d rho_dot| max %.1e m/s)' % (e6, np.abs(drd).max()))
        sk, sp = 2e-7, 0.02
        A = np.vstack([Akbr / sk, Apc / sp, Apd / sp])
        l = np.concatenate([(rd_obs - rd_ref) / sk, (obs_pos['C'] - ref['C'][:, :3]).ravel() / sp,
                            (obs_pos['D'] - ref['D'][:, :3]).ravel() / sp])
        cn = np.linalg.norm(A, axis=0)                     # equilibrate: the columns span ~20 decades
        x = np.linalg.lstsq(A / cn, l, rcond=None)[0] / cn
        ap['c'][est] += x[:npar]
        for k, off in (('C', npar), ('D', npar + 12)):
            ap['y0'][k] += x[off:off + 6]; ap['b'][k] += x[off + 6:off + 9]; ap['s'][k] += x[off + 9:off + 12]
        e7 = np.linalg.norm(ap['c'][est] - c_ref[est] - dc[est]) / np.linalg.norm(dc[est])
        errs.append(e7)
        print('7. closed loop, iteration %d (one day, degree %d, %d SH): |c - c_true| / |dc| = %.1e;  '
              'scale C err %.1e;  bias C err %.1e m/s^2;  kbr empiricals %s'
              % (it + 1, N, npar, e7, np.abs(ap['s']['C'] - s_t['C']).max(),
                 np.abs(ap['b']['C'] - b_t['C']).max(), np.array2string(x[-4:], precision=2)))
    return e6, errs

if __name__ == '__main__':
    check_core()
    check_integrator()
    check_partials()
    closed_loop()
