#!/usr/bin/env python3
# v10review.py -- numerical evidence for the R2 five-expert review (2026-08-16).
#
# Companion to HANDOVER.md section 8 (REVIEW R2). Each test either
# VERIFIES a formula as implemented in the MATLAB tree (transcribed 1:1), or
# DEMONSTRATES a defect found by reading the code, with the corrected formula
# alongside. numpy only, deterministic, exit 0 iff every check behaves as the
# review claims.
#
# Created 2026-08-16 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
# model Claude Fable 5.

import numpy as np

rng = np.random.default_rng(20260816)
FAIL = []


def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    print(f"[{tag}] {name}" + (f"  -- {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)


# ---------------------------------------------------------------------------
# Shared synthetic block-arrow problem:  l_s = A_s x + B_s y_s + sigma_s e_s
# ---------------------------------------------------------------------------
nc = 6
sysdef = [dict(n=400, nl=4, sig=1.0), dict(n=300, nl=3, sig=2.0)]
x_true = rng.normal(size=nc)
for s in sysdef:
    s["A"] = rng.normal(size=(s["n"], nc))
    s["B"] = rng.normal(size=(s["n"], s["nl"]))
    s["y"] = rng.normal(size=s["nl"])
    s["l"] = s["A"] @ x_true + s["B"] @ s["y"] + s["sig"] * rng.normal(size=s["n"])
m_tot = sum(s["n"] for s in sysdef)
u_tot = nc + sum(s["nl"] for s in sysdef)


def assemble(sig):
    """Whitened per-system blocks, reduced common normal, solution, residuals."""
    Ncc = np.zeros((nc, nc))
    ncc = np.zeros(nc)
    Ss = []
    for s, sg in zip(sysdef, sig):
        Aw, Bw, lw = s["A"] / sg, s["B"] / sg, s["l"] / sg
        Nll = Bw.T @ Bw
        Nlc = Bw.T @ Aw
        S = Aw.T @ Aw - Nlc.T @ np.linalg.solve(Nll, Nlc)
        r = Aw.T @ lw - Nlc.T @ np.linalg.solve(Nll, Bw.T @ lw)
        Ncc += S
        ncc += r
        Ss.append((S, Nll, Nlc, Aw, Bw, lw))
    x = np.linalg.solve(Ncc, ncc)
    out = []
    for (S, Nll, Nlc, Aw, Bw, lw), s in zip(Ss, sysdef):
        y = np.linalg.solve(Nll, Bw.T @ lw - Nlc @ x)
        v = lw - Aw @ x - Bw @ y          # whitened residual
        out.append(dict(S=S, Nll=Nll, Nlc=Nlc, Aw=Aw, Bw=Bw, y=y, vtv=v @ v))
    return Ncc, x, out


def leverage_exact(sig):
    """Per-system absorbed dof and redundancy from the FULL joint normal."""
    rows = []
    for s, sg in zip(sysdef, sig):
        Z1 = s["A"] / sg
        rows.append((Z1, s["B"] / sg))
    ntot = sum(z[0].shape[0] for z in rows)
    M = np.zeros((ntot, u_tot))
    r0, c0 = 0, nc
    spans = []
    for Aw, Bw in rows:
        n, nl = Aw.shape[0], Bw.shape[1]
        M[r0:r0 + n, :nc] = Aw
        M[r0:r0 + n, c0:c0 + nl] = Bw
        spans.append((r0, r0 + n))
        r0 += n
        c0 += nl
    Nf = M.T @ M
    H = np.einsum("ij,ij->i", M, M @ np.linalg.inv(Nf))
    return [H[a:b].sum() for a, b in spans], spans, H


# ---------------------------------------------------------------------------
# T1  leverage block identities (lsqvce.leverage / leverageBlocks)
# ---------------------------------------------------------------------------
sig0 = [1.0, 1.0]
Ncc_red, xh, per = assemble(sig0)
Sxx = np.linalg.inv(Ncc_red)
absorb, spans, Hfull = leverage_exact(sig0)
k = 0  # verify system 1 row-wise against the block-inverse formula
S, Nll, Nlc, Aw, Bw = per[k]["S"], per[k]["Nll"], per[k]["Nlc"], per[k]["Aw"], per[k]["Bw"]
G = np.linalg.solve(Nll, Nlc)
h_blk = (np.einsum("ij,ij->i", Aw @ Sxx, Aw)
         + np.einsum("ij,ij->i", Bw @ (np.linalg.inv(Nll) + G @ Sxx @ G.T), Bw)
         - 2 * np.einsum("ij,ij->i", Aw @ (Sxx @ G.T), Bw))
a, b = spans[k]
check("T1a leverage block formula == dense joint inverse",
      np.allclose(h_blk, Hfull[a:b], atol=1e-10),
      f"max|diff| = {np.max(np.abs(h_blk - Hfull[a:b])):.2e}")
check("T1b sum(h) over all systems == u",
      abs(sum(absorb) - u_tot) < 1e-8, f"sum(h) = {sum(absorb):.6f}, u = {u_tot}")

# ---------------------------------------------------------------------------
# T2  FINDING F1: vce 'hutchinson' substitutes tr(N^-1 N_g^cc) for r_g
# ---------------------------------------------------------------------------
print("\n--- T2: vce traceMethod 'hutchinson' (F1) ------------------------------")
r_exact = [s["n"] - a for s, a in zip(sysdef, absorb)]
t_trace = [np.trace(np.linalg.solve(Ncc_red, per[j]["S"])) for j in range(2)]
print(f"  exact r_g   = {np.round(r_exact, 2)}   sum = {sum(r_exact):.2f}"
      f"  (m-u = {m_tot - u_tot})")
print(f"  trace t_g   = {np.round(t_trace, 3)}   sum = {sum(t_trace):.3f}"
      f"  (nc  = {nc})")
check("T2a group sums differ by construction: sum r_g = m-u, sum t_g = nc",
      abs(sum(r_exact) - (m_tot - u_tot)) < 1e-6 and abs(sum(t_trace) - nc) < 1e-6)


def iterate_sigma(denominator, n_it=80):
    sig = np.array([1.0, 1.0])
    for _ in range(n_it):
        Nr, _, out = assemble(sig)
        if denominator == "exact":
            ab, _, _ = leverage_exact(sig)
            den = np.array([s["n"] - a for s, a in zip(sysdef, ab)])
        else:  # the vce.m 'hutchinson' substitution, exact-trace limit (no MC noise)
            den = np.array([np.trace(np.linalg.solve(Nr, o["S"])) for o in out])
        upd = np.sqrt(np.array([o["vtv"] for o in out]) / den)
        sig = sig * upd
        if np.max(np.abs(upd - 1)) < 1e-12:
            break
    return sig


sig_ok = iterate_sigma("exact")
sig_hu = iterate_sigma("trace")
truth = np.array([s["sig"] for s in sysdef])
resol = 0.5 * np.sqrt(2 / np.array(r_exact))   # 1-sigma resolution of a sigma scale
print(f"  truth sigma          = {truth}   (estimator 1-sigma resolution "
      f"{np.round(100 * resol, 1)} %)")
print(f"  exact-denominator    = {np.round(sig_ok, 4)}")
print(f"  trace-denominator    = {np.array2string(sig_hu, precision=3)}"
      f"   (inflation {np.array2string(sig_hu / truth, precision=2)}x)")
check("T2b exact denominator recovers the truth within 3x its own resolution",
      np.all(np.abs(sig_ok / truth - 1) < 3 * resol),
      f"rel. err. {np.round(100 * np.abs(sig_ok / truth - 1), 2)} %")
check("T2c trace denominator inflates without bound or lands >5x off the truth",
      np.all(sig_hu / truth > 5.0),
      "confirms +lsqvce/vce.m lines 173-186/229-236 substitute the wrong quantity")

# --- WP-32 FIX (added 2026-08-18): r_g = n_g - local_g - tr(N^-1 N_g^cc) ---
# The streamed cc-attribution is the SYMMETRIZED Schur split
#   N_g^cc = Mcc_g - 1/2 (Nlc_g' G_s + Nlc_s' G_g),   G_* = inv(Nll_s) Nlc_*
# (what lsqvce.streamContribution probes). The exact per-group remainder from
# the leverage block identities (inv(N)_ll = Nll^-1 + G Sxx G',
# inv(N)_lc = -G Sxx) is
#   local_g = tr((Nll^-1 + G Sxx G') Mll_g) - tr(Nlc_g Sxx G')
# which collapses to nl_s exactly for whole-system groups. Then
#   r_g = n_g - local_g - tr(Sxx N_g^cc)  ==  n_g - sum_{i in g} h_i  exactly.
print("\n--- T2d/T2e: the WP-32 corrected streamed redundancy ---------------------")
ok_split = True
sum_rfix = 0.0
for k in range(2):
    S, Nll, Nlc, Aw, Bw = (per[k]["S"], per[k]["Nll"], per[k]["Nlc"],
                           per[k]["Aw"], per[k]["Bw"])
    G = np.linalg.solve(Nll, Nlc)
    X = Sxx @ G.T                      # Sxx G'  (nc x nl)
    GSG = G @ X                        # G Sxx G'
    iNll = np.linalg.inv(Nll)
    a, b = spans[k]
    n = b - a
    half = n // 2
    for i0, i1 in [(0, half), (half, n)]:
        Ag, Bg = Aw[i0:i1], Bw[i0:i1]
        Mll, Nlcg = Bg.T @ Bg, Bg.T @ Ag
        Ng_cc = Ag.T @ Ag - 0.5 * (Nlcg.T @ G + Nlc.T @ np.linalg.solve(Nll, Nlcg))
        t_g = np.trace(Sxx @ Ng_cc)
        local = np.trace((iNll + GSG) @ Mll) - np.trace(Nlcg @ X)
        r_fix = (i1 - i0) - local - t_g
        r_ex = (i1 - i0) - Hfull[a + i0:a + i1].sum()
        ok_split &= abs(r_fix - r_ex) < 1e-8
        sum_rfix += r_fix
    loc_sys = (np.trace((iNll + GSG) @ (Bw.T @ Bw))
               - np.trace((Bw.T @ Aw) @ X))
    ok_split &= abs(loc_sys - sysdef[k]["nl"]) < 1e-8
check("T2d corrected split == exact h-sums per SUB-system group; whole-system "
      "local_g == nl_s; sum_g r_g = m-u",
      ok_split and abs(sum_rfix - (m_tot - u_tot)) < 1e-6,
      f"sum r_fix = {sum_rfix:.6f} (m-u = {m_tot - u_tot})")


def iterate_sigma_fix(n_it=80):
    sig = np.array([1.0, 1.0])
    for _ in range(n_it):
        Nr, _, out = assemble(sig)
        Sx = np.linalg.inv(Nr)
        den = np.array([s["n"] - s["nl"] - np.trace(Sx @ o["S"])
                        for s, o in zip(sysdef, out)])
        upd = np.sqrt(np.array([o["vtv"] for o in out]) / den)
        sig = sig * upd
        if np.max(np.abs(upd - 1)) < 1e-12:
            break
    return sig


sig_fx = iterate_sigma_fix()
check("T2e corrected denominator recovers the truth within 3x its resolution "
      "(exact-trace limit)",
      np.all(np.abs(sig_fx / truth - 1) < 3 * resol),
      f"sigma_fix = {np.round(sig_fx, 4)} vs truth {truth}")

# ---------------------------------------------------------------------------
# T3  FINDING F3: joint-path residual omits B*y; projector omits locals
# ---------------------------------------------------------------------------
print("\n--- T3: joint path, local parameters (F3) ------------------------------")
s = sysdef[0]
A, B, l = s["A"], s["B"], s["l"]          # sigma = 1, white: W = I
n, nl = s["n"], s["nl"]
M = np.hstack([A, B])
Nf = M.T @ M
xh_f = np.linalg.solve(Nf, M.T @ l)
xh_c, yh = xh_f[:nc], xh_f[nc:]
e_wrong = l - A @ xh_c                      # sweepStreamLags / oneArc formula
e_right = l - A @ xh_c - B @ yh             # the toolbox's own residuals() formula
Ncc1 = A.T @ A - (B.T @ A).T @ np.linalg.solve(B.T @ B, B.T @ A)
s0_code = n - np.trace(np.linalg.solve(Ncc1, A.T @ A))   # Z->inf limit of eq 6.5.45
s0_full = n - (nc + nl)
al_wrong = (e_wrong @ e_wrong) / s0_code
al_right = (e_right @ e_right) / s0_full
print(f"  Omega_0 with B*y omitted = {e_wrong @ e_wrong:10.1f}"
      f"   correct = {e_right @ e_right:10.1f}"
      f"   (x{(e_wrong @ e_wrong) / (e_right @ e_right):.1f})")
print(f"  alpha2_0: code formula = {al_wrong:.3f}   consistent = {al_right:.3f}"
      f"   (truth: 1.0)")
check("T3a consistent residual+projector gives alpha2 ~ 1 at truth (white noise)",
      abs(al_right - 1) < 0.15)
check("T3b omitting B*y biases alpha2 by > 2x on an nl>0 system",
      al_wrong / al_right > 2.0,
      "confirms adjustVCEjoint.m:283 and adjustVCE.m sweepStreamLags residual")

# --- WP-33 FIX (added 2026-08-18): joint projector via block identities ----
# R_bar = Sigma^-1 - Sigma^-1 [A B] inv(N_joint) [A B]' Sigma^-1 realised
# WITHOUT the joint inverse:  [A B] N^-1 [A B]' = Atil Sxx Atil' + B Nll^-1 B'
# with Atil = A - B G, G = Nll^-1 Nlc, Sxx = inv(Ncc_red) -- the same Lc/Ll/G
# blocks lsqvce.leverage uses. The kernel realises the two terms as
# Rt (shared probe on Atil_bar, stochastic) and K = W^-1 W^-T B Ll^-T (EXACT).
print("\n--- T3c/T3d: the WP-33 corrected joint projector -------------------------")
G1 = np.linalg.solve(B.T @ B, B.T @ A)
Atil = A - B @ G1
Sxx1 = np.linalg.inv(Ncc1)
hat_blocks = Atil @ Sxx1 @ Atil.T + B @ np.linalg.solve(B.T @ B, B.T)
hat_joint = M @ np.linalg.solve(Nf, M.T)
check("T3c block realisation == joint-inverse hat (no joint inverse formed)",
      np.allclose(hat_blocks, hat_joint, atol=1e-9),
      f"max|diff| = {np.max(np.abs(hat_blocks - hat_joint)):.2e}")
# s_0 of the fixed kernel (Z->inf limit): n - tr(hat) = n - (nc + nl); and the
# consistent pair (e_right, s0_full) is what the fixed drivers now form.
s0_fix = n - np.trace(hat_blocks)
al_fix = (e_right @ e_right) / s0_fix
check("T3d fixed kernel: s_0 -> n-(nc+nl) and alpha2 -> 1 at truth",
      abs(s0_fix - s0_full) < 1e-8 and abs(al_fix - 1) < 0.15,
      f"s0 = {s0_fix:.6f} (n-u_s = {s0_full}), alpha2 = {al_fix:.4f}")

# ---------------------------------------------------------------------------
# T4  levinsonDurbin + armaACV transcriptions vs. ground truth
# ---------------------------------------------------------------------------
print("\n--- T4: AR machinery ----------------------------------------------------")


def arma_acv(a, b, nlag):                   # armaACV.m transcription
    a = np.asarray(a, float) / a[0]
    b = np.asarray(b, float) / b[0]
    nh = max(4 * nlag, 512)
    for _ in range(8):
        h = np.zeros(nh)
        for t in range(nh):
            acc = b[t] if t < len(b) else 0.0
            for j in range(1, min(len(a), t + 1)):
                acc -= a[j] * h[t - j]
            h[t] = acc
        if abs(h[-1]) < np.finfo(float).eps * np.max(np.abs(h)) or nh >= 2**20:
            break
        nh *= 2
    nfft = 1 << int(np.ceil(np.log2(2 * nh)))
    R = np.fft.ifft(np.abs(np.fft.fft(h, nfft)) ** 2).real
    return R[:nlag]


def levinson(r, p):                          # levinsonDurbin.m transcription
    a = np.array([1.0])
    s2 = r[0]
    for k in range(1, p + 1):
        kappa = -(r[k] + sum(a[j] * r[k - j] for j in range(1, len(a)))) / s2
        if not np.isfinite(kappa) or abs(kappa) >= 1:
            return a, s2
        a = np.concatenate([a, [0.0]]) + kappa * np.concatenate([[0.0], a[::-1]])
        s2 *= (1 - kappa ** 2)
    return a, s2


a_true = np.array([1.0, -0.75, 0.4])         # stationary AR(2), unit innovation
r_th = arma_acv(a_true, [1.0], 12)
a_fit, s2_fit = levinson(r_th, 2)
check("T4a levinsonDurbin recovers the AR from its own theoretical ACV",
      np.allclose(a_fit, a_true, atol=1e-8) and abs(s2_fit - 1) < 1e-8,
      f"a = {np.round(a_fit, 6)}, s2 = {s2_fit:.8f}")
Tm = np.array([[r_th[abs(i - j)] for j in range(2)] for i in range(2)])
yw = np.linalg.solve(Tm, -r_th[1:3])
check("T4b Levinson agrees with the direct Yule-Walker solve",
      np.allclose(a_fit[1:], yw, atol=1e-8))
h_long = np.zeros(20000)
for t in range(20000):
    acc = 1.0 if t == 0 else 0.0
    for j in range(1, min(3, t + 1)):
        acc -= a_true[j] * h_long[t - j]
    h_long[t] = acc
r_dir = np.array([h_long[:20000 - k] @ h_long[k:] for k in range(12)])
check("T4c armaACV FFT autocorrelation == direct impulse autocorrelation",
      np.allclose(r_th, r_dir, rtol=1e-9))

# ---------------------------------------------------------------------------
# T5  dctI norms and psd2acv identity
# ---------------------------------------------------------------------------
print("\n--- T5: spectral primitives ----------------------------------------------")
N = 33
mI, nI = np.meshgrid(np.arange(N), np.arange(N), indexing="ij")
k_o = np.ones(N); k_o[0] = k_o[-1] = 1 / np.sqrt(2)
k_e = np.ones(N); k_e[0] = 1 / np.sqrt(2)
Xo = np.sqrt(2 / (N - 1)) * np.outer(k_o, k_o) * np.cos(mI * nI * np.pi / (N - 1))
Xe = np.sqrt(2 / (N - 1)) * np.outer(k_e, k_e) * np.cos(mI * nI * np.pi / (N - 1))
check("T5a dctI 'ortho' is self-inverse (X X = I)",
      np.allclose(Xo @ Xo, np.eye(N), atol=1e-12))
condXe = np.linalg.cond(Xe)
v = rng.normal(size=N)
check("T5b dctI 'ellmer' roundtrip via explicit inverse; cond ~ 1.5",
      np.allclose(np.linalg.solve(Xe, Xe @ v), v) and condXe < 2.0,
      f"cond = {condXe:.3f}")
f = np.linspace(0, 0.1, 4096)
S = 1 + np.sin(8 * f / 0.1) ** 2
r0 = np.trapezoid(S * np.cos(2 * np.pi * 0 * f), f)
check("T5c psd2acv identity r(0) == trapz(f,S)",
      abs(r0 - np.trapezoid(S, f)) < 1e-12)

# ---------------------------------------------------------------------------
# T6  FINDING F4: psdVCE 'abs' vs 'floor' at a zero-power bin
# ---------------------------------------------------------------------------
print("\n--- T6: psdVCE negative-amplitude policy (F4) ---------------------------")
Sprior, sd, floor_ = 1e-3, 4.0, 1e-12
alpha = rng.normal(0.0, sd, size=200000)     # noisy ratio around a zero-power bin
shat_abs = np.maximum(np.abs(alpha * Sprior), floor_)     # psdVCE.m line 169
shat_flr = np.maximum(alpha * Sprior, floor_)             # what the comment says
ratio = shat_abs.mean() / shat_flr.mean()
print(f"  E[Shat] abs-policy = {shat_abs.mean():.3e}   floor-policy = "
      f"{shat_flr.mean():.3e}   ratio = {ratio:.2f}")
check("T6  abs() injects ~2x the power a floor would at zero-power bins",
      1.8 < ratio < 2.2)

# ---------------------------------------------------------------------------
# T7  redundancyLags kernel: MC estimator vs exact projector
# ---------------------------------------------------------------------------
print("\n--- T7: per-lag kernel (eqs. 6.5.45-6.5.56) ------------------------------")
Nn, P, Z = 60, 4, 20000
c = 0.8 ** np.arange(Nn)                     # AR(1)-like ACV
Sig = np.array([[c[abs(i - j)] for j in range(Nn)] for i in range(Nn)])
W = np.linalg.cholesky(Sig).T                # upper, Sig = W'W
Aa = rng.normal(size=(Nn, P))
Ab = np.linalg.solve(W.T, Aa)
Nmat = Ab.T @ Ab
U = np.linalg.cholesky(Nmat).T
Zb = np.where(rng.random((P, Z)) < 0.5, -1.0, 1.0) / np.sqrt(Z)
Rt = np.linalg.solve(W, Ab @ np.linalg.solve(U, Zb))
X = np.linalg.solve(W, np.eye(Nn))
Rbar_mc = X @ X.T - Rt @ Rt.T
Rbar_ex = np.linalg.inv(Sig) - np.linalg.inv(Sig) @ Aa @ np.linalg.solve(
    Aa.T @ np.linalg.solve(Sig, Aa), Aa.T) @ np.linalg.inv(Sig)


def lag_sum(Rb, nlag):
    out = np.empty(nlag)
    out[0] = np.trace(Rb)
    for q in range(1, nlag):
        out[q] = 2 * np.trace(Rb, offset=q)
    return out


s_mc, s_ex = lag_sum(Rbar_mc, 8), lag_sum(Rbar_ex, 8)
rel = np.max(np.abs(s_mc - s_ex)) / abs(s_ex[0])
check("T7a MC per-lag redundancies match the exact projector (Z=20000)",
      rel < 0.02, f"max dev / s_0 = {rel:.4f}")
e = rng.normal(size=Nn)
et = np.linalg.solve(W, np.linalg.solve(W.T, e))   # Sigma^-1 e via two solves
om3 = 2 * (et[3:] @ et[:-3])
Vn = np.zeros((Nn, Nn)); Vn[np.arange(Nn - 3), np.arange(3, Nn)] = 1
Vn += Vn.T
check("T7b Omega_n unrolled sum == e' Sig^-T V_n Sig^-1 e",
      abs(om3 - et @ (Vn @ et)) < 1e-9)

# ---------------------------------------------------------------------------
print()
if FAIL:
    print("FAILURES:", ", ".join(FAIL))
    raise SystemExit(1)
print("ALL CHECKS BEHAVE AS THE R2 REVIEW CLAIMS (T1-T7; seed 20260816).")
