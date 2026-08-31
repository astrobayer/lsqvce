#!/usr/bin/env python3
"""Content of lsqvce_theory.pdf -- theory, practical considerations,
implementation and ERROR HANDLING of least-squares adjustment with variance
component estimation, tied to the decisions pinned in the lsqvce toolbox.

Source text: the R4 nine-expert review (HANDOVER §11, 2026-08-30) reconciled
against the tree; every claim about the toolbox names the work package that
pinned it. Built by docs_src/build_theory_pdf.py with the master-manual
pipeline (same styles, drawn equations, no external toolchain).

Created 2026-09-01 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Fable 5 -- R4 Level L3, WP-59.
"""
from build_master_pdf import P, B, N, code, callout, table, BODY_W
from pdf_graphics import eq
from reportlab.platypus import PageBreak, Spacer


def _sp(h=4):
    return Spacer(1, h)


def front():
    st = []
    st += [_sp(60)]
    st += [P('lsqvce', 'h1')]
    st += [P('Least-squares adjustment with variance component estimation', 'h2')]
    st += [P('Theory, practical considerations, implementation, error handling',
             'h3')]
    st += [_sp(14)]
    st += [P('This document collects the theory the toolbox rests on and the '
             'reasons for every place where the implementation deviates from '
             'the textbook. It is written for a reader who knows least squares '
             'and wants to know exactly what lsqvce computes, why, and where it '
             'refuses to compute. Each statement about the toolbox names the '
             'work package (WP) that pinned it, so the claim can be checked '
             'against HANDOVER.md rather than believed.', 'body')]
    st += [_sp(6)]
    st += [P('Sources. Koch (1999); Teunissen (2000); Förstner (1979); Koch and '
             'Kusche (2002); Kusche (2003); Teunissen and Amiri-Simkooei (2008); '
             'Amiri-Simkooei (2016); Ellmer (2018); Kvas and Mayer-Gürr (2019); '
             'Kvas et al. (2019); Golub and Van Loan (2013); Higham (2002). '
             'Full list in section 10.', 'body')]
    st += [_sp(6)]
    st += [P('Created 2026-09-01 (CEST, UTC+02:00) with assistance from Claude '
             '(Anthropic), model Claude Fable 5 -- R4 Level L3, WP-59. '
             'Status: supersedes the theory chapter of lsqvce_manual.pdf once '
             'decision R4-D3 is taken; until then both exist and this one is '
             'authoritative for the method choices.', 'bodyc')]
    st += [PageBreak()]
    return st


def s1_model():
    st = []
    st += [P('1&nbsp;&nbsp;The model and the two ways to solve it', 'h1')]
    st += [P('The linear Gauss-Markov model relates <i>n</i> reduced observations '
             '<b>l</b> (observed minus computed) to <i>u</i> unknowns <b>x</b>:',
             'body')]
    st += [eq([('b', 'l'), ('sp', ''), ('r', '='), ('sp', ''), ('b', 'A'), ('b', 'x'),
               ('sp', ''), ('r', '+'), ('sp', ''), ('b', 'e'), ('r', ','), ('sp', ''),
               ('sp', ''), ('r', 'E{'), ('b', 'e'), ('r', '} = 0,'), ('sp', ''),
               ('r', 'D{'), ('b', 'e'), ('r', '} = '), ('g', 'sigma'), ('sub', '0'),
               ('sup', '2'), ('b', 'Q'), ('r', ' = '), ('g', 'sigma'), ('sub', '0'),
               ('sup', '2'), ('b', 'P'), ('sup', '-1')], width=BODY_W, number='1.1')]
    st += [P('The best linear unbiased estimator is the weighted least-squares '
             'solution and its cofactor matrix', 'body')]
    st += [eq([('v', 'x'), ('hat', ''), ('sp', ''), ('r', '='), ('sp', ''),
               ('r', '('), ('b', 'A'), ('sup', 'T'), ('b', 'P'), ('b', 'A'), ('r', ')'),
               ('sup', '-1'), ('b', 'A'), ('sup', 'T'), ('b', 'P'), ('b', 'l'), ('r', ','),
               ('sp', ''), ('sp', ''), ('b', 'Q'), ('sub', 'xx'), ('r', ' = ('),
               ('b', 'A'), ('sup', 'T'), ('b', 'P'), ('b', 'A'), ('r', ')'), ('sup', '-1')],
              width=BODY_W, number='1.2')]
    st += [P('with residuals <b>v</b> = <b>l</b> &minus; <b>A</b><b>x</b>&#770; '
             'and the a-posteriori variance of unit weight', 'body')]
    st += [eq([('g', 'sigma'), ('hat', ''), ('sub', '0'), ('sup', '2'), ('sp', ''),
               ('r', '='), ('sp', ''), ('b', 'v'), ('sup', 'T'), ('b', 'P'), ('b', 'v'),
               ('r', ' / '), ('v', 'r'), ('r', ','), ('sp', ''), ('sp', ''),
               ('v', 'r'), ('r', ' = '), ('v', 'n'), ('r', ' \u2212 '), ('v', 'u')],
              width=BODY_W, number='1.3')]
    st += [P('Everything the toolbox reports about quality descends from (1.3): '
             'the redundancy <i>r</i> is the budget of degrees of freedom that '
             'the variance components will later divide among observation groups '
             '(section 3). Sign convention pinned in the toolbox: <b>l</b> is O '
             '&minus; C and <b>v</b> = <b>l</b> &minus; <b>A</b><b>x</b> &minus; '
             '<b>B</b><b>y</b> (residuals.m; DATA_CONTRACT §3).', 'body')]

    st += [P('1.1&nbsp;&nbsp;Normal equations versus QR', 'h2')]
    st += [P('Two backends solve (1.2). The normal-equation (NEQ) path forms '
             '<b>N</b> = <b>A</b><sup>T</sup><b>P</b><b>A</b> and factors it by '
             'Cholesky; cost about <i>nu</i><sup>2</sup>/2 + <i>u</i><sup>3</sup>/6, '
             'memory O(<i>u</i><sup>2</sup>), and -- decisively for satellite '
             'gravimetry -- the normal matrices of independent arcs and systems '
             'ADD, so accumulation is a reduction and combination is a sum. The '
             'price is the condition number:', 'body')]
    st += [eq([('g', 'kappa'), ('r', '('), ('b', 'A'), ('sup', 'T'), ('b', 'P'), ('b', 'A'),
               ('r', ') = '), ('g', 'kappa'), ('r', '('), ('b', 'P'), ('sup', '1/2'),
               ('b', 'A'), ('r', ')'), ('sup', '2')], width=BODY_W, number='1.4')]
    st += [P('Once &kappa;(<b>A</b>) exceeds roughly 10<sup>7</sup>-10<sup>8</sup>, '
             '<b>N</b> is numerically singular in IEEE double precision (Golub '
             'and Van Loan §5.3; Higham §20). The orthogonal (QR) path factors '
             '<b>P</b><sup>1/2</sup><b>A</b> = <b>Q</b><b>R</b> and solves '
             '<b>R</b><b>x</b> = <b>Q</b><sup>T</sup><b>P</b><sup>1/2</sup><b>l</b>; '
             'multiplication by an orthogonal matrix neither stretches vectors '
             'nor amplifies error, so the forward error grows with &kappa;(<b>A</b>) '
             'plus a residual-weighted &kappa;<sup>2</sup> term instead of '
             '&kappa;<sup>2</sup> outright. It costs about twice the flops and does '
             'not add across blocks.', 'body')]
    st += [callout('Pinned in lsqvce (WP-7, WP-49/50, R4 §11.2 item 2): BOTH backends '
                   'exist behind one interface (QRSystem, NormalEqSystem). The '
                   'engines default to "qr"; the streaming accumulator merges QR '
                   'factors of arcs and days ("mergeRaw"), so the arc-wise '
                   'accumulation the NEQ path is famous for is available under QR '
                   'as well. Changing the default is a decision, not an option '
                   '(Options.Backend documents it).', 'note')]

    st += [P('1.2&nbsp;&nbsp;Local unknowns and the block-arrow structure', 'h2')]
    st += [P('Satellite data carry unknowns that belong to one arc only -- '
             'initial states, accelerometer biases, empirical accelerations. With '
             '<b>B</b> the design of these <i>local</i> unknowns <b>y</b> the '
             'normal system of one observing system is block-arrow: a common '
             'block <b>N</b><sub>cc</sub>, a local block <b>N</b><sub>ll</sub> '
             'and the coupling <b>N</b><sub>lc</sub>. Pre-elimination by the '
             'Schur complement removes the local unknowns exactly:', 'body')]
    st += [eq([('b', 'N'), ('sub', 'cc'), ('sup', 'red'), ('r', ' = '), ('b', 'N'),
               ('sub', 'cc'), ('r', ' \u2212 '), ('b', 'N'), ('sub', 'lc'), ('sup', 'T'),
               ('b', 'N'), ('sub', 'll'), ('sup', '-1'), ('b', 'N'), ('sub', 'lc')],
              width=BODY_W, number='1.5')]
    st += [P('so the global system carries only the common unknowns; the local '
             'ones are recovered afterwards by back-substitution. Two consequences '
             'matter in practice. First, any regularisation of <b>N</b><sub>ll</sub> '
             'biases <b>x</b><sub>c</sub> through (1.5) -- damping the local block '
             'is never free (BlockAccumulator help, "REGULARISATION BIASES THE '
             'COMMON BLOCK"). Second, a system with no local unknowns is a valid '
             'system: the Schur step is the identity (G-07, decided 2026-08-31).',
             'body')]
    st += [P('1.3&nbsp;&nbsp;Nonlinear models', 'h2')]
    st += [P('Nonlinear observation equations are linearised at approximate '
             'values <b>x</b><sub>0</sub> and iterated (Gauss-Newton): the design '
             'becomes the Jacobian at the current iterate and <b>l</b> the '
             'reduced observation <b>l</b> &minus; <i>f</i>(<b>x</b>). The toolbox '
             'iterates with optional Levenberg damping and stops on the step or '
             'the cost (solveNonlinear; through fit as the batched kind with '
             '<i>x0</i>, WP-56b).', 'body')]
    return st


def s2_conditioning():
    st = []
    st += [P('2&nbsp;&nbsp;Ill-conditioning, rank deficiency, regularisation', 'h1')]
    st += [P('Spherical-harmonic gravity normal matrices are ill-conditioned by '
             'nature: signal power decays with degree, polar gaps and uneven '
             'sampling make high orders nearly unobservable, and the parameters '
             'differ by many orders of magnitude. Three tools address this.',
             'body')]
    st += [B('<b>Scaling.</b> Symmetric Jacobi scaling <b>D</b><b>N</b><b>D</b> with '
             '<b>D</b> = diag(<b>N</b>)<sup>-1/2</sup> equalises the diagonal and '
             'lowers &kappa; before any factorisation. Cheap, exact (the solution is '
             'rescaled back), always worth doing.')]
    st += [B('<b>Rank-revealing detection.</b> A Cholesky factorisation that fails, or '
             'a pivoted QR / eigenvalue ratio near machine precision, signals a '
             'datum defect or a broken stochastic model. The toolbox detects '
             'datum defects from the eigenvalue ratios of the reduced system and '
             'reports the redundancy accordingly (datum.m).')]
    st += [B('<b>Regularisation.</b> Tikhonov / Kaula regularisation adds a prior block, '
             '(<b>N</b> + &lambda;<b>K</b>)<b>x</b> = <b>n</b>. The weight &lambda; can be '
             'chosen by the L-curve, by generalised cross-validation (gcvLambda), or '
             '-- most consistently -- by treating the prior as one more observation '
             'group whose variance component the VCE estimates (Koch and Kusche 2002; '
             'Xu et al. 2006). This is WP-24 (CommonPrior, applyRegularisation): '
             'the prior enters as observations and section 3 weights it.')]
    st += [callout('Never add a silent diagonal "nugget" to make a factorisation succeed. '
                   'A regularisation that was not declared is a result that cannot be '
                   'interpreted. lsqvce records regUsed / lambdaUsed in every system '
                   'diagnostic and refuses to hide the choice.', 'warn')]
    st += [P('2.1&nbsp;&nbsp;Cholesky failure is a signal', 'h2')]
    st += [P('MATLAB\'s two-output form [R, p] = chol(N) returns p &gt; 0 instead of '
             'throwing. That is the right primitive: the failure is information '
             '(rank deficiency, an indefinite covariance, a mis-scaled group), to be '
             'reported and, where a fallback exists, handled by an orthogonal '
             'factorisation -- not masked. Section 7 returns to this under error '
             'handling: what lsqvce does when the boundary itself is fuzzy (a to-eps '
             'singular matrix whose Cholesky "succeeds" on one BLAS and fails on '
             'another).', 'body')]
    return st


def s3_vce():
    st = []
    st += [P('3&nbsp;&nbsp;Variance component estimation', 'h1')]
    st += [P('Write the dispersion of the observations as a linear combination of '
             'known cofactor matrices with unknown variance components:', 'body')]
    st += [eq([('r', 'D{'), ('b', 'l'), ('r', '} = '), ('g', 'Sigma'), ('sub', 'k'),
               ('sp', ''), ('g', 'sigma'), ('sub', 'k'), ('sup', '2'), ('sp', ''),
               ('b', 'Q'), ('sub', 'k')], width=BODY_W, number='3.1')]
    st += [P('The components &sigma;<sub>k</sub><sup>2</sup> are what VCE estimates. '
             'Three rigorous families exist; the toolbox implements the first and '
             'its Monte-Carlo form.', 'body')]
    st += [P('3.1&nbsp;&nbsp;Helmert and the redundancy-contribution form', 'h2')]
    st += [P('Helmert\'s (1907) estimator requires the traces of products of '
             'cofactor matrices with the inverse normal matrix. Förstner (1979) and '
             'Koch reduce it, for uncorrelated groups, to the iterated form', 'body')]
    st += [eq([('g', 'sigma'), ('hat', ''), ('sub', 'k'), ('sup', '2'), ('r', ' = '),
               ('b', 'v'), ('sub', 'k'), ('sup', 'T'), ('b', 'P'), ('sub', 'k'), ('b', 'v'),
               ('sub', 'k'), ('r', ' / '), ('v', 'r'), ('sub', 'k'), ('r', ','), ('sp', ''),
               ('sp', ''), ('v', 'r'), ('sub', 'k'), ('r', ' = tr('), ('b', 'P'), ('sub', 'k'),
               ('b', 'Q'), ('sub', 'vv,k'), ('r', '),'), ('sp', ''), ('g', 'Sigma'),
               ('sub', 'k'), ('v', 'r'), ('sub', 'k'), ('r', ' = '), ('v', 'n'),
               ('r', ' \u2212 '), ('v', 'u')], width=BODY_W, number='3.2')]
    st += [P('<i>r</i><sub>k</sub> is the partial redundancy of group <i>k</i>: the '
             'share of the degrees of freedom that group contributes. Iterate until '
             'every group\'s a-posteriori factor is one; iterated Helmert coincides '
             'with the iterated almost-unbiased estimator and, under normality, with '
             'REML. This is lsqvce.engine.vce (WP-12), with the per-observation redundancy '
             'numbers (redundancy.m) supplying <i>r</i><sub>k</sub> exactly when the '
             'problem is small enough to afford them.', 'body')]
    st += [P('3.2&nbsp;&nbsp;MINQUE and LS-VCE', 'h2')]
    st += [P('Rao\'s MINQUE (1971) is non-iterative and distribution-free but may '
             'return negative components and depends on the a-priori values. '
             'LS-VCE (Teunissen and Amiri-Simkooei 2008) casts the estimation of the '
             'components as a least-squares problem of its own: it delivers the '
             'precision of the components, admits covariance components and a-priori '
             'information, and coincides with BIQUE under normality. lsqvce does not '
             'implement LS-VCE; where the components\' own uncertainty is needed it '
             'uses the Fisher information of the Helmert iteration (WP-23). This is a '
             'documented limitation (section 9), not an oversight.', 'body')]
    st += [P('3.3&nbsp;&nbsp;Monte-Carlo trace estimation', 'h2')]
    st += [P('At GRACE scale <i>r</i><sub>k</sub> in (3.2) needs the inverse of the '
             'full normal matrix -- infeasible. Koch and Kusche (2002) and Kusche '
             '(2003) replace the exact trace by Hutchinson\'s stochastic estimator: '
             'for probe vectors <b>z</b> with independent zero-mean unit-variance '
             'entries (Rademacher &plusmn;1),', 'body')]
    st += [eq([('r', 'E{'), ('b', 'z'), ('sup', 'T'), ('b', 'M'), ('b', 'z'), ('r', '} = tr('),
               ('b', 'M'), ('r', '),'), ('sp', ''), ('sp', ''), ('r', 'Var '), ('g', 'prop'),
               ('r', ' 1/'), ('v', 'm'), ('r', ' for '), ('v', 'm'), ('r', ' probes')],
              width=BODY_W, number='3.3')]
    st += [P('Each probe costs one extra solve of the normal system; a few dozen '
             'probes give traces good to a percent. lsqvce implements this as '
             'vce \'hutchinson\' / vceHutchinson with ONE probe generator and an '
             'honest RNG contract (WP-36: the seed is recorded in Results.meta), and '
             'estimates the REDUNDANCY, not the trace of the wrong matrix (WP-32 fixed '
             'exactly that). The probe count and the estimator\'s spread are '
             'user-visible (Options.Probes, Options.Seed).', 'body')]
    st += [P('3.4&nbsp;&nbsp;Degeneracy and non-negativity', 'h2')]
    st += [P('A component driven towards zero means the group carries no '
             'independent information (<i>r</i><sub>k</sub> &rarr; 0); a negative '
             'estimate means sampling noise or a mis-specified model. Options are to '
             'fix the component, merge groups, or impose non-negativity '
             '(Amiri-Simkooei 2016). lsqvce refuses to estimate a group below a '
             'minimum redundancy (Options.MinRedundancy, default 50) and reports the '
             'per-pass history rather than clamping silently; the negative-factor '
             'clamp and the PSD floor of WP-35/40 are declared policies with '
             'warning identifiers, never quiet corrections.', 'body')]
    return st


def s4_stochastic():
    st = []
    st += [P('4&nbsp;&nbsp;Stochastic modelling of correlated satellite data', 'h1')]
    st += [P('GRACE-FO observables -- K-band range-rate, kinematic positions -- are '
             'strongly correlated in time. lsqvce\'s model, per observation group, is '
             'separable (adjustVCE help "THE MODEL"):', 'body')]
    st += [eq([('b', 'C'), ('sub', 'g'), ('r', ' = '), ('g', 'sigma'), ('sub', 'g'),
               ('sup', '2'), ('sp', ''), ('b', 'D'), ('sub', 'g'), ('r', ' ('), ('b', 'R'),
               ('sub', 't,g'), ('sp', ''), ('g', 'otimes'), ('sp', ''), ('b', 'I'),
               ('sub', 'd'), ('r', ') '), ('b', 'D'), ('sub', 'g')], width=BODY_W,
              number='4.1')]
    st += [P('<b>D</b><sub>g</sub> is a FIXED per-epoch scaling supplied by the loader '
             '(the kinematic-orbit epoch covariance, say); <b>R</b><sub>t,g</sub> the '
             'temporal correlation, the same on each of the <i>d</i> axes, ESTIMATED; '
             '&sigma;<sub>g</sub><sup>2</sup> the variance component of section 3. '
             'No cross-axis correlation is what makes the whitener separable and '
             'cheap (separableWhiten).', 'body')]
    st += [P('4.1&nbsp;&nbsp;Whitening by filters, and the one thing to get right', 'h2')]
    st += [P('The temporal law is applied as a whitening filter -- an AR/ARMA filter '
             'fitted to the residual autocovariance (residualAcv, selectArAcv, '
             'levinsonDurbin), or the exact Toeplitz factor built from an estimated '
             'PSD (psd2acv, toeplitzCofactor). A filter restarts its state at every '
             'batch, so the first <i>p</i> rows of a batch are transient. Overlapping '
             'batches by the filter memory <i>p</i> and discarding the overlap makes '
             'the batched result identical to filtering the whole arc; measured '
             '(BlockAccumulator help): 1.6e-15 relative error with the overlap, '
             '1.3e-3 without -- a factor 10<sup>12</sup> that fails silently. The '
             'overlap is therefore automatic (BatchOverlap = p) and not a knob.',
             'body')]
    st += [P('4.2&nbsp;&nbsp;The two loops', 'h2')]
    st += [P('The SEPARATE loop (Options.Model = "vce") alternates: solve under the '
             'current covariance; re-estimate <b>R</b><sub>t,g</sub> from whitened '
             'residuals; estimate &sigma;<sub>g</sub><sup>2</sup> by (3.2); push both '
             'back into the loaders as a covariance override; repeat until every '
             'stream is white (Ljung-Box, Options.WhiteAlpha) and every factor is one '
             '(Options.VceTol, OuterTol). The POOLED loop ("pooled") estimates one '
             'spectrum over all streams per pass. The JOINT loop ("ellmer") is Ellmer\'s '
             '(2018, §6.5) ten-step co-estimation of PSD amplitudes and keyed per-arc '
             'variance factors with Hutchinson probes -- the ITSG-Grace2018 scheme '
             '(Kvas et al. 2019) -- as adjustVCEjoint (WP-13/33/34, through fit since '
             'WP-56c). The two can converge to DIFFERENT models because amplitudes and '
             'scales multiply the same covariance; the toolbox documents this coupling '
             'and the production guidance of WP-40 (open loop or at most three damped '
             'passes) instead of pretending one answer.', 'body')]
    st += [P('4.3&nbsp;&nbsp;Non-stationary and background-model noise', 'h2')]
    st += [P('Ellmer (2018) propagates the star-camera orientation covariance into '
             'the K-band antenna-offset correction, giving a non-stationary noise '
             'model recomputed every month and separating stationary ranging noise '
             'from the non-stationary correction -- 27-46 % lower noise than the '
             'official centres\' solutions at 300 km smoothing (Kvas et al. 2019). '
             'Kvas and Mayer-Gürr (2019) derive stationary covariance models for the '
             'atmosphere-ocean background-model error and weight them by VCE. In '
             'lsqvce the per-epoch scaling <b>D</b><sub>g</sub> is the hook for the '
             'first (epochCovariance in the data contract), a covariance group with '
             'its own component the hook for the second.', 'body')]
    return st


def s5_robust():
    st = []
    st += [P('5&nbsp;&nbsp;Outliers: screening and robust estimation', 'h1')]
    st += [P('Baarda\'s data snooping tests the normalised residual '
             '<i>w</i><sub>i</sub> = <i>v</i><sub>i</sub> / (&sigma;<sub>0</sub> '
             '&radic;<i>q</i><sub>vv,ii</sub>) against N(0,1) and removes or '
             'downweights the largest exceedance; its reliability theory gives the '
             'minimal detectable bias (MDB), and Förstner (1983) extends it to the '
             'Gauss-Markov model. lsqvce\'s diagnose computes <i>w</i>, the redundancy '
             'numbers <i>r</i><sub>i</sub>, MDB and reliability per observation; snoop '
             'screens block-natively (WP-15); report turns them into a five-level '
             'verdict. Two measured traps are documented in README: Pope\'s &tau; goes '
             'blind under contamination, and the cure has its own cliff -- which is '
             'why the verdict never rests on |<i>w</i>| or <i>r</i> alone.', 'body')]
    st += [P('M-estimation replaces rejection by iteratively reweighted least '
             'squares with an influence function -- Huber (quadratic core, linear '
             'tails), IGG3 (keep / downweight / reject zones), bisquare (redescending). '
             'This is lsqvce.engine.robust (WP-25), reached through fit with Options.Robust = '
             '"irls" on the in-memory kinds (WP-56b). Robustness and VCE interact: '
             'outliers inflate the components; screen first, then estimate, or use '
             'the robust arc-factor mean of the joint loop (Options.MeanType).', 'body')]
    return st


def s6_gravity():
    st = []
    st += [P('6&nbsp;&nbsp;Spherical-harmonic gravity recovery', 'h1')]
    st += [P('The disturbing potential is expanded in fully normalised spherical '
             'harmonics up to degree and order <i>L</i><sub>max</sub>; the common '
             'unknowns are the Stokes coefficients. The design matrix holds the partial '
             'derivatives of the observable -- range-rate, positions, potential '
             'differences -- with respect to each coefficient, i.e. a spherical-harmonic '
             'synthesis of the functional. Production processing (ITSG: Kvas et al. '
             '2019) uses variational equations over daily arcs with arc-local nuisance '
             'parameters -- exactly the block-arrow structure of section 1.2 -- and '
             'stacks the arcs.', 'body')]
    st += [P('lsqvce carries a physically structured closed-loop oracle for this: '
             'shDesign (WP-48) generates deterministic SH partials with realistic '
             'column correlation, and demo_gravity (WP-51/53) recovers planted '
             'coefficients from a GRACE-like, ragged, three-stream fixture with '
             'measured numbers pinned (&sigma;<sub>0</sub> 0.9998, stream ratios '
             '0.996 / 1.000 / 1.019, planted anomaly recovered at 1.33 &sigma;). That '
             'oracle, not a Gaussian design, is the reference authority for the '
             'estimator (HANDOVER §6.1.4).', 'body')]
    return st


def s7_errors():
    st = []
    st += [P('7&nbsp;&nbsp;Error handling', 'h1')]
    st += [P('Two kinds of error meet in an adjustment toolbox: numerical -- '
             'floating-point arithmetic destroying information -- and software -- a '
             'caller violating a contract. lsqvce treats both under one principle:',
             'body')]
    st += [callout('Guards are not fixes. A guard that turns a wrong number into a '
                   'plausible one has hidden a defect; a guard that turns it into a '
                   'refusal has reported one. lsqvce refuses.', 'warn')]
    st += [P('7.1&nbsp;&nbsp;Catastrophic cancellation in v<sup>T</sup>Pv', 'h2')]
    st += [P('The quadratic form can be obtained without re-forming residuals from '
             'quantities the accumulation already holds:', 'body')]
    st += [eq([('b', 'v'), ('sup', 'T'), ('b', 'P'), ('b', 'v'), ('r', ' = '), ('b', 'l'),
               ('sup', 'T'), ('b', 'P'), ('b', 'l'), ('r', ' \u2212 '), ('v', 'x'),
               ('hat', ''), ('sup', 'T'), ('b', 'A'), ('sup', 'T'), ('b', 'P'), ('b', 'l')],
              width=BODY_W, number='7.1')]
    st += [P('When the fit is good the two terms are large and nearly equal and '
             'their difference loses most of its digits -- it can even come out '
             'slightly negative (Higham 2002, ch. 1). Three responses exist. Clamping '
             'to zero fabricates a plausible variance; iterative refinement of the '
             'solution masks the conditioning that caused the loss; detecting the '
             'loss and returning NaN reports it. lsqvce does the third on BOTH '
             'backends (WP-49 normal path: the authoritative v<sup>T</sup>Pv comes '
             'from the residual pass, and cancellation yields NaN under ONE identifier '
             'per namespace; WP-50 QR path: the merge-fold cancellation that produced '
             '&sigma;<sub>0</sub> = 0 was found and is guarded by a measured three-anchor '
             'threshold with a rotation backstop). A NaN &sigma;<sub>0</sub> is a '
             'result: it says the problem\'s conditioning exceeded what double '
             'precision can carry through this formula.', 'body')]
    st += [P('7.2&nbsp;&nbsp;Factorisation at the singular boundary', 'h2')]
    st += [P('Whether dpotrf reports failure on a to-eps-singular matrix depends on '
             'the BLAS/LAPACK build: the trailing pivot rounds to &plusmn;tiny. A test '
             'or a code path that asserts a factoriser\'s verdict AT the boundary is '
             'therefore asserting the library, not the mathematics. Pinned after a '
             'live catch (R4 §11.8, 2026-08-30): assert the eigenvalue or rank '
             'statement, keep factoriser assertions to clearly definite regimes, and '
             'treat regime-dependent behaviour as such (the WP-47 REGIME doctrine).',
             'body')]
    st += [P('7.3&nbsp;&nbsp;Contracts fail fast and name the rule', 'h2')]
    st += [P('Every public function validates its inputs at the signature with typed '
             'arguments blocks (WP-39); cross-argument rules are checked in the body '
             'under lsqvce:&lt;unit&gt;:&lt;rule&gt; identifiers, one identifier per rule '
             '(the b1-r3 rule). fit validates the whole data declaration BEFORE any '
             'day is read, so a typo costs a second, not a pass over a month '
             '(DATA_CONTRACT §3.5). Unknown struct fields are errors on purpose. '
             'Errors raised while reading a day carry the loader\'s identifiers with '
             'the user\'s own error attached as the cause.', 'body')]
    st += [P('7.4&nbsp;&nbsp;Warnings are inventory, not noise', 'h2')]
    st += [P('The CI gate refuses a run whose warning inventory gained an identifier '
             'nobody declared (WP-47 rule 7). A test that provokes a warning on '
             'purpose mutes exactly that identifier for its own duration. The '
             'consequence for the toolbox: every warning identifier is a documented '
             'event with a documented cause, or it does not exist.', 'body')]
    st += [P('7.5&nbsp;&nbsp;Determinism and provenance', 'h2')]
    st += [P('The residual pass re-reads the data; a loader that is not deterministic '
             'puts its non-determinism into &sigma;<sub>0</sub>. Every Results carries '
             'an immutable meta record -- toolbox version, MATLAB release, UTC '
             'timestamp, every RNG seed resolved, the resolved option set (WP-44) -- '
             'so a file can answer "which run produced this" a year later.', 'body')]
    return st


def s8_implementation():
    st = []
    st += [P('8&nbsp;&nbsp;Implementation: the workflow in lsqvce', 'h1')]
    st += [P('The user surface is four names (R4, decisions R4-D1/D4): a data '
             'struct describing the problem (DATA_CONTRACT.md), lsqvce.Options for '
             'what to do, lsqvce.fit as the one verb, lsqvce.Results with report(). '
             'Below it the engines and kernels of the tier census '
             '(docs_src/tiers.csv) do the work:', 'body')]
    st += [N('<b>Declare.</b> data.common states the common unknowns ONCE; '
             'data.systems(k) points to buildDay -- the one function the user writes, '
             'returning A, B, l (O &minus; C), t for a day, raw and unwhitened.', 1)]
    st += [N('<b>Build loaders.</b> fit wraps each buildDay through makeLoader: the '
             'meta probe (which arcs a day carries, per day -- ragged schedules are '
             'native, WP-52), arc slicing, a single-day cache, the whitener from the '
             'covariance spec.', 2)]
    st += [N('<b>Accumulate.</b> BlockAccumulator visits system, day, observation type, '
             'arc; whitens each arc with the batch overlap of §4.1; accumulates into '
             'the (nc+nl) system; after a system\'s days pre-eliminates its local '
             'block by (1.5). Parallel axes that pay: days within a system '
             '(DistributeDays, WP-62), arcs in the redundancy pass, systems. Never '
             'within an arc (filter state), never inside the dense solve (BLAS).', 3)]
    st += [N('<b>Solve.</b> Cholesky or QR on the reduced common system; datum defect '
             'detection; declared regularisation.', 4)]
    st += [N('<b>Residual pass.</b> The data is read a second time (determinism!); '
             'residuals, whitened residuals, the authoritative v<sup>T</sup>Pv (§7.1), '
             'redundancy numbers exact or by Hutchinson probes.', 5)]
    st += [N('<b>Iterate.</b> Model "vce" / "pooled" / "ellmer" re-estimate the '
             'stochastic model and the components and push them back into the '
             'loaders as a covariance override (opt.covOverride.bySys, resolved by the '
             'ONE merge applyCovOverride -- G-06), until white and balanced.', 6)]
    st += [N('<b>Judge.</b> diagnose and report: per-observation w, r, MDB; the '
             'five-level verdict; coefficients() with confidence bounds; everything '
             'reachable on the Results for experts (R.acc, R.systems, R.noise).', 7)]
    st += [P('Cost model. Every pass reads each day twice; the VCE loop multiplies by '
             'the pass count. At global SH scale the redundancy / VCE data pass is '
             'where the hours go (about 4 h per iteration at degree 90 serial; 30-40 '
             'min with 16 workers, BlockAccumulator help) -- build the workflow around '
             'parallelising THAT axis.', 'body')]
    st += [table([['Kind of data', 'Fields present', 'Engine', 'Options.Model'],
                  ['Streaming (GRACE-FO)', 'systems, common, days', 'adjustBlock / adjustVCE / adjustVCEjoint', 'gls, vce, pooled, ellmer'],
                  ['In-memory matrices', 'A, l, (p)', 'adjust / robust', 'gls (+Robust)'],
                  ['Batched callback', 'batchFcn, numUnknowns, numBatches, (x0)', 'adjust / robust / solveNonlinear', 'gls (+Robust)']],
                 widths=[95, 130, 150, 75])]
    return st


def s9_limits():
    st = []
    st += [P('9&nbsp;&nbsp;Limitations of the method and of the implementation', 'h1')]
    st += [B('Helmert / redundancy-contribution VCE estimates variance components of '
             'UNCORRELATED groups. Covariance components between groups (LS-VCE) are '
             'not implemented.')]
    st += [B('The separable model (4.1) assumes one temporal law shared by the axes of a '
             'stream and no cross-axis correlation; a vector-AR model (varWhiten, whittle) '
             'exists for the multi-axis case but is not the default.')]
    st += [B('Monte-Carlo redundancy is unbiased but noisy: its spread scales as '
             '1/&radic;m in the probe count. Components stable to a percent need tens of '
             'probes; the seed is recorded so a run can be reproduced exactly.')]
    st += [B('The joint (Ellmer) and the separate loop can converge to different '
             'stochastic models; the alternation can be bistable (F-15). Production '
             'guidance -- open loop or at most three damped passes -- is documented; '
             'the monotone / trust-region cure is a backlog item (WP-42).')]
    st += [B('Cancellation in v<sup>T</sup>Pv is REPORTED (NaN), not repaired; a problem '
             'whose conditioning exceeds double precision through this formula needs '
             'rescaling or the QR path, which the diagnostics recommend but do not '
             'switch to silently.')]
    st += [B('A NaN &sigma;<sub>0</sub>, a group below the redundancy floor, a '
             'non-converged loop are results with warning identifiers -- the toolbox '
             'never returns a plausible substitute in their place.')]
    return st


def s10_refs():
    st = []
    st += [P('10&nbsp;&nbsp;References', 'h1')]
    refs = [
        'Amiri-Simkooei, A. R. (2016). Non-negative least-squares variance component estimation with application to GPS time series. J Geod 90(5):451-466.',
        'Baarda, W. (1968). A testing procedure for use in geodetic networks. Netherlands Geodetic Commission, Publ. on Geodesy 2(5).',
        'Ellmer, M. (2018). Contributions to GRACE gravity field recovery: improvements in dynamic orbit integration, stochastic modelling of the antenna offset correction, and co-estimation of satellite orientations. Dissertation, TU Graz.',
        'Förstner, W. (1979). Ein Verfahren zur Schätzung von Varianz- und Kovarianzkomponenten. Allg. Vermessungsnachrichten 86:446-453.',
        'Förstner, W. (1983). Reliability and discernability of extended Gauss-Markov models. DGK Reihe A 98:79-103.',
        'Golub, G. H., Van Loan, C. F. (2013). Matrix Computations, 4th ed. Johns Hopkins University Press.',
        'Helmert, F. R. (1907). Die Ausgleichungsrechnung nach der Methode der kleinsten Quadrate, 2nd ed. Teubner.',
        'Higham, N. J. (2002). Accuracy and Stability of Numerical Algorithms, 2nd ed. SIAM.',
        'Hutchinson, M. F. (1989). A stochastic estimator of the trace of the influence matrix for Laplacian smoothing splines. Commun. Statist. Simul. 18:1059-1076.',
        'Koch, K.-R. (1999). Parameter Estimation and Hypothesis Testing in Linear Models, 2nd ed. Springer.',
        'Koch, K.-R., Kusche, J. (2002). Regularization of geopotential determination from satellite data by variance components. J Geod 76(5):259-268.',
        'Kusche, J. (2003). A Monte-Carlo technique for weight estimation in satellite geodesy. J Geod 76(11-12):641-652.',
        'Kvas, A. (2019). GRACE gravity field recovery with background model uncertainties. Dissertation, TU Graz; Kvas, A., Mayer-Gürr, T. (2019). J Geod 93(12):2543-2552.',
        'Kvas, A., Behzadpour, S., Ellmer, M., Klinger, B., Strasser, S., Zehentner, N., Mayer-Gürr, T. (2019). ITSG-Grace2018: Overview and evaluation of a new GRACE-only gravity field time series. J Geophys Res Solid Earth 124(8):9332-9344.',
        'Rao, C. R. (1971). Estimation of variance and covariance components -- MINQUE theory. J Multivariate Anal 1:257-275.',
        'Teunissen, P. J. G. (2000). Adjustment Theory: an Introduction. Delft University Press.',
        'Teunissen, P. J. G., Amiri-Simkooei, A. R. (2008). Least-squares variance component estimation. J Geod 82(2):65-82.',
        'Xu, P., Shen, Y., Fukuda, Y., Liu, Y. (2006). Variance component estimation in linear inverse ill-posed models. J Geod 80(1):69-81.',
    ]
    for r in refs:
        st += [P(r, 'bodyc')]
    return st


def theory_story():
    st = []
    st += front()
    st += s1_model()
    st += s2_conditioning()
    st += s3_vce()
    st += s4_stochastic()
    st += s5_robust()
    st += s6_gravity()
    st += s7_errors()
    st += s8_implementation()
    st += s9_limits()
    st += s10_refs()
    return st
