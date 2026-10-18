function out = adjustVCEjoint(accFcn, opt)
%LSQVCE.ENGINE.ADJUSTVCEJOINT  Joint (co-estimated) stochastic-model loop.
%
%   out = LSQVCE.ENGINE.ADJUSTVCEJOINT(accFcn)
%   out = LSQVCE.ENGINE.ADJUSTVCEJOINT(accFcn, 'Name', value, ...)
%
%   Runs Ellmer's (2018, §6.5.7) ten-step algorithm, in which ALL variance
%   components are updated TOGETHER from one residual pass and one shared probe:
%   the per-bin spectral amplitudes S_j (the temporal law), the arc-wise
%   variance factors sigma2_m (the within-group variance factors), and -- when a prior is
%   present -- the geophysical factor sigma2_gpm.
%
%   JOINT, NOT ALTERNATING -- WHY
%     The retired separable driver alternated: fit the temporal whitening, then
%     the between-group variance components, each with its own criterion.
%     Ellmer instead treats every PSD frequency and every arc factor as variance
%     components of ONE Helmert VCE and updates them SIMULTANEOUSLY (step 9).
%     This matters: alternating (block-coordinate) and joint updates share fixed
%     points only if the blocks are SEPARABLE, and here they are NOT -- the arc
%     factor and the frequency amplitudes both scale the same Sigma_ll (eqs.
%     6.5.7, 6.5.18), which is exactly why the gamma normalisation (eq. 6.5.28)
%     exists to prevent double counting. The two schemes can therefore converge
%     to DIFFERENT Sigma. This driver is the joint scheme, and since R6 the
%     only stochastic-model estimator of the toolbox (Estimate = "noise").
%
%   THE TEN STEPS (Ellmer §6.5.7), and where each lives
%     1-3  build Delta_l, A; form Sigma from the current model; decorrelate
%            -> the caller's accFcn does this (see INPUT below)
%     4-5  accumulate N, n; factor N = U'U and solve
%            -> acc.solve() and acc.commonFactor() (the backend's own
%               square root, the prior folded in -- WP-134; chol of the
%               squared normals before)
%     6    draw Zbar and form U^-1 Zbar ONCE for all arcs
%            -> LSQVCE.ENGINE.SHAREDPROBE
%     7    decorrelated residuals ebar = W^-T (l - A x - B y), the FULL
%          Table-6.2 vector (WP-33; the local fit is subtracted per system)
%            -> re-streamed per arc via the store's dataFcn
%     8    Rtilde, etilde, and from them Omega and s per lag
%            -> LSQVCE.KERNEL.REDUNDANCYLAGS
%     9    updated S_j, sigma2_m (and sigma2_gpm)
%            -> LSQVCE.ENGINE.PSDVCE and LSQVCE.ENGINE.ARCVARIANCEFACTORS
%     10   if not converged, rebuild Sigma and go to step 2
%
%   THE LAW: POSITIVE COMPONENTS (R7 WP-111 part 3, user decision 2026-10-11)
%     The temporal law is a NON-NEGATIVE spectrum S on L = NumLags bins,
%     piecewise linear in frequency, with the exact covariance r = K S at
%     every lag of every arc (LSQVCE.KERNEL.HATACV). Each bin is a variance
%     component whose matrix is the ACV of a hat spectrum -- positive
%     semi-definite -- so its Helmert factor is >= 0 by construction and the
%     law needs neither a floor nor a zero-frequency rule. Ellmer's
%     components were the columns of the inverse type-I DCT of a truncated
%     ACF: indefinite, so a bin's factor could go negative, a floored bin made
%     the declared covariance near-singular, and at phi 0.8 / 0.9 the level ran
%     away x5e7 (review R7, R7-02; HANDOVER 12.44). Ellmer's white start (eq.
%     6.5.62) is a regular start here: S = ones is white noise, r = [1; 0; ...]
%     (his degenerate fixed point was a property of the DCT weighting).
%
%   INPUT
%     accFcn  function_handle  @(model) -> a SOLVED-READY lsqvce.engine.BlockAccumulator
%             built under the given stochastic model. It is called once per outer
%             iteration and must apply the model when whitening (steps 1-3).
%             `model` is a struct with fields:
%               .S       L x 1     the current law: the bin spectrum, scaled
%                                  so that r(0) = 1 (the level lives in the
%                                  factors); .Sgroups one per pool group
%               .C       L x 1     its covariance at lags 0..L-1, K S (for the
%                                  record; .Cgroups likewise)
%               .sigma2  M x 1     the current arc variance factors
%               .iter    1x1       the iteration number (1 on the first call)
%               .final   1x1 logical  true on the closing accumulation under the
%                                  returned model (same .iter as the last pass;
%                                  WP-135: .lastSol is then THAT pass's solve;
%                                  WP-144: .lastSol.Ninv, the commons' cofactor
%                                  N^-1, when .s0sq is not finite)
%               .streamSigma2, .arcFactor  the per-stream level and the per-arc
%                                  factor relative to it (R7-13, since 6.1.6);
%                                  LSQVCE.ENGINE.JOINTMODELOPT turns the whole
%                                  model into the loader's covOverride
%             On the FIRST call the caller may ignore the model and use its own
%             a-priori weighting; the driver then estimates from there. On later
%             calls an accFcn that drops .arcFactor weights an arc-factor stream
%             by a level the update does not assume: its factors then compound
%             (the class of finding F-01), and arcFactorsDiverging says so.
%
%   OPTIONS (Name/Value)
%     'NumLags'   1x1 int    L, the number of spectral bins (default 60)
%     'S0'        L x 1 | [] starting spectrum (>= 0, not all zero). Default
%                            ones: white noise (Ellmer's eq. 6.5.62)
%     'Sigma2_0'  M x 1 | [] starting arc factors (default all ones, eq. 6.5.62),
%                            aligned to out.arcKeys; they weight pass 1 per arc
%                            (WP-132), at the price of one extra accumulation
%                            that enumerates the arcs first; finite and > 0
%                            (WP-134: a 0 reached the clamp, a NaN stayed)
%     'MaxIter'   1x1 int    outer iterations (default 10)
%     'Tol'       1x1        the stop rule's tolerance (default 1e-2): every
%                            spectral bin j moved by at most Tol max(1,
%                            sqrt(200 / r_j)) -- Tol for a bin of redundancy
%                            r_j >= 200, growing with the bin's own standard
%                            error sqrt(2 / r_j) below that (at 1e-2 a tenth
%                            of it; WP-147) -- and every arc factor and
%                            PriorDev by at most Tol
%     'Probes'    1x1 int    Z Hutchinson probes (default 100, Ellmer's value)
%     'Seed'      1x1 | []   RNG seed for reproducible probes (default []: one
%                            seed is pinned for the run and returned, out.seed)
%     'Influence' string     "none" (default) | "huber" | "igg3" | "bisquare"
%              R6 WP-105: the ROBUST PASS INSIDE THE LOOP, in the KOCH
%              convention -- the weight is a statement about the
%              OBSERVATION. Every pass's sweep forms ONE statistic per epoch
%              (below: Baarda's w since WP-156), which goes through
%              LSQVCE.KERNEL.ROBUSTWEIGHTS per stream, and the factors
%              enter the NEXT pass as an additive NUGGET on the temporal
%              covariance (model.rowWeights -> UserOpt.epochWeights ->
%              covModel's cm.nugget, q = 1/w - 1), i.e. BEFORE the temporal
%              law: the epoch's variance becomes Sigma_ee / w, its
%              covariances with the neighbours stay the law's (R7 WP-115
%              part 3b, user decision (b) on R7-10: the former FOLD into Cep,
%              sigma / sqrt(w), scaled those too, and a zero-weight epoch
%              stayed a phantom observation under a correlated law --
%              formal/true sd 0.31 at phi 0.9, measured). MEASURED
%              (rounds 4-8): weighting the law-whitened INNOVATION instead
%              (Yang 1994, adjustRobust's convention) neutralises the
%              filter's whole tail around every blunder and the colour of a
%              blunder-laden stream crawled 0.006 -> 0.054 in 16 passes;
%              with the same blunders declared through Cep the plain
%              estimator recovered it. NOT Ellmer (2018): he screens outliers in
%              preprocessing and down-weights whole ARCS through the variance
%              factors; this is Koch's iterated re-weighting inside a VCE.
%              THE LEVEL (R7-12, measured at tier A, HANDOVER 12.43): with
%              the factors inside the component update the level of a CLEAN
%              stream came back 14-19 % LOW (20 seeds, IGG3) -- the influence
%              selects the epochs with large residuals, which no redundancy
%              count undoes. So (user decision 2026-10-11) the component
%              update (the sweep) sees only the REJECTED epochs, factor 0, as
%              nuggets (deleted); partial factors weight the solution only.
%              THE STATISTIC AND THE TRUNCATION (WP-156, design F, user
%              decision 2026-10-18; HANDOVER 12.91). Deleting epoch i removes
%              exactly w_i^2 = (R v)_i^2 / R_ii from v'Pv (R = P Qvv P), Baarda's
%              w for CORRELATED observations; the raw residual rejected on
%              until then correlates with it only by sqrt(c_i / Sigma_ii), so
%              on a coloured stream the kept epochs were truncated in the
%              wrong direction (review 12.84 N3: the level 1-2.5 % high) and a
%              red stream's blunders went unseen (2 % x 2 sigma at phi 0.98:
%              35 % caught, the level 3.6 x). Now: (1) the statistic is w per
%              epoch on the arc's law (LSQVCE.KERNEL.WSNOOP; a d-axis epoch
%              the sum of its axes' w^2, mapped to a normal deviate as R7-08),
%              in units of ONE SCALE PER STREAM, the previous pass's MAD of the
%              final system (never the estimated level: perStreamWeights says
%              why); (2) every arc is snooped -- the largest |w| beyond the
%              rejection point deleted, the arc recomputed, until none is left
%              -- rebuilt from the data's own exclusions every pass (pass 0,
%              whose scale is not known, rejects simultaneously); (3) the kept
%              epochs' truncation enters the REDUNDANCY lag by lag
%              (LSQVCE.KERNEL.TRUNCATIONDEFICIT: s_n - (1 - tau) [R D^-1 R]_n,
%              and the whole law's pair), so the law, its level and sigma0
%              are consistent for any law -- the level is no longer divided
%              by tau, which is exact for the total only (a red stream's
%              variance read +2.4 %); (4) a rejected set that cycles with
%              period 2 is frozen at the union of the two (freezeCycle;
%              out.robust.frozen). tau = F_{d+2}(q_k) / F_d(q_k), P(chi2_d >
%              q_k) = P(|N(0,1)| > kReject), is 0.9733 for IGG3's k1 = 3 at d
%              = 1 and 1 for Huber, which rejects nothing.
%              The RETURNED solution is hard: after the loop one
%              accumulation with the rejected epochs deleted and every
%              partial factor 1, its sigma0 / Cxx corrected for the rejected
%              rows and the truncation -- exact for the total under w
%              (out.robust.truncation; measured: with the partial factors in
%              that solve, sigma*sigma0 stayed 0.83). MEASURED (HANDOVER
%              12.91, 20 seeds paired against "none"): a red stream's variance
%              -0.24 / -0.10 / -0.11 % at phi 0.5 / 0.9 / 0.98 (the raw rule
%              +0.3 / +1.4 / +1.0 %), sigma0^2 -0.08 %, x unchanged; 2 % x 2
%              sigma blunders at phi 0.9: phi_b 0.900 (the raw rule 0.843).
%              Cost: one explicit law inverse per arc (O(m^2) for a Toeplitz
%              law, LSQVCE.KERNEL.LAWINVERSE) and the deficit's lag sums
%              (O(m^2 log m)), against the sweep's O(m^3) factor. A HELD arc
%              keeps its factor and the bins' stop test reads the change of
%              shape (WP-152). Prior streams are never re-weighted. Prefer
%              "igg3" / "bisquare" for gross blunders: Huber keeps c/|u| of
%              their weight. "none" is bitwise the plain estimator.
%     'Tuning'    [] | 1xK   tuning constants of the influence (its defaults)
%     'Verbose'   logical    per-iteration banner (default false)
%     'StepSchedule' 'none' (default) | 'decay'.  WP-40 r2 (F-15): the
%              closed-loop alternation is a stochastic approximation; the
%              plain step is BISTABLE for short arcs (recovered-PSD corr can
%              collapse to negative values on unlucky seeds). 'decay' uses
%              the Robbins-Monro step Gamma_it = 2/(2+it) via PSDVCE's
%              'Gamma'. HONEST CLAIM (psd2acv-faithful mirror, WP-40 r6):
%              damping is DAMAGE LIMITATION -- it prevents the undamped
%              collapse (corr can hit ~0.2) but the closed loop does not
%              beat OPEN-LOOP pass-1 quality on the tested fixtures;
%              production guidance is open loop or <= 3 damped passes.
%              'none' is bitwise the old behaviour.
%     'AvgTail' 1x1 integer >= 0, default 0.  Polyak averaging of the
%              RETURNED spectrum over the last AvgTail passes (log domain) --
%              every pool group's law (out.Sgroups / Cgroups, which the final
%              solve reads), not only out.S (WP-134);
%              (honest r6 claim: helps the damped loop settle; cannot
%              rescue an undamped collapse, and longer tails average over
%              late-pass drift -- keep it <= 4) and does nothing useful
%              without it (mirror-measured). The iteration itself always
%              runs on the raw S.
%     'Pool'     "stream" (default; "all" until WP-157) | "obstype" | "all": which streams SHARE
%                one estimated law (one psdVCE per pool group, out.Sgroups)
%     'Estimate' [] | containers.Map 'satKey|obsType' -> "none" | "scale" |
%                "law": the per-stream policy ([] = "law" for every stream)
%     'ArcFactors' [] | containers.Map 'satKey|obsType' -> logical: whether
%                the stream's arcs carry factors ([] = true for every stream)
%     'OnPass'   [] | @(ctx)  the per-pass inspection hook (lsqvce.util.onPass:
%                return false to stop after this pass)           default []
%     'OnPassContext' struct  fields merged into every ctx (lsqvce.fit
%                passes .data and .processingOptions)            default struct()
%     'SkipThreshold' 1x1 in [0, 1]  the largest fraction of ONE STREAM's arcs a
%                pass may skip (gap / nonPD) before it is refused (default 0.2;
%                WP-134: it was the rate over all arcs, so a small stream could
%                be skipped whole among large ones, silently)
%     'PriorDev' [] | @() [dev, pending]  components the accFcn updates itself
%                -- lsqvce.fit's common prior groups with estimate "scale"
%                (Koch-Kusche, from model.lastSol): dev the largest relative
%                variance change of this pass's update, pending true while a
%                group is in the solve but not yet updated. The loop converges
%                only when dev <= Tol and nothing is pending (WP-134, review 3
%                finding 1: the prior was outside the stop rule, and a loop that
%                settled on pass 1 never estimated it).   default []
%
%   OUTPUT
%     out  1x1 struct
%       .S         L x 1      the converged law: the bin spectrum on
%                             f_j = j / (2(L-1)) cycles per sample, r(0) = 1
%                             (the 'all' group, else the first pool group)
%       .C         L x 1      its covariance at lags 0..L-1 (K S)
%       .Sgroups   containers.Map  pool group -> its S; .Cgroups -> its C
%       .sigma2    M x 1      the converged arc variance factors, ALIGNED TO
%                             .arcKeys (the pass-1 enumeration order; WP-34)
%       .checks    struct     the model checks of the last pass (WP-154):
%                             .stream / .group / .whitePval / .dispersion /
%                             .dispersionBound per
%                             stream (the Ljung-Box p under its law; the
%                             scatter of its arcs' levels in units of their
%                             sampling variance, ~1 when one level fits), and
%                             .lawGroup / .corrTime / .arcMedian per law (the
%                             correlation time -1/log(rho_1) and the median arc,
%                             in epochs); the warnings below read them
%       .seed      1x1        the probe seed used (the given Seed, or the one
%                             pinned for an unseeded run -- R7-37; passing it
%                             back as 'Seed' replays the run)
%       .robust    struct     .influence, .scale (1 x S, the scale of the w
%                             statistic per stream, w units: the MAD of the
%                             last robust sweep's final system; WP-156),
%                             .streams (1 x S 'satKey|obsType', sorted), .factor
%                             (nEpoch x 1, the final factors, ONE PER EPOCH
%                             in sweep order -- a multi-axis epoch's d rows
%                             share it), .keys (A x 4 {sysID, obsType, day,
%                             arc}), .truncation (1 x S, the truncation factor
%                             tau per stream, which the hard sigma0 is divided
%                             by -- the level no longer is, WP-156), .frozen
%                             (the pass a cycling rejected set was frozen at,
%                             0 none; WP-156); empty under "none" (R6 WP-105)
%       .arcKeys   M x 1 cell arc keys s<sys>_<obsType>_d<day>_a<arc> pinned
%                             on pass 1 (LSQVCE.ENGINE.REDUNDANCY's 'arc' format) --
%                             the ONE identity .sigma2 and Sigma2_0 align to
%       .arcStreams M x 1 cell each arc's stream key 'satKey|obsType'
%                             (LSQVCE.UTIL.STREAMKEY), aligned to .arcKeys --
%                             carried, so nobody parses an arc key (R7-43)
%       .skips     S x 3 cell {pass, key, reason} for every arc skipped in
%                             any pass ('gap' | 'nonPD'); a skipped arc KEEPS
%                             its previous factor
%       .acc       1x1        the final accumulator
%       .sol       1x1        the final solution struct
%       .converged 1x1 logical whether every bin settled (binRatio <= 1),
%                             and every factor's change and PriorDev's
%                             change fell below Tol
%       .nIter     1x1        outer iterations actually run
%       .history   struct     .maxDev (nIter x 1, max abs(alpha2 / abar - 1),
%                             the change of the law's shape, psdVCE: a
%                             reading), .binRatio (nIter x 1, the stop
%                             quantity: max over bins of the bin's change
%                             over Tol max(1, sqrt(200 / r_j)); <= 1 settled),
%                             .xStep (nIter x 1, max abs change of x over its
%                             standard error since the pass before; NaN on
%                             pass 1), .S (nIter x L), .levelDev, .priorDev
%                             (nIter x 1, 0 without PriorDev)
%
%   EXAMPLE
%     accFcn = @(model) buildAccumulatorWithModel(cfg, model);
%     r = lsqvce.engine.adjustVCEjoint(accFcn, 'NumLags', 60, 'MaxIter', 8);
%     semilogy(r.S);             % the estimated noise spectrum (the bins)
%     plot(r.sigma2);            % the arc-wise weights
%
%   ACCURACY
%     The update is Ellmer's PSD-VCE law (HANDOVER 6.1.4, the Ellmer (2018)
%     row: test_ellmer's references, and its recorded refutation of the
%     closed-loop alternation, F-15). The per-lag sums it consumes
%     (LSQVCE.KERNEL.REDUNDANCYLAGS) equal the dense trace(Rbar V_n) to 1e-8
%     relative, the nugget route the general one to 1e-9 of s_0
%     (test_ellmer). The levels and laws are statistical estimates: their
%     measured sampling error and bias are stated with the options above
%     (R7-12, HANDOVER 12.43). WP-138.
%
%   ERRORS
%     lsqvce:adjustVCEjoint:needCfcn  the joint mode needs a cross-covariance
%                                     function and none was supplied
%     lsqvce:adjustVCEjoint:badS0     the starting spectrum is malformed
%     lsqvce:adjustVCEjoint:notPD     the reduced common normal matrix (the
%                                     prior included) has no factor: singular to
%                                     the QR solve's own pivot tolerance, or not
%                                     positive definite on the normal backend
%     lsqvce:adjustVCEjoint:skipThreshold  a pass skipped more than
%                                   SkipThreshold of one stream's arcs (R5 M1;
%                                   per stream since WP-134)
%     lsqvce:adjustVCEjoint:badAxisRows  a multi-axis arc whose row count is not a
%                                   multiple of its axis count (R6 WP-90: an arc
%                                   carrying a per-epoch Cep, scalar or multi-axis,
%                                   is descaled per epoch before its lag statistics;
%                                   the former noCep refusal is gone -- the loader
%                                   composes D = I and says so once)
%     lsqvce:adjustVCEjoint:noArcs    no arc produced statistics in a pass: none
%                                     was enumerated, or every one was skipped
%                                     or inactive under SkipThreshold = 1
%     lsqvce:adjustVCEjoint:collapsed the influence function left fewer live
%                                     EPOCHS than parameters + 1: a bodily
%                                     residual offset or near-total
%                                     contamination, NOT 50 % (R6 WP-105;
%                                     measured in HANDOVER 12.34)
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:robust:atBreakdown                  [warning] the robust pass zeroed > 45 % of
%                                                the epochs (once per run; R7 WP-115)
%     lsqvce:adjustVCEjoint:localNotPD           System X
%     lsqvce:adjustVCEjoint:notConverged         [warning] The joint stochastic model did not settle in X outer iterations (last...
%     lsqvce:adjustVCEjoint:arcFactorsDiverging  [warning] arc-factor spread grew (F-01)
%     lsqvce:adjustVCEjoint:poolNotWhite         [warning] a pool group of >= 2 streams with a
%                                                member not white under its law (WP-154, M2)
%     lsqvce:adjustVCEjoint:levelVaries          [warning] a stream without arc factors whose
%                                                arcs' levels scatter beyond their sampling
%                                                error (WP-154, M3)
%     lsqvce:adjustVCEjoint:memoryVsArc          [warning] a law's correlation time above
%                                                half its median arc (WP-154, M4)
%
%   WHAT eq. 6.5.28 MULTIPLIES, AND WHY IT MATTERS (R6 WP-98, finding F-01)
%     gamma2_m = alpha2_m * sigma2_m is correct only when sigma2_m is the factor
%     the WEIGHTING used. Since 6.1.6 (R7-13) that IS the per-arc factor for a
%     stream with arc factors: LSQVCE.ENGINE.JOINTMODELOPT carries it to the
%     loader (covOverride.byArc), so the arcs are weighted the way Ellmer
%     weights them. Until then only the per-STREAM level reached the loader,
%     and multiplying by the per-arc factor anyway made the update a geometric
%     recursion with NO FIXED POINT, because alpha2_m cancels its own arc's
%     factor exactly -- the arc's numerator and denominator carry the same s2_m -- and
%     so never responds to it. MEASURED on a fixture with planted factors
%     1 / 9 / 1: the ratio ran 9.4, 79, 671, 5626, 4.7e4 over five passes, a
%     constant 8.4x per pass; it now settles at 9.4, 8.6, 8.6, 8.6, 8.6.
%     For a stream with arcFactors = false the two numbers are identical (every
%     arc carries the one stream level), so those runs are unchanged.
%     The factors are ABSOLUTE (gamma2_m itself, not gamma2_m / gammaBar^2):
%     the law is a shape here, so the level lives in the factors (R7-09 and
%     WP-116 part 2, HANDOVER 12.40). One arc 10x in sigma among nine:
%     factors 100.0 / 1.00 against planted 100 / 1, and the true scatter of x
%     13-47 % smaller than under the stream-mean weighting (16 seeds).
%
%   See also LSQVCE.ENGINE.PSDVCE, LSQVCE.ENGINE.ARCVARIANCEFACTORS,
%            LSQVCE.KERNEL.REDUNDANCYLAGS, LSQVCE.ENGINE.SHAREDPROBE
%
%   Updated 2026-10-18 (CEST, UTC+02:00), WP-156 (design F, user decision
%   2026-10-18): the robust pass rejects on Baarda's w per epoch, formed and
%   snooped in the sweep (the separate residual pass is gone), the kept
%   epochs' truncation enters the redundancy lag by lag (no division of the
%   level by tau), a cycling rejected set is frozen at its union
%   (history.robustFrozen, out.robust.frozen); N3 retired. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-18 (CEST, UTC+02:00), WP-157: Pool defaults to "stream", as
%   lsqvce.Options does. With assistance from Claude (Anthropic), model Claude
%   Opus 5.5.
%   Updated 2026-10-18 (CEST, UTC+02:00), WP-154 (review 12.84, M2-M4): out.checks
%   and three warnings for the assumptions the data can contradict -- one law
%   for a pool group's spectra, one level for a stream's arcs, a memory short
%   against the arcs (modelChecks); the sweep keeps the arcs' raw lag
%   statistics (projectLags). With assistance from Claude (Anthropic), model
%   Claude Opus 5.5.
%   Updated 2026-10-18 (CEST, UTC+02:00), WP-152 (review 12.84, N1/N2/N4): the
%   bins' stop test reads the change of shape (a robust run's every bin
%   settles at tau, so IGG3 never converged on clean data); a held arc is not
%   divided by tau (it grew 2.7 % per pass); the sweep deletes an arc's
%   EXCLUDED rows by the nugget it deletes rejected epochs with, per axis (a
%   'diag' stream with 10 % excluded rows read its level 5 % low). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-17 (CEST, UTC+02:00), WP-149: a robust hard solution
%   carries its truncation tau and dof also when the backend refused v'Pv, so
%   the residual pass that replaces the NaN (lsqvce.engine.verdict.supersede)
%   divides by tau. With assistance from Claude (Anthropic), model Claude
%   Opus 5.5.
%   Updated 2026-10-17 (CEST, UTC+02:00), WP-144 (A6): model.lastSol carries
%   Ninv (the commons' cofactor) when the backend refused a cancelled v'Pv,
%   so a prior's "scale" step does not divide by NaN. With assistance from
%   Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-16 (CEST, UTC+02:00), WP-147: a spectral bin has settled
%   when it moved by at most Tol max(1, sqrt(200 / r_j)) -- at the default
%   Tol a tenth of its own standard error where that exceeds Tol;
%   .history carries binRatio (the stop quantity) and xStep, and the
%   not-converged warning names the bins that still move and the solution's
%   last step. With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-138 (code review of WP-135): the
%   prior pseudo-streams are recognised by lsqvce.util.priorStream on their
%   obsType (the stream key was parsed for '|prior_'); ACCURACY. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-14 (CEST, UTC+02:00), WP-135 (review 3, lane 4): the closing
%   accumulation's model carries .final = true, so an accFcn that updates its
%   own components (fit's prior groups) knows .lastSol is the last pass's own
%   solve. With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-134 (third review, findings 1, 7,
%   15 and lows): PriorDev puts the accFcn's own components (fit's prior
%   groups) into the stop rule; the shared probe uses acc.commonFactor(),
%   not chol of the squared normals (notPD on cond(A) ~ 3e8 the QR solve
%   accepts); SkipThreshold is per stream and checked before noArcs; the
%   robust sweep leaves the prior pseudo-arcs out (they counted as live
%   epochs); Sigma2_0 values validated; AvgTail averages every pool group.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-133 (review 3 finding 18): the
%   hard solution's rejected rows leave sol.dof through the loader's
%   arc.excluded, and sol.truncation carries tau, so lsqvce.engine.residuals
%   cross-checks against the same number (sigma0Mismatch fired on every
%   igg3 / bisquare noise run). With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-132: Sigma2_0 weights pass 1 per
%   arc (the arcs are enumerated and pinned first); the arc factor's Helmert
%   pair takes the whole law as its component in closed form (redundancyLags
%   .omegaLaw / .sLaw), no longer r'Omega / r's, which assumed it Toeplitz.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-127: a local prior arc before its
%   group's from pass is inactive, not a skip; the per-pass record has one
%   column per pinned stream (a stream absent in one pass broke it). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-111 part 3 (R7-01, R7-02,
%   R7-20; user decision 2026-10-11): the law is a positive spectrum S on L
%   bins (state .S / .Sgroups, start 'S0' = white); every bin's component is
%   PSD and exact at every lag (sweepArcs projects the arc's full lag sums on
%   LSQVCE.KERNEL.HATACV, scaled by the arc's own variance), the arc factor is
%   the Helmert factor of the arc's whole covariance; 'C0' and 'Norm' are
%   retired with the DCT. With assistance from Claude (Anthropic), model
%   Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-115 part 4 (R7-43): stream
%   and arc keys from lsqvce.util.streamKey / arcKey; the stream of an arc is
%   carried (out.arcStreams), the regex that parsed it back is gone. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-115 part 3b (R7-10, R7-12;
%   user decisions (b) and 2026-10-11): the robust factors are an additive
%   nugget (no longer folded into Cep); the component update sees only the
%   rejected epochs and the level is divided by the rejection point's
%   Gaussian truncation factor (out.robust.truncation). With assistance from
%   Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-118 part 2 (R7-38): the per-
%   stream Ljung-Box statistic divides by the pooled pair counts N_k;
%   Noise.MeanType / 'MeanType' retired (user decision 2026-10-10). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-118 part 1 (R7-34, R7-35,
%   R7-37): the separable-driver history rewritten; Pool / Estimate /
%   ArcFactors / OnPass / OnPassContext documented; an unseeded run pins and
%   returns its probe seed (out.seed). With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-116 part 2 (R7-13): the per-arc
%   factors weight the solve (model.arcFactor -> jointModelOpt byArc) and eq.
%   6.5.28 multiplies by them; the factors are the absolute gamma2_m -- the
%   6.1.5 per-stream form counted the level twice on a stream whose arcs
%   differ. With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-116 part 1: gamma-bar per stream with
%   the stream's absolute level (R7-09); a stop rule that also needs the levels
%   and arc factors to settle, and the returned solution accumulated under the
%   returned model (R7-11); 'Smooth' / 'FreezeRedundancy' reach psdVCE (R7-22). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-115: the collapsed refusal
%   and its ERRORS line no longer claim "beyond the breakdown point" (R7-25:
%   the claim WP-110 retired from adjustRobust survived here, line-wrapped
%   past audit check 28) and count EPOCHS, which is what the guard compares;
%   the "level is UNBIASED" claim is withdrawn to an open question (R7-12).
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-09-09 (CEST, UTC+02:00), R6 WP-105: the robust pass inside
%   the loop ('Influence', 'Tuning', model.rowWeights, out.robust, the
%   collapsed refusal) -- the one-loop composition of IRLS and the joint
%   estimator. With assistance from Claude (Anthropic), model Claude Fable 5.1.
%   Updated 2026-09-08 (CEST, UTC+02:00), R6 WP-98 (finding F-01 FIXED): eq.
%   6.5.28 now multiplies alpha2_m by the factor the WEIGHTING used -- the
%   per-stream level -- instead of the per-arc factor, which alpha2 cancels
%   and therefore never responded to. The arc-factor ratio settles instead of
%   growing 8.4x a pass; runs with arcFactors = false are unchanged, because
%   there the two numbers are identical. With assistance from Claude
%   (Anthropic), model Claude Opus 5.
%   Updated 2026-09-08 (CEST, UTC+02:00), R6 WP-97 (finding F-01): the loop
%   records the SPREAD of every arc-factor stream pass by pass and warns
%   (arcFactorsDiverging) when one keeps growing -- the per-arc factors are
%   estimated but never fed back into the weighting, so they compound. The
%   guard reports; it does not correct, because the mechanism is not yet
%   established. With assistance from Claude (Anthropic), model Claude Opus 5.
%   Updated 2026-09-06 (CEST, UTC+02:00), R6 WP-90: the per-axis descale runs
%   for every arc that carries a per-epoch Cep (scalar streams included); the
%   noCep refusal is retired. With assistance from Claude (Anthropic), model
%   Claude Fable 5.1.
%   Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (a) (F-07): typed
%   `arguments` signature replaces inputParser. Option names were already
%   UpperCamelCase and are UNCHANGED. NOTE (premise corrected round a1-2,
%   measured on CI): MATLAB matches comma-syntax Name-Value NAMES
%   case-insensitively even under arguments blocks, so legacy casings
%   keep working; what tightened are the typed constraints, unknown-name
%   rejection, and enum VALUE casing (StepSchedule accepts exactly
%   "none"/"decay" -- strcmpi used to take 'Decay'). With assistance from
%   Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-40 r2 (F-15): optional
%   stabilised alternation -- 'StepSchedule','decay' (Robbins-Monro damped
%   psdVCE update) and 'AvgTail' (Polyak tail average of the returned
%   spectrum); defaults are bitwise the old path. With assistance from
%   Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-34 (fixes F-02, CRITICAL): arc
%   factors are matched by KEY (s<sys>_<obsType>_d<day>_a<arc>, pinned on
%   the pass-1 enumeration), never by position; skipped arcs keep their
%   previous factor and are counted in out.skips; a changed key set raises
%   arcMismatch. New outputs .arcKeys/.skips; Sigma2_0 aligns to .arcKeys.
%   With assistance from Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-18 (CEST, UTC+02:00), WP-33 (fixes F-03, CRITICAL): step
%   7 subtracts the per-system local fit B*y (recovered once per pass as
%   redundancy does) and the per-lag projector absorbs the local block
%   exactly -- shared implementation jointArcParts + redundancyLags
%   'LocalBbar'/'LocalLl'; evidence v10review.py T3. With assistance from
%   Claude (Anthropic), model Claude Fable 5.
%     lsqvce:adjustVCEjoint:arcMismatch  the arc KEY SET changed between
%                                        passes (the loader must be
%                                        deterministic; a covariance failure
%                                        is a SKIP, a vanished enumeration
%                                        entry is an identity error) (WP-34)
%     lsqvce:adjustVCEjoint:badSigma2_0  Sigma2_0 does not match the pass-1
%                                        arc count (align to out.arcKeys), or
%                                        holds a value not finite and > 0
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    % WP-39 stage (a): typed `arguments` signature (rule 5) replaces
    % inputParser. Names were already UpperCamelCase here; NV name
    % matching stays case-insensitive (MATLAB semantics, measured round
    % a1-2) -- the typed constraints and enum-value casing are what
    % tightened. Release note: HANDOVER §9 WP-39.
arguments
    accFcn (1,1) function_handle
    opt.NumLags  (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(opt.NumLags, 2)} = 60
    opt.S0       double = []
    opt.Sigma2_0 double = []
    opt.MaxIter  (1,1) double {mustBeInteger, mustBePositive} = 10
    opt.Tol      (1,1) double {mustBePositive} = 1e-2
    opt.Probes   (1,1) double {mustBeInteger, mustBePositive} = 100
    opt.Seed     double {mustBeScalarOrEmpty} = []
    opt.StepSchedule (1,1) string {mustBeMember(opt.StepSchedule, ["none","decay"])} = "none"
    opt.AvgTail  (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    opt.SkipThreshold (1,1) double {mustBeInRange(opt.SkipThreshold, 0, 1)} = 0.2
    opt.Smooth   (1,1) double {mustBeInteger, mustBeNonnegative} = 0      % R7-22: psdVCE's moving average
    opt.FreezeRedundancy (1,1) double {mustBeNonnegative} = 0             % R7-22: psdVCE's freeze
    % R5 M3 3d: per-stream policy (WP-83, decision R5-D4). Streams are keyed
    % 'satKey|obsType' (matlab.lang.makeValidName of the system id).
    opt.Pool (1,1) string {mustBeMember(opt.Pool, ["all","obsType","stream"])} = "stream"   % WP-157
    opt.Estimate = []            % containers.Map key -> "none"|"scale"|"law"; [] = "law" for all
    opt.ArcFactors = []          % containers.Map key -> logical;              [] = true for all
    % R6 WP-105: the ROBUST pass inside the loop. "none" = the plain estimator
    % (bitwise the old behaviour). Otherwise every pass ends with one residual
    % sweep whose LAW-whitened residuals give the next pass's per-observation
    % factors (lsqvce.kernel.robustWeights); the factors ride inside every
    % arc's whitener (model.rowWeights -> the factory's UserOpt.rowWeights),
    % so the accumulation, the probe redundancies and the component update
    % of the next pass all see ONE re-weighted problem. Prior streams are
    % never re-weighted (a down-weighted datum is a rank defect, not an
    % outlier).
    opt.Influence (1,1) string {mustBeMember(opt.Influence, ["none","huber","igg3","bisquare"])} = "none"
    opt.Tuning double = []       % tuning constants of the influence ([] = its defaults)
    opt.OnPass = []              % R5 M5 5c': function handle fired after EVERY pass
                                 % with the ctx below (see lsqvce.util.onPass); return
                                 % false to stop the loop (stopReason "userStop")
    opt.OnPassContext struct = struct()   % extra ctx fields (fit: data, processingOptions)
    opt.PriorDev = []            % WP-134: @() [dev, pending] of the accFcn's own components
    opt.Verbose  (1,1) logical = false
end
    o = opt;
    % R7-37: an unseeded run PINS one probe seed up front (the package-private
    % probe stream; the global RNG is untouched), as LSQVCE.ENGINE.VCE does. Every
    % pass then probes with the same Z, and out.seed is the number that
    % replays the run -- before, each pass drew a fresh private seed and none
    % was recoverable, against the "seeds in provenance" contract.
    if isempty(o.Seed)
        [~, o.Seed] = probeMatrix(1, 1, 'rademacher', []);
    end
    v = logical(o.Verbose);

    nLag = o.NumLags;
    Kc = lsqvce.kernel.hatACV(nLag, nLag);   % the law's ACV at lags 0..L-1 (the record, C)

    % ---- starting model (step 1 / eq. 6.5.62): white unless S0 says else ---
    S = ones(nLag, 1);
    if ~isempty(o.S0), S = o.S0(:); end
    if numel(S) ~= nLag || any(~isfinite(S)) || any(S < 0) || ~any(S > 0)
        error('lsqvce:adjustVCEjoint:badS0', ...
             'S0 must be %d finite non-negative bin amplitudes, not all zero (got %d entries).', nLag, numel(S));
    end
    S = S / (Kc(1, :) * S);  C = Kc * S;     % r(0) = 1: the law is a shape
    if ~isempty(o.Sigma2_0) && ~all(isfinite(o.Sigma2_0(:)) & o.Sigma2_0(:) > 0)
        % WP-134: a 0 skipped the arc in arcFactorMap but reached the clamp
        % (1e8 times the weight next pass); a NaN stayed NaN in out.sigma2
        error('lsqvce:adjustVCEjoint:badSigma2_0', ...
              'Sigma2_0 must hold finite factors > 0 (got %s).', mat2str(o.Sigma2_0(:).', 3));
    end
    sigma2 = o.Sigma2_0(:);              % may be empty; sized on first sweep

    hist = struct('maxDev', [], 'S', [], 'skips', {cell(0,3)}, ...
                  'streams', {{}}, 'sigma', [], 'white', [], 'whitePval', [], 'levelDev', [], ...
                  'priorDev', [], 'binRatio', [], 'xStep', []);
    xPrev = [];                          % WP-147: x of the pass before, for .xStep
    converged = false;  userStop = false;
    accF = []; solF = [];
    Smap = containers.Map('KeyType', 'char', 'ValueType', 'any');   % R5 M3 3d: the law per pool group
    Cmap = containers.Map('KeyType', 'char', 'ValueType', 'any');   % ... and its covariance (derived)
    pinnedKeys = {};                     % WP-34: pass-1 arc identity
    pinnedStream = {};                   % R6 WP-97: their stream keys, once
    pinnedRow = {};                      % R7-13: their loader keys (rowWeightKey)
    s2Map = [];                          % keyed factors; built on pass 1
    spreadHist = containers.Map('KeyType', 'char', 'ValueType', 'any');   % R6 WP-97 (F-01)

    shatTail = {};                               % WP-40 r2: AvgTail buffer
    gTail = {};                                  % WP-134: ... and every pool group's law
    lastSol = [];                                % R5 M7 7b: the previous pass's solution
    robust = ~strcmp(o.Influence, "none");       % R6 WP-105
    W = containers.Map('KeyType', 'char', 'ValueType', 'any');   % per-epoch factors
    fac = [];  rkeys = cell(0, 4);  rscale = [];  rstreams = {};  tauMap = [];  hardFix = [];
    hist.robustScale = [];  hist.nDown = [];  hist.nZero = [];  hist.robustFrozen = 0;
    rob = robustState(robust, o);                % WP-156: the w statistic's scale, the freeze
    if ~isempty(o.Sigma2_0)
        % WP-132 (review low 22): Sigma2_0 weights pass 1 PER ARC. Pass 1 used to
        % run before anything was pinned, so jointModelOpt saw neither stream
        % levels nor arc factors and weighted every arc with the MEAN of
        % Sigma2_0 -- while eq. 6.5.28 multiplied each arc's own start factor:
        % a non-uniform start computed pass 1 wrong, and a warm start at the
        % fixed point left it. The arcs are now pinned on an enumeration first
        % (the loader's 'meta' calls, which need the stores of one accumulation
        % under the old pass-1 model).
        accE = accFcn(struct('S', S, 'C', C, 'sigma2', sigma2, 'iter', 1, 'streamSigma2', [], ...
                             'Sgroups', Smap, 'Cgroups', Cmap, 'Pool', char(o.Pool), 'Estimate', o.Estimate, ...
                             'lastSol', [], 'rowWeights', W));
        [pinnedKeys, pinnedRow, pinnedStream] = enumerateArcs(accE);
        s2Map = pinFactors(pinnedKeys, o.Sigma2_0);
        accE = [];                                   %#ok<NASGU> release its stores
    end
    if robust
        % PASS 0 (round 10, measured): the first robust sweep runs BEFORE the
        % first spectral update. Started without it, pass 1 fitted the law to
        % the blunders (a near-white, negative-lag shape) and the
        % multiplicative PSD update crawled out of it at ~2 % a pass for 16
        % passes (0.006 -> 0.08); the control that declared the blunders from
        % the start converged. So: one accumulation under the start model,
        % one residual sweep, the factors -- and only then the loop.
        model0 = struct('S', S, 'C', C, 'sigma2', sigma2, 'iter', 1, ...   % iter 1: priors with from = 1 active
                        'streamSigma2', streamMeans(pinnedStream, sigma2), ...
                        'arcFactor', arcFactorMap(pinnedStream, pinnedRow, sigma2, o), ...   % WP-132: Sigma2_0
                        'Sgroups', Smap, 'Cgroups', Cmap, 'Pool', char(o.Pool), 'Estimate', o.Estimate, ...
                        'lastSol', [], 'rowWeights', W);
        [W, tauMap, hardFix, rob] = robustPassZero(accFcn(model0), o, nLag, v, rob, W, model0.streamSigma2);
    end
    for it = 1:o.MaxIter
        % ---- steps 1-5: build under the current model, accumulate, solve --
        % R5 M3 3b: the caller's override needs a LEVEL per stream (the pooled
        % law is a shape): the mean arc factor of every 'satKey|obsType' stream.
        % R5 M7 7b: model.lastSol (x, Cxx, s0sq of the previous pass; [] on
        % pass 1) lets the factory update anything that has its own variance
        % component outside the streams -- the prior groups (Koch-Kusche).
        model = struct('S', S, 'C', C, 'sigma2', sigma2, 'iter', it, ...
                       'streamSigma2', streamMeans(pinnedStream, sigma2), ...
                       'arcFactor', arcFactorMap(pinnedStream, pinnedRow, sigma2, o), ...   % R7-13
                       'Sgroups', Smap, 'Cgroups', Cmap, 'Pool', char(o.Pool), 'Estimate', o.Estimate, ...
                       'lastSol', lastSol, 'rowWeights', W);   % R6 WP-105: the robust factors
        acc   = accFcn(model);
        sol   = acc.solve();
        % WP-134 (third review, finding 1): components the accFcn updates itself
        % -- fit's common prior groups, one Koch-Kusche step per pass from the
        % previous solution -- belong to the stop rule. The loop stopped on the
        % arc factors and the law alone, and a run that settled on pass 1 (the
        % group is first updated on pass 2) never estimated the prior at all:
        % declared sigma 100 returned as 100, |x - x_true| 154 against 5.4.
        priorDev = 0;  priorPending = false;
        if ~isempty(o.PriorDev), [priorDev, priorPending] = o.PriorDev(); end
        lastSol = lastSolOf(acc, sol);       % WP-144 (A6): + N^-1 when s0sq is refused
        accF  = acc;  solF = sol;
        s2Prev = sigma2;                     % R7-11: the levels this pass was weighted with

        U = probeFactor(acc, it);            % WP-134: the backend's own N = U'U (finding 7)

        % ---- step 6: the shared probe, ONCE for every arc -----------------
        sp = lsqvce.engine.sharedProbe(U, 'numProbes', o.Probes, 'seed', o.Seed);

        % ---- steps 7-8: re-stream arcs, per-bin and per-arc statistics -----
        % WP-34 (fixes F-02): arcs are KEYED, the factor map is pinned on the
        % PASS-1 ENUMERATION, and a skipped arc keeps its previous factor --
        % one skip can no longer shift every later arc's factor by position.
        [Qb, Tb, qa, ta, keys, keysAll, skips, wstats, rowAll, streamAll, streamOf, lagRaw, robS] = ...
            sweepArcs(acc, sol, sp.UiZ, nLag, v, robustContext(rob, W, model.streamSigma2));   % WP-156
        % R5 M1: a skip RATE above the threshold is a refusal, not a smaller
        % estimate. The skipped arcs keep their factors and the run continues
        % on the rest only while they are few; beyond SkipThreshold the law
        % estimated from the survivors would describe a different data set.
        % WP-134 (third review, finding 15): the rate is PER STREAM, as the
        % Manual and the theory state -- over all arcs, a small stream skipped
        % whole among large ones (10 of 110 arcs: 9 %) kept its declared level
        % and could report converged; and it is checked BEFORE noArcs, whose
        % "nothing enumerated" hid an all-skipped pass.
        skipGuard(skips, keysAll, streamAll, it, o.SkipThreshold);
        nArc = numel(keys);
        if nArc == 0
            error('lsqvce:adjustVCEjoint:noArcs', ...
                  ['No arc produced statistics in pass %d: none was enumerated, or every ' ...
                   'one was skipped or inactive (SkipThreshold %g).'], it, o.SkipThreshold);
        end
        if isempty(pinnedKeys)
            % Pass 1: pin the arc identity on the ENUMERATION (skipped arcs
            % included -- a pass-1 skip keeps the default factor 1).
            [pinnedKeys, pinnedStream, pinnedRow] = deal(keysAll, streamAll, rowAll);
            s2Map = pinFactors(pinnedKeys, o.Sigma2_0);   % Sigma2_0 empty here: all ones
        else
            sameArcs(keysAll, pinnedKeys);
        end
        if ~isempty(skips)
            for kk = 1:size(skips, 1)
                hist.skips(end+1, :) = {it, skips{kk,1}, skips{kk,2}}; %#ok<AGROW>
            end
            if v
                lsqvce.util.log.msg(['  [joint %2d] SKIPPED %d arc(s) this ' ...
                    'pass (factors kept): %s\n'], it, size(skips,1), ...
                    strjoin(skips(:,1).', ', '));
            end
        end

        % ---- step 9: co-update the PSD AND the arc factors ----------------
        % WP-40 r2 (F-15): the alternation is a stochastic approximation --
        % the plain step (Gamma = 1 every pass) re-amplifies per-bin
        % estimation noise and is BISTABLE for short arcs (corr collapse to
        % negative values on unlucky seeds; wp33_loop_mirror.py grid in the
        % round note). 'StepSchedule','decay' uses the Robbins-Monro step
        % Gamma_it = 2/(2+it): full first steps, averaging later. Measured
        % across the mirror grid: min corr rises from -0.66 to ~0.9 with
        % barely any ceiling loss. Default 'none' is bitwise the old path.
        if strcmpi(o.StepSchedule, 'decay')
            gammaIt = 2 / (2 + it);
        else
            gammaIt = 1;
        end
        % R5 M3 3d: the law is estimated PER GROUP of streams (Pool), from the
        % arcs of streams whose policy is "law"; "scale" streams get arc
        % factors only (their mean is the Helmert factor); "none" streams are
        % untouched. out.S stays the 'all'-group S (the first group's otherwise).
        streamOf = streamOf(:);
        polOf    = cellfun(@(sk) policyOf(o.Estimate, sk), streamOf, 'UniformOutput', false);
        grpOf    = cellfun(@(sk) groupKeyOf(o.Pool, sk), streamOf, 'UniformOutput', false);
        lawCols  = strcmp(polOf, 'law');
        groups   = unique(grpOf(lawCols));
        % R7 WP-111 part 3: the bin statistics are exact for EVERY arc of the
        % group (all lags, the arc's own variance); the update is Helmert's
        % ratio per bin and the new law is rescaled to r(0) = 1 -- the level
        % change it carried moves into the factors through alpha2 below.
        % WP-147: binMax and blocking -- the bins' stop test (binsAgainstBound)
        pvMax = 0;  binMax = 0;  blocking = cell(0, 6);
        for g = 1:numel(groups)
            cols = lawCols & strcmp(grpOf, groups{g});
            Sg = S;  if isKey(Smap, groups{g}), Sg = Smap(groups{g}); end
            pvg = lsqvce.engine.psdVCE(Qb(:, cols), Tb(:, cols), Sg, 'Gamma', gammaIt, ...
                                       'Smooth', o.Smooth, 'FreezeRedundancy', o.FreezeRedundancy);
            Smap(groups{g}) = pvg.Shat / (Kc(1, :) * pvg.Shat);
            Cmap(groups{g}) = Kc * Smap(groups{g});
            pvMax = max(pvMax, pvg.maxDev);
            [binMax, blocking] = binsAgainstBound(pvg, o.Tol, groups{g}, binMax, blocking);
        end
        if ~isempty(groups)
            if isKey(Smap, 'all'), S = Smap('all'); else, S = Smap(groups{1}); end
            C = Kc * S;
        end
        pv = struct('maxDev', pvMax, 'binRatio', binMax);
        s2Cols = cellfun(@(k) s2Map(k), keys(:));    % survived, keyed lookup
        % R6 WP-98 (finding F-01): eq. 6.5.28 multiplies alpha2_m by the factor
        % the WEIGHTING used. What the weighting used is the per-STREAM level --
        % LSQVCE.ENGINE.JOINTMODELOPT can express nothing finer -- so that is
        % what must multiply here. Multiplying by the per-ARC factor instead
        % makes the update a geometric recursion with no fixed point, because
        % alpha2_m cancels its own arc's factor exactly (the arc's numerator and
        % denominator both carry it) and so never responds to it: MEASURED
        % growth of 8.4x per pass over four decades. For a stream with
        % arcFactors = false every arc already carries the one stream level, so
        % these two are the same number and such runs cannot move.
        % R7-13 (WP-116 part 2): an arc-factor stream's factors now REACH the
        % weighting (arcFactorMap -> jointModelOpt byArc -> the loader), so for
        % its arcs the factor the weighting used IS the per-arc factor, and that
        % is what multiplies here -- the same rule, with nothing left to cancel.
        % A stream without arc factors is weighted by its mean, as before.
        s2Applied = s2Cols;
        for sk = unique(streamOf(:)).'
            if arcFactorsOf(o.ArcFactors, sk{1}), continue, end
            cols = strcmp(streamOf, sk{1});
            s2Applied(cols) = mean(s2Cols(cols));
        end
        % R7-09 / WP-116 part 2: an arc-factor stream's factor is the arc's
        % ABSOLUTE Helmert factor gamma2_m = alpha2_m * sigma2_m (eq. 6.5.28),
        % the same number a stream without arc factors averages into its level
        % (fStream below). Ellmer's normalisation by gamma-bar (eqs. 6.5.29-30)
        % exists because there the PSD carries the level; here the law is a
        % SHAPE (normalised to r(0) = 1 above), so the level must
        % stay in the factors. Two wrong turns are on the record. gamma-bar over
        % EVERY column tied a stream's factors to UNRELATED streams' levels
        % (R7-09: 1.01 / 0.79 / 1.15 as only another stream changed). And the
        % 6.1.5 form (gamma2_m / gamma-bar_s^2) * mean(gamma2_s) counted the
        % level twice as soon as a stream's arcs differ: one arc 10x in sigma
        % among nine came back 1176 / 11.65 against 100 / 1, and the stream's
        % level 8x too high against the other streams. MEASURED (16 seeds,
        % mirrors/wp116/perarc_weighting.py): absolute factors 100.0 / 1.00.
        % gamma-bar keeps one job: the scale of the clamp floor (ClampEps), with
        % the kernel's default robust mean -- the 'MeanType' option that chose it
        % is retired (WP-118 part 2, user decision 2026-10-10).
        av = struct('gamma2', nan(nArc, 1), 'sigma2', nan(nArc, 1), 'gammaBar', []);
        lvl = containers.Map('KeyType', 'char', 'ValueType', 'double');   % pooled level per stream
        dispMap = containers.Map('KeyType', 'char', 'ValueType', 'any');      % WP-154 (M3)
        for sk = unique(streamOf(:)).'
            cols = strcmp(streamOf, sk{1});
            % R7 WP-111 part 3: the arc's component is its WHOLE covariance T, so
            % alpha2_m = etilde' T etilde / trace(Rbar T) (qa / ta from sweepArcs:
            % redundancyLags .omegaLaw / .sLaw, WP-132)
            avk = lsqvce.engine.arcVarianceFactors(qa(cols), ta(cols), s2Applied(cols));   % 1 x M rows
            % WP-156 (design F): the sweep deleted the REJECTED epochs and took
            % the truncation of the kept ones into the redundancy, lag by lag and
            % in the whole law's pair (lsqvce.kernel.truncationDeficit) -- the
            % level is no longer divided by tau (R7 WP-115 part 3b did, which is
            % exact for the total only: a red stream's variance read +2.4 %)
            [g2, s2, lvl(sk{1})] = arcLevels(avk, s2Applied(cols));   % WP-125/126: lvl pooled
            dispMap(sk{1}) = levelDispersion(avk, ta(cols));          % WP-154 (M3)
            av.gamma2(cols) = g2;  av.sigma2(cols) = s2;
            av.gammaBar(end+1) = avk.gammaBar;
        end
        % WP-40 r2 (F-15): Polyak tail averaging of the RETURNED spectrum
        % (log domain, last AvgTail passes). The loop keeps iterating on the
        % raw S; only the final S/C are averaged -- measured to add ~0.02
        % corr on top of the decay schedule (and to do NOTHING without it:
        % averaging cannot rescue an undamped collapse; mirror grid).
        if o.AvgTail > 0
            shatTail{end+1} = log(max(S, realmin));          %#ok<AGROW>
            gTail{end+1} = logLaws(Smap);                    %#ok<AGROW>  WP-134: every group
        end
        % R5 M4 4b: a stream WITHOUT arc factors still needs its LEVEL estimated
        % (the PSD shape is normalised; the level rode on the arc factors).
        % The level is the Helmert factor of the stream: the mean over its arcs
        % of gamma2_m = alpha2_m * sigma2_m (eq. 6.5.28, the arc's NEW absolute
        % factor) WITHOUT the unit-mean normalisation of the ellmer mode (that
        % normalisation exists because the PSD carries the level there; here
        % the shape is normalised, so the level must be absolute). One common
        % factor per stream. Found by the gravity demo: planted 1 / 2 / 0.5
        % came back 1 / 1 / 1 without it.
        % WP-125 (user decision 2026-10-12): the level is the POOLED Helmert
        % factor of the stream, sigma2 * sum_m q_m / sum_m t_m (arcVarianceFactors'
        % .pooled since WP-126, where it is unit-tested) -- the Helmert
        % estimator of ONE component over all its arcs -- not the mean of the
        % arcs' ratios, which weights a 30-epoch arc like a 3000-epoch one and,
        % with probe redundancies, is biased by E[1/t] > 1/E[t]. MEASURED
        % (dev/validation/wp125_level_mirror.py, 20 arcs, 6 locals each): on
        % mixed arcs (10 x 30, 10 x 2000) sd 0.047 -> 0.010; with 10 % probe
        % noise on t the bias +0.0105 -> +0.0011. Equal arcs: no change.
        fStream = containers.Map('KeyType', 'char', 'ValueType', 'double');
        for sk = unique(streamOf(:)).'
            cols = strcmp(streamOf, sk{1}) & ~strcmp(polOf, 'none') & ~arcFactorsOf(o.ArcFactors, sk{1});
            if ~any(cols), continue, end
            if isfinite(lvl(sk{1})) && lvl(sk{1}) > 0, fStream(sk{1}) = lvl(sk{1}); end   % .pooled
        end
        for kk = 1:nArc
            if strcmp(polOf{kk}, 'none'), continue, end
            if arcFactorsOf(o.ArcFactors, streamOf{kk})
                s2Map(keys{kk}) = av.sigma2(kk);     % keyed write-back, per arc
            elseif isKey(fStream, streamOf{kk})
                s2Map(keys{kk}) = fStream(streamOf{kk});    % one level per stream
            end
        end
        sigma2 = cellfun(@(k) s2Map(k), pinnedKeys(:));   % pinned order

        % R6 WP-97 (finding F-01): record the SPREAD of every arc-factor
        % stream, pass by pass. An arc factor is a RELATIVE weight inside its
        % stream, so its spread must settle; one that keeps growing is not an
        % estimate. arcSpreadGuard reads this after the loop.
        for sk = unique(streamOf(:)).'
            if ~arcFactorsOf(o.ArcFactors, sk{1}), continue, end
            f = sigma2(strcmp(pinnedStream, sk{1}));
            f = f(isfinite(f) & f > 0);
            if numel(f) < 2, continue, end
            if isKey(spreadHist, sk{1}), r0 = spreadHist(sk{1}); else, r0 = []; end
            spreadHist(sk{1}) = [r0; max(f)/min(f)];
        end

        % ---- step 9c (R6 WP-105): the robust pass -- next pass's factors ----
        % WP-156 (design F): THIS pass's sweep formed Baarda's w per epoch and
        % snooped every arc (robS); the influence turns the statistics into
        % the factors, written for the NEXT pass, one pass behind the model --
        % the same lag adjustRobust has. A period-2 cycle of the rejected set
        % is frozen at its union (freezeCycle, user decision 2026-10-18).
        if robust && rob.frozen == 0
            [W, fac, rkeys, rscale, rstreams, nDown, nZero, tauMap, hardFix, rob] = ...
                robustSweep(robS, o, size(sol.x, 1), rob);
            [W, fac, nZero, hardFix, rob] = freezeCycle(W, fac, rkeys, robS, hardFix, rob, it, v);
            hist.robustFrozen = rob.frozen;
        end
        if robust
            hist.robustScale(it,1:numel(rscale)) = rscale;       %#ok<AGROW>  one scale per stream (w units)
            hist.nDown(it,1) = nDown;  hist.nZero(it,1) = nZero;  %#ok<AGROW>
            if v
                lsqvce.util.log.msg('  [joint %2d] robust (%s): scales %s, %d down-weighted, %d at zero\n', ...
                    it, o.Influence, mat2str(rscale, 3), nDown, nZero);
            end
        end

        [hist, xPrev] = recordStep(hist, it, pv, sol, xPrev);   % maxDev, WP-147 binRatio / xStep
        % R7-11: the LEVELS and ARC FACTORS are free components too. A run
        % whose only free components are levels (no "law" stream: pvMax = 0)
        % stopped "converged" on pass 1 under the DECLARED weights, and prior
        % groups were never updated. The stop rule now also needs every arc's
        % factor to have settled: max |sigma2_new / sigma2_old - 1| <= Tol.
        if numel(s2Prev) ~= numel(sigma2)        % pass 1 without Sigma2_0: the start is all ones
            s2Prev = ones(size(sigma2));
        end
        okL = isfinite(sigma2) & isfinite(s2Prev) & s2Prev > 0;
        levelDev = 0;
        if any(okL), levelDev = max(abs(sigma2(okL) ./ s2Prev(okL) - 1)); end
        hist.levelDev(it,1) = levelDev;              %#ok<AGROW>
        hist.priorDev(it,1) = priorDev;              %#ok<AGROW>  WP-134
        hist.S(it,:)       = S(:).';                 %#ok<AGROW>
        % R5 M4 4a: the per-stream pass record Results.convergence reads --
        % the stream's level factor and its Ljung-Box whiteness on the
        % law-whitened residuals of THIS pass (K = 20 lags, alpha 0.05).
        [hist.streams, sigRow, whiteRow, pRow] = whitenessRow(wstats, pinnedStream, sigma2);
        hist.sigma(it, :)     = sigRow;              %#ok<AGROW>
        hist.white(it, :)     = whiteRow;            %#ok<AGROW>
        hist.whitePval(it, :) = pRow;                %#ok<AGROW>

        if v
            lsqvce.util.log.msg(['  [joint %2d] shape change %.4e   bins %.2f of ' ...
                'their bound   arcs = %d   gammaBar (per stream) = %s\n'], it, pv.maxDev, ...
                pv.binRatio, nArc, mat2str(av.gammaBar, 4));
        end

        % ---- step 9b: the user's per-pass inspection hook (R5 M5 5c') -----
        ctx = struct('driver', 'adjustVCEjoint', 'pass', it, 'x', sol.x(:), ...
                     'Cxx', sol.Cxx, 'sigma0', sol.sigma0, 's0sq', sol.s0sq, 'S', S, 'C', C, ...
                     'Sgroups', Smap, 'Cgroups', Cmap, ...
                     'sigma2', sigma2, 'arcKeys', {pinnedKeys}, 'streams', {hist.streams}, ...
                     'sigma', hist.sigma(it, :), 'white', hist.white(it, :), ...
                     'whitePval', hist.whitePval(it, :), 'maxDev', pv.maxDev, ...
                     'binRatio', pv.binRatio, 'xStep', hist.xStep(it), ...   % WP-147
                     'skips', {skips}, 'acc', acc, ...
                     'factor', fac, 'keys', {rkeys}, 'scale', rscale);   % R6 WP-105 ([] when Influence "none")
        userStop = ~lsqvce.util.onPass(opt.OnPass, opt.OnPassContext, ctx);

        % ---- step 10: converged? -----------------------------------------
        if pv.binRatio <= 1 && levelDev <= o.Tol && priorDev <= o.Tol && ~priorPending
            converged = true;
            break
        end
        if userStop, break, end
    end

    if ~converged && ~userStop           % WP-147: the text says what still moves
        lsqvce.util.log.warn('lsqvce:adjustVCEjoint:notConverged', '%s', ...
            notSettledMessage(o.MaxIter, o.Tol, hist, blocking));
    end

    arcSpreadGuard(spreadHist);          % R6 WP-97 (F-01): report, never silence
    checks = modelChecks(hist, Cmap, dispMap, lagRaw, streamOf, o);   % WP-154 (M2-M4)

    if o.AvgTail > 0 && numel(shatTail) >= o.AvgTail
        S = exp(mean(cat(2, shatTail{end-o.AvgTail+1:end}), 2));
        S = S / (Kc(1, :) * S);  C = Kc * S;
        % WP-134: every pool group's law too -- out.Sgroups kept the raw last
        % pass (out.Sgroups('all') ~= out.S), and under Pool ~= "all" the final
        % solve reads the groups, so the average never reached it
        [Smap, Cmap] = tailAverage(Smap, Cmap, gTail(end-o.AvgTail+1:end), Kc);
    end

    % R7-11: the returned SOLUTION under the returned MODEL. The loop's last
    % accumulation was weighted with the model BEFORE that pass's update, so
    % x, Cxx and sigma0 were one pass behind the C / levels / factors handed
    % back next to them. One more accumulation under the final model (same
    % pass index, so the prior activation is that of the last pass; lastSol
    % lets the prior groups take their last update). Not after a user stop:
    % the hook asked for the loop to end where it was.
    % R7 WP-115 part 3b (user decision 2026-10-11, measured on CI): a robust
    % run's final solution is HARD -- the rejected epochs deleted, partial
    % factors back to 1 -- and its sigma0 / Cxx are corrected for the
    % truncation (hardSolution). With the partial factors in this solve,
    % v'C^-1 v carried the selection bias and sigma0, recomputed under the
    % corrected levels, gave it all back: sigma*sigma0 stayed 0.83.
    if ~userStop
        modelF = struct('S', S, 'C', C, 'sigma2', sigma2, 'iter', numel(hist.maxDev), ...
                        'streamSigma2', streamMeans(pinnedStream, sigma2), ...
                        'arcFactor', arcFactorMap(pinnedStream, pinnedRow, sigma2, o), ...   % R7-13
                        'Sgroups', Smap, 'Cgroups', Cmap, 'Pool', char(o.Pool), 'Estimate', o.Estimate, ...
                        'lastSol', lastSol, 'rowWeights', hardWeights(W), 'final', true);
        accF = accFcn(modelF);
        solF = hardSolution(accF.solve(), hardFix, accF.excludedRows());
    end
    out = struct('S', S, 'C', C, 'sigma2', sigma2, 'acc', accF, ...
                 'userStop', userStop, 'Sgroups', Smap, 'Cgroups', Cmap, ...
                 'sol', solF, 'converged', converged, 'nIter', numel(hist.maxDev), ...
                 'history', hist, ...
                 'arcKeys', {pinnedKeys(:)}, 'arcStreams', {pinnedStream(:)}, ...
                 'skips', {hist.skips}, ...
                 'seed', o.Seed, 'checks', checks, ...                 % R7-37; WP-154
                 'robust', struct('influence', char(o.Influence), 'scale', rscale, ...
                                  'streams', {rstreams}, 'factor', fac, 'keys', {rkeys}, ...      % R6 WP-105
                                  'truncation', tauRow(tauMap, rstreams), ...                    % R7 WP-115
                                  'frozen', rob.frozen));                                        % WP-156
end

%==========================================================================
function [W, fac, rkeys, rscale, rstreams, nDown, nZero, tauMap, hardFix, rob] = robustSweep(robS, o, nc, rob)
%ROBUSTSWEEP  The robust pass (R6 WP-105; WP-156 design F): the factor map
%   W for the NEXT pass (one per epoch) from the statistics THIS pass's
%   sweep formed -- robS from sweepArcs: per epoch Baarda's w for correlated
%   observations in units of the stream's scale (lsqvce.kernel.wSnoop), +-Inf
%   on an epoch the arc's snooping deleted, NaN on one excluded on every
%   axis. The influence acts per stream (perStreamWeights); its MAD of the
%   statistic, rscale, carries the stream's scale to the next pass (rob.scale
%   x rscale, level-free -- robustContext -- so the snooping threshold settles
%   where the MAD of the final system is one). tauMap, the truncation factor per stream at the
%   influence's rejection point (averaged over the stream's epochs by their
%   axis count), and hardFix (.tau the row-weighted factor, .nRej the
%   rejected ROWS) correct the final hard solution's sigma0 (hardSolution):
%   the total v'Pv of a sample truncated in its w directions is exactly tau
%   times the untruncated one. The level is no longer corrected here (the
%   deficit is in the sweep's redundancy). The collapse refusal and the
%   breakdown warning as before.
    uEp = robS.u;  rkeys = robS.keys;  sizes = robS.sizes;  dEp = robS.d;
    [fac, rscale, rstreams, kReject, rowsOf] = perStreamWeights(uEp, rkeys, sizes, o.Influence, o.Tuning);
    tauMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for si = 1:numel(rstreams)
        [ud, ~, j] = unique(dEp(rowsOf{si}));          % one factor per axis count
        tv = arrayfun(@(d) truncationFactor(kReject, d), ud);
        tauMap(rstreams{si}) = mean(tv(j));
        s0 = 1;  if isKey(rob.scale, rstreams{si}), s0 = rob.scale(rstreams{si}); end
        if isfinite(rscale(si)) && rscale(si) > 0, rob.scale(rstreams{si}) = s0 * rscale(si); end
    end
    rob.active = isfinite(kReject);                   % the deficit applies from the next sweep on
    hardFix = hardFixOf(fac, dEp, kReject);
    nAlive = sum(fac > 0);
    if nAlive < nc + 1                                              % live EPOCHS
        error('lsqvce:adjustVCEjoint:collapsed', ...
             ['The influence function left %d epochs with a non-zero weight, ' ...
              'fewer than the %d parameters + 1: almost every epoch was rejected. ' ...
              'Two things do that, and neither is an ordinary outlier problem. ' ...
              '(1) The residuals are offset BODILY from zero by several MAD -- a ' ...
              'bias the model did not absorb -- because the influence standardises ' ...
              'u/s, not (u - med)/s. (2) Contamination so near total that fewer ' ...
              'than %d clean epochs are left; contamination merely past 50 %% does ' ...
              'NOT reach here (HANDOVER 12.34). Check the model, not the data.'], ...
              nAlive, nc, nc + 1);
    end
    warnAtBreakdown(fac, 'adjustVCEjoint');        % engine-private: > 45 % of the epochs zeroed (R7 WP-115)
    W = factorMap(fac, rkeys, sizes);
    rscale = cellfun(@(k) rob.scale(k), rstreams);
    nDown = sum(fac < 0.5 & fac > 0);
    nZero = sum(fac == 0);
end

function [W, tauMap, hardFix, rob] = robustPassZero(acc0, o, nLag, v, rob, W, lvl)
%ROBUSTPASSZERO  PASS 0 of a robust run (round 10): the first robust pass
%   runs BEFORE the first spectral update, on ACC0, the accumulation under
%   the start model. WP-156 (design F): the statistic is Baarda's w, formed
%   per epoch by the sweep, so pass 0 sweeps once for it alone. It does not
%   snoop -- the stream's scale is not known yet -- and rejects
%   simultaneously under the MAD of the undeleted system (a warm start);
%   from pass 1 the sweep snoops with the previous pass's scale.
    sol0 = acc0.solve();
    sp0 = lsqvce.engine.sharedProbe(probeFactor(acc0, 0), 'numProbes', o.Probes, 'seed', o.Seed);
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, robS] = sweepArcs(acc0, sol0, sp0.UiZ, nLag, v, robustContext(rob, W, lvl));
    [W, ~, ~, rscale0, ~, nDown0, nZero0, tauMap, hardFix, rob] = robustSweep(robS, o, size(sol0.x, 1), rob);
    if v
        lsqvce.util.log.msg('  [joint  0] robust (%s): scales %s, %d down-weighted, %d at zero\n', ...
            o.Influence, mat2str(rscale0, 3), nDown0, nZero0);
    end
end

function hardFix = hardFixOf(fac, dEp, kReject)
%HARDFIXOF  The final hard solution's correction: .tau the truncation factor
%   weighted by the rows of every epoch, .nRej the rejected rows.
    ud = unique(dEp);  tAll = arrayfun(@(d) truncationFactor(kReject, d), ud);
    [~, j] = ismember(dEp, ud);
    hardFix = struct('tau', sum(dEp .* tAll(j)) / sum(dEp), 'nRej', sum(dEp(fac == 0)));
end

function W = factorMap(fac, rkeys, sizes)
%FACTORMAP  The per-epoch factors keyed by the loader's arc key.
    W = containers.Map('KeyType', 'char', 'ValueType', 'any');
    off = 0;
    for kk = 1:size(rkeys, 1)
        n = sizes(kk);                                              % EPOCHS of the arc
        W(lsqvce.util.rowWeightKey(rkeys{kk,1}, rkeys{kk,3}, rkeys{kk,2}, rkeys{kk,4})) = fac(off + (1:n));
        off = off + n;
    end
end

function rob = robustState(robust, o)
%ROBUSTSTATE  WP-156: the robust pass's state across passes -- .on (a robust
%   run), .kRej the influence's rejection point (Inf for Huber), .tau(d) its
%   truncation factor for a d-axis epoch, .scale the stream's scale of the
%   w statistic (empty before pass 0's MAD), .active (a robust pass has
%   rejected: the sweep subtracts the truncation deficit), .frozen (the pass
%   a period-2 cycle was frozen at, 0 none) and .zeroSets (the last two
%   rejected sets, freezeCycle).
    kRej = Inf;
    if robust
        [~, ~, ~, ~, kRej] = lsqvce.kernel.robustWeights(linspace(-2, 2, 41).', o.Influence, o.Tuning);
    end
    rob = struct('on', robust, 'kRej', kRej, 'tau', @(d) truncationFactor(kRej, d), ...
                 'scale', containers.Map('KeyType', 'char', 'ValueType', 'double'), ...
                 'active', false, 'frozen', 0, 'zeroSets', {{}});
end

function rc = robustContext(rob, W, lvl)
%ROBUSTCONTEXT  What the sweep needs from the robust pass (WP-156): .on (form
%   the statistic: a robust run whose set is not frozen), .thr the snooping
%   threshold (Inf until a scale exists -- pass 0 does not snoop), .scale
%   (stream -> s, level-free: the scale of w sqrt(L)), .level (stream -> L,
%   the stream's level factor the arcs were weighted with, model.
%   streamSigma2; 1 where none), .W the factors the arcs were weighted with
%   (which epochs the robust pass deleted, as opposed to the data's own
%   exclusions), and .tauBar(d) = 1 - tau(d) once a robust pass has rejected
%   (else []). WHY level-free: w is standardised by the law, which carries the
%   estimated level, so w scales as 1 / sqrt(L); a scale carried in w units
%   from the last pass lagged every change of the level (pass 0 -> 1 moves it
%   from the declared to the estimated one: 4 x for a stream declared 1 and
%   truly 0.5) and one pass snooped 13 % of the epochs -- measured, the level
%   then read 6-11 % low and nearly every run froze a cycle.
    if isempty(lvl), lvl = containers.Map('KeyType', 'char', 'ValueType', 'double'); end
    rc = struct('on', rob.on && rob.frozen == 0, 'thr', Inf, 'scale', rob.scale, 'level', lvl, ...
                'W', W, 'tauBar', []);
    if rob.scale.Count > 0, rc.thr = rob.kRej; end
    if rob.on && rob.active && isfinite(rob.kRej), rc.tauBar = @(d) 1 - rob.tau(d); end
end

function [W, fac, nZero, hardFix, rob] = freezeCycle(W, fac, rkeys, robS, hardFix, rob, it, v)
%FREEZECYCLE  WP-156 (user decision 2026-10-18): a rejected set that returns
%   to the set of two passes ago, having differed from the last one, is a
%   period-2 cycle -- an epoch or two at the threshold that the law's
%   feedback on w flips each pass (measured: 2-5 of 20 seeds per scenario
%   without the freeze). It is frozen at the UNION of the two sets (the
%   borderline epochs deleted) and the loop settles the law under it; the
%   sweep keeps subtracting the deficit and no longer forms the statistic.
%   Reported in the log, history.robustFrozen and out.robust.frozen.
    rowKeys = arcRowKeys(rkeys);
    z = zeroSetOf(fac, rowKeys, robS.sizes);
    H = rob.zeroSets;
    if numel(H) >= 2 && isequal(z, H{end-1}) && ~isequal(z, H{end})
        add = setdiff(H{end}, z);
        off = [0, cumsum(robS.sizes(:).')];
        for a = 1:numel(add)
            cut = find(add{a} == '#', 1, 'last');
            kk = find(strcmp(rowKeys, add{a}(1:cut-1)), 1);
            if ~isempty(kk), fac(off(kk) + str2double(add{a}(cut+1:end))) = 0; end
        end
        W = factorMap(fac, rkeys, robS.sizes);
        hardFix = hardFixOf(fac, robS.d, rob.kRej);
        rob.frozen = it;
        if v
            lsqvce.util.log.msg('  [joint %2d] robust: the rejected set cycles; frozen at the union (%d epochs)\n', ...
                it, sum(fac == 0));
        end
    end
    rob.zeroSets = [H(max(1, end-1):end), {z}];
    nZero = sum(fac == 0);
end

function z = zeroSetOf(fac, rowKeys, sizes)
%ZEROSETOF  The rejected epochs as sorted 'arcKey#epoch' labels (the loader's
%   arc key, so a label names the same epoch every pass, skips or not).
    z = {};
    off = 0;
    for kk = 1:numel(sizes)
        e = find(fac(off + (1:sizes(kk))) == 0);
        z = [z, arrayfun(@(x) sprintf('%s#%d', rowKeys{kk}, x), e(:).', 'UniformOutput', false)]; %#ok<AGROW>
        off = off + sizes(kk);
    end
    z = sort(z);
end

function k = arcRowKeys(rkeys)
%ARCROWKEYS  The loader's key of every arc of rkeys {sysID, obsType, day, arc}.
    k = cell(1, size(rkeys, 1));
    for kk = 1:size(rkeys, 1)
        k{kk} = lsqvce.util.rowWeightKey(rkeys{kk,1}, rkeys{kk,3}, rkeys{kk,2}, rkeys{kk,4});
    end
end

%==========================================================================
function [fac, scale, streams, kReject, rowsOf] = perStreamWeights(vw, keys, sizes, influence, tuning)
%PERSTREAMWEIGHTS  Robust factors PER STREAM (R6 WP-105, round 2).
%   The influence function standardises residuals by a MAD scale, and that
%   scale must be ONE PER STREAM and must NOT be the estimated level: the
%   first design used a global MAD and it drove the level of a clean stream
%   into a spiral (0.377 measured against 0.80 expected) -- a global scale
%   mis-standardises every stream whose level differs from the mix, the
%   next level is estimated from the mis-weighted residuals, and the next
%   scale sees that. With a MAD per stream the standardisation is
%   v / MAD(v), invariant to the level the whitener applied, so the weights
%   do not feed on the level they will shape: the level converges to the
%   retained-observation ratio, not to zero. Prior streams get factor 1.
%   A MULTI-AXIS arc (axes(kk) = d > 1, rows epoch-major) is weighted PER
%   EPOCH (round 4): the influence acts on ONE statistic per epoch, and the d
%   rows of an epoch share the factor. Since R7-08 that statistic is the
%   chi-calibrated deviate (it was the RMS, which zeroed a quarter of clean
%   3-axis epochs under IGG3); since WP-156 of the epoch's w (lsqvce.kernel.
%   wSnoop), and an epoch excluded on every axis (NaN) keeps factor 1 and stays
%   out of the MAD. That is what the sweep needs: its temporal covariance is
%   per axis and per epoch, so a factor must be per epoch to be expressible
%   there.
%   OUTPUT fac n x 1; scale 1 x S (the MAD per stream); streams 1 x S cell
%          'satKey|obsType', sorted; kReject 1x1 the influence's rejection
%          point (Inf for Huber); rowsOf 1 x S cell, each stream's epochs.
    n = numel(vw);
    fac = ones(n, 1);
    skey = cell(size(keys, 1), 1);
    for kk = 1:size(keys, 1)
        skey{kk} = lsqvce.util.streamKey(keys{kk,1}, keys{kk,2});
    end
    % prior pseudo-streams keep factor 1 -- judged on the obsType each key was
    % built from (WP-138: the key was parsed for '|prior_')
    isPri = cellfun(@(t) lsqvce.util.priorStream.is(char(string(t))), keys(:, 2));
    streams = unique(skey(~isPri)).';
    scale = nan(1, numel(streams));
    rowsOf = cell(1, numel(streams));
    kReject = Inf;
    ends = cumsum(sizes(:));  starts = ends - sizes(:) + 1;
    for si = 1:numel(streams)
        arcs = find(strcmp(skey, streams{si}));
        rows = [];
        for kk = arcs(:).'
            rows = [rows, starts(kk):ends(kk)]; %#ok<AGROW>
        end
        live = rows(~isnan(vw(rows)));                 % WP-156: NaN = excluded on every axis, factor 1
        [fac(live), scale(si), ~, ~, kReject] = lsqvce.kernel.robustWeights(vw(live), influence, tuning);
        rowsOf{si} = rows;
    end
end

function tau = truncationFactor(k, d)
%TRUNCATIONFACTOR  E[q | q < q_k] / d for q ~ chi2_d, the rejection point k
%   mapped through the two-sided normal tail: P(chi2_d > q_k) = P(|N| > k)
%   (the chi-calibrated statistic of a d-axis epoch is that N, R7-08). The
%   mean square of a Gaussian sample truncated there is tau times the
%   untruncated one, tau = F_{d+2}(q_k) / F_d(q_k); d = 1 gives
%   1 - 2 k phi(k) / (2 Phi(k) - 1) (0.97334 at k = 3). k = Inf: 1.
%   R7 WP-115 part 3b, user decision 2026-10-11.
    if ~isfinite(k), tau = 1; return, end
    p = erfc(k / sqrt(2));
    if d == 1
        qk = k^2;
    else
        lo = 0;  hi = d + 50 * sqrt(2 * d) + k^2;
        for it = 1:200
            qk = (lo + hi) / 2;
            if gammainc(qk / 2, d / 2, 'upper') > p, lo = qk; else, hi = qk; end
            if hi - lo <= 4 * eps(qk), break, end
        end
    end
    tau = gammainc(qk / 2, d / 2 + 1) / gammainc(qk / 2, d / 2);
end

function H = hardWeights(W)
%HARDWEIGHTS  The factor map of the final hard solution: 0 stays 0 (rejected,
%   a deletion-sized nugget), every partial factor becomes 1. [] / empty in,
%   the same out (no robust pass).
    H = W;
    if isempty(W) || W.Count == 0, return, end
    H = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for k = keys(W)
        H(k{1}) = double(W(k{1}) > 0);
    end
end

function sol = hardSolution(sol, hardFix, nExcl)
%HARDSOLUTION  sigma0, s0sq and Cxx of the hard solution, corrected for the
%   truncation: the kept rows are a Gaussian sample truncated at the
%   rejection point, E = tau. The rejected rows leave sol.dof already when
%   the loader marked them (WP-133: makeLoader's arc.excluded, nExcl of them
%   in the accumulator); any the loader did not mark are taken out here, as
%   before. So s0sq * dof / ((dof - (nRej - nExcl)) tau), and sol.dof is the
%   live count. sol.truncation = tau, which lsqvce.engine.residuals applies
%   to ITS s0^2 (review 3 finding 18: it cross-checked against the corrected
%   number and warned sigma0Mismatch on every igg3 / bisquare noise run).
%   hardFix = [] (no robust pass): unchanged. A refused v'Pv (s0sq NaN) stays
%   NaN but still carries dof and tau: the verdict's residual pass, which
%   replaces it (lsqvce.engine.verdict.supersede), divides by sol.truncation --
%   skipped, a normal-backend igg3 run returned sigma0 1.3 % off QR's (WP-149).
    if isempty(hardFix), return, end
    nLeft = max(hardFix.nRej - nExcl, 0);
    f = sol.dof / max(sol.dof - nLeft, 1) / hardFix.tau;
    sol.dof = sol.dof - nLeft;
    sol.s0sq = sol.s0sq * f;
    sol.sigma0 = sqrt(sol.s0sq);
    sol.Cxx = sol.Cxx * f;
    sol.truncation = hardFix.tau;
end

function [gamma2, sigma2, pooled] = arcLevels(avk, s2Applied)
%ARCLEVELS  One stream's arc factors and pooled level from arcVarianceFactors'
%   AVK. WP-156: no division by a truncation factor any more -- a robust
%   run's truncation is in the sweep's redundancy (lsqvce.kernel.
%   truncationDeficit); until then (R7 WP-115 part 3b, WP-152 N2) every
%   factor but a held arc's was divided by tau here.
%   INPUT   avk struct, s2Applied M x 1 (the factors as weighted)
%   OUTPUT  gamma2 / sigma2 M x 1 (sigma2: the clamped gamma2_m, or s2Applied
%           when the mean was unusable -- badMean: unchanged this pass), pooled 1x1
    gamma2 = avk.gamma2;
    pooled = avk.pooled;
    if isfinite(avk.gammaBar) && avk.gammaBar > 0
        sigma2 = avk.sigma2 * avk.gammaBar^2;
    else
        sigma2 = s2Applied;
    end
end

function t = tauRow(tauMap, streams)
%TAUROW  The truncation factors of STREAMS (1 x S, aligned): the level's
%   divisor per stream, 1 where none applies (no robust pass, Huber, a
%   prior stream, or tauMap = [] before the first robust sweep).
    t = ones(1, numel(streams));
    if isempty(tauMap), return, end
    for si = 1:numel(streams)
        if isKey(tauMap, streams{si}), t(si) = tauMap(streams{si}); end
    end
end

%==========================================================================
function [Qb, Tb, qa, ta, keys, keysAll, skips, wstats, rowAll, streamAll, keyStream, lagRaw, robS] = ...
    sweepArcs(acc, sol, UiZ, L, v, rc)
%SWEEPARCS  Steps 7-8: re-stream every arc and return its Helmert statistics,
%   one COLUMN per arc (R7 WP-111 part 3):
%     Qb(:, m) = lambda_m K' Omega^m,  Tb(:, m) = lambda_m K' s^m   (L x M)
%       the per-BIN numerators / denominators psdVCE sums per pool group, with
%       Omega^m / s^m the arc's lag sums over ALL its lags
%       (LSQVCE.KERNEL.REDUNDANCYLAGS), K = LSQVCE.KERNEL.HATACV(L, arc length)
%       and lambda_m = r_m(0) the arc's own law variance (level x factor, as
%       weighted): the bin component of the arc is lambda_m S_j toeplitz(K_j);
%       K' is applied in hatACV's projection form (fold + one FFT, no m x L K)
%     qa(m) = etilde' T etilde,  ta(m) = trace(Rbar T)              (1 x M)
%       the arc's numerator / denominator with its WHOLE law covariance T as
%       the component -- the arc factor (eq. 6.5.27), in closed form
%       (redundancyLags .omegaLaw / .sLaw, WP-132; r_m' Omega^m / r_m' s^m
%       before, which assumed T Toeplitz)
%   Mirrors the arc enumeration used by LSQVCE.ENGINE.REDUNDANCY and
%   LSQVCE.ENGINE.STREAMCONTRIBUTION, so grouping and whitening stay consistent.
%   WP-34 (fixes F-02): every arc is KEYED s<sys>_<obsType>_d<day>_a<arc>
%   (LSQVCE.ENGINE.REDUNDANCY's 'arc' key format, verbatim). OUTPUTS:
%     keys     1 x M   cell  keys of the arcs that produced columns, ALIGNED
%                            to the columns of Qb/Tb/qa/ta
%     keysAll  1 x Ma  cell  keys of every ENUMERATED arc (the identity the
%                            caller pins on pass 1)
%     skips    S x 2   cell  {key, reason} for every arc skipped THIS pass
%     rowAll   Ma x 1  cell  the LOADER's key of every enumerated arc,
%                            LSQVCE.UTIL.ROWWEIGHTKEY(sysID, day, obsType, arc),
%                            aligned to keysAll (R7-13: how a per-arc factor
%                            reaches the arc it belongs to)
%     streamAll Ma x 1 cell  the STREAM key of every enumerated arc
%                            (lsqvce.util.streamKey), aligned to keysAll
%     keyStream 1 x M  cell  the same, aligned to keys
%   R7-43: the stream of an arc is CARRIED, never parsed back out of its key.
%     robS     struct  WP-156: the robust pass's statistics (rc.on), per epoch
%                      of every arc in sweep order -- .u (lsqvce.kernel.wSnoop's
%                      z), .d (axes), .keys (A x 4 {sysID, obsType, day, arc}),
%                      .sizes (1 x A epochs); rc (robustContext) also carries the
%                      truncation deficit oneArc subtracts from the redundancy
%
%   NO SHARED WEIGHT ANY MORE (R6 WP-100 is superseded). Until part 3 one
%   pooled C weighted the lag statistics of every arc in eqs. 6.5.22-6.5.23
%   while each pool group's prior was its own law -- a compromise that was
%   measured better than threading each group's freshly fitted law into its
%   arcs (planted lag-1 0.50 came back 1.331). With positive components the
%   weight is not a choice: an arc's bin component is lambda_m S_j
%   toeplitz(K_j) under the covariance it WAS weighted with, so lambda_m is
%   read from that covariance and K does not depend on the law at all.
%   The behaviour is pinned by test_vcePolicy>testPoolStreamGivesEachStreamItsOwnLaw.
    ids    = acc.systemIDs();
    x      = sol.x(:);
    qa = zeros(1, 0);  ta = zeros(1, 0);
    lagRaw = struct('om', {{}}, 'sn', {{}}, 'lam', zeros(1, 0));   % WP-154: per arc, before the bins
    keys   = {};
    keysAll = {};
    rowAll  = {};
    streamAll = {};  keyStream = {};
    skips  = cell(0, 2);
    nArc   = 0;
    wstats = containers.Map('KeyType', 'char', 'ValueType', 'any');   % stream -> [lag sums; pair counts] (R7-38)
    robS = struct('u', zeros(0, 1), 'd', zeros(0, 1), 'keys', {cell(0, 4)}, 'sizes', zeros(1, 0));   % WP-156

    for k = 1:numel(ids)
        st = acc.getStore(ids(k));
        uo = st.uo;
        % WP-33: recover this system's LOCAL solution and Schur pieces ONCE
        % per pass, exactly as LSQVCE.ENGINE.REDUNDANCY does -- the residual is the
        % full-vector one (Ellmer step 7, Table 6.2) and the projector gains
        % the exact local block via jointArcParts + redundancyLags.
        nlS = size(st.Nlc, 1);
        if nlS > 0
            [ySys, Gsys, LlS, plS] = lsqvce.engine.BlockAccumulator.localSolveOf(st, x);   % WP-128
            if plS > 0
                error('lsqvce:adjustVCEjoint:localNotPD', ...
                     ['System %s: the local normal matrix is not positive ' ...
                      'definite -- the local block cannot be absorbed into ' ...
                      'the projector.'], string(st.sysID));
            end
        else
            ySys = [];  Gsys = [];  LlS = [];
        end
        for d = st.dayList(:).'
            meta = st.dataFcn('meta', d, [], uo);
            for ot = 1:numel(meta.obsTypes)
                obsType = meta.obsTypes{ot};
                for a = meta.arcs{ot}(:).'
                    key = lsqvce.util.arcKey(st.sysID, obsType, d, a);
                    sk  = lsqvce.util.streamKey(st.sysID, obsType);
                    keysAll{end+1} = key;                    %#ok<AGROW>
                    streamAll{end+1, 1} = sk;                %#ok<AGROW>
                    rowAll{end+1, 1} = lsqvce.util.rowWeightKey(st.sysID, d, obsType, a); %#ok<AGROW>
                    [om, sn, r, why, wsum, wm, qw, tw, zE, dE] = oneArc(st, obsType, a, d, uo, UiZ, ...
                                                                x, ySys, Gsys, LlS, rc, sk);
                    if isempty(om)
                        % an inactive prior arc is not a skip (WP-127): it counted
                        % against SkipThreshold before its group's 'from' pass
                        if ~strcmp(why, 'inactive'), skips(end+1, :) = {key, why}; end   %#ok<AGROW>
                        continue
                    end
                    if isKey(wstats, sk), w0 = wstats(sk); else, w0 = zeros(2*(WHITE_LAGS() + 1), 1); end
                    wstats(sk) = w0 + [wsum; wm];
                    nArc = nArc + 1;
                    keys{end+1}     = key;                   %#ok<AGROW>
                    keyStream{end+1} = sk;                   %#ok<AGROW>
                    lagRaw.om{nArc} = om;  lagRaw.sn{nArc} = sn;  lagRaw.lam(nArc) = r(1);   % lambda_m
                    qa(nArc) = qw;                           %#ok<AGROW>  WP-132: the whole law,
                    ta(nArc) = tw;                           %#ok<AGROW>  no Toeplitz assumption
                    if ~isempty(zE)                          % WP-156: the robust statistic
                        robS.u = [robS.u; zE];  robS.d = [robS.d; dE * ones(numel(zE), 1)];
                        robS.keys(end+1, :) = {st.sysID, obsType, d, a};
                        robS.sizes(end+1) = numel(zE);
                    end
                end
            end
        end
        if v
            lsqvce.util.log.msg('%sswept system %s\n', lsqvce.util.verbosity.banner('Joint estimation'), string(st.sysID));
        end
    end
    [Qb, Tb] = projectLags(lagRaw, L);                % WP-154: the arcs' raw statistics are kept
end

function [Qb, Tb] = projectLags(lagRaw, L)
%PROJECTLAGS  The bin statistics of every arc, lambda_m K' [Omega^m s^m]
%   (LSQVCE.KERNEL.HATACV's projection form, no m x L K), one COLUMN per arc.
%   WP-154: the sweep keeps each arc's raw lag statistics (lagRaw), whose
%   lengths the model checks read (the arcs' median, M4).
    M = numel(lagRaw.om);
    Qb = zeros(L, M);  Tb = zeros(L, M);
    for a = 1:M
        qt = lagRaw.lam(a) * lsqvce.kernel.hatACV(L, numel(lagRaw.om{a}), [lagRaw.om{a} lagRaw.sn{a}]);
        Qb(:, a) = qt(:, 1);  Tb(:, a) = qt(:, 2);
    end
end

function sameArcs(keysAll, pinnedKeys)
%SAMEARCS  Refuse a pass whose arc key set differs from pass 1's (WP-34).
%   ERRORS  lsqvce:adjustVCEjoint:arcMismatch
    if isequal(sort(keysAll), sort(pinnedKeys)), return, end
    missing = setdiff(pinnedKeys, keysAll);
    newk    = setdiff(keysAll, pinnedKeys);
    error('lsqvce:adjustVCEjoint:arcMismatch', ...
        ['The arc key set changed between passes -- the loader must ' ...
         'be deterministic across passes (an arc that merely fails ' ...
         'its covariance is SKIPPED and kept, but one that vanishes ' ...
         'from the enumeration is an identity error). Missing: %s. ' ...
         'New: %s.'], ...
        strjoin([missing, {'(none)'}], ', '), ...
        strjoin([newk,    {'(none)'}], ', '));
end

function Dt = levelDispersion(avk, ta)
%LEVELDISPERSION  WP-154 (M3): the scatter of a stream's per-arc Helmert
%   ratios alpha2_m about their redundancy-weighted mean, in units of the
%   sampling variance 2 / r_m of a ratio of redundancy r_m (Gaussian):
%   D = sum_m (alpha2_m / abar - 1)^2 r_m / 2 / (M - 1), about 1 for a stream
%   whose level is the same on every arc (chi^2_{M-1} / (M - 1) for white
%   noise). Its bound: the 99.9 % point of that, at least 4 -- a coloured
%   stream's ratios overlap and scatter more (MEASURED on 12 arcs of 256
%   epochs, AR(1) 0.5: D 0.35-1.42 homogeneous, 3.2-5.0 for levels 1 / 1.1 /
%   1.2 by day, ~120 for 0.5 / 1 / 2).
%   OUTPUT  Dt  1x2 [D, bound]; [NaN NaN] below 3 usable arcs
    a = avk.alpha2(:);  r = ta(:);
    ok = ~avk.held(:) & isfinite(a) & isfinite(r) & r > 0;
    Dt = [NaN NaN];
    M = nnz(ok);
    if M < 3, return, end
    abar = sum(a(ok) .* r(ok)) / sum(r(ok));
    D = sum((a(ok) / abar - 1).^2 .* r(ok) / 2) / (M - 1);
    Dt = [D, max(4, 2 * gammaincinv(0.999, (M - 1) / 2) / (M - 1))];
end

function checks = modelChecks(hist, Cmap, dispMap, lagRaw, streamOf, o)
%MODELCHECKS  The model assumptions the data can contradict (WP-154, review
%   12.84 M2-M4), read after the loop; each warns once and is reported in
%   out.checks (one row per stream / group):
%     M2 a pool group of two or more "law" streams whose member is not white
%        under the group's law (Ljung-Box p < 1e-3, the last pass): one law
%        for different spectra -- measured eff 36, levels +96 % / -22 %
%     M3 a stream without arc factors whose arcs' levels scatter beyond their
%        sampling error (levelDispersion above its bound): a level that changes
%        from arc to arc -- measured eff 4.5, the level 1.33 x the truth (0.5 /
%        1 / 2 by day); levels 1 / 1.1 / 1.2 by day are flagged too (D 3.3-4.9
%        against the bound 4), a homogeneous stream reaches 1.4. Only a stream
%        white under its law is judged: the bound is the sampling variance of
%        whitened residuals, and a misfit law (a pooled stream, M2) scatters
%        the ratios with no change of level
%     M4 a law whose correlation time -1 / log(rho_1) exceeds half the median
%        arc of its group: the arcs' locals absorb the long memory. MEASURED
%        on 256-epoch arcs (Grid 60, Smooth 1): phi 0.99 (true 100 epochs,
%        estimated 59-103) level 0.976 -- not flagged; phi 0.995 (200,
%        estimated 138-202) level 0.88 -- flagged
%   OUTPUT  checks  struct .stream, .group, .whitePval, .dispersion, .dispersionBound (S x 1),
%                   .lawGroup, .corrTime, .arcMedian (G x 1)
    streams = hist.streams(:);
    S = numel(streams);
    grp = cellfun(@(sk) groupKeyOf(o.Pool, sk), streams, 'UniformOutput', false);
    pol = cellfun(@(sk) policyOf(o.Estimate, sk), streams, 'UniformOutput', false);
    p = nan(S, 1);  if ~isempty(hist.whitePval), p = hist.whitePval(end, :).'; end
    D = nan(S, 1);  Dmax = nan(S, 1);
    for s = 1:S
        if isKey(dispMap, streams{s}), Dt = dispMap(streams{s});  D(s) = Dt(1);  Dmax(s) = Dt(2); end
    end
    law = strcmp(pol, 'law');
    for g = unique(grp(law)).'
        in = law & strcmp(grp, g{1});
        bad = find(in & p < 1e-3);
        if nnz(in) >= 2 && ~isempty(bad)
            lsqvce.util.log.warnOnce('lsqvce:adjustVCEjoint:poolNotWhite', ...
                ['Pool group %s shares one law among %d streams, and %s is not white under it ' ...
                 '(Ljung-Box p = %.2g): the law does not fit that stream -- streams of different ' ...
                 'spectra in one group (which biases both levels and the formal errors), or ' ...
                 'something unmodelled in it (an event, a drift). Compare Pool = "stream" ' ...
                 '(R.noise.checks).'], g{1}, nnz(in), streams{bad(1)}, p(bad(1)));
        end
    end
    noAF = ~strcmp(pol, 'none') & ~cellfun(@(sk) arcFactorsOf(o.ArcFactors, sk), streams);
    % only a WHITE stream: the bound is the sampling variance of whitened
    % residuals, and a law that does not fit (a pooled group's stream, M2)
    % scatters its arcs' ratios without any change of level (D 4.6-6.9 seen)
    bad = find(noAF & D > Dmax & p >= 1e-3);
    if ~isempty(bad)
        lsqvce.util.log.warnOnce('lsqvce:adjustVCEjoint:levelVaries', ...
            ['Stream %s carries one level for all its arcs, but the arcs'' levels scatter %.3g ' ...
             'times their sampling variance (bound %.3g): the level changes from arc to arc (or day to ' ...
             'day). Give the stream arc factors (ArcFactors = true, or arcFactors on its ' ...
             'covariance spec) (R.noise.checks).'], streams{bad(1)}, D(bad(1)), Dmax(bad(1)));
    end
    gk = Cmap.keys;  G = numel(gk);
    tau = nan(G, 1);  mMed = nan(G, 1);
    arcGrp = cellfun(@(sk) groupKeyOf(o.Pool, sk), streamOf(:), 'UniformOutput', false);
    for gi = 1:G
        Cg = Cmap(gk{gi});
        m = cellfun(@numel, lagRaw.om(strcmp(arcGrp, gk{gi})));
        if numel(Cg) < 2 || isempty(m) || ~(Cg(2) > 0 && Cg(2) < Cg(1)), continue, end
        tau(gi) = -1 / log(Cg(2) / Cg(1));  mMed(gi) = median(m);
    end
    bad = find(tau > mMed / 2);
    if ~isempty(bad)
        lsqvce.util.log.warnOnce('lsqvce:adjustVCEjoint:memoryVsArc', ...
            ['The law of %s has a correlation time of %.3g epochs, more than half of its ' ...
             'median arc (%d epochs): the arcs'' local parameters absorb the long memory and ' ...
             'its level comes out low (measured: -12 %% at 200 against 256 epochs). Longer ' ...
             'arcs, or a law that models the low frequencies, read it (R.noise.checks).'], ...
             gk{bad(1)}, tau(bad(1)), round(mMed(bad(1))));
    end
    checks = struct('stream', {streams}, 'group', {grp}, 'whitePval', p, 'dispersion', D, ...
                    'dispersionBound', Dmax, ...
                    'lawGroup', {gk(:)}, 'corrTime', tau, 'arcMedian', mMed);
end

%==========================================================================
function [om, sn, r, why, wsum, wm, qw, tw, zE, dE] = oneArc(st, obsType, arc, day, uo, UiZ, x, ySys, Gsys, LlS, rc, sk)
%ONEARC  Decorrelate one arc, form its residual, and return the per-lag Omega
%   and s from LSQVCE.KERNEL.REDUNDANCYLAGS over ALL the arc's lags, and r, the
%   first column of the arc's law covariance (without the robust nugget; R7
%   WP-111 part 3). WP-34: a skipped arc REPORTS why --
%   `why` is '' on success, 'gap' (no data) or 'nonPD' (covariance not
%   positive definite this pass); silence was finding F-02's accomplice.
%   WP-33: ySys/Gsys/LlS are the system's local solution, Schur factor and
%   local Cholesky (computed once per pass by sweepArcs); the residual is
%   ebar = W^-T(l - A x - B y) and the projector absorbs the local block
%   exactly (jointArcParts + redundancyLags 'LocalBbar'/'LocalLl').
%   qw / tw: the arc factor's Helmert pair with the arc's WHOLE law as the
%   component (redundancyLags .omegaLaw / .sLaw, WP-132), summed over the axes.
%   WP-156: zE (mE x 1, [] when rc does not ask for it) the robust pass's
%   statistic per epoch and dE the arc's axis count (arcRobust); sn and tw
%   carry the truncation deficit once a robust pass has rejected.
    uo.obsType = obsType;
    uo.day     = day;
    uo.arc     = arc;
    why = '';
    wsum = zeros(WHITE_LAGS + 1, 1);  wm = zeros(WHITE_LAGS + 1, 1);   % R5 M4 4a / R7-38: sums and pair counts
    qw = 0;  tw = 0;  zE = [];  dE = 1;
    a = st.dataFcn(day, obsType, arc, uo);
    if isempty(a) || ~isfield(a, 'l') || isempty(a.l)
        % WP-127: a local prior group's pseudo-arc is listed on every pass (the
        % arc identity is pinned on pass 1) but carries rows only from its
        % 'from' pass on -- empty before that is INACTIVE, not a data gap
        why = 'gap';  if lsqvce.util.priorStream.is(obsType), why = 'inactive'; end
        om = [];  sn = [];  r = [];  return
    end

    A = a.A;  l = a.l(:);
    m = size(A, 1);
    B = zeros(m, 0);
    if isfield(a, 'B') && ~isempty(a.B), B = a.B; end
    nAx = lsqvce.util.getfielddef(a, 'nAxis', 1);
    hasCep = isfield(a, 'Cep') && ~isempty(a.Cep);
    if nAx > 1 || hasCep
        % R5 M3 (finding H-15): a d-axis stream is d single-axis series that
        % share the TEMPORAL law once the FIXED per-epoch covariance has been
        % removed (the separable model). So: descale every epoch by Cep^-1/2
        % (rows stay epoch-major), then treat axis k = rows k:d:m as one arc of
        % m/d epochs against the epoch-indexed temporal covariance a.Cfcn, and
        % ADD the per-lag statistics over the axes -- the lag is a lag in
        % epochs, the dense factor is (m/d)^2 instead of m^2, and one arc
        % factor covers the arc. R6 WP-90: the branch is entered whenever the
        % arc CARRIES per-epoch error information -- a scalar stream with a
        % per-epoch sigma (d = 1) is descaled exactly like a 3-axis one; a
        % multi-axis arc without Cep (a raw loader) is D = I.
        mE = m / nAx;
        if mE ~= round(mE)
            error('lsqvce:adjustVCEjoint:badAxisRows', ...
                  'Arc (sys %s, %s, day %s, arc %s): %d rows are not a multiple of %d axes.', ...
                  string(st.sysID), string(obsType), string(day), string(arc), m, nAx);
        end
        Cep = ones(mE, nAx);
        if hasCep, Cep = a.Cep; end
        % pure descale: AR(0) (a = 1, s = 1, r = 1) through the separable kernel
        [Ad, Bd, ld] = lsqvce.kernel.separableWhiten(1, 1, 1, Cep, nAx, A, B, l, true);
        Sig = a.Cfcn(1:mE, uo);                       % TEMPORAL, over epochs
        if isvector(Sig), Sig = diag(Sig(:)); end
        Sig = (Sig + Sig.')/2;
        [W, pd] = chol(Sig);
        if pd > 0
            om = [];  sn = [];  r = [];  why = 'nonPD';  return
        end
        % R7 WP-115 part 3b (R7-12, user decision 2026-10-11): the component
        % update sees the REJECTED epochs only -- factor 0, floored to 1e-12 by
        % the loader, i.e. q = 1/f - 1 ~ 1e12 -- as nuggets, which deletes
        % them to rounding; partial factors weighted the solution (the
        % accumulator's full nugget) and are NOT repeated here, because the
        % influence selects exactly the epochs with large residuals and the
        % level of a clean stream came back 14-19 % low with them inside
        % (HANDOVER 12.43). The truncation this leaves is corrected per
        % stream (truncationFactor).
        % WP-152 (N4): the EXCLUDED rows (arc.excluded, weight 0) are deleted
        % the same way, per axis: under a 'diag' or 'var' law the loader leaves
        % them as zero rows without a nugget, and the sweep counted each as an
        % observation with a zero residual and full redundancy -- 10 % excluded
        % came back as a level 0.926 against 1.018 for its clean twin
        WL = W;
        q = lsqvce.util.getfielddef(a, 'nugget', []);
        rej0 = false(mE, 1);
        if ~isempty(q), rej0 = q(:) > 1e11; end
        exAx = excludedAxes(a, mE, nAx);
        om = zeros(mE, 1);  sn = zeros(mE, 1);  r = Sig(:, 1);
        rejW = [];
        for k = 1:nAx
            rej = rej0 | exAx(:, k);
            if ~isequal(rej, rejW)                    % one factor per deletion pattern
                [W, cN, pd] = deletedFactor(Sig, WL, rej, q);
                if pd > 0
                    om = [];  sn = [];  r = [];  why = 'nonPD';  return
                end
                rejW = rej;
            end
            rows = k:nAx:m;
            [AtilBar, Bbar, ebar] = jointArcParts(Ad(rows, :), Bd(rows, :), ld(rows), W, ...
                                                  x(1:size(A,2)), ySys, Gsys);
            rl = lsqvce.kernel.redundancyLags(AtilBar, W, UiZ, ebar, mE, ...
                                       'LocalBbar', Bbar, 'LocalLl', LlS, 'Nugget', cN, 'LawW', WL);
            om = om + finiteOf(rl.omega);
            sn = sn + finiteOf(rl.s);
            qw = qw + rl.omegaLaw;  tw = tw + rl.sLaw;      % WP-132: the whole law
            wsum = wsum + whiteLagSums(ebar);  wm = wm + lagPairs(numel(ebar));
        end
        [zE, du, dt] = arcRobust(rc, sk, st, obsType, day, arc, Sig, Ad, Bd, ld, x, ySys, exAx, numel(sn));
        sn = sn - du;  tw = tw - dt;  dE = nAx;
        return
    end
    if isfield(a, 'Cfcn') && ~isempty(a.Cfcn)
        Sig = a.Cfcn(1:m, uo);
        if isvector(Sig), Sig = diag(Sig(:)); end
        Sig = (Sig + Sig.')/2;
        [WL, pd] = chol(Sig);
        if pd > 0
            om = [];  sn = [];  r = [];  why = 'nonPD';  return   % skipped, REPORTED
        end
        [W, cN, pd] = deletedFactor(Sig, WL, excludedAxes(a, m, 1), []);   % WP-152 (N4)
        if pd > 0
            om = [];  sn = [];  r = [];  why = 'nonPD';  return
        end
    else
        % A Wfcn-only loader hands back a whitening OPERATOR, not the explicit
        % triangular W that redundancyLags needs for Rtilde (eq. 6.5.48) and
        % etilde (eq. 6.5.54). Rather than silently drop the arc -- which would
        % bias the estimate by an invisible amount -- say so.
        error('lsqvce:adjustVCEjoint:needCfcn', ...
             ['Arc (sys %s, %s, day %s, arc %s) supplies only a Wfcn whitening ' ...
              'operator. The joint spectral estimator needs the EXPLICIT arc ' ...
              'covariance to form W with Sigma = W''W, so the loader must also ' ...
              'provide Cfcn. Most LSQVCE.ENGINE.COVMODEL types do; a hand-written ' ...
              'loader may need one added.'], string(st.sysID), string(obsType), ...
              string(day), string(arc));
    end

    % Step 7 with the FULL Table-6.2 vector: subtract B*y, reduce the common
    % design (Atil = A - B G) so the shared probe matches the factored N, and
    % hand the exact local block to the projector (WP-33, fixes F-03).
    [AtilBar, Bbar, ebar] = jointArcParts(A, B, l, W, x(1:size(A,2)), ...
                                          ySys, Gsys);
    rl = lsqvce.kernel.redundancyLags(AtilBar, W, UiZ, ebar, m, ...
                               'LocalBbar', Bbar, 'LocalLl', LlS, 'Nugget', cN, 'LawW', WL);
    om = finiteOf(rl.omega);
    sn = finiteOf(rl.s);
    r  = Sig(:, 1);
    qw = rl.omegaLaw;  tw = rl.sLaw;                    % WP-132: the whole law
    wsum = whiteLagSums(ebar);  wm = lagPairs(numel(ebar));
    [zE, du, dt] = arcRobust(rc, sk, st, obsType, day, arc, Sig, A, B, l, x, ySys, excludedAxes(a, m, 1), numel(sn));
    sn = sn - du;  tw = tw - dt;
end

function [zE, du, dt] = arcRobust(rc, sk, st, obsType, day, arc, Sig, A, B, l, x, y, exAx, nL)
%ARCROBUST  WP-156 (design F): what the robust pass needs from one arc, on the
%   law Sig the arc's sweep uses (temporal, over epochs) and its explicit
%   inverse (lsqvce.kernel.lawInverse, O(m^2) for a Toeplitz law):
%     zE  m_E x 1  Baarda's w per epoch with iterative data snooping
%                  (lsqvce.kernel.wSnoop) in units of the stream's scale
%                  (rc.scale / sqrt(rc.level): level-free), when rc.on -- the residual of the
%                  full fit l - A x - B y per axis (rows epoch-major, A / B / l
%                  descaled by Cep), the arc's local columns as the parameters;
%                  the set rebuilt from the data's own exclusions every pass
%     du / dt      the truncation deficit of the redundancy per lag and of the
%                  whole law's pair (lsqvce.kernel.truncationDeficit), once a
%                  robust pass has rejected (rc.tauBar), on the epochs the
%                  arcs were weighted with: rc.W's zeros deleted
%   The robust zeros and the data's exclusions both reach the arc as
%   arc.excluded (WP-133); rc.W tells them apart. A local prior group's
%   pseudo-arc is never re-weighted (WP-134): zE = [], no deficit.
    zE = [];  du = zeros(nL, 1);  dt = 0;
    if ~(rc.on || ~isempty(rc.tauBar)) || lsqvce.util.priorStream.is(char(string(obsType))), return, end
    [mE, nAx] = size(exAx);
    key = lsqvce.util.rowWeightKey(st.sysID, day, obsType, arc);
    rr = false(mE, 1);
    if isKey(rc.W, key), f = rc.W(key);  rr = f(:) <= 0; end
    dataEx = logical(exAx);  dataEx(rr, :) = false;
    act = any(B ~= 0, 1);
    E = zeros(mE, nAx);  X = cell(1, nAx);
    nc = size(A, 2);
    for k = 1:nAx
        rows = k:nAx:size(A, 1);
        e = l(rows) - A(rows, :) * x(1:nc);
        if any(act), e = e - B(rows, act) * y(act); end
        E(:, k) = e;
        X{k} = full(B(rows, act));
    end
    P = lsqvce.kernel.lawInverse(Sig);
    if rc.on
        s = 1;  if isKey(rc.scale, sk), s = rc.scale(sk); end
        L = 1;  if isKey(rc.level, sk), L = rc.level(sk); end
        zE = lsqvce.kernel.wSnoop(P, E, X, dataEx, rc.thr, s / sqrt(L));   % w units: w ~ 1 / sqrt(L)
    end
    if ~isempty(rc.tauBar)
        [du, dt] = lsqvce.kernel.truncationDeficit(P, X, dataEx, rr, rc.tauBar(nAx), nL);
    end
end

function ex = excludedAxes(a, mE, nAx)
%EXCLUDEDAXES  arc.excluded (m x 1, epoch-major) as an mE x nAx logical: the
%   rows of weight 0 per axis; all false when the arc marks none (WP-152).
    ex = false(mE, nAx);
    if isfield(a, 'excluded') && ~isempty(a.excluded)
        ex = reshape(logical(a.excluded(:)), nAx, mE).';
    end
end

function [W, cN, pd] = deletedFactor(Sig, WL, rej, q)
%DELETEDFACTOR  The Cholesky factor of Sig with the epochs REJ deleted by a
%   nugget (R7 WP-115 part 3b; WP-152 for the excluded rows): the loader's q
%   where it set one (a robust factor 0, or an excluded epoch under a
%   correlated law: q ~ 1e12), 1e12 elsewhere, times the epoch's variance --
%   deletion to rounding, the lag structure kept. No epoch deleted: W = WL
%   (the law's own factor) and cN = [], bitwise the undeleted arc.
%   INPUT   Sig mE x mE, WL its upper factor, rej mE x 1 logical, q [] or mE x 1
%   OUTPUT  W upper factor of Sig + diag(cN), cN [] or mE x 1, pd chol's flag
    pd = 0;  W = WL;  cN = [];
    if ~any(rej), return, end
    qq = 1e12 * ones(numel(rej), 1);
    if ~isempty(q), has = q(:) > 1e11;  qq(has) = q(has); end
    cN = zeros(numel(rej), 1);
    dSig = diag(Sig);
    cN(rej) = qq(rej) .* dSig(rej);
    [W, pd] = chol(Sig + diag(cN));
end

function K = WHITE_LAGS()
%WHITE_LAGS  Lags of the Ljung-Box whiteness statistic per stream (20).
    K = 20;
end

function c = lagPairs(n)
%LAGPAIRS  Number of products e(i) e(i+k) in one arc's lag-k sum, k = 0..K:
%   max(n - k, 0). Pooled over a stream's arcs, these are the N_k the
%   Ljung-Box statistic must divide by (R7-38).
    c = max(n - (0:WHITE_LAGS()).', 0);
end

function s = whiteLagSums(e)
%WHITELAGSUMS  sum_i e(i+n) e(i), n = 0..K, of one arc's law-whitened residual
%   (zero-padded FFT, exact to rounding); aggregated per stream by sweepArcs.
    e = e(:);  K = WHITE_LAGS();
    n = numel(e);
    s = zeros(K + 1, 1);
    if n < 2, return, end
    nfft = 2^nextpow2(2*n - 1);
    F = fft(e, nfft);
    r = ifft(real(F .* conj(F)), 'symmetric');
    k = min(K + 1, n);
    s(1:k) = r(1:k);
end

%==========================================================================
function y = finiteOf(v)
%FINITEOF  A per-lag vector as a column with non-finite entries set to zero.
    y = v(:);
    y(~isfinite(y)) = 0;
end

function s = skipSummary(skips)
%SKIPSUMMARY  'gap: 3, nonPD: 12' from the {key, reason} skip list.
    if isempty(skips), s = 'none'; return, end
    [why, ~, j] = unique(skips(:, 2));
    n = accumarray(j(:), 1);
    parts = cell(1, numel(why));
    for k = 1:numel(why), parts{k} = sprintf('%s: %d', why{k}, n(k)); end
    s = strjoin(parts, ', ');
end

function U = probeFactor(acc, it)
%PROBEFACTOR  The square root N = U'U the shared probe uses (WP-134, review 3
%   finding 7): acc.commonFactor(), the backend's own factor -- QR's R1 with
%   the prior folded in by QR, or the chol the normal backend's solve() used.
%   The probe needs ANY square root of N; chol of the SQUARED reduced normals
%   refused runs of cond(A) ~ 3e8 that the QR solve accepts (up to ~4e12),
%   and where it succeeded the probe term carried its error (4-25 %).
%   ERRORS  lsqvce:adjustVCEjoint:notPD  no factor, or a pivot the solve refuses
    why = '';
    try
        U = acc.commonFactor();
    catch ME
        U = NaN;  why = sprintf(' (%s)', ME.message);    % the normal backend's chol refused
    end
    dU = abs(diag(U));
    % tol: 1e3*eps relative pivot -- QRSystem.RankTol, the pivot the solve
    % itself refuses, so the loop refuses nothing the solve accepted
    if isempty(dU) || ~all(isfinite(U(:))) || min(dU) <= 1e3 * eps * max(dU)
        error('lsqvce:adjustVCEjoint:notPD', ...
             ['The reduced common normal matrix (the prior included) has no factor at ' ...
              'outer iteration %d%s: it is singular to the solve''s own pivot tolerance. ' ...
              'Add constraints or a prior before estimating the stochastic model.'], it, why);
    end
end

function skipGuard(skips, keysAll, streamAll, it, thr)
%SKIPGUARD  Refuse a pass that skipped more than thr of ONE stream's arcs
%   (WP-134, review 3 finding 15: the rate was over all arcs). skips S x 2
%   {key, reason}; keysAll 1 x M and streamAll M x 1 aligned (sweepArcs).
    if isempty(skips) || isempty(keysAll), return, end
    [~, at] = ismember(skips(:, 1), keysAll);
    skStream = streamAll(at);
    for sk = unique(skStream(:)).'
        nS = sum(strcmp(streamAll, sk{1}));
        mine = strcmp(skStream, sk{1});
        if nnz(mine) / nS > thr
            error('lsqvce:adjustVCEjoint:skipThreshold', ...
                 ['Pass %d skipped %d of the %d arcs of stream %s (%.0f %% > SkipThreshold ' ...
                  '%.0f %%): %s. A noise law or level estimated from its remaining arcs would ' ...
                  'not describe the stream; fix its covariance model or raise ' ...
                  'Options.SkipThreshold deliberately.'], it, nnz(mine), nS, sk{1}, ...
                  100 * nnz(mine) / nS, 100 * thr, skipSummary(skips(mine, :)));
        end
    end
end

function L = logLaws(Smap)
%LOGLAWS  A snapshot of every pool group's law, log domain (AvgTail, WP-134).
    L = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for g = keys(Smap)
        L(g{1}) = log(max(Smap(g{1}), realmin));
    end
end

function [Smap, Cmap] = tailAverage(Smap, Cmap, tail, Kc)
%TAILAVERAGE  Polyak average of every pool group's law over the tail passes
%   (log domain, r(0) = 1 again); a group not estimated on every tail pass
%   keeps its last law (WP-134).
    for g = keys(Smap)
        if ~all(cellfun(@(t) isKey(t, g{1}), tail)), continue, end
        Lg = cellfun(@(t) t(g{1}), tail, 'UniformOutput', false);
        Sg = exp(mean(cat(2, Lg{:}), 2));
        Smap(g{1}) = Sg / (Kc(1, :) * Sg);
        Cmap(g{1}) = Kc * Smap(g{1});
    end
end

function [keysAll, rowAll, streamAll] = enumerateArcs(acc)
%ENUMERATEARCS  The arc identity sweepArcs pins -- arc keys, loader keys and
%   stream keys of every enumerated arc, in sweepArcs' order -- from the
%   loader's 'meta' calls alone (WP-132: lets Sigma2_0 weight pass 1). Any
%   disagreement with sweepArcs' own enumeration is refused by pass 1's
%   arcMismatch check. OUTPUT keysAll 1 x M cell, rowAll / streamAll M x 1 cell.
    keysAll = {};  rowAll = {};  streamAll = {};
    ids = acc.systemIDs();
    for k = 1:numel(ids)
        st = acc.getStore(ids(k));
        for d = st.dayList(:).'
            meta = st.dataFcn('meta', d, [], st.uo);
            for ot = 1:numel(meta.obsTypes)
                obsType = meta.obsTypes{ot};
                for a = meta.arcs{ot}(:).'
                    keysAll{end+1} = lsqvce.util.arcKey(st.sysID, obsType, d, a);        %#ok<AGROW>
                    streamAll{end+1, 1} = lsqvce.util.streamKey(st.sysID, obsType);       %#ok<AGROW>
                    rowAll{end+1, 1} = lsqvce.util.rowWeightKey(st.sysID, d, obsType, a); %#ok<AGROW>
                end
            end
        end
    end
end

function s2Map = pinFactors(pinnedKeys, sigma2_0)
%PINFACTORS  The keyed start factors: Sigma2_0 aligned to the pinned arcs, or
%   all ones without it. Raises badSigma2_0 on a count mismatch.
    s2Map = containers.Map('KeyType', 'char', 'ValueType', 'double');
    if isempty(sigma2_0)
        for kk = 1:numel(pinnedKeys), s2Map(pinnedKeys{kk}) = 1; end
        return
    end
    if numel(sigma2_0) ~= numel(pinnedKeys)
        error('lsqvce:adjustVCEjoint:badSigma2_0', ...
            ['Sigma2_0 has %d entries but pass 1 enumerated %d ' ...
             'arcs (see out.arcKeys for the pinned order).'], ...
             numel(sigma2_0), numel(pinnedKeys));
    end
    for kk = 1:numel(pinnedKeys)
        s2Map(pinnedKeys{kk}) = sigma2_0(kk);
    end
end

function M = streamMeans(arcStreams, sigma2)
%STREAMMEANS  containers.Map 'satKey|obsType' -> mean arc factor of the stream.
%   arcStreams are the pinned arcs' STREAM keys, aligned to sigma2 (R7-43:
%   carried from sweepArcs, no longer parsed out of the arc keys); empty
%   before pass 1 pins them (then every stream's level is its declared
%   sigma, factor 1).
    M = containers.Map('KeyType', 'char', 'ValueType', 'double');
    if isempty(arcStreams) || isempty(sigma2) || numel(arcStreams) ~= numel(sigma2), return, end
    acc = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for k = 1:numel(arcStreams)
        sk = arcStreams{k};
        if ~isKey(acc, sk), acc(sk) = []; end
        v = sigma2(k);
        if isfinite(v) && v > 0, acc(sk) = [acc(sk); v]; end
    end
    ks = keys(acc);
    for k = 1:numel(ks)
        v = acc(ks{k});
        if ~isempty(v), M(ks{k}) = mean(v); end
    end
end

function F = arcFactorMap(arcStreams, rowKeys, sigma2, o)
%ARCFACTORMAP  containers.Map rowWeightKey -> the arc's variance factor
%   RELATIVE to its stream's level (sigma2_m / mean over the stream), for the
%   arcs of every stream that carries arc factors and is estimated (R7-13).
%   LSQVCE.ENGINE.JOINTMODELOPT turns it into covOverride.byArc; the stream's
%   level rides on bySys as before, so level x factor = sigma2_m is what the
%   loader weights the arc with. Empty before pass 1 pins the arcs.
    F = containers.Map('KeyType', 'char', 'ValueType', 'double');
    if isempty(arcStreams) || numel(rowKeys) ~= numel(arcStreams) || numel(sigma2) ~= numel(arcStreams), return, end
    lvl = streamMeans(arcStreams, sigma2);
    for k = 1:numel(arcStreams)
        sk = arcStreams{k};
        if ~arcFactorsOf(o.ArcFactors, sk) || strcmp(policyOf(o.Estimate, sk), 'none'), continue, end
        if ~isKey(lvl, sk) || ~(isfinite(sigma2(k)) && sigma2(k) > 0), continue, end
        F(rowKeys{k}) = sigma2(k) / lvl(sk);
    end
end

function arcSpreadGuard(spreadHist)
%ARCSPREADGUARD  Refuse to return a DIVERGING set of arc factors in silence.
%   R6 WP-97, kept as a RATCHET by WP-98 (HANDOVER 12.21, 12.22). An arc
%   variance factor is a RELATIVE weight inside its stream, so over the passes
%   its spread must settle. Finding F-01 was that it did not: eq. 6.5.28
%   multiplied by a factor the weighting never used, which made the update a
%   geometric recursion with no fixed point (measured 8.4x per pass). WP-98
%   fixed that, so on this tree nothing provokes this guard any more.
%
%   IT STAYS ANYWAY, and deliberately: it is what would catch this class
%   returning on a user's own data, where no test of ours runs. Its FIRING
%   direction was demonstrated against the pre-fix code and is on the record
%   with numbers -- 3/3 at spreads 4670 / 5600 / 6109 with a degraded arc,
%   silent 3/3 without one. Do not delete it because it is quiet: quiet is now
%   the correct reading, and dev/tests pins that reading.
%
%   THE RULE, and where its numbers come from: fire only when the spread BOTH
%   grew at every one of the last three passes AND ended at least 4x its first
%   pass. Measured silent on the shipped demos even before the fix:
%   demo_grace 4.13 -> 3.27 over six passes.
    ks = keys(spreadHist);
    for i = 1:numel(ks)
        r = spreadHist(ks{i});
        if numel(r) < 3, continue, end
        if all(diff(r(end-2:end)) > 0) && r(end) >= 4 * r(1)
            lsqvce.util.log.warnOnce('lsqvce:adjustVCEjoint:arcFactorsDiverging', ...
               ['Stream %s: the spread of its per-arc variance factors grew from %.3g to ' ...
                '%.3g over %d passes and was still growing at the last one. Per-arc factors ' ...
                'that compound instead of settling are not an estimate (the class of finding ' ...
                'F-01). Take the ONE-PASS estimate (MaxPasses = 1) for per-arc factors, or set ' ...
                'arcFactors = false on this stream and read its level instead. A stream with ' ...
                'arcFactors = false is unaffected; the level of THIS stream is the mean of ' ...
                'these factors, so it drifts with them.'], ...
                ks{i}, r(1), r(end), numel(r));
        end
    end
end

function p = policyOf(E, sk)
%POLICYOF  'none' | 'scale' | 'law' for a stream; an empty map = 'law' for everybody.
    p = 'law';
    if isempty(E), return, end
    if isKey(E, sk), p = lower(char(E(sk))); end
end

function tf = arcFactorsOf(A, sk)
%ARCFACTORSOF  whether the stream's arcs get factors; an empty map = true.
    tf = true;
    if isempty(A), return, end
    if isKey(A, sk), tf = logical(A(sk)); end
end

function g = groupKeyOf(pool, sk)
%GROUPKEYOF  the pool group a stream's law belongs to.
    switch lower(char(pool))
        case 'all'
            g = 'all';
        case 'obstype'
            g = extractAfter(sk, '|');
        otherwise
            g = sk;                                  % 'stream'
    end
end

function [streams, sig, white, pval] = whitenessRow(wstats, pinnedStream, sigma2)
%WHITENESSROW  One convergence-table row: per stream, the level factor
%   sqrt(mean arc factor) and the Ljung-Box test on the pass's law-whitened
%   residuals, POOLED over the stream's arcs: rho_k = s_k / s_0 from the
%   summed lag products, Q = N (N+2) sum_k rho_k^2 / N_k ~ chi2_K' under
%   whiteness, with N the residual count, N_k = sum over arcs of max(m_a - k,
%   0) the products actually in s_k, and K' the lags with N_k > 0; white =
%   p > 0.05. For ONE arc this is Ljung-Box exactly. R7-38: it divided by
%   (N - k), as if the stream were one series, while each arc loses k
%   products -- on 64 arcs of 30 epochs a white stream had size 0.001 (E[Q]
%   13 instead of 20) and on 200 x 12 epochs 0.000; pooled, 0.047-0.058 on
%   every layout (mirrors/wp118/ljungbox_pooled.py, 1000-2000 trials).
%   WP-127: the columns are the PINNED streams (pass 1, sorted), not the
%   streams that delivered residuals THIS pass -- a stream with no surviving
%   arc in one pass (a prior group before its 'from', a whole stream skipped
%   nonPD) changed the row width and the history assignment failed. Such a
%   stream gets its level, white = false and p = NaN.
    streams = unique(pinnedStream(:)).';         % sorted
    sMeans  = streamMeans(pinnedStream, sigma2);
    S = numel(streams);
    sig = ones(1, S);  white = false(1, S);  pval = nan(1, S);
    K = WHITE_LAGS();
    for j = 1:S
        if isKey(sMeans, streams{j}), sig(j) = sqrt(max(sMeans(streams{j}), eps)); end
        if ~isKey(wstats, streams{j}), continue, end
        w = wstats(streams{j});
        sk = w(1:K+1);  Nk = w(K+2:end);              % lag sums; pair counts (R7-38)
        m = Nk(1);  r0 = sk(1);
        ok = Nk(2:end) > 0;  Ke = nnz(ok);
        if Ke == 0 || r0 <= 0, continue, end
        rho = sk(2:end) / r0;
        Q = m * (m + 2) * sum(rho(ok).^2 ./ Nk([false; ok]));
        pval(j)  = gammainc(Q/2, Ke/2, 'upper');     % 1 - chi2cdf(Q, K'), base MATLAB
        white(j) = pval(j) > 0.05;
    end
end

%==========================================================================
function [binMax, blocking] = binsAgainstBound(pvg, tol, group, binMax, blocking)
%BINSAGAINSTBOUND  The bins' stop test of one pool group (WP-147).
%   A bin has SETTLED when it moved by at most tol max(1, sqrt(200 / r_j)):
%   tol for a bin of redundancy >= 200, below that a bound that grows with
%   the bin's own standard error sqrt(2 / r_j) -- at the default tol 1e-2 a
%   tenth of it. Scaled BY tol, not floored at a tenth, so a tighter tol
%   still tightens every bin (a floor would have made tol 1e-3 mean nothing
%   for them). The plain max|alpha2 - 1| <= tol (WP-152: now the shape change
%   max|alpha2 / abar - 1|, psdVCE) asked bins the data barely
%   determine -- the zero-frequency bin the per-arc locals absorb, 60 bins on
%   60-sample arcs -- for a precision a hundred times finer than their
%   estimate has: 14-168 passes on the fixtures of HANDOVER 12.79, while x
%   was within 0.1-0.2 of its standard error of the fixed point after 10. A
%   run never stops later than under the plain rule (the bound is >= tol).
%
%   INPUT   pvg  1x1 struct  LSQVCE.ENGINE.PSDVCE's output (.dev, .redundancy)
%           tol  1x1 double  the loop's Tol;  group  1xK char  the pool group
%           binMax 1x1, blocking B x 6 cell  the running values of this pass
%   OUTPUT  binMax  1x1  the largest dev_j / bound_j so far (<= 1: settled)
%           blocking (B+b) x 6 cell  {group, bin, L, r_j, median r, dev_j} for
%                   every bin of this group above its bound, worst first
    bound = tol * max(1, sqrt(200 ./ max(pvg.redundancy, realmin)));
    ratio = pvg.dev ./ bound;
    binMax = max([binMax; ratio]);
    [rs, ord] = sort(ratio, 'descend');
    rMed = median(pvg.redundancy(pvg.redundancy > 0));
    for j = reshape(ord(rs > 1), 1, [])
        blocking(end+1, :) = {group, j, numel(ratio), pvg.redundancy(j), rMed, pvg.dev(j)}; %#ok<AGROW>
    end
end

function [hist, xPrev] = recordStep(hist, it, pv, sol, xPrev)
%RECORDSTEP  The pass's stop quantities into the history (WP-147).
%   .maxDev (max abs(alpha2 / abar - 1), a reading), .binRatio (the bins' stop
%   quantity) and .xStep: the solution's largest change since the pass
%   before in units of its own standard error sqrt(diag(Cxx)) sigma0 -- the
%   number a user needs to judge a stop on the pass budget; NaN on pass 1
%   or when the sizes do not match.
%   INPUT   hist struct, it 1x1, pv 1x1 struct (.maxDev .binRatio), sol the
%           pass's solution (.x .Cxx .sigma0), xPrev the previous x (or [])
%   OUTPUT  hist with row it filled; xPrev = this pass's x
    hist.maxDev(it, 1) = pv.maxDev;
    hist.binRatio(it, 1) = pv.binRatio;
    sx = sqrt(max(diag(sol.Cxx), 0)) * sol.sigma0;
    step = NaN;
    if numel(xPrev) == numel(sol.x) && numel(sx) == numel(sol.x) && any(sx > 0)
        ok = sx > 0;
        step = max(abs(sol.x(ok) - xPrev(ok)) ./ sx(ok));
    end
    hist.xStep(it, 1) = step;
    xPrev = sol.x(:);
end

function msg = notSettledMessage(maxIter, tol, hist, blocking)
%NOTSETTLEDMESSAGE  The not-converged warning's text (WP-147): what still
%   moves -- the bins above their bound with their redundancy against their
%   group's median, at most three named -- and the solution's last step.
%   INPUT   maxIter, tol 1x1;  hist the history;  blocking B x 6 cell (above)
%   OUTPUT  msg 1xK char
    if isempty(blocking)
        moving = 'no spectral bin';
    else
        nb = size(blocking, 1);
        shown = arrayfun(@(q) sprintf('%s bin %d of %d (redundancy %.3g, median %.3g, moved %.2e)', ...
                    blocking{q, :}), 1:min(3, nb), 'UniformOutput', false);
        moving = strjoin(shown, '; ');
        if nb > 3, moving = sprintf('%s; and %d more', moving, nb - 3); end
    end
    msg = sprintf(['The joint stochastic model did not settle in %d outer iterations ' ...
        '(last shape change %.3e, last max level change %.3e, last prior change ' ...
        '%.3e, tol %.1e). Still moving by more than tol max(1, sqrt(200 / redundancy)) ' ...
        '(at tol 1e-2 a tenth of their standard error): %s. The solution moved by %.3g ' ...
        'of its standard error in the last pass. Inspect out.history.S, .binRatio, ' ...
        '.levelDev and .priorDev: a spectrum that keeps moving usually means too few ' ...
        'arcs per frequency, a prior that keeps moving a declared sigma far from the data.'], ...
        maxIter, hist.maxDev(end), hist.levelDev(end), hist.priorDev(end), tol, ...
        moving, hist.xStep(end));
end

function last = lastSolOf(acc, sol)
%LASTSOLOF  The previous pass's solution as the accFcn's model.lastSol sees it.
%   .x, .Cxx and .s0sq as solved; when the backend refused a cancelled v'Pv
%   (s0sq NaN, and with it every entry of Cxx) also .Ninv = N^-1, the commons'
%   cofactor (BlockAccumulator.commonCofactor), which is finite. WP-144 (A6):
%   fit's prior "scale" step took N^-1 = Cxx / s0sq = NaN, held the group at
%   its declared sigma and reported PriorDev 0 -- the loop stopped
%   "converged" with the prior never estimated (reproduced: 1e6 held, the QR
%   backend found 4.9e4). Computed only in that case: Cxx already carries it.
%   INPUT   acc  the pass's accumulator;  sol  its solution struct
%   OUTPUT  last 1x1 struct  .x (nc x 1) .Cxx .s0sq [.Ninv nc x nc]
    last = struct('x', sol.x(:), 'Cxx', sol.Cxx, 's0sq', sol.s0sq);
    if ~isfinite(sol.s0sq)
        last.Ninv = acc.commonCofactor();
    end
end
