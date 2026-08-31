% LSQVCE  Batch least squares: accumulation, variance components, and judgment.
% Version 1.0  15-Jul-2026
%
% One flat package. Everything is a batch-least-squares tool at the same level;
% the two accumulators differ only in problem structure:
%   Accumulator       flat parameters      (general adjustments)
%   BlockAccumulator  common + local       (block-arrow, e.g. gravity recovery)
%
% Getting started
%   lsqvce_setup   - Add to path and register the Help browser docs
%   doc lsqvce     - The toolbox landing page
%   demo_workflow  - How to READ the results (start here)
%   demo_lsqvce    - Six failure modes, one per scenario
%
% Solving
%   Accumulator      - Streaming normal-equation / QR accumulator (flat)
%   BlockAccumulator - Streaming block-arrow accumulator: common parameters
%                      shared across systems + per-system locals pre-eliminated
%                      by Schur complement (gravity recovery, sensor networks)
%   adjust           - One-call driver: accumulate, solve, diagnose, judge
%   solveNonlinear   - Re-linearisation loop for nonlinear problems
%   diagnoseAccum        - Diagnose a BlockAccumulator adjustment with the full
%                      five-level verdict (diagnoses the FULL parameter set,
%                      so r_i counts the DoF the local parameters absorb)
%   QRSystem, NormalEqSystem, ReducedContribution - accumulator backends
%
% Access layer (one-call front door over the streaming pipeline)
%   Options        - THE one settings object for fit (R4 WP-55): every knob
%                    named, typed, defaulted to the engines' own defaults;
%                    engineArgs(target) maps it onto one existing driver.
%   adjustBlock    - One line: accumulate -> solve -> diagnoseAccum -> diagnose,
%                    returning a Results object (the base of adjustVCE). The
%                    streaming block-arrow solve; pairs with ADJUST (flat) and
%                    ADJUSTVCE (iterative). (Was recover; that name was REMOVED
%                    2026-07-18.)
%   Results        - Handle object carrying x / sigma0 / diag / acc plus
%                    report() / redundancy() / residuals() / disp()
%   makeLoader     - Turn a day-builder buildDay(day,obsType) into an
%                    accumulator-ready Loader (meta / per-arc / one-day cache)
%   adjustVCE      - Iterative separable-GLS engine: stochastic model (temporal
%                    whitening) + variance components (VCE), fluent builder over adjustBlock.
%                    Pairs by name with the in-memory shortcut ADJUST.
%                    P.empiricalStochastic() runs empirical stochastic
%                    modelling end to end (homogenise -> loop -> validate);
%                    P.validate(R) is the formal-vs-empirical acceptance test.
%                    (Was SeparableProblem; that name was REMOVED 2026-07-18.)
%   resolveNc      - Derive the common-parameter count nc from the systems
%                    themselves (a Loader's 'commonParams'/'numCommon', or a
%                    raw spec's .nc), so nc is stated exactly ONCE. Used by
%                    adjustBlock and adjustVCE; call it directly only to
%                    inspect where each system's nc came from.
%
% Diagnostics
%   diagnose       - Per-observation v, r_i, w, MDB, reliability, contamination
%   report         - Five-level verdict with actionable messages
%   systematic     - Ljung-Box, runs test, covariate binning
%   residuals      - Second-pass v, whitened w, sigma0 (BlockAccumulator path)
%   redundancy     - Per-group r_i and variance components (block-arrow exact)
%   contribution   - Per-group information share of the solution: Level-1 scalar
%                    share r_g/sum(r_g) (exact) + Level-2 per-parameter diagonal
%                    diag(N^-1 N_g^cc) (streamed). "Which system determines what."
%   streamContribution - Matrix-free per-group N^-1 N_g^cc probing over ONE data
%                    sweep (never forms a group normal); the shared primitive
%                    behind contribution() Level 2 and vce's hutchinson trace
%   leverage, leverageBlocks, localDiagnostics - leverage / conditioning
%
% Outliers
%   snoop          - Iterative rejection via Cholesky downdate
%   robust         - Huber / bisquare / IGG3 reweighting
%   outliers       - Data snooping via partial redundancy (block-arrow path)
%
% Stochastic model / variance components
%   varComp        - Helmert variance components (flat), with a resolution limit
%   vce            - Helmert variance components per arc (block-arrow path)
%   vceHutchinson  - Helmert variance components via the Ellmer Hutchinson trace
%                    estimator: a SPEED-oriented sibling of varComp for LARGE
%                    common blocks (no dense inverse). Strict Ellmer: per-group
%                    components + estimated r_g ONLY -- no per-obs r_i, snooping,
%                    or sum(r_i)=m-u guard. Prefer varComp below a few thousand u.
%   poolGroups     - Suggest/build a COARSER VCE grouping when a group's
%                    redundancy is too low to estimate a reliable component
%                    (roadmap item 4). Returns a remap; the caller re-runs VCE.
%   traceEstimate  - Monte-Carlo (Hutchinson) trace estimator (Ellmer 2018
%                    6.5.6): tr(N^-1 N_g) without forming N^-1. OPT-IN
%                    accelerator for VCE at high nc; exact leverage is the
%                    default (see the accuracy note in lsqvce.leverage).
%   iterateNoiseModel     - Iterate the noise model until sigma0 ~ 1 and residuals white
%
% Non-parametric (Ellmer 2018) stochastic model — THREE LAYERS
%   Everything here is public and documented, but they are NOT all things you
%   call day to day. Read top-down: the DRIVER is the front door, the ADAPTER
%   only matters if you drive it yourself, and the PRIMITIVES are the machinery
%   — exposed because they are independently useful and because this is meant to
%   be a readable reference implementation, not because you must touch them.
%
%   -- DRIVER (start here) --
%   adjustVCEjoint - The JOINT stochastic model: Ellmer's 10-step loop (§6.5.7),
%                    co-estimating the PSD, the arc factors and the group scale
%                    in ONE Helmert VCE per pass. Contrast lsqvce.adjustVCE,
%                    which ALTERNATES shape and scale; the two can converge to
%                    DIFFERENT models because the components are coupled.
%                    Usually you want adjustVCE(...,'vceMode','joint') instead,
%                    which wires all of this up for you.
%
%   -- ADAPTER (only if you drive adjustVCEjoint yourself) --
%   jointModelOpt  - Turn the driver's `model` struct into loader UserOpt.
%                    Encodes three conventions that otherwise fail SILENTLY (the
%                    loop diverges rather than erroring): sigma2 is a VARIANCE
%                    but a covariance spec wants a STANDARD DEVIATION; the
%                    override must land in covOverride.bySys.(sat).(obsType);
%                    and psdMode is read from the UserOpt, not the spec.
%
%   -- PRIMITIVES (the machinery; useful standalone) --
%   psdVCE         - Per-frequency PSD amplitude by VCE (eqs. 6.5.19-6.5.26).
%                    The non-parametric alternative to an AR/ARMA shape fit.
%   arcVarianceFactors - Relative per-arc weights within a group, with the gamma
%                    normalisation and robust mean (eqs. 6.5.27-6.5.30, 6.5.57).
%   redundancyLags - Per-lag redundancies s_n = tr(Rbar V_n) and residual sums
%                    Omega_n for one arc (eqs. 6.5.45-6.5.56).
%   sharedProbe    - The once-per-iteration shared block U^-1*Zbar (eq. 6.5.48),
%                    reused across every arc. Ellmer's central optimisation.
%   toeplitzCofactor - Symmetric Toeplitz cofactor from a lag vector (eq. 6.5.6),
%                    dense or matrix-free.
%   covComponent   - Banded lag-component V_n apply (eq. 6.5.16), never forming V_n.
%   dctI           - Type-I DCT, the covariance <-> PSD pairing (eq. 6.5.12).
%                    Orthonormal by default; 'Norm','ellmer' for bit-exact ITSG.
%
% Covariance and whitening
%   covModel - covariance spec -> {Cfcn, Wfcn}: diag/acf/psd/arma/var
%   armaWhiten      - exact FIR whitener (= chol(C)^-1)
%   separableWhiten - per-epoch D x temporal R_t (multi-axis positions)
%   varWhiten, whittle, varACV - vector-AR path (cross-axis temporal correlation)
%   decorrFactor    - dense W with W'W = inv(C)
%   crossCorr       - does the separable assumption hold? ('needs-var' verdict)
%   checkPSD        - verify the AR reduction reproduces your PSD
%
% Signal processing
%   psd2acv, levinsonDurbin, armaACV, ar2arma, selectArOrder
%   selectArAcv    - automatic AR order from a sample ACV (no grid); the
%                    ACV-domain sibling of selectArOrder
%
% Separable GLS helpers (the GRACE-FO worked case lives in examples/)
%   residualAcv - temporal ACV of the D^-1-descaled residuals
%   collectResiduals    - stack a residual field over the arcs of a stream
%   whitenessCorr          - largest |autocorrelation| over a lag range (whiteness)
%   resolveSatSpecs     - per-satellite config: overrides over defaults
%   normalizeSatList    - normalise a satellite list to a cell of ids
%   satKeyOf            - valid struct field name for a satellite id
%
% Your own case (not part of the package; call the toolbox via lsqvce.)
%   template/lsqvce_template          - START HERE. Returns the Loaders; you fill
%                    two documented stubs (myBuildDay, myEpochCov). The
%                    observationGroups loop is real code; the stubs error until
%                    written. See `doc lsqvce` -> "Your own case (template)".
%
% Examples (not part of the package; call the toolbox via the lsqvce. prefix)
%   examples/demo_grace_dataLoader      - the GRACE-FO reference setup: C/D pair
%                    (merged, range-rate couples both arc states) + independent
%                    satellite E on the same gravity field. A WORKED case to read
%                    alongside the template.
%   examples/demo_grace           - the driver that runs it end to end
%   examples/demo_gravity_design        - fully-normalised Legendre + potential
%                    partials along a sampled orbit. The ONE physically
%                    structured design generator; the demo and the oracle
%                    test suite both call it.
%   examples/demo_gravity_dataLoader    - GRACE-like fixture on that design:
%                    two satellites plus the energy-balance KBR difference,
%                    a DIFFERENT local count per system, a ragged arc
%                    schedule and a different length for every arc.
%   examples/demo_gravity         - the driver: weights, coefficient recovery
%                    in units of the formal error, error degree variances
%   examples/demo_messaging           - the SAME run at every verbosity level,
%                    side by side, so you can see what a long job will say
%
% Regularisation
%   CommonPrior, applyRegularisation, gcvLambda
%
% Model testing / free networks
%   ftest          - Nested-model F-test
%   ellipse        - Confidence ellipses
%   datum          - Null-space basis; minimum-constraint solutions
%
% Parallelism / IO
%   mapreduce      - serial | parfor | spmd strategy layer
%   PrefetchLoader - background one-day-ahead prefetch
%
% Support / multi-dataset
%   stat           - Distribution kernel; no Statistics Toolbox required
%   util.log       - Toolbox-wide logging: one file identifier (logfid) + one
%                    screen flag (logflag) route ALL status text and warnings
%   util.showExport - Show a figure briefly and export it as PNG (the engine
%                    side of the adjustVCE per-obsType plot hooks)
%   util.verbosity - Message LEVELS (0 silent | 1 progress | 2 stage |
%                    3 diagnostic) and the shared banner formatter, so every
%                    log line aligns and each level includes the ones below it
%   util.getfielddef - Field value with a default (struct access that tolerates
%                    a missing field, used across the option plumbing)
%   util.mergeStruct - Overlay the fields of one struct onto another (how the
%                    fluent configurators merge Name/Value into the opt structs)
%   util.hasBackgroundPool - True when backgroundPool can actually START
%                    (probe-by-USE; the ONE pool-availability decision)
%   Combine        - multi-dataset combination
%   BatchIndex     - batch bookkeeping
%   synthLoader    - synthetic ground truth and the LOADER TEMPLATE
%
% THE FIVE LEVELS OF THE REPORT
%   0 INTEGRITY   sum(r_i) == n-u.  Tests the CODE, not the data.
%   1 GLOBAL      variance factor vs chi-square.  READ THE RATIO, NOT THE VERDICT.
%   2 LOCALISE    data snooping, Bonferroni, contamination gauge.
%   3 RELIABILITY redundancy numbers, MDB, (w,r) quadrants.
%                 The level that says what you CANNOT check.
%   4 PARAMETERS  correlations, conditioning, precision.
%
% TWO HAZARDS THIS TOOLBOX EXISTS TO CATCH
%   1. Pope's tau goes BLIND at ~10% contamination: the blunders inflate s0,
%      tau divides by s0, and they hide one another. Iterating does not help.
%      Use 'test','mad' and heed the contamination gauge.
%   2. A small residual is NOT evidence of a good observation. A residual
%      reveals only the fraction r_i of any error: at r_i = 0 a 40-sigma
%      blunder gives |w| = 0.000. Read the QUADRANTS, not w alone.
%
% Provenance: merged from lsqAccumToolbox (accumulation) and lsqdiag
% (diagnostics) into one flat framework. Generated with assistance from
% Claude (Anthropic), Claude Opus 4.8.
