function tests = test_fit
%TEST_FIT  Tests for lsqvce.fit, the one user verb (R4 WP-56, decision R4-D4).
%
%     lsqvce.fit   data struct (DATA_CONTRACT.md) + lsqvce.Options -> Results;
%                  builds the loaders through makeLoader and dispatches on
%                  opts.Model to adjustBlock / adjustVCE
%
%   Claims under test: (1) EXACTNESS of the dispatch -- fit reproduces the
%   direct driver call to solver precision; (2) every rule of the v3 data
%   contract fires its documented lsqvce:fit id BEFORE any data is read;
%   (3) the declaration is order-free (fields, systems); (4) nl = 0 is the
%   silent default (neither localsPerArc nor localParams); (5) the derived
%   quantities (nc from commonParams rows, obsTypes from covariance
%   fieldnames, u from x0) drive the run.
%
%   Run: runtests('tests/test_fit')
%
%   Updated 2026-08-31 (CEST, UTC+02:00), v3 data contract (Level L10,
%   WP-72): counts derived, observationTypes/systemId/numLocal/numCommon/
%   loader retired from the surface (unknown fields now), localParams RxC,
%   id with index default, uCount on the batched nonlinear form. With
%   assistance from Claude (Anthropic), model Claude Fable 5.
%   Created 2026-08-31 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5 -- R4 Level L2, WP-56.
    tests = functiontests(localfunctions);
end

%% ===================================================== reference: dispatch
function testFitGlsReproducesAdjustBlock(tc)
% categories: reference, precision
    % Same declaration, two roads: fit(data) and the makeLoader + adjustBlock
    % call it stands for. Same code path, so agreement is solver precision.
    [data, cfg] = declaration();
    opts = lsqvce.Options();                          % Model "gls"
    Rf = lsqvce.fit(data, opts);
    L  = loaderOf(cfg, data.systems(1), data);
    a  = opts.engineArgs("adjustBlock");
    Rd = lsqvce.engine.adjustBlock({L}, data.days, a.nv{:});
    tc.verifyClass(Rf, 'lsqvce.Results');
    % tol: 1e-12 -- identical code path; guards multithreaded-BLAS reduction order only (§10.8 WP-47 r1e)
    tc.verifyEqual(Rf.x, Rd.x, 'RelTol', 1e-12);
    tc.verifyEqual(Rf.sigma0, Rd.sigma0, 'RelTol', 1e-12);
    tc.verifyEqual(Rf.sol.Cxx, Rd.sol.Cxx, 'RelTol', 1e-12);
    tc.verifyLessThan(max(abs(Rf.x - cfg.x_c)), 0.15);            % and it is right
    tc.verifyTrue(isstruct(Rf.systems) && isfield(Rf.systems, 'dataFcn') && ...
        isa(Rf.systems(1).dataFcn, 'function_handle'), ...
        'R.systems carries the normalised system specs incl. the loader closure');
end

function testFitVceReproducesTheEngine(tc)
% categories: reference, precision | tags: diagnostics
    % WP-58 parity: fit(data, Options(Model="vce", ...)) versus the fluent
    % configurator chain it replaces. Exact trace mode (vce default) makes
    % the separate loop deterministic; two outer passes keep the CI cost low.
    % The fixture's 'pos' stream is a SCALAR ARMA stream (buildDayImpl), but
    % adjustVCE's default treats 'pos' as 3 components per epoch and then
    % needs a per-epoch scaling (lsqvce:covModel:noCep -- the r1 CI catch).
    % Both runs therefore declare pos -> 1, which also exercises the
    % ComponentsPerEpoch mapping of Options.engineArgs.
    % WP-47 gate rule 7: MaxOuter = 2 stops BEFORE convergence by design, so
    % the two convergence warnings are provoked deliberately -- muted with
    % the ONE helper (cleanup-restored), never left in the inventory.
    muteWarn(tc, 'lsqvce:adjustVCE:notConverged');
    muteWarn(tc, 'lsqvce:vce:notConverged');
    [data, cfg] = declaration();
    axes = struct('pos', 1);
    opts = lsqvce.Options(Model="vce", MaxOuter=2, VceMaxIter=3, Verbose=false, ...
                          ComponentsPerEpoch=axes);
    Rf = lsqvce.fit(data, opts);
    L  = loaderOf(cfg, data.systems(1), data);
    Rd = lsqvce.engine.adjustVCE({L}, data.days, 'MaxOuter', 2, 'Verbose', false, 'Scale', struct('MaxIter', 3), 'Axes', axes).run();
    % tol: 1e-12 -- identical code path (see the gls test)
    tc.verifyEqual(Rf.x, Rd.x, 'RelTol', 1e-12);
    tc.verifyEqual(Rf.sigma0, Rd.sigma0, 'RelTol', 1e-12);
    tc.verifyEqual(Rf.history.sigma, Rd.history.sigma, 'RelTol', 1e-12, ...
        'the per-pass variance components must match the fluent run');
    tc.verifyEqual(numel(Rf.history.streams), numel(Rd.history.streams));
end

%% ============================================================== contract
function testDeclarationRulesFireBeforeAnyRead(tc)
% categories: contract
    % Every contract rule has ONE lsqvce:fit id, and all of them fire before
    % buildDay is called: the builder here THROWS if reached.
    [data, ~] = declaration();
    bomb = @(varargin) error('test:reached', 'buildDay must not be called during validation');
    data.systems(1).buildDay = bomb;

    d = data;  d.nonsense = 1;
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:unknownField');
    d = data;  d.systems(1).arcz = 2;
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:unknownField');
    % the retired v2 fields are unknown fields now -- loud, named failures
    d = data;  d.common.nc = 6;                    % the v2 sub-struct itself
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:unknownField');
    d = data;  d.systems(1).numLocal = 4;
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:unknownField');
    d = data;  d.systems(1).observationTypes = {'pos', 'kbr'};
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:unknownField');
    d = data;  d.systems(1).systemId = 1;
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:unknownField');
    d = data;  d.systems(1).loader = bomb;
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:unknownField');
    d = rmfield(data, 'systems');
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:noSystems');
    d = data;  d.systems = struct([]);
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:noSystems');
    d = data;  d.systems = rmfield(d.systems, 'covariance');
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:missingField');
    d = data;  d.systems(1).covariance = struct();               % no obs types
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:missingField');
    d = data;  d.systems(1).id = 'CD';
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:badId');
    d = data;  d.systems(1).localsPerArc = 2;                    % TWO ways
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:localCount');
    d = data;  d.systems(1).arcs = 'fromDay';                    % ENGINE-tier only
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:badArcs');
    d = rmfield(data, 'commonParams');
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:noCommon');
    d = data;  d.commonParams = 6;                               % a count, not a table
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:noCommon');
    d = rmfield(data, 'days');
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:noDays');
    d = data;  d.days = [1 2.5];
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:noDays');
    d = data;  d.systems(2) = d.systems(1);                      % same explicit id
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:duplicateSystem');
    d = data;  d.systems(1).id = 2;  d.systems(2) = d.systems(1);
    d.systems(2).id = [];                                        % explicit 2 vs default 2
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:duplicateSystem');
    d = data;  d.systems(1).buildDay = 'notAHandle';
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:builderKind');
    % a makeLoader CONSTRUCTION error is relabelled under lsqvce:fit with the
    % original attached as the cause
    d = data;  d.systems(1).covariance = rmfield(d.systems(1).covariance, 'kbr');
    try
        lsqvce.fit(d);
        tc.verifyFail('a covariance spec missing an obsType must be refused');
    catch e
        tc.verifyEqual(e.identifier, 'lsqvce:fit:noCov');
        tc.verifyEqual(numel(e.cause), 1);
        tc.verifyEqual(e.cause{1}.identifier, 'lsqvce:makeLoader:noCov');
    end
end

function testNotYetPathsAreRefusedExplicitly(tc)
% categories: contract
    % Robust IRLS is the in-memory kinds' path (WP-56b); on the streaming
    % kind fit refuses explicitly -- never a silent fallback. (Model "ellmer"
    % was refused here until WP-56c wired it; its reference test is below.)
    [data, ~] = declaration();
    tc.verifyEqual(idOf(@() lsqvce.fit(data, lsqvce.Options(Robust="irls"))), 'lsqvce:fit:notYet');
end

function testSignatureRejectsNonStructAndNonOptions(tc)
% categories: contract
    % Positional typing at the signature (MATLAB ids): the in-memory
    % fit(A, l, ...) form is WP-56b, so a matrix is refused today.
    ok = false;
    try
        lsqvce.fit(randn(5, 2));
    catch e
        ok = startsWith(e.identifier, 'MATLAB:');
    end
    tc.verifyTrue(ok, 'a matrix first argument must be rejected by the signature');
    [data, ~] = declaration();
    ok = false;
    try
        lsqvce.fit(data, struct('Model', 'gls'));
    catch e
        ok = startsWith(e.identifier, 'MATLAB:');
    end
    tc.verifyTrue(ok, 'opts must be an lsqvce.Options, not a struct');
end

%% ============================================================ robustness
function testDeclarationIsOrderFree(tc)
% categories: robustness
    % Field order in the system struct is irrelevant (same loader, same run).
    [data, ~] = declaration();
    s  = data.systems(1);
    f  = fieldnames(s);
    s2 = orderfields(s, f(end:-1:1));
    d2 = data;  d2.systems = s2;
    R1 = lsqvce.fit(data);
    R2 = lsqvce.fit(d2);
    % tol: 1e-12 -- identical code path (see the reference tests)
    tc.verifyEqual(R2.x, R1.x, 'RelTol', 1e-12);
    tc.verifyEqual(R2.sigma0, R1.sigma0, 'RelTol', 1e-12);
end

function testSystemWithoutLocalsThroughFit(tc)
% categories: contract, robustness
    % G-07 via the front door, v3 form: nl = 0 is the DEFAULT -- neither
    % localsPerArc nor localParams declared, B = zeros(m, 0).
    [data, cfg] = declaration();
    cfg.nl = 0;
    m = cfg.arcsPerDay * cfg.mArc;
    data.systems(1).buildDay = @(d, ot) setfield(buildDayImpl(cfg, d, ot), 'B', zeros(m, 0)); %#ok<SFLD>
    data.systems(1).localParams = [];                            % explicit 'not used'
    R = lsqvce.fit(data);
    tc.verifyClass(R, 'lsqvce.Results');
    tc.verifyEqual(numel(R.x), cfg.nc);
    tc.verifyLessThan(max(abs(R.x - cfg.x_c)), 0.15);
end

function testDerivedQuantitiesDriveTheRun(tc)
% categories: contract
    % v3: the derivations ARE the declaration -- nc = commonParams rows
    % (checked against A on the first read: arcBadNc relabelled is the proof
    % the derived value is enforced), obsTypes = covariance fieldnames in
    % declaration order, id defaults to the index.
    [data, cfg] = declaration();
    R = lsqvce.fit(data);
    tc.verifyEqual(numel(R.x), size(data.commonParams, 1), ...
        'nc must be the commonParams row count');
    tc.verifyEqual(R.systems(1).sysID, 1, 'id must default to the system index');
    d = data;  d.commonParams = num2cell((1:cfg.nc + 1).');      % one row too many
    ok = false;
    try
        lsqvce.fit(d);
    catch e
        ok = contains(e.identifier, 'arcBadNc') || contains(e.identifier, 'BadNc');
    end
    tc.verifyTrue(ok, 'a wrong commonParams row count must be enforced against A');
end

%% ============================================== WP-57: Results.coefficients
function testCoefficientsTable(tc)
% categories: contract, precision
    % R4 WP-57: estimates + formal errors + confidence bounds as ONE table
    % (the fitlm Coefficients/coefCI pair folded together), z from erfcinv.
    [data, cfg] = declaration();
    R = lsqvce.fit(data, lsqvce.Options(Verbose=false));
    T = R.coefficients();
    tc.verifyClass(T, 'table');
    tc.verifyEqual(size(T), [cfg.nc 5]);
    tc.verifyEqual(T.Properties.VariableNames, {'index', 'estimate', 'se', 'lower', 'upper'});
    % tol: exact -- the same arithmetic on the same numbers
    tc.verifyEqual(T.estimate, R.x);
    tc.verifyEqual(T.se, sqrt(max(diag(R.sol.Cxx), 0)));
    % tol: 1e-12 -- z(0.05) = 1.959963984540054 in closed form (16 digits)
    tc.verifyEqual((T.upper - T.estimate) ./ T.se, 1.959963984540054 * ones(cfg.nc, 1), 'RelTol', 1e-12);
    T1 = R.coefficients(Alpha = 1 - erf(1/sqrt(2)));     % P(|Z| > 1): one-sigma bounds
    % tol: 1e-12 -- erfc / erfcinv round trip, z = 1 exactly in exact arithmetic
    tc.verifyEqual(T1.upper - T1.estimate, T1.se, 'RelTol', 1e-12);
    tc.verifyTrue(all(T.lower < T.estimate & T.estimate < T.upper));
    ok = false;
    try
        R.coefficients(Alpha = 1.5);
    catch e
        ok = startsWith(e.identifier, 'MATLAB:');
    end
    tc.verifyTrue(ok, 'Alpha outside (0,1) must be rejected by the signature');
    sol = rmfield(R.sol, 'Cxx');                          % a Results without a covariance
    R2 = lsqvce.Results(sol, R.diag, R.acc, R.source, R.batchFcn, R.nBatch, R.systems);
    tc.verifyEqual(idOf(@() R2.coefficients()), 'lsqvce:Results:noCovariance');
end

%% ================================================================ helpers
%% ================================================= WP-56c: Model "ellmer"
function testEllmerThroughFitReproducesDirect(tc)
% categories: reference, precision | tags: diagnostics, provenance
    % R4 WP-56c: fit(Model="ellmer") is lsqvce.engine.adjustVCEjoint over an accFcn
    % that fit builds from the loaders (BlockAccumulator + jointModelOpt for
    % every 'satKey|obsType' stream), bridged into a Results. The direct route
    % below hand-builds the same accFcn; same Seed -> same probes -> identical.
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');   % JointMaxIter 2 stops early by design
    [data, cfg] = declaration();
    % jointModelOpt maps the model through an 'acf' spec that is ZERO beyond
    % C's last lag; with NumLags < arc length the Toeplitz matrix is a
    % truncated correlation, indefinite in general, every arc is skipped as
    % nonPD and the driver reports noArcs (the L3 r1/r2 catches, together
    % with G-08's lag-unit defect). NumLags = arc length avoids truncation.
    nLag = cfg.mArc;
    opts = lsqvce.Options(Model="ellmer", NumLags=nLag, JointMaxIter=2, Probes=20, ...
                          Seed=7, Verbose=false);
    R = lsqvce.fit(data, opts);

    L = loaderOf(cfg, data.systems(1), data);
    streams = {'x1|pos', 'x1|kbr'};                       % makeValidName('1') = 'x1'
    accFcn = @(model) directJointAcc(cfg.nc, L, data.days, ...
                 lsqvce.engine.jointModelOpt(model, streams, 'PsdMode', 'toeplitz'));
    out = lsqvce.engine.adjustVCEjoint(accFcn, 'NumLags', nLag, 'MaxIter', 2, 'Probes', 20, ...
                                'Seed', 7, 'Verbose', false);
    % tol: 1e-12 -- identical code path; the seed fixes the Hutchinson probes
    tc.verifyEqual(R.x, out.sol.x(:), 'RelTol', 1e-12);
    tc.verifyEqual(R.noise.C, out.C, 'RelTol', 1e-12);
    tc.verifyEqual(R.noise.sigma2, out.sigma2, 'RelTol', 1e-12);
    tc.verifyEqual(R.history.kind, 'ellmer');
    tc.verifyEqual(R.history.nIter, out.nIter);
    tc.verifyEqual(R.converged, out.converged);
    tc.verifyClass(R.acc, 'lsqvce.engine.BlockAccumulator');
    tc.verifyEqual(R.meta.seeds.jointSeed, 7, 'the probe seed is recorded in the provenance');
    tc.verifyEqual(R.meta.options.Model, "ellmer");
    % adjustVCE's history readers refuse cleanly; the verdict works
    tc.verifyEqual(idOf(@() R.summary()), 'lsqvce:Results:noHistory');
    tc.verifyEqual(idOf(@() R.report('Verbose', false)), '');
end

function acc = directJointAcc(nc, L, days, uo)
    acc = lsqvce.engine.BlockAccumulator(nc, 'KeepLocal', true, 'Verbose', false, 'UserOpt', uo);
    acc.addSystem(L.sysID, L.nl, L.closure, days);
end

%% ================================================== WP-62: parallel policy
function testDistributeDaysThroughFitEqualsFused(tc)
% categories: precision, robustness
    % R4 WP-62: Options.DistributeDays / ParallelThreshold reach the
    % accumulator through fit -> adjustBlock. Below the pool threshold (and
    % without PCT) mapreduce runs the DISTRIBUTED code path serially -- per-day
    % raw partials merged in day order -- so the numerical claim is the same
    % one test_parallel pins for the engine: equal to the fused loop up to
    % the summation-order gap, no pool needed on CI.
    [data, ~] = declaration();
    R1 = lsqvce.fit(data, lsqvce.Options(Verbose=false));
    R2 = lsqvce.fit(data, lsqvce.Options(Verbose=false, Parallel="parfor", ...
                                         DistributeDays=true, ParallelThreshold=2));
    % tol: 1e-11 -- serial-vs-distributed reduction order (test_parallel
    % testDistributeDaysEqualsFusedSerial's measured bound); bitwise NOT promised
    tc.verifyEqual(R2.x,      R1.x,      'RelTol', 1e-11);
    tc.verifyEqual(R2.sigma0, R1.sigma0, 'RelTol', 1e-11);
    % adjustBlock records its RESOLVED option struct (lowerCamel, WP-44)
    tc.verifyEqual(R2.meta.options.distributeDays, true, 'the policy is recorded in the provenance');
    tc.verifyEqual(R2.meta.options.parallelThreshold, 2);
    tc.verifyEqual(R1.meta.options.distributeDays, false);
    % the same knobs are accepted by the VCE engine's constructor (signature)
    L = loaderOf([], data.systems(1), data);
    P = lsqvce.engine.adjustVCE({L}, data.days, 'DistributeDays', true, 'ParallelThreshold', 2, 'Verbose', false);
    tc.verifyEqual(P.blockOpt.DistributeDays, true);
    tc.verifyEqual(P.blockOpt.ParallelThreshold, 2);
end

%% ================================================== WP-67: the user-path demo
function testDemoFitRuns(tc)
% categories: contract, robustness | tags: diagnostics
    % R4 Level L7: examples/demo_fit.m is the 2.0 user path on all three data
    % kinds. It is a template users copy, so it must keep running: every
    % Results it returns is checked here (the demo mutes the two convergence
    % warnings its two-pass "vce" run provokes and restores them on exit).
    here = fileparts(mfilename('fullpath'));
    exDir = fullfile(here, '..', '..', 'examples');    % dev/tests -> root (L8)
    addpath(exDir);
    tc.addTeardown(@() rmpath(exDir));
    out = demo_fit(Verbose=false);
    for f = {'matrix', 'robust', 'batch', 'nonlinear', 'gls', 'vce'}
        R = out.(f{1});
        tc.verifyClass(R, 'lsqvce.Results', f{1});
        tc.verifyTrue(isfinite(R.sigma0) && R.sigma0 > 0, [f{1} ': sigma0 finite and positive']);
    end
    % tol: 3 formal sigma per unknown -- 30 observations with noise ~ 1/sqrt(p),
    % p in [0.5, 1.5], give formal errors of ~0.3-0.4; a literal bound would
    % pin the fixture's noise, not the estimator (the L7 r1 CI catch)
    z = (out.matrix.x - out.truth.x) ./ sqrt(diag(out.matrix.sol.Cxx));
    tc.verifyLessThan(max(abs(z)), 3, 'the in-memory truth is recovered within 3 formal sigma');
    % tol: 1e-12 -- the same rows accumulated in two batches vs one; summation order only
    tc.verifyEqual(out.batch.x, out.matrix.x, 'RelTol', 1e-12, 'two batches == one batch (same rows)');
    tc.verifyLessThan(out.robust.source.factor(7), 0.5, 'the planted gross error is downweighted');
    % tol: 5e-3 -- 40 points at 1e-3 noise (the WP-56b pin)
    tc.verifyEqual(out.nonlinear.x, out.truth.ab, 'AbsTol', 5e-3);
    tc.verifyEqual(size(out.vce.history.sigma, 1), 2, 'two outer passes recorded');
    tc.verifyClass(out.gls.acc, 'lsqvce.engine.BlockAccumulator');
    tc.verifyEqual(numel(out.vce.systems), 3, 'fit built the three systems from the data struct');
    % tol: 5 formal sigma per coefficient -- test_physOracle's own bound
    z = (out.gls.x - out.truth.xGravity) ./ (sqrt(diag(out.gls.sol.Cxx)) * out.gls.sigma0);
    tc.verifyLessThan(max(abs(z)), 5, 'the physical oracle is recovered within 5 formal sigma');
end

%% ================================================ WP-56b: kinds B and C
function testMatrixKindReproducesAdjust(tc)
% categories: reference, contract
    % R4 WP-56b: fit(data) with data.A / .l / .p is exactly lsqvce.engine.adjust on
    % one batch -- same driver, same options -- so x, sigma0, Cxx and the
    % per-observation diagnostics are byte-identical.
    [A, l, p] = matrixFixture();
    data = struct('A', A, 'l', l, 'p', p);
    R = lsqvce.fit(data, lsqvce.Options(Verbose=false));
    [sol, D, V] = lsqvce.engine.adjust(size(A, 2), @(k) deal(A, l, p), 1, 'Verbose', false);
    tc.verifyClass(R, 'lsqvce.Results');
    % tol: exact -- identical code path and inputs
    tc.verifyEqual(R.x, sol.x(:));
    tc.verifyEqual(R.diag.v, D.v);
    % the bridge restores the STREAMING Results contract: sigma0 a posteriori,
    % Cxx = s0sq * Qxx (the batch driver's sol carries the a-priori Sigma0 and
    % no Cxx -- the r1 CI catch of this level)
    % tol: exact -- same arithmetic as lsqvce.engine.diagnose, same numbers
    tc.verifyEqual(R.sigma0, sqrt(sol.s0sq));
    tc.verifyEqual(R.sol.sigma0Prior, sol.sigma0);
    tc.verifyEqual(R.sol.Cxx, sol.s0sq * D.Qxx);
    tc.verifyEqual(R.source.verdict, V.verdict, 'the verdict rides on R.source');
    tc.verifyEmpty(R.acc);
    tc.verifyEqual(R.meta.options.kind, 'matrix');
    % the user-facing conveniences work without an accumulator ...
    T = R.coefficients();
    tc.verifyEqual(T.estimate, R.x);
    tc.verifyEqual(idOf(@() R.report('Verbose', false)), '', 'report() must run on a kind-B Results');
    % ... and the streaming-only ones refuse with ONE clear id
    tc.verifyEqual(idOf(@() R.residuals()), 'lsqvce:Results:noAccumulator');
    tc.verifyEqual(idOf(@() R.contribution()), 'lsqvce:Results:noAccumulator');
end

function testMatrixKindRobustReproducesRobust(tc)
% categories: reference, contract
    % Robust="irls" on kind B is lsqvce.engine.robust with the Options' Influence /
    % Tuning / RobustMaxIter / RobustTol; one planted outlier.
    [A, l, p] = matrixFixture();
    l(7) = l(7) + 25;                                   % gross error
    data = struct('A', A, 'l', l, 'p', p);
    opts = lsqvce.Options(Robust="irls", Influence="huber", RobustMaxIter=20, Verbose=false);
    R = lsqvce.fit(data, opts);
    [sol, D, ROB] = lsqvce.engine.robust(size(A, 2), @(k) deal(A, l, p), 1, ...
        'Influence', 'huber', 'MaxIter', 20, 'Verbose', false);
    % tol: exact -- identical code path and inputs
    tc.verifyEqual(R.x, sol.x(:));
    tc.verifyEqual(R.sigma0, sqrt(sol.s0sq));
    tc.verifyEqual(R.diag.v, D.v);
    tc.verifyEqual(fieldnames(R.source), fieldnames(ROB), 'the IRLS record rides on R.source');
end

function testBatchKindNonlinearReproducesSolveNonlinear(tc)
% categories: reference, contract
    % Kind C with x0: [A, w, p] = batchFcn(k, x) is lsqvce.engine.solveNonlinear.
    % Model: l = exp(a t) + b, unknowns [a; b], from a known truth.
    rs = RandStream('mt19937ar', 'Seed', 5);
    t  = linspace(0, 1, 40).';
    xt = [0.7; -0.2];
    lobs = exp(xt(1)*t) + xt(2) + 1e-3*randn(rs, 40, 1);
    idx = {1:20, 21:40};
    bf = @(k, x) deal([t(idx{k}).*exp(x(1)*t(idx{k})), ones(20, 1)], ...
                      lobs(idx{k}) - (exp(x(1)*t(idx{k})) + x(2)), []);
    data = struct('batchFcn', bf, 'numBatches', 2, 'x0', [0.5; 0]);   % u = numel(x0), v3
    opts = lsqvce.Options(NonlinearMaxIter=15, TolX=1e-9, Verbose=false);
    R = lsqvce.fit(data, opts);
    [sol, D, ~, IT] = lsqvce.engine.solveNonlinear(2, bf, 2, [0.5; 0], 'MaxIter', 15, 'TolX', 1e-9, 'Verbose', false);
    % tol: exact -- identical code path and inputs
    tc.verifyEqual(R.x, sol.x(:));
    tc.verifyEqual(R.diag.v, D.v);
    tc.verifyEqual(fieldnames(R.source), fieldnames(IT), 'the iteration record rides on R.source');
    % tol: 5e-3 -- 40 points at 1e-3 noise recover the truth to a few 1e-3
    tc.verifyEqual(R.x, xt, 'AbsTol', 5e-3);
    tc.verifyEqual(R.meta.options.kind, 'batch');
end

function testBatchKindLinearReproducesAdjust(tc)
% categories: reference
    [A, l, p] = matrixFixture();
    idx = {1:15, 16:30};
    bf = @(k) deal(A(idx{k}, :), l(idx{k}), p(idx{k}));
    data = struct('batchFcn', bf, 'numUnknowns', size(A, 2), 'numBatches', 2);
    R = lsqvce.fit(data, lsqvce.Options(Verbose=false));
    [sol, ~, ~] = lsqvce.engine.adjust(size(A, 2), bf, 2, 'Verbose', false);
    % tol: exact -- identical code path and inputs
    tc.verifyEqual(R.x, sol.x(:));
    tc.verifyEqual(R.nBatch, 2);
end

function testKindsBAndCRulesFireBeforeAnyRead(tc)
% categories: contract
    [A, l, p] = matrixFixture();
    tc.verifyEqual(idOf(@() lsqvce.fit(struct('A', A, 'l', l(1:end-1)))), 'lsqvce:fit:badMatrix');
    tc.verifyEqual(idOf(@() lsqvce.fit(struct('A', A, 'l', l, 'p', p(1:3)))), 'lsqvce:fit:badMatrix');
    tc.verifyEqual(idOf(@() lsqvce.fit(struct('A', A))), 'lsqvce:fit:missingField');
    tc.verifyEqual(idOf(@() lsqvce.fit(struct('A', A, 'l', l, 'weights', p))), 'lsqvce:fit:unknownField');
    tc.verifyEqual(idOf(@() lsqvce.fit(struct('A', A, 'l', l), lsqvce.Options(Model="vce"))), 'lsqvce:fit:modelKind');
    tc.verifyEqual(idOf(@() lsqvce.fit(struct('batchFcn', @(k) deal(A, l, p)))), 'lsqvce:fit:missingField');
    tc.verifyEqual(idOf(@() lsqvce.fit(struct('batchFcn', @(k) deal(A, l, p), ...
        'numBatches', 1))), 'lsqvce:fit:missingField');          % linear form needs u
    tc.verifyEqual(idOf(@() lsqvce.fit(struct('batchFcn', @(k, x) deal(A, l, p), ...
        'numBatches', 1, 'x0', zeros(3, 1), 'numUnknowns', 3))), 'lsqvce:fit:uCount');
    tc.verifyEqual(idOf(@() lsqvce.fit(struct('batchFcn', @(k) deal(A, l, p), 'numUnknowns', 3, ...
        'numBatches', 1), lsqvce.Options(Model="pooled"))), 'lsqvce:fit:modelKind');
end

function [A, l, p] = matrixFixture()
% A small well-conditioned in-memory problem with a known truth.
    rs = RandStream('mt19937ar', 'Seed', 17);
    A  = [ones(30, 1), randn(rs, 30, 2)];
    xt = [1; -2; 0.5];
    p  = 0.5 + rand(rs, 30, 1);
    l  = A*xt + randn(rs, 30, 1) ./ sqrt(p);
end

function id = idOf(f)
% Identifier of the error f raises ('' if it does not raise).
    id = '';
    try
        f();
    catch e
        id = e.identifier;
    end
end

function [data, cfg] = declaration()
% The canonical single-system declaration over the makeLoader-path fixture
% (v3 contract: the tables declare, the counts are derived).
    cfg = buildCfg();
    data = struct();
    data.commonParams = num2cell((1:cfg.nc).');       % nc = rows (RxC, C = 1)
    data.days = 1:3;
    s = struct();
    s.buildDay = @(d, ot) buildDayImpl(cfg, d, ot);
    s.covariance = cfg.cov;                           % fieldnames = pos, kbr
    s.localParams = num2cell((1:cfg.nl).');           % nl = rows
    s.arcs = cfg.arcsPerDay;
    data.systems = s;                                 % id defaults to 1
end

function L = loaderOf(~, s, data)
% What fit builds for one system -- written out, so the tests compare two
% independent spellings of the same declaration (v3: obsTypes from the
% covariance fieldnames, the tables carry the counts, id -> SystemId).
    L = lsqvce.engine.makeLoader(s.buildDay, ...
            'ObservationTypes', fieldnames(s.covariance).', ...
            'Covariance', s.covariance, 'SystemId', 1, ...
            'LocalParams', s.localParams, ...
            'Arcs', s.arcs, 'Days', data.days, 'CommonParams', data.commonParams);
end

function cfg = buildCfg()
% cfg for the makeLoader buildDay path (mirrors tests/test_adjustBlock.m).
    cfg = struct('nc',6, 'nl',4, 'arcsPerDay',2, 'mArc',300, 'dt',5, 'seed',7);
    rs = RandStream('mt19937ar','Seed',11);
    cfg.x_c = randn(rs, 6, 1);
    cfg.y1  = randn(rs, 4, 1);
    cfg.cov.pos = struct('type','arma','ar',[1 -1.6 0.9 -0.25],'ma',1,'sigma',0.7);
    cfg.cov.kbr = struct('type','diag','sigma',0.2);
end

function D = buildDayImpl(cfg, day, obsType)
    M  = cfg.arcsPerDay * cfg.mArc;
    sd = mod(cfg.seed + 7*day + 31*sum(double(obsType)), 2^31-1);
    rs = RandStream('mt19937ar', 'Seed', sd);
    D.t = (0:M-1).' * cfg.dt;
    D.A = randn(rs, M, cfg.nc);
    D.B = randn(rs, M, cfg.nl);
    if strcmp(obsType, 'kbr'), D.A(:,1:3) = 0.1 * D.A(:,1:3); end
    e   = synthE(rs, cfg.cov.(obsType), M);
    D.l = D.A*cfg.x_c + D.B*cfg.y1(1:cfg.nl) + e;
    D.arcIdx = {1:cfg.mArc, cfg.mArc + (1:cfg.mArc)};
end

function e = synthE(rs, cm, M)
    nBurn = 2000;
    switch lower(cm.type)
        case 'diag'
            e = cm.sigma * randn(rs, M, 1);
        case 'arma'
            w = cm.sigma * randn(rs, M + nBurn, 1);
            e = filter(cm.ma, cm.ar, w);
            e = e(nBurn+1:end);
    end
end
