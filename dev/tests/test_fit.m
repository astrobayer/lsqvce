function tests = test_fit
%TEST_FIT  Tests for lsqvce.fit, the one user verb (R4 WP-56, decision R4-D4).
%
%     lsqvce.fit   data struct (the Manual, chapter 3) + lsqvce.Options -> Results;
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
%   Updated 2026-09-09 (CEST, UTC+02:00), v6 (R6 WP-102): THE CELLS ARE THE
%   ARCS -- every fixture day is cut by the test-private arcCells into the
%   1 x nArc cell contract; NEW points: the schedule is optional on the fit
%   surface (discovered == declared, bitwise), an empty cell is the gap of a
%   scheduled arc (== the resolver that omits it), 'day' rows without a
%   schedule live on every calendar day and an empty day is refused
%   (emptyDay), 'arc' rows need the schedule (localPlan, relabelled). With
%   assistance from Claude (Anthropic), model Claude Fable 5.1.
%   Updated 2026-08-31 (CEST, UTC+02:00), v4 SCOPED tables (Level L11,
%   WP-76): fixtures on the scope column; badScope cases (direct + relabel);
%   localsPerArc -> unknownField; NEW reference parity scoped == hand-numbered
%   to 1e-12 and the day-without-arcs gap rules (common + local side). With
%   assistance from Claude (Anthropic), model Claude Fable 5.
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
% categories: reference, contract
    % R5 M4 4b (decision R5-D4): Model "vce" is a POLICY of the one joint
    % estimator -- one law per stream (Pool "stream"), no arc factors -- and
    % Model "pooled" is Pool "all" without arc factors. Pinned as EXACT
    % equalities against the explicit spelling with Model "ellmer" (same
    % seed => same probes => same arithmetic), so the mapping can never drift.
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');
    [data, ~] = declaration();
    % the same budget on both spellings: vce reads MaxOuter, ellmer JointMaxIter
    base = {'MaxPasses', 2, 'Noise', struct('Grid', 8), 'Seed', 7, 'Verbose', false};
    Rv = lsqvce.fit(data, lsqvce.Options(base{:}, Estimate="noise", Pool="stream", ArcFactors=false));
    dataE = data;
    ots = fieldnames(dataE.systems(1).covariance);
    for k = 1:numel(ots), dataE.systems(1).covariance.(ots{k}).arcFactors = false; end
    Re = lsqvce.fit(dataE, lsqvce.Options(base{:}, Estimate="noise", Pool="stream"));
    % tol: exact -- one code path, one seed
    tc.verifyEqual(Rv.x, Re.x);
    tc.verifyEqual(Rv.noise.C, Re.noise.C);
    % vce: no PER-ARC factors -- every arc of a stream carries the stream's
    % ONE level (the Helmert factor), so within a stream the factors are equal
    for s = unique(cellfun(@(k) regexprep(k, '_d.*$', ''), Rv.noise.arcKeys, 'UniformOutput', false)).'
        f = Rv.noise.sigma2(startsWith(Rv.noise.arcKeys, s{1}));
        % tol: exact -- one number written to every arc of the stream
        tc.verifyEqual(f, repmat(f(1), size(f)), sprintf('%s: one level', s{1}));
    end
    tc.verifyEqual(Rv.history.kind, 'ellmer');
    tc.verifyEqual(Rv.meta.options.Pool, "stream", 'the provenance keeps the resolved policy');
    tc.verifyFalse(Rv.meta.options.ArcFactors);
    % the per-stream summary every model shares
    tc.verifyEqual(numel(Rv.noise.stream), numel(ots));
    tc.verifyTrue(all(cellfun(@(m) isfield(m, 'sigma') && m.sigma > 0, Rv.noise.model)));
    % "pooled" = Pool "all" without arc factors
    Rp = lsqvce.fit(data, lsqvce.Options(base{:}, Estimate="noise", ArcFactors=false));
    Rpe = lsqvce.fit(dataE, lsqvce.Options(base{:}, Estimate="noise", Pool="all"));
    % tol: exact -- as above
    tc.verifyEqual(Rp.x, Rpe.x);
    tc.verifyEqual(keys(Rp.noise.Cgroups), {'all'});
end
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
    d = data;  d.systems(1).localsPerArc = 2;                    % retired in v4
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:unknownField');
    d = data;  d.commonParams{2, 1} = 'arc';                     % not a common scope
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:badScope');
    d = data;  d.commonParams{1, 1} = 42;                        % not a token at all
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:badScope');
    d = data;  d.systems(1).arcs = 'fromDay';                    % v6: the string form is gone
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:badArcs');
    d = rmfield(data, 'commonParams');
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:noCommon');
    d = data;  d.commonParams = 6;                               % a count, not a table
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:noCommon');
    d = rmfield(data, 'days');
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:noDays');
    d = data;  d.days = [1 2.5];
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:noDays');
    d = data;  d.systems(1).id = 7;  d.systems(2) = d.systems(1);   % same explicit id
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:duplicateSystem');
    d = data;  d.systems(1).id = 2;  d.systems(2) = d.systems(1);
    d.systems(2).id = [];                                        % explicit 2 vs default 2
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:duplicateSystem');
    d = data;  d.systems(1).buildDay = 'notAHandle';
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:builderKind');
    % a makeLoader CONSTRUCTION error is relabelled under lsqvce:fit with the
    % original attached as the cause. (noCov is UNREACHABLE through fit since
    % v3 -- the type list IS fieldnames(covariance) -- so the mechanism is
    % proven on paramType: a localParams cell that is not a 2-D matrix. A
    % non-cell would be stopped by makeLoader's SIGNATURE with a MATLAB
    % validation id, before the body runs -- the b1-r3 split.)
    d = data;  d.systems(1).localParams = cell(2, 2, 2);
    try
        lsqvce.fit(d);
        tc.verifyFail('a non-matrix localParams cell must be refused by makeLoader via fit');
    catch e
        tc.verifyEqual(e.identifier, 'lsqvce:fit:paramType');
        tc.verifyEqual(numel(e.cause), 1);
        tc.verifyEqual(e.cause{1}.identifier, 'lsqvce:makeLoader:paramType');
    end
    % v4: a bad LOCAL scope token is makeLoader's badScope, relabelled
    d = data;  d.systems(1).localParams{1, 1} = 'daily';
    try
        lsqvce.fit(d);
        tc.verifyFail('a bad localParams scope token must be refused');
    catch e
        tc.verifyEqual(e.identifier, 'lsqvce:fit:badScope');
        tc.verifyEqual(e.cause{1}.identifier, 'lsqvce:makeLoader:badScope');
    end
end

function testNotYetPathsAreRefusedExplicitly(tc)
% categories: contract
    % R5 M2 2b: Robust IRLS now RUNS on the streaming kind (a policy on the
    % block path); what is still refused explicitly is its composition with a
    % VCE model (the one-loop feature of M4) -- never a silent fallback.
    [data, ~] = declaration();
    tc.verifyEqual(idOf(@() lsqvce.fit(data, lsqvce.Options(Robust="huber", Estimate="noise", Pool="stream", ArcFactors=false))), 'lsqvce:fit:notYet');
end

function testRobustRunsOnTheStreamingKind(tc)
% categories: contract, robustness
    % R5 M2 2b: the streaming declaration with Robust="huber" converges, keeps
    % every observation's factor, and agrees with the GLS solve to the
    % influence's tolerance when there is no outlier (all factors ~ 1)
    [data, ~] = declaration();
    R  = lsqvce.fit(data, lsqvce.Options(Robust="huber", Verbose=false));
    R0 = lsqvce.fit(data, lsqvce.Options(Verbose=false));
    tc.verifyEqual(R.stopReason, "converged");
    fac = R.history(end).factor;
    tc.verifyGreaterThan(mean(fac), 0.9, 'clean data: almost every factor is 1');
    % tol: 3e-3 absolute -- Huber leaves ~5 % of Gaussian rows below 1, so the
    % robust and the GLS estimates differ by a fraction of those rows' share
    % (parameters are O(1); the smallest, O(0.03), moved by 1.3e-3 on r1)
    tc.verifyEqual(R.x, R0.x, 'AbsTol', 3e-3);
    tc.verifyEqual(R.meta.options.kind, 'stream');
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
    data.systems(1).buildDay = @(d, ot) setfield(buildDayImpl(cfg, d, ot), 'B', ...   %#ok<SFLD>
                                        {zeros(cfg.mArc, 0), zeros(cfg.mArc, 0)});   % v6: one cell per arc
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
    d = data;  d.commonParams = [repmat({'common'}, cfg.nc + 1, 1), ...
                                 num2cell((1:cfg.nc + 1).')];    % one row too many
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
    % R4 WP-56c: fit(Estimate="noise") is lsqvce.engine.adjustVCEjoint over an accFcn
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
    opts = lsqvce.Options(Estimate="noise", Noise=struct('Grid', nLag), MaxPasses=2, Probes=20, ...
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
    tc.verifyEqual(R.meta.options.Estimate, "noise");
    % R5 M4 4a: the joint sweep collects whiteness per stream per pass, so
    % the history readers WORK on an ellmer Results (before 4a they refused
    % with lsqvce:Results:noHistory); the verdict works too
    tc.verifyEqual(idOf(@() R.summary()), '');
    tc.verifyEqual(idOf(@() R.convergence()), '');
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
    R2 = lsqvce.fit(data, lsqvce.Options(Verbose=false, Parallel=struct('Mode', "parfor"), ...
                                         Parallel=struct('DistributeDays', true, 'Threshold', 2)));
    % tol: 1e-11 -- serial-vs-distributed reduction order (test_parallel
    % testDistributeDaysEqualsFusedSerial's measured bound); bitwise NOT promised
    tc.verifyEqual(R2.x,      R1.x,      'RelTol', 1e-11);
    tc.verifyEqual(R2.sigma0, R1.sigma0, 'RelTol', 1e-11);
    % adjustBlock records its RESOLVED option struct (lowerCamel, WP-44)
    % R5 M5 5b: every kind records the USER's resolved Options as provenance
    tc.verifyEqual(R2.meta.options.Parallel.DistributeDays, true, 'the policy is recorded in the provenance');
    tc.verifyEqual(R2.meta.options.Parallel.Threshold, 2);
    tc.verifyEqual(R1.meta.options.Parallel.DistributeDays, false);
    % the same knobs are accepted by the VCE engine's constructor (signature)
end

%% ================================================== WP-67: the user-path demo
function testDemoFitRuns(tc)
% categories: contract, robustness | tags: diagnostics
    % R5 M4 4b: the demo's two-pass vce run is the joint estimator now and
    % says honestly that two passes did not settle (gate rule 7: muted here)
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');
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
    tc.verifyLessThan(out.robust.history(end).factor(7), 0.5, 'the planted gross error is downweighted');
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
    % R5 M2 (WP-81 2a): fit(data) with data.A / .l / .p rides the BLOCK path
    % as one local-free system (finding H-02). The flat driver
    % lsqvce.engine.adjust is the REFERENCE until round 2b retires it: x,
    % sigma0 and Cxx agree to rounding; the per-observation residual of the
    % block path is the WHITENED one, sqrt(p).*v (the rows were pre-whitened
    % so the declared law is diag(1)) -- pinned as that exact relation.
    [A, l, p] = matrixFixture();
    data = struct('A', A, 'l', l, 'p', p);
    R = lsqvce.fit(data, lsqvce.Options(Verbose=false));
    % R5 M5 5b': the reference is the DENSE closed form (O1), not a second
    % driver -- the flat one is retired
    P = diag(p);  N = A.'*P*A;  xd = N \ (A.'*P*l);  vd = A*xd - l;
    dof = numel(l) - size(A, 2);  s0sq = (vd.'*P*vd) / dof;  Qxx = inv(N);
    tc.verifyClass(R, 'lsqvce.Results');
    % tol: 1e-12 relative -- block QR vs dense normal equations on one
    % well-conditioned 30 x 3 problem, cond ~ 1e1
    tc.verifyEqual(R.x, xd, 'RelTol', 1e-12);
    tc.verifyEqual(R.sigma0, sqrt(s0sq), 'RelTol', 1e-12);
    % tol: 1e-10 relative / 1e-14 absolute -- Cxx and the residual carry the
    % squared conditioning; 100x the x bound
    tc.verifyEqual(R.sol.Cxx, s0sq * Qxx, 'RelTol', 1e-10, 'AbsTol', 1e-14);
    % tol: as above (the whitened residual sqrt(p).*v)
    tc.verifyEqual(R.diag.v, sqrt(p(:)) .* vd, 'RelTol', 1e-10, 'AbsTol', 1e-14, ...
        'the block path reports the whitened residual sqrt(p).*v');
    tc.verifyNotEmpty(R.acc, 'kind B now carries the accumulator');
    tc.verifyEqual(R.meta.options.kind, 'matrix');
    tc.verifyEqual(R.meta.options.Estimate, "none", 'the resolved USER options are the provenance');
    % the user-facing conveniences ...
    T = R.coefficients();
    tc.verifyEqual(T.estimate, R.x);
    tc.verifyEqual(idOf(@() R.report('Verbose', false)), '', 'report() must run on a kind-B Results');
    % ... and the streaming-only ones now WORK on kind B (the point of 2a)
    res = R.residuals();
    tc.verifyTrue(isstruct(res) || iscell(res) || istable(res), 'residuals() runs on kind B');
    tc.verifyEqual(idOf(@() R.contribution()), '', 'contribution() runs on kind B');
end

function testMatrixKindRobustReproducesRobust(tc)
% categories: reference, contract
    % R5 M5 5b': the IRLS fixed point is checked against its own DEFINITION
    % (the flat reference driver is retired): at convergence the factors are
    % the influence function of the standardised whitened residuals, and the
    % estimate is the weighted GLS solution with those factors -- a closed
    % form on the final weights. The planted error is down-weighted.
    [A, l, p] = matrixFixture();
    l(7) = l(7) + 25;                                   % gross error
    data = struct('A', A, 'l', l, 'p', p);
    opts = lsqvce.Options(Robust="huber", RobustOpts=struct('MaxIter', 40, 'Tol', 1e-12), Verbose=false);
    R = lsqvce.fit(data, opts);
    fac = R.history(end).factor;
    tc.verifyEqual(numel(fac), numel(l));
    tc.verifyLessThan(fac(7), 0.2, 'the planted 25-sigma error is down-weighted');
    tc.verifyEqual(R.stopReason, "converged");
    % the fixed point: x = argmin sum fac_i p_i v_i^2 with fac from the LAST pass
    Pw = diag(p(:) .* fac(:));  xd = (A.'*Pw*A) \ (A.'*Pw*l);
    % tol: 1e-6 relative -- the last factors were computed from the PREVIOUS
    % estimate; at Tol 1e-12 on x the two differ by the influence's Lipschitz
    % constant times the step, well below 1e-6
    tc.verifyEqual(R.x, xd, 'RelTol', 1e-6);
    % and the factors are Huber's on the standardised whitened residuals
    vw = sqrt(p(:)) .* (A*R.x - l);
    facDef = lsqvce.kernel.robustWeights(vw, 'huber', []);
    % tol: 1e-6 absolute -- as above
    tc.verifyEqual(fac, facDef, 'AbsTol', 1e-6);
    tc.verifyEqual(R.meta.options.kind, 'matrix');
    tc.verifyNotEmpty(R.acc, 'kind B robust carries the accumulator too');
end
function testBatchKindNonlinearReproducesSolveNonlinear(tc)
% categories: reference, contract
    % R5 M2 2c': kind C with x0 is Gauss-Newton on the BLOCK path
    % (lsqvce.engine.adjustNonlinear); the flat driver solveNonlinear stays
    % the reference until the flat family retires. Same step control, same
    % fixed point: x agrees to rounding; the truth is recovered.
    % Model: l = exp(a t) + b, unknowns [a; b], from a known truth.
    rs = RandStream('mt19937ar', 'Seed', 5);
    t  = linspace(0, 1, 40).';
    xt = [0.7; -0.2];
    lobs = exp(xt(1)*t) + xt(2) + 1e-3*randn(rs, 40, 1);
    idx = {1:20, 21:40};
    bf = @(k, x) deal([t(idx{k}).*exp(x(1)*t(idx{k})), ones(20, 1)], ...
                      lobs(idx{k}) - (exp(x(1)*t(idx{k})) + x(2)), []);
    data = struct('batchFcn', bf, 'numBatches', 2, 'x0', [0.5; 0]);   % u = numel(x0), v3
    opts = lsqvce.Options(MaxPasses=15, Tolerance=1e-9, Verbose=false);
    R = lsqvce.fit(data, opts);
    % R5 M5 5b': the reference is the Gauss-Newton FIXED POINT itself (the
    % flat driver is retired): at x the normal equations of the linearised
    % problem have a zero right-hand side, J' (l - f(x)) = 0
    J = [t.*exp(R.x(1)*t), ones(40, 1)];  w = lobs - (exp(R.x(1)*t) + R.x(2));
    % tol: 1e-8 absolute -- the step tolerance was 1e-9 on x; the gradient
    % scales with J (O(1)) times the residual (1e-3)
    tc.verifyEqual(J.'*w, [0; 0], 'AbsTol', 1e-8, 'the Gauss-Newton stationarity condition holds');
    tc.verifyEqual(R.stopReason, "converged");
    tc.verifyEqual(R.noise.kind, 'nonlinear');
    tc.verifyNotEmpty(R.history, 'the pass record is kept');
    % tol: 5e-3 -- 40 points at 1e-3 noise recover the truth to a few 1e-3
    tc.verifyEqual(R.x, xt, 'AbsTol', 5e-3);
    tc.verifyEqual(R.meta.options.kind, 'batch');
    tc.verifyNotEmpty(R.acc, 'kind C nonlinear carries the accumulator');
end

function testDatumIsADeclaredPrior(tc)
% categories: contract, reference, robustness
    % R5 M5 5b: a free network's DATUM is a declared prior on the null space
    % (data.priors, L = B'), solved on the BLOCK path -- the same accumulator
    % as everything else, no flat driver. With B spanning the null space the
    % solution is the exact minimum-constraint one for ANY prior weight (the
    % data's range and the null space are orthogonal): pinned against the
    % flat driver's DatumB solution to rounding, and B'x = 0 exactly.
    rs = RandStream('mt19937ar', 'Seed', 5);
    u = 4;  pairs = [1 2; 2 3; 3 4; 4 1; 1 3];
    A = zeros(size(pairs, 1), u);
    for k = 1:size(pairs, 1), A(k, pairs(k, 1)) = -1;  A(k, pairs(k, 2)) = 1; end
    hTrue = [0; 1.5; -0.7; 2.2];
    l = A*hTrue + 0.01*randn(rs, size(A, 1), 1);
    B = ones(u, 1);                                       % the null space: a common level
    data = struct('A', A, 'l', l, 'priors', struct('L', B.', 'sigma', 1, 'name', 'datum'));
    R = lsqvce.fit(data, Verbose=false);
    % tol: 1e-12 absolute -- the constraint holds algebraically on O(1) numbers
    tc.verifyEqual(B.'*R.x, 0, 'AbsTol', 1e-12, 'the minimum constraint holds');
    tc.verifyNotEmpty(R.acc, 'the block path carries the accumulator');
    % the minimum-constraint solution in closed form: the pseudo-inverse
    % solution of the free network (B spans the null space)
    xd = pinv(A.'*A) * (A.'*l);
    % tol: 1e-10 relative -- QR block with the prior vs the dense pseudo-inverse
    tc.verifyEqual(R.x, xd, 'RelTol', 1e-10, 'AbsTol', 1e-12);
    % the weight does not matter when B spans the null space
    data.priors.sigma = 1e3;
    R2 = lsqvce.fit(data, Verbose=false);
    % tol: 1e-8 relative -- a 1e6 weight change moves the solution only through rounding
    tc.verifyEqual(R2.x, R.x, 'RelTol', 1e-8, 'AbsTol', 1e-12);
    % no datum: refused on the block path, cleanly
    lastwarn('');
    tc.verifyEqual(idOf(@() lsqvce.fit(struct('A', A, 'l', l), Verbose=false)), 'lsqvce:qr:singularR');
    [~, wid] = lastwarn;  tc.verifyEmpty(wid);
    % a malformed prior is refused by name
    tc.verifyEqual(idOf(@() lsqvce.fit(struct('A', A, 'l', l, 'priors', struct('L', ones(1, 3))), Verbose=false)), 'lsqvce:fit:badPrior');
end
function testBatchKindLinearReproducesAdjust(tc)
% categories: reference
    [A, l, p] = matrixFixture();
    idx = {1:15, 16:30};
    bf = @(k) deal(A(idx{k}, :), l(idx{k}), p(idx{k}));
    data = struct('batchFcn', bf, 'numUnknowns', size(A, 2), 'numBatches', 2);
    R = lsqvce.fit(data, lsqvce.Options(Verbose=false));
    P = diag(p);  xd = (A.'*P*A) \ (A.'*P*l);        % the dense closed form (O1)
    % tol: 1e-12 relative -- block path (batch k = day k) vs dense normal equations
    tc.verifyEqual(R.x, xd, 'RelTol', 1e-12);
    tc.verifyEqual(R.meta.options.kind, 'batch');
    tc.verifyEqual(R.meta.options.Noise.Grid, 60, 'resolved options carried');
end

function testKindsBAndCRulesFireBeforeAnyRead(tc)
% categories: contract
    [A, l, p] = matrixFixture();
    tc.verifyEqual(idOf(@() lsqvce.fit(struct('A', A, 'l', l(1:end-1)))), 'lsqvce:fit:badMatrix');
    tc.verifyEqual(idOf(@() lsqvce.fit(struct('A', A, 'l', l, 'p', p(1:3)))), 'lsqvce:fit:badMatrix');
    tc.verifyEqual(idOf(@() lsqvce.fit(struct('A', A))), 'lsqvce:fit:missingField');
    tc.verifyEqual(idOf(@() lsqvce.fit(struct('A', A, 'l', l, 'weights', p))), 'lsqvce:fit:unknownField');
    tc.verifyEqual(idOf(@() lsqvce.fit(struct('A', A, 'l', l), lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false))), 'lsqvce:fit:modelKind');
    tc.verifyEqual(idOf(@() lsqvce.fit(struct('batchFcn', @(k) deal(A, l, p)))), 'lsqvce:fit:missingField');
    tc.verifyEqual(idOf(@() lsqvce.fit(struct('batchFcn', @(k) deal(A, l, p), ...
        'numBatches', 1))), 'lsqvce:fit:missingField');          % linear form needs u
    tc.verifyEqual(idOf(@() lsqvce.fit(struct('batchFcn', @(k, x) deal(A, l, p), ...
        'numBatches', 1, 'x0', zeros(3, 1), 'numUnknowns', 3))), 'lsqvce:fit:uCount');
    tc.verifyEqual(idOf(@() lsqvce.fit(struct('batchFcn', @(k) deal(A, l, p), 'numUnknowns', 3, ...
        'numBatches', 1), lsqvce.Options(Estimate="noise", ArcFactors=false))), 'lsqvce:fit:modelKind');
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

function testScopedExpansionReproducesHandNumbering(tc)
% categories: reference, contract, precision
    % THE v4 CLAIM: the scoped tables + info-driven buildDay give EXACTLY the
    % adjustment a hand-numbered ENGINE declaration of the same problem gives.
    % One system, 2 days x 2 arcs; common: 2 'common' + 1 'day' -> nc = 2+2=4;
    % local: 1 'arc' + 1 'day' + 1 'system' -> nl = 2*(2*1+1)+1 = 7.
    cfg = scopedCfg();
    data = struct();
    data.commonParams = [{'common', 'g1'}; {'common', 'g2'}; {'day', 'accScale'}];
    data.days = 1:2;
    s = struct();
    s.buildDay   = @(d, ot, common, local, po, info) scopedBuildDay(cfg, d, ot, info);
    s.covariance = cfg.cov;
    s.localParams = [{'arc', 'ic'}; {'day', 'bias'}; {'system', 'scale'}];
    s.arcs = 2;
    data.systems = s;
    Rf = lsqvce.fit(data, lsqvce.Options(Verbose=false));
    tc.verifyEqual(numel(Rf.x), 4, 'nc = commonCount + nDays*dayRows');
    tc.verifyEqual(Rf.systems(1).nl, 7, 'nl from the scope expansion');

    % the hand-numbered road: ENGINE loader, free 7-row local table, the SAME
    % deterministic day content placed at the SAME (hand-computed) offsets
    L = lsqvce.engine.makeLoader(@(d, ot) handBuildDay(cfg, d, ot), ...
            'ObservationTypes', fieldnames(cfg.cov).', 'Covariance', cfg.cov, ...
            'SystemId', 1, 'LocalParams', num2cell((1:7).'), ...
            'CommonParams', num2cell((1:4).'), 'Arcs', 2, 'Days', 1:2);
    a  = lsqvce.Options(Verbose=false).engineArgs("adjustBlock");
    Rd = lsqvce.engine.adjustBlock({L}, 1:2, a.nv{:});
    % tol: 1e-12 -- identical numbers on two spellings (BLAS order, HANDOVER 10.8)
    tc.verifyEqual(Rf.x, Rd.x, 'RelTol', 1e-12);
    tc.verifyEqual(Rf.sigma0, Rd.sigma0, 'RelTol', 1e-12);
    tc.verifyEqual(Rf.sol.Cxx, Rd.sol.Cxx, 'RelTol', 1e-12);
end

function testDayWithoutArcsGetsNoBlocks(tc)
% categories: contract, robustness
    % The gap rules, both sides: a calendar day whose SCHEDULE carries no arc
    % gets NO common day block (nc shrinks) and NO local blocks (nl shrinks),
    % and every later offset stays gap-safe -- the run goes through.
    cfg = scopedCfg();
    data = struct();
    data.commonParams = [{'common', 'g1'}; {'common', 'g2'}; {'day', 'accScale'}];
    data.days = 1:3;
    s = struct();
    s.buildDay   = @(d, ot, common, local, po, info) scopedBuildDay(cfg, d, ot, info);
    s.covariance = cfg.cov;
    s.localParams = [{'arc', 'ic'}; {'day', 'bias'}; {'system', 'scale'}];
    s.arcs = @(day, obsType) arcListWithGap(day);    % day 2: manoeuvre, no arcs
    data.systems = s;
    R = lsqvce.fit(data, lsqvce.Options(Verbose=false));
    tc.verifyEqual(numel(R.x), 2 + 2*1, 'day 2 must carry no common day block');
    tc.verifyEqual(R.systems(1).nl, 2*(2*1 + 1) + 1, ...
        'day 2 must carry no local blocks; days 1 and 3 unaffected');
end

function ids = arcListWithGap(day)
    if day == 2, ids = zeros(1, 0); else, ids = [1 2]; end
end

%% ============================================ v6 (R6 WP-102): the cells are the arcs
function testScheduleIsOptionalOnTheFitSurface(tc)
% categories: contract, metamorphic
    % v6: `arcs` is optional. Without it the cells buildDay returns are the
    % schedule (ids 1:nArc), and the run is BITWISE the run that declares the
    % count it would have discovered -- the declaration states nothing the
    % data does not already state. 'system'-scoped locals: no plan depends
    % on the arc count.
    [data, ~] = declaration();
    Rs = lsqvce.fit(data, lsqvce.Options(Verbose=false));           % arcs = 2 declared
    d  = data;  d.systems = rmfield(d.systems, 'arcs');
    Rd = lsqvce.fit(d, lsqvce.Options(Verbose=false));              % discovered
    tc.verifyEqual(Rd.x, Rs.x, 'discovered == declared, bitwise');
    tc.verifyEqual(Rd.sigma0, Rs.sigma0);
    tc.verifyEqual(Rd.sol.Cxx, Rs.sol.Cxx);
    tc.verifyEqual(Rd.redundancy(), Rs.redundancy());
end

function testEmptyCellIsTheGapOfAScheduledArc(tc)
% categories: contract, metamorphic
    % v6: under a declared schedule the day returns one cell per scheduled id
    % and an EMPTY cell is that id's data gap. Scheduling a third arc that
    % every day returns empty must give BITWISE the two-arc run: a gap adds
    % nothing, and the redundancy count says so.
    [data, cfg] = declaration();
    R2 = lsqvce.fit(data, lsqvce.Options(Verbose=false));
    d  = data;
    d.systems(1).arcs = [1 2 3];                                     % arc 3 is scheduled ...
    d.systems(1).buildDay = @(dd, ot) padGap(buildDayImpl(cfg, dd, ot), cfg);  % ... and empty
    R3 = lsqvce.fit(d, lsqvce.Options(Verbose=false));
    tc.verifyEqual(R3.x, R2.x, 'a gap cell contributes nothing');
    tc.verifyEqual(R3.sigma0, R2.sigma0);
    tc.verifyEqual(R3.redundancy(), R2.redundancy());
    % the count is enforced: a schedule of 3 and two cells is refused by name
    d.systems(1).buildDay = @(dd, ot) buildDayImpl(cfg, dd, ot);
    tc.verifyEqual(idOf(@() lsqvce.fit(d, lsqvce.Options(Verbose=false))), 'lsqvce:makeLoader:arcCount');
end

function D = padGap(D, cfg)
% Append the empty cell of a scheduled-but-absent third arc.
    D.A{3} = zeros(0, cfg.nc);  D.B{3} = zeros(0, cfg.nl);  D.l{3} = zeros(0, 1);  D.t{3} = zeros(0, 1);
end

function testDayRowsWithoutScheduleLiveOnEveryCalendarDay(tc)
% categories: contract, robustness
    % v6: without `arcs`, the calendar is the declaration -- every calendar
    % day owns its 'day' blocks (common and local) and nc / nl count them --
    % and a calendar day that then returns no arc in any stream is refused
    % at its 'meta' probe (emptyDay), before any solve could go rank-deficient
    % on an unobserved block.
    cfg = scopedCfg();
    data = struct();
    data.commonParams = [{'common', 'g1'}; {'common', 'g2'}; {'day', 'accScale'}];
    data.days = 1:3;
    s = struct();
    s.buildDay   = @(d, ot, common, local, po, info) scopedBuildDay(cfg, d, ot, info, 2);
    s.covariance = cfg.cov;
    s.localParams = [{'day', 'bias'}; {'system', 'scale'}];         % no 'arc' rows: no schedule needed
    data.systems = s;
    R = lsqvce.fit(data, lsqvce.Options(Verbose=false));
    tc.verifyEqual(numel(R.x), 2 + 3, 'every calendar day owns a common day block');
    tc.verifyEqual(R.systems(1).nl, 3 + 1, 'every calendar day owns its local day block');
    % day 2 returns no cells: refused by name, at the read
    d = data;
    d.systems(1).buildDay = @(dd, ot, common, local, po, info) emptyOnDay2(cfg, dd, ot, info);
    tc.verifyEqual(idOf(@() lsqvce.fit(d, lsqvce.Options(Verbose=false))), 'lsqvce:makeLoader:emptyDay');
    % the same day taken off the calendar is the declared way to say so
    d.systems(1).days = [1 3];
    R2 = lsqvce.fit(d, lsqvce.Options(Verbose=false));
    tc.verifyEqual(numel(R2.x), 2 + 2, 'day 2 off this system''s calendar: no common day block');
    tc.verifyEqual(R2.systems(1).nl, 2 + 1);
    % 'arc' rows without a schedule: refused at the declaration, relabelled
    d = data;  d.systems(1).localParams = [{'arc', 'ic'}; {'day', 'bias'}];
    tc.verifyEqual(idOf(@() lsqvce.fit(d, lsqvce.Options(Verbose=false))), 'lsqvce:fit:localPlan');
end

function D = emptyOnDay2(cfg, day, obsType, info)
    if day == 2
        D = struct('A', {{}}, 'B', {{}}, 'l', {{}}, 't', {{}});
    else
        D = scopedBuildDay(cfg, day, obsType, info, 2);
    end
end

function cfg = scopedCfg()
% Deterministic content generator shared by BOTH roads of the parity test.
    cfg = struct('mArc', 120, 'dt', 5, 'seed', 31);
    cfg.cov.pos = struct('type', 'diag', 'sigma', 0.5);
end

function [Ael, Bel, l, tt] = scopedDayContent(cfg, day, obsType, nArc)
% The day's raw content, keyed ONLY on (day, obsType): per-arc design columns
% for 1 IC each, one day column, one system column, deterministic rhs.
    M  = nArc * cfg.mArc;
    sd = mod(cfg.seed + 7*day + 31*sum(double(obsType)), 2^31-1);
    rs = RandStream('mt19937ar', 'Seed', sd);
    Ael = randn(rs, M, 4);                       % nc = 4 columns of partials
    Bel = randn(rs, M, nArc + 2);                % [arc cols | day col | sys col]
    l   = randn(rs, M, 1);
    tt  = (0:M-1).' * cfg.dt;
end

function D = scopedBuildDay(cfg, day, obsType, info, nArc)
% The v4 road: every placement comes from INFO, no literal offset anywhere.
% v6: without a declared schedule info.arcId is empty and the caller says
% how many arcs (cells) the day carries.
    if nargin < 5, nArc = numel(info.arcId); end
    [Ael, Bel, D.l, D.t] = scopedDayContent(cfg, day, obsType, nArc);
    M = size(Ael, 1);
    D.A = zeros(M, info.nc);
    D.A(:, 1:info.commonCount) = Ael(:, 1:info.commonCount);
    if info.commonDayCount > 0
        D.A(:, info.commonDayOffset + (1:info.commonDayCount)) = ...
            Ael(:, info.commonCount + (1:info.commonDayCount));
    end
    D.B = zeros(M, info.nl);
    D.arcIdx = cell(1, nArc);
    for a = 1:nArc
        D.arcIdx{a} = (a-1)*cfg.mArc + (1:cfg.mArc);
        if ~isempty(info.arcId)                      % 'arc' rows exist only with a schedule
            D.B(D.arcIdx{a}, info.offsetOf(info.arcId(a)) + 1) = Bel(D.arcIdx{a}, a);
        end
    end
    if info.dayCount > 0
        D.B(:, info.dayOffset + 1) = Bel(:, nArc + 1);
    end
    D.B(:, info.sysOffset + 1) = Bel(:, nArc + 2);
    D = arcCells(D);                                 % v6: the cells are the arcs
end

function D = handBuildDay(cfg, day, obsType)
% The hand-numbered road: SAME content, literal v4-layout offsets for the
% 2-day, 2-arc, gap-free case (nc: [2 common | day1 | day2]; nl per day:
% [arc1 | arc2 | day], system block last).
    nArc = 2;
    [Ael, Bel, D.l, D.t] = scopedDayContent(cfg, day, obsType, nArc);
    M = size(Ael, 1);
    D.A = zeros(M, 4);
    D.A(:, 1:2) = Ael(:, 1:2);
    D.A(:, 2 + day) = Ael(:, 3);                 % this day's common day column
    D.B = zeros(M, 7);
    base = (day - 1) * 3;                        % 3 local cols per day
    D.arcIdx = cell(1, nArc);
    for a = 1:nArc
        D.arcIdx{a} = (a-1)*cfg.mArc + (1:cfg.mArc);
        D.B(D.arcIdx{a}, base + a) = Bel(D.arcIdx{a}, a);
    end
    D.B(:, base + 3) = Bel(:, nArc + 1);         % the day block
    D.B(:, 7)        = Bel(:, nArc + 2);         % the system block
    D = arcCells(D);                             % v6: the cells are the arcs
end

function [data, cfg] = declaration()
% The canonical single-system declaration over the makeLoader-path fixture
% (v3 contract: the tables declare, the counts are derived).
    cfg = buildCfg();
    data = struct();
    data.commonParams = [repmat({'common'}, cfg.nc, 1), num2cell((1:cfg.nc).')];
    data.days = 1:3;                                  % v4: all-'common' -> nc = rows
    s = struct();
    s.buildDay = @(d, ot) buildDayImpl(cfg, d, ot);
    s.covariance = cfg.cov;                           % fieldnames = pos, kbr
    s.localParams = [repmat({'system'}, cfg.nl, 1), num2cell((1:cfg.nl).')];
    s.arcs = cfg.arcsPerDay;                          % 'system' rows -> nl = rows,
    data.systems = s;                                 % B layout as in v3; id = 1
end

function L = loaderOf(~, s, data)
% What fit builds for one system -- written out, so the tests compare two
% independent spellings of the same declaration (v4: obsTypes from the
% covariance fieldnames, SCOPED tables + the common plan, id -> SystemId).
    D = numel(data.days);
    cp = struct('nc', size(data.commonParams, 1), ...
                'commonCount', size(data.commonParams, 1), ...
                'days', data.days, 'dayOffset', zeros(1, D), ...
                'dayCount', zeros(1, D));       % all-'common' table: trivial plan
    L = lsqvce.engine.makeLoader(s.buildDay, ...
            'ObservationTypes', fieldnames(s.covariance).', ...
            'Covariance', s.covariance, 'SystemId', 1, ...
            'LocalParams', s.localParams, 'Scoped', true, 'CommonPlan', cp, ...
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
    D = arcCells(D);                                  % v6: the cells are the arcs
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

function testLocalPriorsAreTheCounterpartOfDataPriors(tc)
% categories: contract, precision
    % R5-D8 round 5: data.systems(k).priors(l) mirrors data.priors(j) -- one
    % struct per GROUP with L in the columns of B (lsqvce.layout announces
    % them), x0, sigma, from, name, estimate. Each group rides as its own
    % stream prior_<name>. A tight group pins the local to x0; several
    % structs are several groups with their own sigma; a wrong width, a row
    % coupling two days and a bad policy are refused by name; a loose group
    % leaves the solution alone; from gates a group out of a single solve.
    [data, ~] = declaration();
    R0 = lsqvce.fit(data, Verbose=false);
    y0 = R0.residuals('Verbose', false).sys(1).y;
    lay = lsqvce.layout(data);
    nl = lay.systems(1).nl;
    tc.verifyEqual(nl, numel(y0), 'the layout announces B''s width');
    j = 2;  target = 0.123;
    cj = lay.systems(1).col(j, [], []);                 % a 'system' row: one column
    tc.verifyGreaterThan(cj, 0);
    e = @(c) double((1:nl) == c);
    data.systems(1).priors = [struct('L', e(cj), 'x0', target, 'sigma', 1e-6, 'from', 1, 'name', 'pin', 'estimate', "none"), ...
                              struct('L', e(lay.systems(1).col(1, [], [])), 'x0', 0, 'sigma', 1e3, 'from', 1, 'name', 'soft', 'estimate', "none")];
    R1 = lsqvce.fit(data, Verbose=false);
    y1 = R1.residuals('Verbose', false).sys(1).y;
    % tol: 1e-5 absolute -- the group's sigma (1e-6) against the data's O(1)
    % leverage on that parameter; the pin holds to a few sigma
    tc.verifyEqual(y1(j), target, 'AbsTol', 1e-5, 'the pinned local equals x0');
    tc.verifyEqual(R1.diag.n - R0.diag.n, 2, 'two pseudo-observations, one per group row');
    tc.verifyNotEqual(y1(j), y0(j));
    % a loose group alone leaves the solution where it was
    data.systems(1).priors = struct('L', e(lay.systems(1).col(1, [], [])), 'x0', 0, 'sigma', 1e6, 'name', 'loose');
    Rl = lsqvce.fit(data, Verbose=false);
    % tol: 1e-8 relative -- a 1e6-sigma prior on O(1) parameters perturbs at 1e-12
    tc.verifyEqual(Rl.x, R0.x, 'RelTol', 1e-8);
    % from gates the group out of a single solve
    data.systems(1).priors = struct('L', e(cj), 'x0', target, 'sigma', 1e-6, 'from', 2, 'name', 'late');
    R2 = lsqvce.fit(data, Verbose=false);
    tc.verifyEqual(R2.diag.n, R0.diag.n);
    % the same group by NAME PATTERN (column 2 of localParams): fit builds L
    names = data.systems(1).localParams(:, 2);
    if all(cellfun(@(v) ischar(v) || isstring(v), names))
        pat = ['^' regexptranslate('escape', char(string(names{j}))) '$'];
        data.systems(1).priors = struct('group', pat, 'x0', target, 'sigma', 1e-6);
        Rp = lsqvce.fit(data, Verbose=false);
        yp = Rp.residuals('Verbose', false).sys(1).y;
        % tol: 1e-5 absolute -- as above
        tc.verifyEqual(yp(j), target, 'AbsTol', 1e-5, 'the pattern form pins the same local');
    end
    % refusals by name (the loader's ids re-homed by fit at build time)
    data.systems(1).priors = struct('L', ones(1, nl + 3), 'sigma', 1);
    tc.verifyEqual(idOf(@() lsqvce.fit(data, Verbose=false)), 'lsqvce:fit:badPrior');
    data.systems(1).priors = struct('L', e(cj), 'sigma', 1, 'estimate', "law");
    tc.verifyEqual(idOf(@() lsqvce.fit(data, Verbose=false)), 'lsqvce:fit:badPrior');
end

function testLocalPriorGroupIsEstimatedLikeACommonOne(tc)
% categories: contract, reference
    % R5-D8 round 5 (Koch-Kusche on locals): a local group with estimate
    % "scale" is a stream of its own, so under Estimate="noise" its sigma is
    % the stream's level -- the Helmert factor the joint estimator computes
    % for every stream -- exactly the common groups' rule. Fixture: the demo
    % declaration; a damping group on all locals of system 1 started at a
    % wrong (tight) sigma must be loosened by the data, and its level appears
    % on R.noise.model.
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(fullfile(root, 'examples'));
    [~, data] = demo_grace_dataLoader();
    lay = lsqvce.layout(data);
    S = lay.systems(1);  nl = S.nl;
    L = eye(nl);                                   % all locals of system 1 (its system block)
    data.systems(1).priors = struct('L', L, 'x0', zeros(size(L, 1), 1), 'sigma', 1e-4, 'name', 'damp', 'estimate', "scale");
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');
    R = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, MaxPasses=5, ...
                     Noise=struct('Grid', 16), Seed=3, Verbose=false));
    k = find(strcmp(R.noise.stream, 'x1|prior_damp'), 1);
    tc.verifyNotEmpty(k, 'the group is a stream of its own');
    tc.verifyGreaterThan(R.noise.model{k}.sigma, 1, 'the tight damping (1e-4 on O(1) locals) is loosened by the data');
    tc.verifyTrue(any(strcmp(R.history.streams, 'x1|prior_damp')));
end

function testDemoPriorsIsTheWorkedLayoutExample(tc)
% categories: contract, precision
    % R6 WP-96: examples/demo_priors.m is the ONE worked example of
    % lsqvce.layout -- a rank-deficient calibration case whose datum can only
    % be written in the columns the layout announces. A demo nobody executes
    % rots, and this one carries measured claims, so both are checked: it must
    % RUN, and the claims of its ACCURACY block must hold.
    here  = fileparts(mfilename('fullpath'));
    exDir = fullfile(here, '..', '..', 'examples');    % dev/tests -> root (L8)
    addpath(exDir);
    tc.addTeardown(@() rmpath(exDir));
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');
    out = evalc('[R, s] = demo_priors();');

    tc.verifyClass(R, 'lsqvce.Results');
    tc.verifyTrue(contains(out, 'lsqvce.layout says'), 'it prints the layout census');
    tc.verifyEqual([s.nc, s.nl], [7, 13], 'nc and nl are derived from the two tables');
    tc.verifyEqual(s.nullDim, 3, 'one null direction per day');
    tc.verifyEqual(s.unregularised, 'lsqvce:qr:rankDeficient', ...
                   'without priors the singular case is REFUSED, not solved');
    tc.verifyEqual(s.nPseudo, 9, '1 datum + 2 smoothness + 6 drift rows');
    % the two ways of writing ONE group must give ONE solution: same operator,
    % same normals, so the difference is exactly zero, not merely small
    tc.verifyEqual(s.formGap, 0, 'the L form and the group pattern agree bitwise');
    % tol: 1e-6 absolute -- the observable part is datum-INVARIANT by algebra;
    % what is left is the round-off of two independent solves (measured 1.6e-11)
    tc.verifyLessThan(s.datumGap, 1e-6, 'moving the datum does not move the signal');
    tc.verifyLessThan(max(abs(s.zSignal)), 3, 'the signal is recovered within 3 sigma');
    tc.verifyGreaterThan(s.priorLevel, 0, 'the prior group is a stream with a level');
end
