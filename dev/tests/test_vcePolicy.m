function tests = test_vcePolicy
%TEST_VCEPOLICY  The iterative Models as POLICIES of the one estimator (R5 M4 4c).
%   tests = test_vcePolicy
%
%   Re-homes the user-facing capabilities the retired separate mode
%   (lsqvce.engine.adjustVCE, R5 finding H-03, decision R5-D4) used to be
%   tested for, on the USER surface through lsqvce.fit with Model "vce" /
%   "pooled" / "ellmer": the convergence table and the per-stream summary,
%   the contribution shares, the pooled mode end to end, the honest gate, one
%   law per pool group, the independent-satellite stream, verbosity, and the
%   recovery of planted relative weights (a reference against a known truth
%   -- the kind of oracle that caught three defects in round 4b where a
%   parity test could not).
%
%   See also LSQVCE.FIT, LSQVCE.ENGINE.ADJUSTVCEJOINT, TEST_ELLMER, TEST_PHYSORACLE
%
%   Created 2026-09-03 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M4, WP-85 round 4c.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    pathBefore = path;   % R7-45: a standalone run (run_lsqvce / runtests) leaves the path as found
    tc.addTeardown(@() path(pathBefore));
    here = fileparts(mfilename('fullpath'));
    root = fileparts(fileparts(here));
    addpath(root);  addpath(fullfile(root, 'examples'));
    w = warning('off', 'lsqvce:adjustVCEjoint:notConverged');
    tc.addTeardown(@() warning(w));
    [~, data] = demo_grace_dataLoader();
    tc.TestData.data = data;
    % ONE shared run per model (the rebuild-per-test lesson of the old suite)
    tc.TestData.Rv = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, MaxPasses=3, Noise=struct('Grid', 16), Seed=3, Verbose=false));
    tc.TestData.Rp = lsqvce.fit(data, lsqvce.Options(Estimate="noise", ArcFactors=false, MaxPasses=3, Noise=struct('Grid', 16), Seed=3, Verbose=false));
end

%% ---------------------------------------------------------------- history
function testConvergenceTableShape(tc)
% categories: contract | tags: diagnostics
    R = tc.TestData.Rv;
    T = R.convergence();
    tc.verifyEqual(height(T), R.history.nIter, 'one row per pass');
    tc.verifyEqual(size(R.history.sigma, 2), numel(R.history.streams));
    tc.verifyEqual(size(R.history.white), size(R.history.sigma));
    tc.verifyEqual(size(R.history.whitePval), size(R.history.sigma));
end

function testSummaryIsPerStream(tc)
% categories: contract | tags: diagnostics
    R = tc.TestData.Rv;
    S = R.summary();
    tc.verifyEqual(height(S), numel(R.history.streams), 'one row per stream');
end

function testConvergenceErrorsWithoutHistory(tc)
% categories: contract, robustness | tags: diagnostics
    rs = RandStream('mt19937ar', 'Seed', 2);
    A = randn(rs, 20, 2);  l = A*[1; 2] + 0.01*randn(rs, 20, 1);
    R = lsqvce.fit(struct('A', A, 'l', l), lsqvce.Options(Verbose=false));
    tc.verifyError(@() R.convergence(), 'lsqvce:Results:noHistory');
    tc.verifyError(@() R.summary(),     'lsqvce:Results:noHistory');
end

function testStreamsIncludeTheIndependentSatellite(tc)
% categories: contract
    % the demo's system 2 ('E') is its own stream key, keyed satKey|obsType
    R = tc.TestData.Rv;
    tc.verifyEqual(numel(R.history.streams), 4);
    tc.verifyTrue(any(endsWith(R.history.streams, '|posE')), 'the independent satellite is a stream');
    tc.verifyEqual(numel(R.noise.stream), 4);
end

%% ------------------------------------------------------------ contribution
function testContributionRowsSumToOne(tc)
% categories: precision, robustness
    R = tc.TestData.Rv;
    con = R.contribution();
    shares = arrayfun(@(g) g.share, con.group);
    % tol: 1e-10 absolute -- the shares are normalised to one over the groups;
    % the normalisation is exact, the rounding is not
    tc.verifyEqual(sum(shares), 1, 'AbsTol', 1e-10);
    tc.verifyTrue(all(shares >= 0));
end

%% ------------------------------------------------------------ pooled mode
function testPooledRunsEndToEnd(tc)
% categories: contract, reference, robustness
    R = tc.TestData.Rp;
    tc.verifyClass(R, 'lsqvce.Results');
    tc.verifyEqual(keys(R.noise.Cgroups), {'all'}, 'pooled: ONE law');
    tc.verifyTrue(all(isfinite(R.x)));
    tc.verifyTrue(all(isfinite(R.noise.C)));
end

function testGateRequiresSettledSpectrum(tc)
% categories: robustness
    % three passes on the demo do not settle the PSD: converged is FALSE and
    % says so, the last pass is returned (convergence failure is a RESULT)
    R = tc.TestData.Rp;
    tc.verifyFalse(R.converged);
    tc.verifyEqual(R.stopReason, "maxPasses");
    tc.verifyTrue(all(isfinite(R.x)));
end

function testPolicyUpdatesEveryStream(tc)
% categories: contract
    % every stream's level moved away from 1 on a demo declared white
    R = tc.TestData.Rv;
    tc.verifyTrue(all(abs(R.history.sigma(end, :) - 1) > 1e-6), 'every stream was updated');
    tc.verifyEqual(numel(keys(R.noise.Cgroups)), 4, 'vce: one law per stream');
end

function testPoolByObsTypeCoarsens(tc)
% categories: contract
    % Pool "obsType" under Model "ellmer": one law per observation type --
    % the three position streams (posC, posD, posE) are three types on the
    % demo, so 4 groups collapse to... the number of distinct types
    data = tc.TestData.data;
    R = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="obsType", MaxPasses=2, Noise=struct('Grid', 16), Seed=3, Verbose=false));
    types = unique(cellfun(@(s) extractAfter(s, '|'), R.history.streams, 'UniformOutput', false));
    tc.verifyEqual(sort(keys(R.noise.Cgroups)), sort(types));
end

function testStreamModelCarriesOnlyItsOwnEstimatedLaw(tc)
% categories: contract
    % R6 WP-99: R.noise.model{k}.C is documented as "the stream's law". Until
    % this round it fell back to the GLOBAL out.C for a stream whose law was
    % not estimated, so a "none" or "scale" stream reported a law fitted from
    % somebody else's arcs -- measured on a four-stream run with ONE "law"
    % stream, all four .C came back identical. A plausible wrong value returned
    % in silence is what R4 forbids, so the field is now EMPTY there: the law
    % that stream used is the one the user declared, and the toolbox did not
    % touch it. Fixture: the demo declaration with ONE stream demoted to
    % "scale", so both branches are exercised in one run.
    data = tc.TestData.data;
    ots = fieldnames(data.systems(1).covariance);
    demoted = ots{1};
    data.systems(1).covariance.(demoted).estimate = "scale";
    R = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", MaxPasses=2, ...
                     Noise=struct('Grid', 16), Seed=3, Verbose=false));
    kDem = find(strcmp(extractAfter(string(R.noise.stream), '|'), demoted), 1);
    tc.verifyNotEmpty(kDem, 'the demoted stream is still reported');
    tc.verifyEmpty(R.noise.model{kDem}.C, ...
                   'a stream whose law was NOT estimated reports no law');
    tc.verifyGreaterThan(R.noise.model{kDem}.sigma, 0, ...
                         'its LEVEL is still estimated -- that is what "scale" means');
    others = setdiff(1:numel(R.noise.stream), kDem);
    for k = others
        tc.verifyNotEmpty(R.noise.model{k}.C, 'a "law" stream reports the law it fitted');
        % and it is ITS pool group's law, not some other group's
        tc.verifyEqual(R.noise.model{k}.C, R.noise.Cgroups(R.noise.stream{k}), ...
                       'the law reported is the one keyed on this stream');
    end
end

function testReportSaysWhetherTheStochasticModelSettled(tc)
% categories: contract
    seedFor(tc, 162);   % R7-45: drew from the global stream unseeded (suite-order dependent)
    % R6 WP-99: lsqvce.report read NOTHING of the iteration state, so a run of
    % Estimate="noise" that stopped on maxPasses -- whose weights were still
    % moving -- got a verdict that said nothing about it. The engine warns, and
    % warnings are routinely muted (this suite mutes that one in three places),
    % so the verdict has to carry it. Three cases, because the interesting part
    % is what it must NOT say: a single solve has no loop and must stay silent.
    A = randn(200, 5);  xt = (1:5).';
    R0 = lsqvce.fit(struct('A', A, 'l', A*xt + 0.1*randn(200, 1)), ...
                    lsqvce.Options(Verbose=false));
    V0 = R0.report('Verbose', false);
    tc.verifyEmpty(codesOf(V0, 'NOISE_'), 'a single solve says nothing about a loop it never ran');

    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');
    data = tc.TestData.data;
    Ru = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", MaxPasses=1, ...
                     Noise=struct('Grid', 16), Seed=3, Verbose=false));
    tc.assumeFalse(Ru.converged, 'the one-pass run is the unsettled case');
    Vu = Ru.report('Verbose', false);
    tc.verifyEqual(codesOf(Vu, 'NOISE_'), {'NOISE_UNSETTLED'}, 'it names the unsettled model');
    f = Vu.findings{find(strcmp(cellfun(@(q) q.code, Vu.findings, 'UniformOutput', false), ...
                               'NOISE_UNSETTLED'), 1)};
    tc.verifyEqual(f.severity, 'WARN', 'a WARN: it qualifies every level below, it does not void them');
    tc.verifyEqual(f.level, 1, 'level 1 -- the weights are what level 1 judges');
    tc.verifyTrue(contains(f.msg, 'maxPasses'), 'and it says WHY the loop stopped');
end

function codes = codesOf(V, prefix)
%CODESOF  The finding codes of V that start with prefix, in order.
    all = cellfun(@(f) f.code, V.findings, 'UniformOutput', false);
    codes = all(startsWith(all, prefix));
end

%% ------------------------------------------------------------- messaging
function testVerboseZeroIsSilent(tc)
% categories: contract, robustness | tags: diagnostics
    data = tc.TestData.data;
    o = lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, MaxPasses=1, Noise=struct('Grid', 16), Seed=3, Verbose=0);
    out = evalc('lsqvce.fit(data, o);');
    tc.verifyEmpty(strtrim(out), 'Verbose 0 prints nothing');
end

function testVerboseOnePrintsThePassRecord(tc)
% categories: contract | tags: diagnostics
    data = tc.TestData.data;
    o = lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, MaxPasses=1, Noise=struct('Grid', 16), Seed=3, Verbose=1);
    out = evalc('lsqvce.fit(data, o);');
    tc.verifyNotEmpty(strtrim(out), 'Verbose 1 prints the run''s progress');
end

%% ---------------------------------------------------- reference: planted weights
function testRecoversPlantedRelativeWeights(tc)
% categories: reference, precision
    % The gravity demo plants stream sigmas 1 / 2 / 0.5 (potC, potD, kbr)
    % behind a declaration of 1 everywhere. Model "vce" must recover the
    % RELATIVE weights (a common factor is sigma0's). This is the oracle
    % that found the missing level estimate and the ignored psd level in
    % round 4b -- promoted to a counted reference point (R5 M7 item, early).
    [cfg, data] = demo_gravity_dataLoader();
    R = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, MaxPasses=6, Seed=1, Verbose=false));
    got = cellfun(@(m) m.sigma, R.noise.model);
    want = cellfun(@(s) cfg.trueSigma.(extractAfter(s, '|')), R.noise.stream);
    gotR = got / got(1);  wantR = want / want(1);
    % tol: 25 % relative -- the coloured potC stream carries the shape/level
    % split's bias on 60-lag laws from few arcs (measured 1.18 / 2.10 / 0.51
    % vs 1 / 2 / 0.5 on 2026-09-02); the same bound test_physOracle uses
    tc.verifyEqual(gotR, wantR, 'RelTol', 0.25, 'the planted relative weights are recovered');
end

%% ---------------------------------------------------------- re-reads (M4 4d)
function testCacheDaysCountsTheBuilds(tc)
% categories: contract, speed
    % R5 M4 4d: the number a sweep level is judged by. A two-stage pass builds
    % every (day, obsType) twice (accumulate + sweep) and a VCE run of n
    % passes 2n+1 times (+ the final diagnostics); with CacheDays = Inf the
    % calendar stays resident and every day is built ONCE per run. Counted
    % through a wrapped buildDay on the demo declaration (2 systems, 3 days).
    data = tc.TestData.data;
    for k = 1:numel(data.systems)
        cnt = containers.Map('KeyType', 'char', 'ValueType', 'double');
        inner = data.systems(k).buildDay;
        data.systems(k).buildDay = @(varargin) countBuild(cnt, inner, varargin{:});
        data.systems(k).processingOptions.counter = cnt;      % keep the map reachable
    end
    n = 2;
    cnts = cellfun(@(s) s.counter, {data.systems.processingOptions}, 'UniformOutput', false);
    lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, MaxPasses=n, Noise=struct('Grid', 16), Seed=3, Verbose=false, CacheDays=1));
    v1 = cell2mat(cellfun(@(c) cell2mat(values(c)), cnts, 'UniformOutput', false));
    tc.verifyEqual(unique(v1), 2*n + 1, 'CacheDays 1: accumulate + sweep per pass, + the final diagnostics');
    for c = cnts, remove(c{1}, keys(c{1})); end
    lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, MaxPasses=n, Noise=struct('Grid', 16), Seed=3, Verbose=false, CacheDays=Inf));
    vInf = cell2mat(cellfun(@(c) cell2mat(values(c)), cnts, 'UniformOutput', false));
    tc.verifyEqual(unique(vInf), 1, 'CacheDays Inf: every day is built once per run');
end

function D = countBuild(cnt, inner, varargin)
    key = sprintf('%d|%s', varargin{1}, char(varargin{2}));
    if isKey(cnt, key), cnt(key) = cnt(key) + 1; else, cnt(key) = 1; end
    n = nargin(inner);                          % pass what the builder declares
    if n < 0, n = numel(varargin); end
    D = inner(varargin{1:min(n, numel(varargin))});
end

%% ------------------------------------------------------------- plot (M5 5c)
function testPlotReturnsFiguresForEveryKind(tc)
% categories: contract | tags: diagnostics
    % R5 M5 5c: R.plot(kind) replaces the retired plot hooks -- one figure per
    % question, nothing exported by the toolbox. Headless-safe (Visible=false).
    R = tc.TestData.Rv;
    fhs = gobjects(0);
    for kind = ["residuals", "contribution", "psd", "convergence"]
        fh = R.plot(kind, Visible=false);
        tc.verifyClass(fh, 'matlab.ui.Figure', sprintf('%s returns a figure', kind));
        tc.verifyFalse(isempty(findall(fh, 'Type', 'axes')), sprintf('%s draws axes', kind));
        fhs(end+1) = fh; %#ok<AGROW>
    end
    close(fhs);
    fh = R.plot("residuals", Stream=R.history.streams{1}, Visible=false);
    tc.verifyEqual(numel(findall(fh, 'Type', 'axes')), 1, 'Stream= selects one panel');
    close(fh);
    tc.verifyError(@() R.plot("hooks", Visible=false), 'lsqvce:Results:badPlotKind');
    tc.verifyError(@() R.plot("residuals", Stream="x9|nothing", Visible=false), 'lsqvce:Results:badPlotKind');
end

function testPlotRefusesWhatTheResultsLacks(tc)
% categories: contract, robustness | tags: diagnostics
    rs = RandStream('mt19937ar', 'Seed', 2);
    A = randn(rs, 20, 2);  l = A*[1; 2] + 0.01*randn(rs, 20, 1);
    R = lsqvce.fit(struct('A', A, 'l', l), lsqvce.Options(Verbose=false));
    tc.verifyError(@() R.plot("psd", Visible=false), 'lsqvce:Results:noHistory');
    tc.verifyError(@() R.plot("convergence", Visible=false), 'lsqvce:Results:noHistory');
    fh = R.plot("residuals", Visible=false);            % a single solve has residuals
    tc.verifyClass(fh, 'matlab.ui.Figure');  close(fh);
end

%% ------------------------------------------------- OnPass (M5 5c'): the inspection seam
function testOnPassFiresEveryPassWithTheContext(tc)
% categories: contract | tags: diagnostics
    % The project's opening requirement: after EVERY iteration, a hook for
    % inspection. Options.OnPass fires at the end of every pass of the joint
    % estimator with the declaration, the processing options, the live
    % accumulator and the pass state; the hook may plot, save, or stop.
    data = tc.TestData.data;
    seen = containers.Map('KeyType', 'double', 'ValueType', 'any');
    hook = @(ctx) recordPass(seen, ctx);
    R = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, ...
                     MaxPasses=3, Noise=struct('Grid', 16), Seed=3, Verbose=false, OnPass=hook));
    tc.verifyEqual(double(seen.Count), 3, 'one call per pass');
    c = seen(2);
    tc.verifyEqual(c.pass, 2);
    tc.verifyEqual(c.driver, 'adjustVCEjoint');
    tc.verifyEqual(numel(c.x), numel(R.x));
    tc.verifyTrue(isfield(c, 'data') && isequal(c.data.days, data.days), 'the declaration as given');
    tc.verifyTrue(isfield(c, 'processingOptions') && iscell(c.processingOptions) && ...
                  numel(c.processingOptions) == numel(data.systems), 'one processingOptions per system');
    tc.verifyEqual(c.processingOptions{1}, data.systems(1).processingOptions);
    tc.verifyClass(c.acc, 'lsqvce.engine.BlockAccumulator', 'the live accumulator');
    tc.verifyEqual(numel(c.streams), numel(c.sigma));
    tc.verifyTrue(isfield(c, 'Cgroups') && isfield(c, 'sigma2') && isfield(c, 'whitePval'));
    % the accumulator inside the hook is usable: a residual sweep at pass 2
    res = lsqvce.engine.residuals(c.acc, struct('x', c.x, 'sigma0', c.sigma0, 's0sq', c.s0sq, 'Cxx', c.Cxx), 'Verbose', false);
    tc.verifyTrue(isstruct(res) && isfield(res, 'sys'));
end

function testOnPassCanStopTheLoop(tc)
% categories: contract, robustness | tags: diagnostics
    data = tc.TestData.data;
    stopAt2 = @(ctx) ctx.pass < 2;                % false after the second pass
    R = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, ...
                     MaxPasses=6, Noise=struct('Grid', 16), Seed=3, Verbose=false, OnPass=stopAt2));
    tc.verifyEqual(R.history.nIter, 2);
    tc.verifyEqual(R.stopReason, "userStop");
    tc.verifyFalse(R.converged);
    tc.verifyTrue(all(isfinite(R.x)), 'the last pass is returned with its solution');
    % IRLS and Gauss-Newton fire the same hook
    rs = RandStream('mt19937ar', 'Seed', 6);
    A = randn(rs, 40, 3);  l = A*[1; -2; 0.5] + 0.01*randn(rs, 40, 1);  l(5) = l(5) + 10;
    seen = containers.Map('KeyType', 'double', 'ValueType', 'any');
    Rr = lsqvce.fit(struct('A', A, 'l', l), Robust="huber", Verbose=false, OnPass=@(ctx) recordPass(seen, ctx));
    tc.verifyEqual(double(seen.Count), numel(Rr.history), 'IRLS: one call per pass');
    tc.verifyEqual(seen(1).driver, 'adjustRobust');
    tc.verifyTrue(isfield(seen(1), 'factor') && isfield(seen(1), 'data'));
    tc.verifyEqual(seen(1).data.A, A, 'kind B: the declaration is the ctx.data');
    seen2 = containers.Map('KeyType', 'double', 'ValueType', 'any');
    t = linspace(0, 1, 30).';  lobs = exp(0.7*t) - 0.2 + 1e-3*randn(rs, 30, 1);
    bf = @(k, x) deal([t.*exp(x(1)*t), ones(30, 1)], lobs - (exp(x(1)*t) + x(2)), []);
    Rn = lsqvce.fit(struct('batchFcn', bf, 'numBatches', 1, 'x0', [0.5; 0]), Verbose=false, ...
                    OnPass=@(ctx) recordPass(seen2, ctx));
    tc.verifyEqual(double(seen2.Count), numel(Rn.history), 'Gauss-Newton: one call per pass');
    tc.verifyEqual(seen2(1).driver, 'adjustNonlinear');
    % a throwing hook is not swallowed: it arrives with the cause
    try
        lsqvce.fit(data, lsqvce.Options(Estimate="noise", MaxPasses=2, Noise=struct('Grid', 16), ...
                   Seed=3, Verbose=false, OnPass=@(ctx) error('mine:hook', 'plot failed')));
        tc.verifyFail('must throw');
    catch ME
        tc.verifyEqual(ME.identifier, 'lsqvce:callUser:threw');
        tc.verifyEqual(ME.cause{1}.identifier, 'mine:hook');
    end
    tc.verifyError(@() lsqvce.Options(OnPass=42), 'lsqvce:Options:badOnPass');
end

function recordPass(seen, ctx)
    c = ctx;  c.acc = ctx.acc;                    % keep the handle, drop nothing
    seen(ctx.pass) = c;
end

%% --------------------------------------------- the prior's own component (M7 7b)
function testPriorScaleReachesTheHelmertFixedPoint(tc)
% categories: contract, precision | tags: diagnostics
    % R5 M7 7b (Koch-Kusche): a prior group with estimate "scale" gets one
    % Helmert step per pass -- sigma_k^2 <- sigma_k^2 * Omega_k / r_k -- and
    % at the fixed point the prior's weighted residual sum equals its
    % redundancy, Omega_k = r_k. Fixture: the demo declaration with a plain
    % damping prior on all common parameters, started at a WRONG sigma.
    data = tc.TestData.data;
    R0 = lsqvce.fit(data, Verbose=false);
    nc = numel(R0.x);
    data.priors = struct('L', eye(nc), 'x0', zeros(nc, 1), 'sigma', 0.01, 'name', 'damping', 'estimate', "scale");
    R = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, ...
                     MaxPasses=8, Noise=struct('Grid', 16), Seed=3, Verbose=false));
    tc.verifyTrue(isfield(R.noise, 'priors') && numel(R.noise.priors) == 1);
    g = R.noise.priors(1);
    tc.verifyGreaterThan(g.sigma, 0.01, 'the deliberately tight prior was loosened by the data');
    recs = R.noise.priorHistory;
    tc.verifyGreaterThanOrEqual(numel(recs), 2, 'the group was updated on the passes after the first');
    last = recs{end};
    % tol: 15 % -- the fixed-point ratio Omega/r after a handful of passes on a
    % 3-day demo (each pass is one multiplicative step; measured 1.00-1.10)
    tc.verifyEqual(last.omega / last.r, 1, 'RelTol', 0.15, 'Helmert condition: Omega = r at the fixed point');
    tc.verifyTrue(last.updated);
    % estimate "none": the sigma stays where it was declared
    data.priors.estimate = "none";
    Rn = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, ...
                      MaxPasses=2, Noise=struct('Grid', 16), Seed=3, Verbose=false));
    tc.verifyEqual(Rn.noise.priors(1).sigma, 0.01);
    % the policy is validated by name
    data.priors.estimate = "law";
    tc.verifyError(@() lsqvce.fit(data, Verbose=false), 'lsqvce:fit:badPrior');
end

%% ------------------------------------------------ compare (M7 7c): ftest onto Results
function testCompareIsTheNestedModelFTest(tc)
% categories: contract, reference
    % R5 M7 7c: R.compare(Rred) re-homes lsqvce.engine.ftest. Reference: the
    % F statistic in closed form from the two residual sums; a needed
    % parameter is kept (reject H0), an unneeded one is not.
    rs = RandStream('mt19937ar', 'Seed', 17);
    n = 80;  t = linspace(0, 1, n).';
    A3 = [ones(n, 1), t, t.^2];  xt = [1; -2; 3];
    l = A3*xt + 0.05*randn(rs, n, 1);
    Rfull = lsqvce.fit(struct('A', A3, 'l', l), Verbose=false);
    Rred  = lsqvce.fit(struct('A', A3(:, 1:2), 'l', l), Verbose=false);   % drops the needed t^2
    F = Rfull.compare(Rred);
    Fref = ((Rred.sol.vtv - Rfull.sol.vtv) / 1) / (Rfull.sol.vtv / (n - 3));
    % tol: 1e-10 relative -- the same two residual sums in a closed form
    tc.verifyEqual(F.F, Fref, 'RelTol', 1e-10);
    tc.verifyTrue(F.significant, 'the quadratic term is needed');
    % an UNNEEDED extra parameter: full = cubic, reduced = the true quadratic
    A4 = [A3, t.^3];
    Rfull4 = lsqvce.fit(struct('A', A4, 'l', l), Verbose=false);
    F4 = Rfull4.compare(Rfull, Alpha=0.01);
    tc.verifyFalse(F4.significant, 'the cubic term is not needed at alpha 0.01');
    % nesting is checked by the engine (not a subset: refused by name)
    Rother = lsqvce.fit(struct('A', [ones(n, 1), sin(t)], 'l', l), Verbose=false);
    tc.verifyError(@() Rred.compare(Rother), 'lsqvce:ftest:notNested');
end

%% --------------------------------------- R5-D8 round 1: from + group label
function testPriorActivationAndGroupPattern(tc)
% categories: contract
    % R5-D8 round 1/6: a group is declared EITHER by a name pattern (matched
    % against the name column of commonParams; fit builds L, one row per
    % matched row per instance) OR by L (used as is); never both. `from`
    % activates the group at a pass; a single-pass run reports and skips a
    % from > 1 group.
    data = tc.TestData.data;
    R0 = lsqvce.fit(data, Verbose=false);  nc = numel(R0.x);
    names = data.commonParams(:, 2);
    tc.assumeTrue(all(cellfun(@(v) ischar(v) || isstring(v), names)), 'the fixture must carry names in column 2');
    pat = ['^' regexptranslate('escape', char(string(names{1}))) '$'];
    lay = lsqvce.layout(data);
    % the columns of A whose template row carries that name (names need not
    % be unique: every matching row, every instance)
    hit = cellfun(@(v) ~isempty(regexp(char(string(v)), pat, 'once')), names);
    colsHit = find(arrayfun(@(c) hit(c.row), lay.common));
    Lwant = zeros(numel(colsHit), nc);
    for q = 1:numel(colsHit), Lwant(q, colsHit(q)) = 1; end
    e = @(j) double((1:nc) == j);
    % (a MATLAB struct array needs one field set: both forms name all fields,
    % the unused one empty -- the toolbox reads "empty" as "not given")
    data.priors = [struct('group', pat, 'L', [], 'x0', 0, 'sigma', 0.5, 'name', '', 'estimate', "none", 'from', 1), ...
                   struct('group', '', 'L', e(3), 'x0', 0, 'sigma', 1.0, 'name', 'late', 'estimate', "none", 'from', 3)];
    R = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, ...
                     MaxPasses=4, Noise=struct('Grid', 16), Seed=3, Verbose=false));
    g = R.noise.priors;
    tc.verifyEqual(numel(g), 2);
    tc.verifyEqual(g(1).name, pat, 'the pattern is the default name');
    tc.verifyEqual(size(g(1).L, 1), numel(colsHit), 'one row per matched column');
    % tol: exact -- the pattern group's L is the unit rows on the matched columns
    tc.verifyEqual(g(1).L, Lwant);
    tc.verifyEqual(g(2).from, 3);
    % the same declaration with an explicit L gives the same group
    data.priors(1) = struct('group', '', 'L', Lwant, 'x0', 0, 'sigma', 0.5, 'name', pat, 'estimate', "none", 'from', 1);
    R2 = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, ...
                     MaxPasses=4, Noise=struct('Grid', 16), Seed=3, Verbose=false));
    % tol: exact -- identical declarations
    tc.verifyEqual(R2.x, R.x);
    % refusals by name: both forms, neither form, a pattern that matches nothing
    data.priors = struct('group', pat, 'L', Lwant, 'sigma', 1);
    tc.verifyError(@() lsqvce.fit(data, Verbose=false), 'lsqvce:fit:badPrior');
    data.priors = struct('sigma', 1);
    tc.verifyError(@() lsqvce.fit(data, Verbose=false), 'lsqvce:fit:badPrior');
    data.priors = struct('group', '^no_such_name$', 'sigma', 1);
    tc.verifyError(@() lsqvce.fit(data, Verbose=false), 'lsqvce:fit:badPrior');
    % a single solve with a from > 1 group reports and skips it (provoked on
    % purpose here: gate rule 7)
    muteWarn(tc, 'lsqvce:fit:priorNotActive');
    data.priors = struct('L', e(3), 'sigma', 1, 'from', 3);
    lastwarn('');
    lsqvce.fit(data, Verbose=false);
    [~, wid] = lastwarn;
    tc.verifyEqual(wid, 'lsqvce:fit:priorNotActive');
end

%% ------------------------------------------------- R5-D8 round 4: the write-back
function testToTableWritesEstimatesIntoTheExpandedTables(tc)
% categories: contract, precision
    seedFor(tc, 511);   % R7-45: drew from the global stream unseeded (suite-order dependent)
    % R.toTable(data): the expanded n x W cells (one row per instantiated
    % parameter) with the estimate and its standard deviation in the named
    % columns, every other column copied from the template row, and a map
    % naming the instance (row, day, arc) of every line.
    data = tc.TestData.data;
    R = lsqvce.fit(data, Verbose=false);
    [C, L, map] = R.toTable(data);
    nc = numel(R.x);
    tc.verifyEqual(size(C, 1), nc, 'one common row per column of A');
    tc.verifyEqual(size(C, 2), max(size(data.commonParams, 2), 10));
    % tol: exact -- the same doubles copied
    tc.verifyEqual(cell2mat(C(:, 9)), R.x);
    tc.verifyEqual(cell2mat(C(:, 10)), sqrt(diag(R.sol.Cxx)));
    tc.verifyEqual(C(:, 1), data.commonParams(cellfun(@(m) m.row, num2cell(map.common)), 1), 'the token column is the template row''s');
    tc.verifyEqual(numel(L), numel(data.systems));
    res = R.residuals('Verbose', false);
    for k = 1:numel(L)
        tc.verifyEqual(size(L{k}, 1), numel(res.sys(k).y), 'one local row per local unknown');
        % tol: exact
        tc.verifyEqual(cell2mat(L{k}(:, 9)), res.sys(k).y(:));
        tc.verifyEqual(cell2mat(L{k}(:, 10)), sqrt(diag(res.sys(k).Cyy)));
        m = map.local{k};
        tc.verifyTrue(all(ismember([m.row], 1:size(data.systems(k).localParams, 1))));
        dayRows = m(strcmp(data.systems(k).localParams([m.row], 1), 'day'));
        tc.verifyTrue(all(ismember([dayRows.day], data.days)), 'day instances carry their day');
    end
    % a Results from another problem: refused by name (the declaration must
    % be the one the Results came from)
    Rb = lsqvce.fit(struct('A', randn(10, 2), 'l', randn(10, 1)), Verbose=false);
    tc.verifyError(@() Rb.toTable(data), 'lsqvce:Results:tableMismatch');
end

%% ------------------------------ R6 WP-90: the per-epoch sigma reaches the estimator
function testPerEpochSigmaIsUndoneByTheEstimator(tc)
% categories: reference, precision
    % A scalar stream whose noise is WHITE with a per-epoch sigma spanning
    % two orders of magnitude (the error information the preprocessing
    % knows), returned as D.Cep; declared: level 1, temporal law white.
    % Estimate="noise" must find the level ~ 1 (the descale reaches the
    % joint sweep for d = 1 -- the v5.2 branch ran only for nAxis > 1); the
    % SAME data declared WITHOUT Cep (the v5.2 situation, E-03) cannot: its
    % level is the RMS of the sigmas, far from 1.
    nc = 4;  nE = 300;  nArc = 2;  nl = 2;
    xTrue = (1:nc).';
    data.commonParams = [repmat({'common'}, nc, 1), num2cell((1:nc).')];
    data.days = 1:3;
    s.buildDay = @(day, ot) perEpochDay(day, nc, nl, nE, nArc, xTrue, true);
    s.covariance = struct('rr', struct('type', 'diag', 'sigma', 1));
    s.arcs = nArc;
    s.localParams = {'system', 'b1'; 'system', 'b2'};
    data.systems = s;
    opts = lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, MaxPasses=3, ...
                          Noise=struct('Grid', 16), Seed=5, Verbose=false);
    R = lsqvce.fit(data, opts);
    tc.verifyEqual(R.noise.stream{1}, 'x1|rr');
    % (6 arcs x 300 epochs, r ~ 1790: the Helmert sigma scatter is ~ 1.7 %;
    % MEASURED 1.0047 on this fixture 2026-09-06, R2026a; without Cep 3.80)
    % tol: 5 % = three times that scatter, the shape/level split inside it
    tc.verifyEqual(R.noise.model{1}.sigma, 1, 'RelTol', 0.05, 'the per-epoch sigma is undone');
    z = (R.x - xTrue) ./ (sqrt(diag(R.sol.Cxx)) * R.sigma0);
    tc.verifyLessThan(max(abs(z)), 5, 'the field is recovered within its formal errors');
    dataNo = data;
    dataNo.systems.buildDay = @(day, ot) perEpochDay(day, nc, nl, nE, nArc, xTrue, false);
    R0 = lsqvce.fit(dataNo, opts);
    % the RMS of a log-uniform sigma on [0.1, 10] is ~ 3.3: without the
    % per-epoch information the level cannot be 1 (a bar, not a value)
    tc.verifyGreaterThan(R0.noise.model{1}.sigma, 2, 'without Cep the level is the RMS of the sigmas');
end

function D = perEpochDay(day, nc, nl, nE, nArc, xTrue, withCep)
%PEREPOCHDAY  A deterministic day: white noise with a per-epoch sigma 0.1 .. 10.
    rs = RandStream('mt19937ar', 'Seed', 1000 + day);
    m = nE * nArc;
    cep = 0.1 * 100 .^ (mod((1:m).' * 7 + day, 13) / 12);      % 0.1 .. 10, deterministic
    A = randn(rs, m, nc);  A(:, 1) = 1;
    B = zeros(m, nl);
    for a = 1:nArc, B((a-1)*nE + (1:nE), :) = 0.3 * randn(rs, nE, nl); end
    y = 0.1 * (1:nl).';
    D.A = A;  D.B = B;
    D.l = A * xTrue + B * y + cep .* randn(rs, m, 1);         % white noise, per-epoch sigma
    D.t = (0:m-1).' * 5;
    D.arcIdx = arrayfun(@(a) (a-1)*nE + (1:nE), 1:nArc, 'UniformOutput', false);
    if withCep, D.Cep = cep; end
    D = arcCells(D);                                 % v6: the cells are the arcs
end

function testDemoNoisePolicyIsTheWorkedPolicyExample(tc)
% categories: contract, reference
    % R6 WP-96: examples/demo_noisePolicy.m is the ONE worked example of the
    % per-stream policy fields and of a 'psd' law. Its claims are a reference
    % against a planted truth, so they are asserted here, not merely printed:
    % "none" must not move the declared sigma, "scale" must find the planted
    % level of a given psd shape, "law" must find the planted AR(1) colour,
    % and arcFactors must isolate the one degraded arc -- all in ONE run since
    % WP-98 fixed F-01, which is itself part of what this test now pins.
    out = evalc('[R, s] = demo_noisePolicy();');
    tc.verifyClass(R, 'lsqvce.Results');
    tc.verifyTrue(contains(out, 'arc factors'), 'the arc-factor section is printed');
    % R6 WP-98: the demo used to need TWO runs because a "law" stream and an
    % arcFactors stream could not settle together (F-01). Now all four policies
    % run at once, so the level for mag comes from the SAME Results as the rest
    tc.verifyNumElements(s.passes, 1, 'one run, one pass count');
    tc.verifyEqual(numel(s.stream), 4, 'all four streams in one run');

    tc.verifyTrue(s.noneFixed, 'estimate "none" leaves the declared sigma exactly alone');
    % tol: 10 % relative -- one level per stream from 9 arcs of 256 epochs has
    % a relative standard error of about 0.5*sqrt(2/9) = 0.24; 10 % is well
    % inside that and still catches a policy that did nothing (a factor 2 here)
    tc.verifyEqual(s.level(2), s.truth(2), 'RelTol', 0.10, ...
                   '"scale" on the psd shape finds the planted level');
    % tol: the same 10 % relative band, for the same reason -- the level of a
    % stream whose shape was fitted too (measured 5.1 % high on this fixture)
    tc.verifyEqual(s.level(3), s.truth(3), 'RelTol', 0.10, ...
                   '"law" finds the level of the colour it fitted');
    % tol: the fitted lag-1 correlation of an AR(1) from 9 arcs; +-0.15 keeps
    % white (0) and the planted 0.5 apart by three times the band
    tc.verifyEqual(s.rho1(1), s.rho1(2), 'AbsTol', 0.15, ...
                   '"law" recovers the planted AR(1) colour');
    % the arc factor is a VARIANCE ratio: planted 9, and it must be nowhere
    % near 1 (no split at all) -- a factor 2 band on a 9 is the honest bound
    tc.verifyGreaterThan(s.arcRatio(1), 4.5, 'the degraded arc is isolated');
    tc.verifyLessThan(s.arcRatio(1), 18, 'and not blown up');
    tc.verifyEqual(s.badEstimate, "lsqvce:fit:badEstimate", ...
                   'a policy that is not one of the three names is refused by name');
end

%% ------------------------------------------- arc factors (R6 WP-97, F-01)
function testArcFactorsRecoverThePlantedRatioInOnePass(tc)
% categories: reference, precision
    % R6 WP-97: the oracle dev/tests never had. test_ellmer pins the MECHANISM
    % of the arc factors ("none" leaves them alone, "scale" estimates them,
    % sigmaScale == sqrt(mean sigma2)); nothing asserted the NUMBER against a
    % planted truth, which is why finding F-01 survived every green run.
    % Planted here: arc 2 of every day is 3x worse in sigma, so the arc
    % VARIANCE factors are 1 / 9 / 1.
    R = lsqvce.fit(arcFactorFixture(3), ...
                   lsqvce.Options(Estimate="noise", Pool="stream", MaxPasses=1, ...
                                  Verbose=false));
    f2 = mean(reshape(R.noise.model{1}.sigma2, 3, []), 2).';
    ratio = f2(2) / mean(f2([1 3]));
    % tol: an EMPIRICAL band. Six seeds gave 7.83 8.49 8.84 9.37 9.62 9.85
    % against the planted 9 (measured 2026-09-08), so [6, 13] clears seed
    % noise by a wide margin while staying far from 1 (no split at all) and
    % from the compounded values the next test pins (> 3700).
    tc.verifyGreaterThan(ratio, 6, 'the degraded arc is found in ONE pass');
    tc.verifyLessThan(ratio, 13, 'and not exaggerated');
    % tol: a factor 2 -- the two undegraded arcs are one population; more than
    % that between them would be structure, not sampling
    tc.verifyLessThan(max(f2([1 3])) / min(f2([1 3])), 2, 'the two good arcs agree');
end

function testArcFactorsStaySettledAcrossPasses(tc)
% categories: reference, robustness
    % R6 WP-98. This test WAS the characterisation of finding F-01 -- it
    % asserted that the ratio exceeds 100 after four passes, and its comment
    % said INVERT ME when F-01 is fixed. This is that inversion.
    % The mechanism (HANDOVER 12.22): alpha2_m cancels its own arc's factor
    % exactly, because sweepArcs scales omegaT and sT by the same s2_m and
    % alpha2 is their ratio through two linear operators; nothing else applied
    % a per-arc factor to the data, so eq. 6.5.28 multiplied by a factor the
    % weighting never used and the update was a geometric recursion with no
    % fixed point -- MEASURED per-pass growth 8.48, 8.44, 8.38, 8.34. It now
    % multiplies by the factor the weighting DID use, and the ratio settles:
    %   passes 1 2 3 4 5  ->  9.37  8.56  8.56  8.56  8.56   (planted 9)
    r1 = arcRatio(lsqvce.fit(arcFactorFixture(3), ...
             lsqvce.Options(Estimate="noise", Pool="stream", MaxPasses=1, Verbose=false)));
    r4 = arcRatio(lsqvce.fit(arcFactorFixture(3), ...
             lsqvce.Options(Estimate="noise", Pool="stream", MaxPasses=4, Verbose=false)));
    % tol: the SAME band the one-pass oracle asserts. Four passes need not
    % improve the estimate, but they must not destroy it
    tc.verifyGreaterThan(r4, 6, 'four passes keep the planted ratio');
    tc.verifyLessThan(r4, 13, 'and do not inflate it');
    % tol: 25 % against the one-pass value -- passes 2..5 measured 8.562,
    % 8.562, 8.563, 8.562, a fixed point; the band leaves room for probe noise
    % and none for a trend (the defect grew by 8.4x per pass)
    tc.verifyEqual(r4, r1, 'RelTol', 0.25, 'the ratio SETTLES, it does not compound');
end

function testArcFactorGuardIsSilentOnASettledRun(tc)
% categories: contract, robustness
    % R6 WP-98: the guard WP-97 added stays in as a ratchet -- it is what would
    % catch this class of defect returning on a user's own data -- and this
    % test pins the half of it that can still be exercised: it must be SILENT
    % on the very fixture that used to make it fire.
    % ITS FIRING DIRECTION WAS DEMONSTRATED AGAINST THE PRE-FIX CODE and is on
    % the record with numbers (WP-97, HANDOVER 12.21): with a degraded arc it
    % fired 3/3 at spreads 4670 / 5600 / 6109, and stayed silent 3/3 without
    % one. It cannot be provoked now, which is the point of the fix; the
    % primary instrument against a regression is the settling test above.
    id = 'lsqvce:adjustVCEjoint:arcFactorsDiverging';
    opts = lsqvce.Options(Estimate="noise", Pool="stream", MaxPasses=4, Verbose=false);
    lsqvce.fit(arcFactorFixture(3), opts);          % fit resets the inventory
    tc.verifyFalse(any(strcmp(lsqvce.util.log.inventory(), id)), ...
                   'a settled run does not warn');
    lsqvce.fit(arcFactorFixture(1), opts);
    tc.verifyFalse(any(strcmp(lsqvce.util.log.inventory(), id)), ...
                   'nor does one with no degraded arc');
end

function q = arcRatio(R)
%ARCRATIO  The bad arc's mean variance factor against the good arcs' (arc 2).
    f2 = mean(reshape(R.noise.model{1}.sigma2, 3, []), 2).';
    q = f2(2) / mean(f2([1 3]));
end

function testArcFactorsReconstructThePlantedPerArcSigma(tc)
% categories: reference, precision
    % R6 WP-100: the COMPLETE statement, which the ratio test above does not
    % make. A stream's per-arc standard deviation is not its level and not its
    % factor: it is
    %       declared sigma  x  sqrt(arc factor)  x  sigma0
    % and the last term is where a least-squares adjustment always keeps the
    % overall scale. Reading the first two alone is what made HANDOVER 12.22
    % record a level of 0.678 "where the arcs were planted at 0.5" as an open
    % question; it was an incomplete reconstruction, not a defect. MEASURED
    % here on the same fixture: 0.517 / 1.499 / 0.507 against the planted
    % 0.5 / 1.5 / 0.5.
    R = lsqvce.fit(arcFactorFixture(3), ...
                   lsqvce.Options(Estimate="noise", Pool="stream", MaxPasses=4, ...
                                  Verbose=false));
    f2 = mean(reshape(R.noise.model{1}.sigma2, 3, []), 2).';
    implied = 0.5 * sqrt(f2) * R.sigma0;              % 0.5 is the fixture's declared sigma
    planted = [0.5 1.5 0.5];
    % A per-arc variance component from 3 x 256 epochs is worth about
    % 0.5*sqrt(2/3) = 41 % before pooling over days, and the INCOMPLETE
    % reconstruction (omitting sigma0) is 29 % low on this fixture.
    % tol: 15 % relative -- inside the sampling bound above, and far outside
    % the error the incomplete reading would show
    tc.verifyEqual(implied, planted, 'RelTol', 0.15, ...
        'declared x sqrt(factor) x sigma0 reconstructs the planted per-arc sigma');
end

function testPoolStreamGivesEachStreamItsOwnLaw(tc)
% categories: reference, robustness
    % R6 WP-100: what Pool="stream" promises, pinned -- and it is pinned
    % because a plausible "fix" would have destroyed it. sweepArcs weights
    % EVERY arc by one C while Cgroups holds a law per pool group, which reads
    % like a conflation; threading each group's own law into the weighting was
    % implemented and MEASURED, and it made the second group diverge every
    % time (planted 0.50 -> 1.331; planted 0.40 -> 1.365; sigma0 down to 0.71).
    % Weighting a group by its own freshly fitted law is self-referential; the
    % shared C is what breaks that loop. The shipped code recovers both
    % colours: 0.198 and 0.477 against 0.20 and 0.50, sigma0 0.997.
    R = lsqvce.fit(colourFixture(0.20, 0.50), ...
                   lsqvce.Options(Estimate="noise", Pool="stream", MaxPasses=16, ...
                                  Verbose=false));
    rho = @(k) localRho(R.noise.Cgroups(k));
    % tol: 0.08 absolute on a lag-1 correlation fitted from 9 arcs of 256
    % epochs -- wide enough for seed noise, narrow enough that the two planted
    % colours (0.20 and 0.50) cannot be confused with each other
    tc.verifyEqual(rho('x1|a'), 0.20, 'AbsTol', 0.08, 'stream a keeps its own colour');
    % tol: 0.08 absolute, the same bound and the same reason as above
    tc.verifyEqual(rho('x1|b'), 0.50, 'AbsTol', 0.08, 'stream b keeps its own colour');
    % tol: 10 % on sigma0 -- the whole point is that the weighting stays
    % consistent; the refuted change drove it to 0.71
    tc.verifyEqual(R.sigma0, 1, 'RelTol', 0.10, 'and the run stays consistently weighted');
end

%% ------------------------------- R6 WP-105: the robust pass inside the joint loop
function testRobustPassInsideTheJointLoopRecoversTheLawUnderOutliers(tc)
% categories: reference, robustness, contract
    % R6 WP-105: Robust together with Estimate "noise" runs -- one residual
    % sweep per pass, the factors inside every arc's whitener on the next
    % pass, so accumulation, probe redundancies and the component update see
    % ONE re-weighted problem. THE ORACLE: the colour fixture (two AR(1)
    % streams, planted 0.20 / 0.50, level 1) with 2 % of stream a's rows
    % replaced by 25-sigma blunders. Plain "noise" must be pulled off the
    % planted law by them; the robust pass must recover it and put the
    % blunders at ZERO weight -- which is why the oracle uses IGG3, a
    % REDESCENDING influence: Huber is monotone and keeps c/|u| of a
    % blunder's weight, so 2 % at 25 sigma retain 0.02*1.345*25 = 0.67 of
    % white variance and dilute the colour to ~0.12 (round 1, measured with
    % a global scale: 0.003). The manual says so.
    plantA = 0.20;  plantB = 0.50;
    data = colourFixture(plantA, plantB);
    inner = data.systems.buildDay;
    data.systems.buildDay = @(day, ot) blunders(inner(day, ot), day, ot, 0.02, 25);   % 2 %, 25 sigma
    base = {'Estimate', "noise", 'Pool', "stream", 'MaxPasses', 16, 'Verbose', false};
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');   % the PLAIN run under blunders may not settle -- that is its point
    Rp = lsqvce.fit(data, lsqvce.Options(base{:}));                      % plain
    trail = containers.Map('KeyType', 'double', 'ValueType', 'any');    % round 6: the pass trajectory
    watch = @(ctx) trailOf(trail, ctx);
    Rr = lsqvce.fit(data, lsqvce.Options(base{:}, 'Robust', "igg3", 'OnPass', watch));    % robust pass inside
    % THE CONTROL (round 7, the user's "measure"): the same blunders declared
    % as data of huge variance through Cep -- the toolbox composes Cep BEFORE
    % the temporal filter, i.e. the KOCH convention (weight, then filter) --
    % and NO robust pass. If the plain estimator recovers the colour here in
    % a few passes, the crawl above is the Yang convention's (weight AFTER the
    % filter: a zero-weight innovation neutralises the filter's whole tail).
    dk = data;
    dk.systems.buildDay = @(day, ot) declareBlunders(blunders(inner(day, ot), day, ot, 0.02, 25), day, ot, 0.02);
    trailK = containers.Map('KeyType', 'double', 'ValueType', 'any');
    Rk = lsqvce.fit(dk, lsqvce.Options(base{:}, 'OnPass', @(ctx) trailOf(trailK, ctx)));
    rho = @(R, k) R.noise.Cgroups(k);
    rhoK = rho(Rk, 'x1|a');  rhoK = rhoK(2) / rhoK(1);
    rhoP = rho(Rp, 'x1|a');  rhoR = rho(Rr, 'x1|a');
    rhoP = rhoP(2) / rhoP(1);  rhoR = rhoR(2) / rhoR(1);
    Ca = Rr.noise.Cgroups('x1|a');
    ia = find(strcmp(Rr.noise.stream, 'x1|a'), 1);
    % tol: 0.08 absolute on a lag-1 correlation from 9 arcs of 256 epochs --
    % the bound testPoolStreamGivesEachStreamItsOwnLaw uses on clean data
    tc.verifyEqual(rhoR, plantA, 'AbsTol', 0.08, sprintf(['the robust pass recovers stream a''s colour under ' ...
        'blunders (plain %.3f; robust C(1:4) = %s, level a %.3f, sigma0 %.3f, passes %d, converged %d; ' ...
        'per pass [rho_a nZero maxDev]: %s)'], ...
        rhoP, mat2str(Ca(1:4).', 3), Rr.noise.model{ia}.sigma, Rr.sigma0, Rr.history.nIter, Rr.converged, ...
        mat2str(cell2mat(values(trail, num2cell(sort(cell2mat(keys(trail)))))'), 3)));
    % tol: 0.08 absolute, the same clean-data bound -- the CONTROL's claim is
    % that Koch-style declaration of the blunders lets the plain estimator
    % recover the colour
    tc.verifyEqual(rhoK, plantA, 'AbsTol', 0.08, sprintf(['CONTROL: blunders declared through Cep (weight ' ...
        'BEFORE the filter), no robust pass -- passes %d, converged %d, sigma0 %.3f; per pass [rho_a nZero maxDev]: %s'], ...
        Rk.history.nIter, Rk.converged, Rk.sigma0, ...
        mat2str(cell2mat(values(trailK, num2cell(sort(cell2mat(keys(trailK)))))'), 3)));
    tc.verifyGreaterThan(abs(rhoP - plantA), abs(rhoR - plantA), ...
        'the plain estimator is further from the planted colour than the robust one');
    % the blunders, located exactly (the same RNG sequence the injector used):
    % every planted row must be at factor 0 (IGG3 rejects beyond k1 = 3
    % sigma; the standardisation is the stream's own MAD, so a 25-sigma
    % blunder is 25 sigma whatever level the whitener applied); the rest of
    % stream a keeps unit weight in the bulk, and the rows the law filter
    % contaminates from a blunder next door may go too -- but no more than
    % that (the 41 extra zeros of round 4 are that neighbourhood)
    f = Rr.noise.robust.factor;  k = Rr.noise.robust.keys;
    planted = false(size(f));  isA = false(size(f));  off = 0;
    for j = 1:size(k, 1)
        n = 256;                                                    % equal arcs
        if strcmp(k{j, 2}, 'a')
            isA(off + (1:n)) = true;
            planted(off + blunderRows(k{j, 3}, k{j, 4}, 0.02, n)) = true;
        end
        off = off + n;
    end
    tc.verifyEqual(sum(planted), 45, 'the injector planted round(0.02*256) = 5 rows in each of 9 arcs');
    tc.verifyEqual(f(planted), zeros(sum(planted), 1), 'every planted blunder is rejected');
    fa = f(isA & ~planted);
    tc.verifyEqual(median(fa), 1, 'the bulk keeps unit weight');
    % tol: the rows a blunder contaminates through the whitening filter --
    % at most its two neighbours -- plus the 0.27 % of N(0,1) beyond 3 sigma
    tc.verifyLessThanOrEqual(sum(fa == 0), 2 * 45 + 10, 'no rejection beyond the blunders'' neighbourhood');
    tc.verifyEqual(Rr.noise.robust.influence, 'igg3');
    tc.verifyEqual(Rr.noise.robust.streams, {'x1|a', 'x1|b'}, 'one MAD scale per stream');
    tc.verifyEqual(size(Rr.history.robustScale), [Rr.history.nIter, 2]);
    % THE LEVEL. WP-105 expected the retained-observation ratio E[w u^2]/E[w]
    % (0.73 IGG3, 0.80 Huber); round 4 read 1.007 on the CLEAN stream b. The
    % "weights cancel in alpha^2" explanation is disputed by R7-12 (s_0 holds
    % tr(Sigma_t^-1) of the unweighted law; a mirror reads 0.83-0.93), and
    % this 15 % band cannot tell 0.86 from 1 -- it pins "no gross bias",
    % not "unbiased". The discriminating multi-seed measurement is WP-115's
    % open tier-A item.
    ib = find(strcmp(Rr.noise.stream, 'x1|b'), 1);
    lvlB = Rr.noise.model{ib}.sigma * Rr.sigma0;
    % tol: 15 % relative on an effective level from 9 arcs of 256 epochs
    % (the clean-data bound of testRecoversPlantedRelativeWeights is 25 %)
    tc.verifyEqual(lvlB, 1, 'RelTol', 0.15, 'the clean stream''s level shows no gross bias under the robust pass');
    % Robust "none" carries no factors at all
    tc.verifyEmpty(Rp.noise.robust.factor);
end

function trailOf(trail, ctx)
%TRAILOF  Per-pass record for the diagnostic: stream a's lag-1 colour, the
%   zero-weight count, the convergence measure.
    r = NaN;
    if isKey(ctx.Cgroups, 'x1|a'), C = ctx.Cgroups('x1|a');  r = C(2) / C(1); end
    nz = NaN;  if ~isempty(ctx.factor), nz = sum(ctx.factor == 0); end
    trail(ctx.pass) = [r, nz, ctx.maxDev];   %#ok<NASGU>
end

function D = blunders(D, day, obsType, fraction, nSigma)
%BLUNDERS  Replace FRACTION of stream a's rows by NSIGMA-sigma blunders (v6 cells).
    if ~strcmp(obsType, 'a'), return, end
    for a = 1:numel(D.l)
        [idx, sgn] = blunderRows(day, a, fraction, numel(D.l{a}));
        D.l{a}(idx) = D.l{a}(idx) + nSigma * sgn;
    end
end

function D = declareBlunders(D, day, obsType, fraction)
%DECLAREBLUNDERS  The control's declaration: the planted rows of stream a get
%   a per-epoch sigma of 1e6 through Cep (composed BEFORE the temporal law:
%   the Koch convention), every other row sigma 1.
    if ~strcmp(obsType, 'a'), return, end
    D.Cep = cell(size(D.l));
    for a = 1:numel(D.l)
        c = ones(numel(D.l{a}), 1);
        c(blunderRows(day, a, fraction, numel(D.l{a}))) = 1e6;
        D.Cep{a} = c;
    end
end

function [idx, sgn] = blunderRows(day, arc, fraction, m)
%BLUNDERROWS  The rows the injector hits in (day, arc), and their signs: ONE
%   RNG sequence per day, consumed arc by arc, so the test can locate every
%   planted row without the injector telling it.
    rs = RandStream('mt19937ar', 'Seed', 9900 + day);
    for a = 1:arc
        idx = randperm(rs, m, round(fraction * m));
        sgn = 2 * (rand(rs, numel(idx), 1) > 0.5) - 1;
    end
    idx = idx(:).';
end

function r = localRho(C)
%LOCALRHO  Lag-1 autocorrelation of a fitted lag covariance.
    r = C(2) / C(1);
end

function data = colourFixture(phiA, phiB)
%COLOURFIXTURE  TWO streams, each AR(1) with its OWN phi, both "law", no
%   locals, so Pool="stream" has two distinct laws to find and nothing else
%   can absorb the difference. Deterministic: seeded on (day, stream).
    cfg = struct('nc', 4, 'days', 1:3, 'arcs', 3, 'mArc', 256, 'dt', 5, ...
                 'xTrue', [0.8 -0.45 0.3 0.15]);
    cfg.phi = struct('a', phiA, 'b', phiB);
    data = struct();
    data.commonParams = [repmat({'common'}, cfg.nc, 1), ...
                         arrayfun(@(j) sprintf('c%d', j), (1:cfg.nc).', ...
                                  'UniformOutput', false)];
    data.days = cfg.days;
    spec = @(nm) struct('type', 'diag', 'sigma', 1, 'estimate', "law", 'arcFactors', false);
    data.systems = struct( ...
        'buildDay', @(day, ot) colourDay(cfg, day, ot), ...
        'covariance', struct('a', spec('a'), 'b', spec('b')), ...
        'arcs', cfg.arcs, 'id', 1, 'name', 'colour');
end

function D = colourDay(cfg, day, obsType)
%COLOURDAY  One (day, stream) of the fixture above.
    m = cfg.arcs * cfg.mArc;
    t = (day - 1)*86400 + (0:m-1).' * cfg.dt;
    w = 2*pi / (0.5 * 86400);
    D.A = [sin(w*t), cos(w*t), sin(2*w*t), cos(2*w*t)];
    D.B = zeros(m, 0);
    rs = RandStream('mt19937ar', 'Seed', 7700 + 10*day + double(obsType(1)));
    e  = filter(1, [1 -cfg.phi.(obsType)], randn(rs, m + 2000, 1));
    D.l = D.A * cfg.xTrue(:) + e(2001:end);
    D.t = t;
    D.arcIdx = arrayfun(@(a) (a-1)*cfg.mArc + (1:cfg.mArc), 1:cfg.arcs, ...
                        'UniformOutput', false);
    D = arcCells(D);                                 % v6: the cells are the arcs
end

function data = arcFactorFixture(badScale)
%ARCFACTORFIXTURE  ONE white stream, three arcs a day, arc 2 badScale x worse
%   in SIGMA -- so the planted arc VARIANCE factors are [1 badScale^2 1] and
%   the stream's level is 0.5. The law is estimated ("law") so the loop keeps
%   iterating and the arc factors have somewhere to drift; no local unknowns,
%   so nothing else can absorb the difference. Deterministic: the noise comes
%   from a stream seeded on the day.
    cfg = struct('nc', 4, 'days', 1:3, 'arcs', 3, 'mArc', 256, 'dt', 5, ...
                 'badArc', 2, 'badScale', badScale, 'sigma', 0.5, ...
                 'xTrue', [0.8 -0.45 0.3 0.15]);
    data = struct();
    data.commonParams = [repmat({'common'}, cfg.nc, 1), ...
                         arrayfun(@(j) sprintf('c%d', j), (1:cfg.nc).', ...
                                  'UniformOutput', false)];
    data.days = cfg.days;
    data.systems = struct( ...
        'buildDay', @(day, ot) arcFactorDay(cfg, day), ...
        'covariance', struct('a', struct('type', 'diag', 'sigma', cfg.sigma, ...
                                         'estimate', "law", 'arcFactors', true)), ...
        'arcs', cfg.arcs, 'id', 1, 'name', 'arcFactorFixture');
end

function D = arcFactorDay(cfg, day)
%ARCFACTORDAY  One day of the fixture above.
    m = cfg.arcs * cfg.mArc;
    t = (day - 1)*86400 + (0:m-1).' * cfg.dt;
    w = 2*pi / (0.5 * 86400);
    D.A = [sin(w*t), cos(w*t), sin(2*w*t), cos(2*w*t)];
    D.B = zeros(m, 0);
    rs = RandStream('mt19937ar', 'Seed', 9100 + day);
    e  = cfg.sigma * randn(rs, m, 1);
    rows = (cfg.badArc - 1)*cfg.mArc + (1:cfg.mArc);
    e(rows) = e(rows) * cfg.badScale;
    D.l = D.A * cfg.xTrue(:) + e;
    D.t = t;
    D.arcIdx = arrayfun(@(a) (a-1)*cfg.mArc + (1:cfg.mArc), 1:cfg.arcs, ...
                        'UniformOutput', false);
    D = arcCells(D);                                 % v6: the cells are the arcs
end
