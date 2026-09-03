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
    tc.verifyTrue(F.reject, 'the quadratic term is needed');
    % an UNNEEDED extra parameter: full = cubic, reduced = the true quadratic
    A4 = [A3, t.^3];
    Rfull4 = lsqvce.fit(struct('A', A4, 'l', l), Verbose=false);
    F4 = Rfull4.compare(Rfull, Alpha=0.01);
    tc.verifyFalse(F4.reject, 'the cubic term is not needed at alpha 0.01');
    % nesting is checked by the engine (not a subset: refused by name)
    Rother = lsqvce.fit(struct('A', [ones(n, 1), sin(t)], 'l', l), Verbose=false);
    tc.verifyError(@() Rred.compare(Rother), 'lsqvce:ftest:notNested');
end
