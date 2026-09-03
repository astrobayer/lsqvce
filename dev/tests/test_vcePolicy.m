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
    tc.TestData.Rv = lsqvce.fit(data, lsqvce.Options(Model="vce", MaxOuter=3, NumLags=16, Seed=3, Verbose=false));
    tc.TestData.Rp = lsqvce.fit(data, lsqvce.Options(Model="pooled", MaxOuter=3, NumLags=16, Seed=3, Verbose=false));
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
    R = lsqvce.fit(data, lsqvce.Options(Model="ellmer", Pool="obsType", JointMaxIter=2, NumLags=16, Seed=3, Verbose=false));
    types = unique(cellfun(@(s) extractAfter(s, '|'), R.history.streams, 'UniformOutput', false));
    tc.verifyEqual(sort(keys(R.noise.Cgroups)), sort(types));
end

%% ------------------------------------------------------------- messaging
function testVerboseZeroIsSilent(tc)
% categories: contract, robustness | tags: diagnostics
    data = tc.TestData.data;
    out = evalc('lsqvce.fit(data, lsqvce.Options(Model="vce", MaxOuter=1, NumLags=16, Seed=3, Verbose=0));');
    tc.verifyEmpty(strtrim(out), 'Verbose 0 prints nothing');
end

function testVerboseOnePrintsThePassRecord(tc)
% categories: contract | tags: diagnostics
    data = tc.TestData.data;
    out = evalc('lsqvce.fit(data, lsqvce.Options(Model="vce", MaxOuter=1, NumLags=16, Seed=3, Verbose=1));');
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
    R = lsqvce.fit(data, lsqvce.Options(Model="vce", MaxOuter=6, Seed=1, Verbose=false));
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
    lsqvce.fit(data, lsqvce.Options(Model="vce", MaxOuter=n, NumLags=16, Seed=3, Verbose=false, CacheDays=1));
    v1 = cell2mat(cellfun(@(c) cell2mat(values(c)), cnts, 'UniformOutput', false));
    tc.verifyEqual(unique(v1), 2*n + 1, 'CacheDays 1: accumulate + sweep per pass, + the final diagnostics');
    for c = cnts, remove(c{1}, keys(c{1})); end
    lsqvce.fit(data, lsqvce.Options(Model="vce", MaxOuter=n, NumLags=16, Seed=3, Verbose=false, CacheDays=Inf));
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
