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
%   Updated 2026-10-17 (CEST, UTC+02:00), WP-144: plot("psd") on a log axis
%   for every pool group (A3); the prior's "scale" update survives a refused
%   v'Pv (A6). With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-16 (CEST, UTC+02:00), WP-147: the bins' stop rule -- a
%   low-information bin no longer blocks the stop, a tighter Tolerance
%   tightens every bin, the not-converged warning says what still moves.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-16 (CEST, UTC+02:00), WP-145: demo_noisePolicy must
%   converge; the captures echo their warnings. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-14 (CEST, UTC+02:00), WP-135 (third review, lane 4): the
%   converged column reads R.converged; a prior from the last pass is scaled
%   in the closing solve. With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-134 (third review): the common prior group is in the stop rule.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-131: a late prior beyond the last
%   pass is reported; a non-robust noise run is not called robust. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
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
    % the POOLED run (Pool "all", the default until WP-157) puts the demo's four
    % streams of different physics in one group, which WP-154's check reports
    % -- off for this run only (a class-wide mute silenced
    % testPoolGroupOfDifferentSpectraWarns too)
    wp = warning('off', 'lsqvce:adjustVCEjoint:poolNotWhite');
    tc.TestData.Rp = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="all", ArcFactors=false, MaxPasses=3, Noise=struct('Grid', 16), Seed=3, Verbose=false));
    warning(wp);
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

function testConvergenceTableReadsTheConvergedFlag(tc)
% categories: contract | tags: diagnostics
    % WP-135 (review 3, lane 4): the converged column read h.converged, which
    % no producer writes -- false on every row of a converged run. With every
    % stream at estimate "none" the loop settles on pass 1.
    data = tc.TestData.data;
    for k = 1:numel(data.systems)
        for f = fieldnames(data.systems(k).covariance).'
            data.systems(k).covariance.(f{1}).estimate = "none";
        end
    end
    R = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, MaxPasses=4, ...
                     Noise=struct('Grid', 16), Seed=3, Verbose=false));
    tc.assumeTrue(R.converged, 'a model with nothing to estimate settles');
    T = R.convergence();
    tc.verifyTrue(T.converged(end), 'the last row of a converged run');
    tc.verifyFalse(any(T.converged(1:end-1)), 'only the last row can be true');
    tc.verifyEqual(tc.TestData.Rv.convergence().converged(end), tc.TestData.Rv.converged);
end

function testPriorFromTheLastPassIsScaledInTheClosingSolve(tc)
% categories: contract | tags: diagnostics
    % WP-135 (review 3, lane 4): the closing solve repeats the last pass's
    % index and updated only groups with from < pass, so a "scale" group that
    % entered on the last pass came back with its declared sigma. MaxPasses = 2
    % and from = 2: the group is in pass 2 (pending, so no convergence there)
    % and must take its first update in the closing solve (model.final).
    data = tc.TestData.data;
    nc = numel(tc.TestData.Rv.x);
    data.priors = struct('L', eye(nc), 'x0', zeros(nc, 1), 'sigma', 0.01, 'from', 2, ...
                         'name', 'late', 'estimate', "scale");
    R = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, MaxPasses=2, ...
                     Noise=struct('Grid', 16), Seed=3, Verbose=false));
    tc.verifyEqual(R.history.nIter, 2);
    tc.verifyNotEqual(R.noise.priors(1).sigma, 0.01, 'the late group is scaled in the closing solve');
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
    % WP-131: the joint loop records influence "none" for a non-robust run;
    % the report called every such run "a robust (none) fit"
    tc.verifyEmpty(codesOf(Vu, 'ROBUST_'), 'a non-robust noise run is not called robust');
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
    echoWarnings(out);   % WP-145: its warnings reach the gate
    tc.verifyEmpty(strtrim(out), 'Verbose 0 prints nothing');
end

function testVerboseOnePrintsThePassRecord(tc)
% categories: contract | tags: diagnostics
    data = tc.TestData.data;
    o = lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, MaxPasses=1, Noise=struct('Grid', 16), Seed=3, Verbose=1);
    out = evalc('lsqvce.fit(data, o);');
    echoWarnings(out);   % WP-145: its warnings reach the gate
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
    % passes 2n+2 times: + the accumulation under the RETURNED model (R7-11,
    % WP-116: the solution used to be one pass behind it) + diagnose's ONE
    % sweep, which RE-READS the data since WP-86 part 1 instead of replaying
    % stored rows (the price of not keeping n x u_full doubles). Since WP-121
    % the robust check is formed in that sweep (2n+3 while it swept on its
    % own; only Test "mad" still sweeps first). The
    % diagnostic source's own enumeration pass (part 1: 2n+4) is gone since
    % part 2: its index comes from the loaders' memoised meta answers and its
    % row count from the accumulator. With CacheDays = Inf the calendar stays
    % resident and every day is built ONCE per run. Counted through a wrapped
    % buildDay on the demo declaration (2 systems, 3 days).
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
    tc.verifyEqual(unique(v1), 2*n + 2, ['CacheDays 1: accumulate + sweep per pass, + the ' ...
                   'solve under the returned model (R7-11), + diagnose''s one re-reading ' ...
                   'sweep (WP-86 part 1, WP-121); the source reads no data itself (part 2)']);
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
function testPsdPlotIsLogarithmicForEveryPoolGroup(tc)
% categories: contract | tags: diagnostics
    % WP-144 (A3): with several pool groups the spectra were drawn after
    % hold on, and semilogy under hold keeps the axes' scale -- they came out
    % on a LINEAR axis (reproduced: 4 groups, YScale "linear"). Spectra span
    % decades; the axis must be logarithmic whatever the number of groups.
    R = tc.TestData.Rv;                              % Pool "stream": one group per stream
    tc.assertGreaterThan(R.noise.Sgroups.Count, 1, 'the fixture has several pool groups');
    fh = R.plot("psd", Visible=false);
    tc.addTeardown(@() close(fh));
    ax = findall(fh, 'Type', 'axes');
    tc.verifyEqual(get(ax(1), 'YScale'), 'log');
    tc.verifyNumElements(findall(ax(1), 'Type', 'line'), R.noise.Sgroups.Count, 'one line per group');
end

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
function testCompareRefusesEstimatedWeights(tc)
% categories: contract, robustness
    % WP-150 (F1): the nested F test assumes ONE fixed weighting. After
    % Estimate "noise" each model's v'Pv is normalised to its own redundancy,
    % so F = 1 whatever the data -- measured: F = 1.000, p 0.32 ("not
    % significant") for a sine term the declared weights test at F = 2705.
    % A robust influence re-weights each model differently. Both are refused.
    rs = RandStream('mt19937ar', 'Seed', 150);
    n = 600;  t = (1:n).' / n;
    l = 1 + 2*t + 0.3*sin(2*pi*5*t) + 0.1*randn(rs, n, 1);
    Af = [ones(n, 1), t, sin(2*pi*5*t)];  Ar = Af(:, 1:2);
    mk = @(A) struct('days', 1, 'commonParams', {[repmat({'common'}, size(A, 2), 1), ...
                        arrayfun(@(j) sprintf('c%d', j), (1:size(A, 2)).', 'UniformOutput', false)]}, ...
        'systems', struct('buildDay', @(d, ot) struct('A', {{A}}, 'B', {{zeros(n, 0)}}, 'l', {{l}}, 't', {{t * 86400}}), ...
                          'covariance', struct('s', struct('type', 'diag', 'sigma', 1, 'estimate', "scale")), ...
                          'arcs', 1, 'id', 1, 'name', 's'));
    o = lsqvce.Options(Estimate="noise", Verbose=0);
    Rf = lsqvce.fit(mk(Af), o);  Rr = lsqvce.fit(mk(Ar), o);
    tc.verifyError(@() Rf.compare(Rr), 'lsqvce:Results:compareEstimatedWeights');
    Rh = lsqvce.fit(struct('A', Af, 'l', l), Robust="huber", Verbose=false);
    Rd = lsqvce.fit(struct('A', Ar, 'l', l), Verbose=false);
    tc.verifyError(@() Rh.compare(Rd), 'lsqvce:Results:compareEstimatedWeights');
    tc.verifyError(@() Rd.compare(Rh), 'lsqvce:Results:compareEstimatedWeights');
    % the same data under the declared weights: the term is overwhelmingly needed
    F = lsqvce.fit(struct('A', Af, 'l', l), Verbose=false).compare(Rd);
    tc.verifyTrue(F.significant && F.F > 100, sprintf('F = %.4g with declared weights', F.F));
end

function testCompareIsTheNestedModelFTest(tc)
% categories: contract, reference
    % R5 M7 7c: R.compare(Rred) re-homes lsqvce.engine.verdict.ftest. Reference: the
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

function testLatePriorBeyondTheLastPassIsReported(tc)
% categories: contract
    % WP-131: a prior group activates at pass `from`; a joint loop that ended
    % before it (MaxPasses < from, early convergence, a user stop) never
    % applied the group -- R.noise.priors still listed it with its sigma and
    % nothing was said (reviewer's mirror: x pinned to 0 +- 1e-6 came out
    % 0.4995). priorNotActive now reports it; a loop that reaches the pass
    % applies it.
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');   % two passes on purpose (gate rule 7)
    muteWarn(tc, 'lsqvce:fit:priorNotActive');            % provoked on purpose; read from the inventory
    data = tc.TestData.data;
    nc = numel(lsqvce.fit(data, Verbose=false).x);
    e3 = double((1:nc) == 3);
    data.priors = struct('L', e3, 'x0', 0, 'sigma', 1e-6, 'name', 'pin', 'from', 3);
    o = @(np) lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, ...
                             MaxPasses=np, Noise=struct('Grid', 16), Seed=3, Verbose=false);
    R2 = lsqvce.fit(data, o(2));
    tc.verifyTrue(any(strcmp(R2.meta.warnings, 'lsqvce:fit:priorNotActive')), ...
        'a group beyond the last pass is reported');
    R5 = lsqvce.fit(data, o(5));
    if R5.history.nIter >= 3
        % tol: 1e-4 -- the pinning prior (sigma 1e-6) against an unpinned x(3) of O(1)
        tc.verifyLessThan(abs(R5.x(3)), 1e-4, 'a loop that reaches pass 3 applies the group');
        tc.verifyFalse(any(strcmp(R5.meta.warnings, 'lsqvce:fit:priorNotActive')));
    else
        tc.verifyTrue(any(strcmp(R5.meta.warnings, 'lsqvce:fit:priorNotActive')), ...
            'a loop that converged before pass 3 says so');
    end
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
    echoWarnings(out);   % WP-145: its warnings reach the gate
    tc.verifyClass(R, 'lsqvce.Results');
    tc.verifyTrue(contains(out, 'arc factors'), 'the arc-factor section is printed');
    % R6 WP-98: the demo used to need TWO runs because a "law" stream and an
    % arcFactors stream could not settle together (F-01). Now all four policies
    % run at once, so the level for mag comes from the SAME Results as the rest
    tc.verifyNumElements(s.passes, 1, 'one run, one pass count');
    tc.verifyEqual(numel(s.stream), 4, 'all four streams in one run');
    % WP-145: the figures below are a fixed point's only if the loop reached
    % one -- it used to stop on its budget (12 passes; the plain law needed 14)
    tc.verifyTrue(R.converged, sprintf('the shipped run must converge (stopReason %s)', R.stopReason));

    tc.verifyTrue(s.noneFixed, 'estimate "none" leaves the declared sigma exactly alone');
    % tol: 10 % relative -- one level per stream from 9 arcs of 256 epochs has
    % a relative standard error of about 0.5*sqrt(2/9) = 0.24; 10 % is well
    % inside that and still catches a policy that did nothing (a factor 2 here)
    tc.verifyEqual(s.level(2), s.truth(2), 'RelTol', 0.10, ...
                   '"scale" on the psd shape finds the planted level');
    % tol: the same 10 % relative band, for the same reason -- the level of a
    % stream whose shape was fitted too (measured 3.5 % high on this fixture)
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
    % here on the same fixture before 6.1.6: 0.517 / 1.499 / 0.507 against the
    % planted 0.5 / 1.5 / 0.5. Since 6.1.6 (R7-13) the factors are absolute
    % and weight the solve, so sigma0 carries less of it; the band is unchanged.
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

function testScaleOnlyRunIteratesAndRecoversTheLevels(tc)
% categories: reference, robustness | tags: convergence
    % R7-11: with no "law" stream the loop's only criterion, max|alpha2 - 1|,
    % was 0 by construction, so the run stopped "converged" on pass 1 and
    % returned the solve under the DECLARED weights (demo_lsqvce case 5:
    % sigma0 2.31 printed under "variance ratio back at ~1"). The stop rule
    % now needs the levels to settle too, and the returned solution is the
    % one under the returned model. Stream b is truly 3x its declared sigma.
    R = lsqvce.fit(scaleOnlyFixture(1, 3, false), lsqvce.Options(Estimate="noise", Pool="stream", ...
                   MaxPasses=8, Verbose=false));
    tc.verifyGreaterThan(R.history.nIter, 1, 'a level-only run does not stop on pass 1');
    tc.verifyTrue(R.converged);
    lvl = @(s) R.noise.model{find(strcmp(R.noise.stream, s), 1)}.sigma * R.sigma0;
    % tol: 10 % relative -- a level from 3 days x 3 arcs x 256 white epochs
    % (sampling ~3 %); the pass-1 answer this replaces was off by ~2x
    tc.verifyEqual([lvl('x1|a'), lvl('x1|b')], [1 3], 'RelTol', 0.10, ...
        'declared x sqrt(level) x sigma0 reconstructs both planted sigmas');
end

function testArcFactorStreamIgnoresOtherStreamsLevels(tc)
% categories: metamorphic, reference
    % R7-09: an arc-factor stream's factors were normalised by gamma-bar over
    % EVERY column, so stream a's reconstructed sigma moved with stream b's
    % truth (review mirror: 1.01 / 0.79 / 1.15 for b = 1 / 4 / 0.25). The
    % review's remedy -- gamma-bar per stream alone -- loses a's own level
    % (mirrors/wp116/gammabar_modes.py: a truly 3x came back 1.93); per-stream
    % gamma-bar times the stream's Helmert level recovers 3.02 for any b.
    opts = lsqvce.Options(Estimate="noise", Pool="stream", MaxPasses=8, Verbose=false);
    implied = @(R, s) sqrt(mean(R.noise.model{find(strcmp(R.noise.stream, s), 1)}.sigma2)) * R.sigma0;
    R1 = lsqvce.fit(scaleOnlyFixture(3, 1, true), opts);
    R4 = lsqvce.fit(scaleOnlyFixture(3, 4, true), opts);
    % tol: 10 % relative -- stream a's level from 9 arcs of 256 white epochs
    % (sampling ~3 %); the shipped normalisation was off by 20-35 % here
    tc.verifyEqual(implied(R1, 'x1|a'), 3, 'RelTol', 0.10, 'a recovers its own level (b = 1)');
    % tol: as above
    tc.verifyEqual(implied(R4, 'x1|a'), 3, 'RelTol', 0.10, 'and the same when only b changes (b = 4)');
    % tol: 5 % relative -- the metamorphic statement itself: the same a-data
    % and a 4x change of b must not move a's level beyond the shared sigma0's
    % own sampling
    tc.verifyEqual(implied(R4, 'x1|a'), implied(R1, 'x1|a'), 'RelTol', 0.05);
end

function testArcFactorsWeightTheSolve(tc)
% categories: contract | tags: convergence
    % R7-13 (WP-116 part 2, user decision 2026-10-10: Ellmer's convention). The
    % per-arc factors estimated on pass 1 must be what weights pass 2's arcs.
    % End to end, through every seam it crosses: adjustVCEjoint keys the factor
    % by the loader's key (sweepArcs, rowWeightKey), jointModelOpt emits it as
    % covOverride.byArc, makeLoader>batchOf hands day and arc to
    % applyCovOverride. Read back from pass 2's LIVE accumulator: the
    % covariance of arc 2 over that of arc 1 is the ratio of their pass-1
    % factors. Before 6.1.6 only the stream mean reached the loader, and this
    % ratio was 1 whatever the factors said.
    store = containers.Map('KeyType', 'double', 'ValueType', 'any');
    opts = lsqvce.Options(Estimate="noise", Pool="stream", MaxPasses=2, Verbose=false, ...
                          OnPass=@(ctx) stashPass(store, ctx));
    lsqvce.fit(arcFactorFixture(3), opts);
    tc.assertTrue(isKey(store, 1) && isKey(store, 2), 'both passes reached the hook');
    p1 = store(1);  acc = store(2).acc;
    ids = acc.systemIDs();
    st = acc.getStore(ids(1));
    c = zeros(1, 2);
    for a = 1:2
        uo = st.uo;  uo.obsType = 'a';  uo.day = 1;  uo.arc = a;
        arc = st.dataFcn(1, 'a', a, uo);
        S = arc.Cfcn(1:numel(arc.l), uo);
        c(a) = S(1);
    end
    f = @(a) p1.sigma2(strcmp(p1.arcKeys, sprintf('s%s_a_d1_a%d', string(ids(1)), a)));
    % tol: 1e-10 relative -- both covariances are one law scaled by a stored
    % factor; only rounding separates them. Before R7-13 the ratio was 1.
    tc.verifyEqual(c(2) / c(1), f(2) / f(1), 'RelTol', 1e-10, ...
        'pass 2 weights each arc with its own pass-1 factor');
    tc.verifyGreaterThan(c(2) / c(1), 4, 'and the degraded arc IS weighted down');
end

function go = stashPass(store, ctx)
%STASHPASS  OnPass hook: keep each pass's factors, keys and live accumulator.
    store(ctx.pass) = struct('sigma2', ctx.sigma2, 'arcKeys', {ctx.arcKeys}, 'acc', ctx.acc);
    go = true;
end

function testArcFactorsAreAbsoluteOnAHeterogeneousStream(tc)
% categories: reference, precision
    % WP-116 part 2: the 6.1.5 normalisation, (gamma2_m / gamma-bar_s^2) *
    % mean(gamma2_s), counted a stream's level twice as soon as its arcs
    % differ -- its homogeneous fixtures could not see it. Mirror
    % (mirrors/wp116/perarc_weighting.py, one arc 10x among nine, 16 seeds):
    % 1176 / 11.65 against the planted variances 100 / 1, and the OTHER
    % stream's level pushed off too (2.8 against 4). Since 6.1.6 a factor is
    % the arc's absolute gamma2_m. Planted here: stream a's arc 2 of every day
    % 10x in sigma, stream b 2x its declared sigma.
    R = lsqvce.fit(scaleOnlyFixture(1, 2, true, 10), lsqvce.Options(Estimate="noise", ...
                   Pool="stream", MaxPasses=8, Verbose=false));
    ka = find(strcmp(R.noise.stream, 'x1|a'), 1);
    f = reshape(R.noise.model{ka}.sigma2, 3, []) * R.sigma0^2;    % declared sigma 1
    % tol: 15 % relative -- a per-arc variance from 3 x 256 white epochs
    % (~5 % sampling); the double count put the good arcs at ~4x and the bad
    % at ~4x its truth on this layout
    tc.verifyEqual(mean(f(2, :)), 100, 'RelTol', 0.15, 'the degraded arcs: planted variance 100');
    % tol: as above
    tc.verifyEqual(mean(mean(f([1 3], :))), 1, 'RelTol', 0.15, 'the good arcs: planted variance 1');
    lvlB = R.noise.model{find(strcmp(R.noise.stream, 'x1|b'), 1)}.sigma * R.sigma0;
    % tol: 10 % relative -- stream b's level from 9 arcs of 256 epochs (~3 %);
    % the double count moved it by ~15 % in the mirror
    tc.verifyEqual(lvlB, 2, 'RelTol', 0.10, 'and the other stream keeps its own level');
end

function testLjungBoxIsSizedOnShortArcs(tc)
% categories: reference, robustness
    % R7-38: the per-stream Ljung-Box statistic pooled per-arc lag sums but
    % divided by (N - k), as if the stream were ONE series; every arc loses k
    % products, so on short arcs Q came out small and a white stream was
    % never flagged -- and never close to it (64 x 30 epochs: size 0.001,
    % E[Q] 13 instead of 20; mirrors/wp118/ljungbox_pooled.py). Pooled with
    % the pair counts N_k it is sized (0.048 there). Under a white null the
    % p-value is uniform: over 16 seeds its mean must sit near 0.5.
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');   % one pass on purpose: the test reads pass 1
    p = zeros(1, 16);
    for sd = 1:16
        R = lsqvce.fit(shortArcFixture(sd), lsqvce.Options(Estimate="noise", Pool="stream", ...
                       ArcFactors=false, MaxPasses=1, Verbose=false));
        p(sd) = R.history.whitePval(end, 1);
    end
    % tol: the mean of 16 U(0,1) draws has sd 0.072; [0.29, 0.71] is three
    % sd around 0.5. The mis-sized statistic put the p-values near 1 (E[Q]
    % 13 against 20 degrees of freedom gives p ~ 0.87 at the mean).
    tc.verifyGreaterThan(mean(p), 0.29, 'p-values are not pushed towards 0');
    tc.verifyLessThan(mean(p), 0.71, 'nor piled up towards 1 (the R7-38 symptom)');
end

%% ------------------------------- WP-144 (A6): the prior's scale after a refused v'Pv
function testPriorScaleSurvivesARefusedVtPv(tc)
% categories: robustness, reference
    % A tight fit (commons of size 1e5 against unit noise) makes the normal
    % backend refuse v'Pv as cancelled: s0sq and every entry of Cxx are NaN.
    % The prior's "scale" step took N^-1 = Cxx / s0sq, so r was NaN, the
    % group was never updated, PriorDev stayed 0 and the loop reported
    % "converged" with the declared sigma 1e6 (reproduced; the QR backend
    % found 4.9e4). N^-1 now comes from the cofactor, which is finite.
    muteWarn(tc, 'lsqvce:accum:cancellation');      % it fires on every solve here, by design
    data = tightPriorFixture();
    fitWith = @(backend) lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ...
        ArcFactors=false, Backend=backend, MaxPasses=6, Noise=struct('Grid', 16), Seed=3, Verbose=false));
    Rn = fitWith("normal");  Rq = fitWith("qr");
    upd = cellfun(@(r) r(1).updated, Rn.noise.priorHistory);
    tc.verifyTrue(all(upd), sprintf('updated on every pass (%s)', mat2str(upd)));
    tc.verifyTrue(all(cellfun(@(r) isfinite(r(1).r), Rn.noise.priorHistory)), 'its redundancy is finite');
    % tol: 1e-6 relative -- the same Helmert step on the same data, the two
    % backends' N^-1 agreeing far below the statistics (measured equal to
    % four printed digits, 4.886e4 both)
    tc.verifyEqual(Rn.noise.priors(1).sigma, Rq.noise.priors(1).sigma, 'RelTol', 1e-6, ...
        'the normal backend estimates the prior as the QR backend does');
    tc.verifyLessThan(Rn.noise.priors(1).sigma, 1e5, 'and moves it away from the declared 1e6');
end

%% ------------------------------- WP-149: the fit's own Results after a refused v'Pv
function testRefusedVtPvIsSupersededInResults(tc)
% categories: robustness, reference
    % The same tight fit: the normal backend refuses v'Pv on every solve, and
    % the Results fit returned carried R.sigma0 = NaN and NaN formal errors
    % although the warning says the residual pass supersedes it. The verdict
    % ran that pass; Results now takes vtv, sigma0 and Cxx from it. The QR
    % backend on the same data is the independent route, the residual pass
    % on R.sol the cross-check (it must no longer have anything to supersede).
    muteWarn(tc, 'lsqvce:accum:cancellation');      % fires on every solve here, by design
    muteWarn(tc, 'lsqvce:normalEq:cancellation');   % (the commons' and the backend's guard)
    muteWarn(tc, 'lsqvce:fit:priorScaleHeld');      % the fixture's "scale" prior, held by the Estimate "none" half
    data = tightPriorFixture();
    noiseNV = {'Pool', "stream", 'ArcFactors', false, 'MaxPasses', 6, ...
               'Noise', struct('Grid', 16), 'Seed', 3};
    for est = ["noise", "none"]
        nv = {'Estimate', est, 'Verbose', false};
        if est == "noise", nv = [nv, noiseNV]; end %#ok<AGROW>
        Rn = lsqvce.fit(data, lsqvce.Options(nv{:}, 'Backend', "normal"));
        Rq = lsqvce.fit(data, lsqvce.Options(nv{:}, 'Backend', "qr"));
        tc.verifyTrue(isfinite(Rn.sigma0) && all(isfinite(Rn.sol.Cxx(:))), ...
            sprintf('%s: sigma0 and Cxx are finite (%g)', est, Rn.sigma0));
        % tol: 1e-9 relative -- the direct residual sums of two backends at
        % the same optimum; measured <= 6e-13 (sigma0) and 6e-13 (formal
        % errors) on both estimates, with and without the prior
        tc.verifyEqual(Rn.sigma0, Rq.sigma0, 'RelTol', 1e-9, ...
            sprintf('%s: sigma0 as the QR backend finds it', est));
        % tol: 1e-9 relative, the same bound -- the formal errors are s0 sqrt(diag(Q))
        tc.verifyEqual(sqrt(diag(Rn.sol.Cxx)), sqrt(diag(Rq.sol.Cxx)), 'RelTol', 1e-9, ...
            sprintf('%s: the formal errors as QR''s', est));
        res = verifyWarningFree(tc, @() lsqvce.engine.residuals(Rn.acc, Rn.sol, ...
                  'Collect', false, 'Verbose', false));
        tc.verifyFalse(res.supersedesSol, sprintf('%s: nothing left to supersede', est));
        % tol: 1e-12 relative -- the same residuals summed in two passes;
        % measured <= 2.3e-13
        tc.verifyEqual(Rn.sigma0, res.sigma0, 'RelTol', 1e-12, ...
            sprintf('%s: R.sigma0 is the residual pass''s', est));
    end
end

function testRefusedVtPvKeepsTheRobustTruncation(tc)
% categories: robustness, reference
    % A robust joint run's hard solution divides s0^2 by the truncation tau
    % (adjustVCEjoint>hardSolution). It skipped a refused v'Pv entirely, so
    % sol.truncation was never set and the superseding residual pass divided
    % by 1: igg3 on the normal backend gave sigma0 1.3 % off QR's (measured,
    % tau 0.973). tau and dof are now carried with the NaN.
    muteWarn(tc, 'lsqvce:accum:cancellation');
    muteWarn(tc, 'lsqvce:normalEq:cancellation');
    muteWarn(tc, 'lsqvce:residuals:supersedesSol');  % the robust sweeps inside the loop
    data = tightPriorFixture();
    fitWith = @(backend) lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ...
        ArcFactors=false, Backend=backend, Robust="igg3", MaxPasses=6, Noise=struct('Grid', 16), ...
        Seed=3, Verbose=false));
    Rn = fitWith("normal");  Rq = fitWith("qr");
    tc.assumeLessThan(Rq.sol.truncation, 1, 'the fixture rejects epochs (tau < 1)');
    % tol: 1e-12 -- tau is a mean of closed-form truncation factors over the same
    % rejected epochs (the two backends' residuals agree to ~1e-12)
    tc.verifyEqual(Rn.sol.truncation, Rq.sol.truncation, 'AbsTol', 1e-12, 'the same tau');
    % tol: 1e-9 relative -- as above; measured 4.3e-13
    tc.verifyEqual(Rn.sigma0, Rq.sigma0, 'RelTol', 1e-9, 'sigma0 corrected for the truncation');
end

function data = tightPriorFixture()
%TIGHTPRIORFIXTURE  One white stream, 3 days x 4 arcs x 64 epochs, 4 commons of
%   size 1e5 against unit noise -- l'Pl / v'Pv ~ 1e10, past the backends'
%   8-digit cancellation guard -- and a loose common prior (sigma 1e6,
%   estimate "scale"). Deterministic (WP-144, A6).
    cfg = struct('nc', 4, 'days', 1:3, 'arcs', 4, 'mArc', 64, 'dt', 5, ...
                 'xTrue', 1e5 * [0.8 -0.45 0.3 0.15]);
    data = struct();
    data.commonParams = [repmat({'common'}, cfg.nc, 1), ...
                         arrayfun(@(j) sprintf('c%d', j), (1:cfg.nc).', 'UniformOutput', false)];
    data.days = cfg.days;
    spec = struct('type', 'diag', 'sigma', 1, 'estimate', "scale", 'arcFactors', false);
    data.systems = struct('buildDay', @(day, ot) tightPriorDay(cfg, day), ...
                          'covariance', struct('a', spec), 'arcs', cfg.arcs, 'id', 1, 'name', 'tight');
    data.priors = struct('L', eye(cfg.nc), 'x0', zeros(cfg.nc, 1), 'sigma', 1e6, ...
                         'name', 'loose', 'estimate', "scale");
end

function D = tightPriorDay(cfg, day)
%TIGHTPRIORDAY  One day of the fixture above.
    m = cfg.arcs * cfg.mArc;
    t = (day - 1)*86400 + (0:m-1).' * cfg.dt;
    w = 2*pi / 3000;
    D.A = [sin(w*t), cos(w*t), sin(2*w*t), cos(2*w*t)];
    D.B = zeros(m, 0);
    rs = RandStream('mt19937ar', 'Seed', 5500 + day);
    D.l = D.A * cfg.xTrue(:) + randn(rs, m, 1);
    D.t = t;
    D.arcIdx = arrayfun(@(a) (a-1)*cfg.mArc + (1:cfg.mArc), 1:cfg.arcs, 'UniformOutput', false);
    D = arcCells(D);                                 % v6: the cells are the arcs
end

%% ------------------------------- WP-147: the stop rule and its low-information bins
function testLowInformationBinsDoNotBlockTheStop(tc)
% categories: contract, reference | tags: diagnostics
    % WP-147: the loop stopped only when EVERY spectral bin moved by at most
    % Tolerance. On short arcs a 60-bin law has bins the data barely
    % determine; held to a precision far finer than their estimate has, they
    % kept the loop going for 14-168 passes on the fixtures of HANDOVER 12.79
    % while the solution had settled. A bin is now held to Tolerance max(1,
    % sqrt(200 / r)): at the default a tenth of its own standard error where
    % that is larger. On this fixture (30-epoch arcs, white, the default 60
    % bins) the loop settles in 8 passes, measured, and the history shows the
    % plain rule could not have stopped there.
    R = lsqvce.fit(lawShortArcFixture(), lsqvce.Options(Estimate="noise", Pool="stream", ...
                   ArcFactors=false, MaxPasses=20, Verbose=false, Seed=5));
    h = R.history;
    tc.verifyTrue(R.converged, sprintf('settles within the budget (binRatio %s)', mat2str(h.binRatio.', 3)));
    tc.verifyGreaterThan(h.maxDev(end), 1e-2, ...
        'the stop is the relaxed bins'' doing: max abs(alpha2 - 1) is still above Tolerance');
    tc.verifyLessThanOrEqual(h.binRatio(end), 1, 'every bin within its bound');
    % tol: the solution's last step in its own standard errors -- the claim
    % the rule makes is that the estimate has settled (measured 0.007)
    tc.verifyLessThan(h.xStep(end), 0.05, 'the solution has settled');
    T = R.convergence();
    tc.verifyEqual(T.binRatio, h.binRatio);
    tc.verifyEqual(T.xStep(2:end), h.xStep(2:end));
    tc.verifyTrue(isnan(T.xStep(1)), 'pass 1 has no step');
end

function testATighterToleranceTightensEveryBin(tc)
% categories: contract, metamorphic
    % The bound SCALES with Tolerance, it is not floored at a tenth of a
    % standard error: a floor would leave the weak bins' bound where it was,
    % so Tolerance 1e-3 could not ask for more digits (demo_grace does).
    % Measured on the fixture: maxDev at the stop 0.0276 at 1e-2, 0.0035 at
    % 1e-3 (49 passes) -- every bin about ten times tighter.
    data = lawShortArcFixture();
    run = @(tol) lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ...
                   ArcFactors=false, MaxPasses=80, Tolerance=tol, Verbose=false, Seed=5));
    R2 = run(1e-2);  R3 = run(1e-3);
    tc.verifyTrue(R2.converged && R3.converged);
    tc.verifyGreaterThan(R3.history.nIter, R2.history.nIter, 'a tighter tolerance takes longer');
    % tol: a factor 2 between the two stops; the scaled bound gives ~10, a
    % floored one ~1 -- the band separates the two designs, not two runs
    tc.verifyLessThan(R3.history.maxDev(end), 0.5 * R2.history.maxDev(end), ...
        'the weak bins tighten with Tolerance too');
end

function testNotConvergedWarningSaysWhatStillMoves(tc)
% categories: contract | tags: diagnostics
    % WP-147 (option B): a stop on the pass budget names the bins that block
    % the rule, with their redundancy against their group's median, and the
    % solution's last step in its own standard errors -- what a user needs
    % to judge whether the budget stop matters. setupOnce mutes the id for
    % the whole suite, so this test switches it back on for itself.
    w = warning('on', 'lsqvce:adjustVCEjoint:notConverged');
    tc.addTeardown(@() warning(w));
    lastwarn('');
    tc.verifyWarning(@() lsqvce.fit(lawShortArcFixture(), lsqvce.Options(Estimate="noise", ...
        Pool="stream", ArcFactors=false, MaxPasses=2, Verbose=false, Seed=5)), ...
        'lsqvce:adjustVCEjoint:notConverged');
    msg = lastwarn;
    tc.verifySubstring(msg, 'Still moving by more than');
    tc.verifyTrue(~isempty(regexp(msg, 'x1\|a bin \d+ of 60 \(redundancy [\d.]+, median [\d.]+, moved', 'once')), msg);
    tc.verifyTrue(~isempty(regexp(msg, 'The solution moved by [\d.e+-]+ of its standard error', 'once')), msg);
end

function data = lawShortArcFixture()
%LAWSHORTARCFIXTURE  ONE white stream whose LAW is estimated on the default
%   60 bins from short arcs: 3 days x 20 arcs x 30 epochs, 4 commons, the
%   level 0.7 against a declared 1 (WP-147). Deterministic.
    cfg = struct('nc', 4, 'days', 1:3, 'arcs', 20, 'mArc', 30, 'dt', 5, ...
                 'xTrue', [0.8 -0.45 0.3 0.15]);
    data = struct();
    data.commonParams = [repmat({'common'}, cfg.nc, 1), ...
                         arrayfun(@(j) sprintf('c%d', j), (1:cfg.nc).', 'UniformOutput', false)];
    data.days = cfg.days;
    spec = struct('type', 'diag', 'sigma', 1, 'estimate', "law", 'arcFactors', false);
    data.systems = struct('buildDay', @(day, ot) lawShortArcDay(cfg, day), ...
                          'covariance', struct('a', spec), 'arcs', cfg.arcs, 'id', 1, 'name', 'lowInfo');
end

function D = lawShortArcDay(cfg, day)
%LAWSHORTARCDAY  One day of the fixture above.
    m = cfg.arcs * cfg.mArc;
    t = (day - 1)*86400 + (0:m-1).' * cfg.dt;
    w = 2*pi / 3000;
    D.A = [sin(w*t), cos(w*t), sin(2*w*t), cos(2*w*t)];
    D.B = zeros(m, 0);
    rs = RandStream('mt19937ar', 'Seed', 4700 + day);
    D.l = D.A * cfg.xTrue(:) + 0.7 * randn(rs, m, 1);
    D.t = t;
    D.arcIdx = arrayfun(@(a) (a-1)*cfg.mArc + (1:cfg.mArc), 1:cfg.arcs, 'UniformOutput', false);
    D = arcCells(D);                                 % v6: the cells are the arcs
end

function data = shortArcFixture(seed)
%SHORTARCFIXTURE  ONE white stream, "scale" on a declared diag law of the
%   right level, 3 days x 10 arcs x 30 epochs, 4 commons. Deterministic per seed.
    cfg = struct('nc', 4, 'days', 1:3, 'arcs', 10, 'mArc', 30, 'dt', 5, ...
                 'xTrue', [0.8 -0.45 0.3 0.15], 'seed', seed);
    data = struct();
    data.commonParams = [repmat({'common'}, cfg.nc, 1), ...
                         arrayfun(@(j) sprintf('c%d', j), (1:cfg.nc).', 'UniformOutput', false)];
    data.days = cfg.days;
    spec = struct('type', 'diag', 'sigma', 1, 'estimate', "scale", 'arcFactors', false);
    data.systems = struct('buildDay', @(day, ot) shortArcDay(cfg, day), ...
                          'covariance', struct('a', spec), 'arcs', cfg.arcs, 'id', 1, 'name', 'shortArcs');
end

function D = shortArcDay(cfg, day)
%SHORTARCDAY  One day of the fixture above.
    m = cfg.arcs * cfg.mArc;
    t = (day - 1)*86400 + (0:m-1).' * cfg.dt;
    w = 2*pi / 3000;
    D.A = [sin(w*t), cos(w*t), sin(2*w*t), cos(2*w*t)];
    D.B = zeros(m, 0);
    rs = RandStream('mt19937ar', 'Seed', 7700 + 100*cfg.seed + day);
    D.l = D.A * cfg.xTrue(:) + randn(rs, m, 1);
    D.t = t;
    D.arcIdx = arrayfun(@(a) (a-1)*cfg.mArc + (1:cfg.mArc), 1:cfg.arcs, 'UniformOutput', false);
    D = arcCells(D);                                 % v6: the cells are the arcs
end

function data = scaleOnlyFixture(scaleA, scaleB, arcFactorsA, badA)
%SCALEONLYFIXTURE  TWO white streams, both "scale" (no law), declared sigma 1;
%   stream a truly scaleA x, stream b scaleB x its declared sigma; stream a
%   carries arc factors when arcFactorsA, and its arc 2 of every day is badA x
%   worse in sigma (default 1: no bad arc). No locals. Deterministic.
    if nargin < 4, badA = 1; end
    cfg = struct('nc', 4, 'days', 1:3, 'arcs', 3, 'mArc', 256, 'dt', 5, ...
                 'xTrue', [0.8 -0.45 0.3 0.15], 'scaleA', scaleA, 'scaleB', scaleB, 'badA', badA);
    data = struct();
    data.commonParams = [repmat({'common'}, cfg.nc, 1), ...
                         arrayfun(@(j) sprintf('c%d', j), (1:cfg.nc).', 'UniformOutput', false)];
    data.days = cfg.days;
    spec = struct('type', 'diag', 'sigma', 1, 'estimate', "scale", 'arcFactors', false);
    specA = spec;  specA.arcFactors = logical(arcFactorsA);
    data.systems = struct('buildDay', @(day, ot) scaleOnlyDay(cfg, day, ot), ...
                          'covariance', struct('a', specA, 'b', spec), ...
                          'arcs', cfg.arcs, 'id', 1, 'name', 'scaleOnly');
end

function D = scaleOnlyDay(cfg, day, obsType)
%SCALEONLYDAY  One (day, stream) of the fixture above.
    m = cfg.arcs * cfg.mArc;
    t = (day - 1)*86400 + (0:m-1).' * cfg.dt;
    w = 2*pi / (0.5 * 86400);
    D.A = [sin(w*t), cos(w*t), sin(2*w*t), cos(2*w*t)];
    D.B = zeros(m, 0);
    rs = RandStream('mt19937ar', 'Seed', 6600 + 10*day + double(obsType(1)));
    sc = cfg.scaleA * ones(m, 1);
    sc(cfg.mArc + (1:cfg.mArc)) = cfg.scaleA * cfg.badA;          % arc 2 of stream a
    if strcmp(obsType, 'b'), sc(:) = cfg.scaleB; end
    D.l = D.A * cfg.xTrue(:) + sc .* randn(rs, m, 1);
    D.t = t;
    D.arcIdx = arrayfun(@(a) (a-1)*cfg.mArc + (1:cfg.mArc), 1:cfg.arcs, 'UniformOutput', false);
    D = arcCells(D);                                 % v6: the cells are the arcs
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
    % THE LEVEL. With the law estimated, sigma*sigma0 is the PROCESS sd of the
    % stream: 1/sqrt(1 - 0.5^2) = 1.155 for the clean AR(1) stream b, not 1.
    % Round 4's 1.007, read against 1 inside a 15 % band, was 0.87 of the
    % truth -- R7-12's bias, measured at tier A in WP-115 part 3a (HANDOVER
    % 12.43: robust/plain 0.83 over 20 seeds) and corrected in part 3b (only
    % rejected epochs leave the component update; the truncation factor).
    % The discriminating test is testRobustLevelOfACleanStreamIsUnbiased.
    ib = find(strcmp(Rr.noise.stream, 'x1|b'), 1);
    lvlB = Rr.noise.model{ib}.sigma * Rr.sigma0;
    % tol: 10 % relative on an effective level from 9 arcs of 256 epochs; the
    % plain run of part 3a read 1.134 against 1.155 (law estimated, 20 seeds)
    tc.verifyEqual(lvlB, 1 / sqrt(1 - plantB^2), 'RelTol', 0.10, ...
        'the clean stream''s level is its process sd under the robust pass');
    % Robust "none" carries no factors at all
    tc.verifyEmpty(Rp.noise.robust.factor);
end

function testRobustLevelOfACleanStreamIsUnbiased(tc)
% categories: reference, robustness
    % R7-12 (WP-115 part 3b, user decision 2026-10-11). Two CLEAN AR(1)
    % streams with their TRUE law declared ('arma', estimate "scale"), so the
    % level is 1 by construction and the plain loop returns 0.997 (part 3a,
    % 20 seeds). With IGG3 inside the component update the level came back
    % 0.81 / 0.86 (robust/plain, se 0.004): the influence selects the epochs
    % with large residuals. Since part 3b only REJECTED epochs leave the
    % update, and the level is divided by the rejection point's Gaussian
    % truncation factor (0.9733 at k1 = 3, R.noise.robust.truncation).
    data = levelFixture(struct('a', 0.2, 'b', 0.5), 3);
    base = {'Estimate', "noise", 'Pool', "stream", 'MaxPasses', 16, 'Verbose', false};
    Rp = lsqvce.fit(data, lsqvce.Options(base{:}));
    Rr = lsqvce.fit(data, lsqvce.Options(base{:}, 'Robust', "igg3"));
    for st = {'x1|a', 'x1|b'}
        lv = @(R) R.noise.model{find(strcmp(R.noise.stream, st{1}), 1)}.sigma * R.sigma0;
        % tol: 6 % -- the per-seed sd of the robust/plain ratio is 0.012-0.016
        % (part 3a, 20 seeds); 6 % is four of them and a third of the bias
        % being guarded (0.81-0.86)
        tc.verifyEqual(lv(Rr) / lv(Rp), 1, 'AbsTol', 0.06, sprintf(['%s: the robust level ' ...
            'equals the plain one on clean data (plain %.4f, robust %.4f)'], st{1}, lv(Rp), lv(Rr)));
        % tol: 6 % -- the plain level is 1 by construction (0.997 +- 0.011)
        tc.verifyEqual(lv(Rr), 1, 'AbsTol', 0.06, sprintf('%s: the robust level is 1', st{1}));
    end
    % tol: closed form 1 - 2 k phi(k) / (2 Phi(k) - 1) at k = 3, to 1e-10
    tc.verifyEqual(Rr.noise.robust.truncation, [1 1] * (1 - 6 * exp(-4.5) / sqrt(2*pi) / erf(3 / sqrt(2))), ...
        'RelTol', 1e-10, 'the truncation factor of IGG3''s k1 = 3, per stream');
end

function testRobustJointLoopConvergesOnCleanData(tc)
% categories: robustness, reference
    % WP-152 (review 12.84, N1). Under a robust influence every bin's Helmert
    % ratio settles at the truncation factor tau, not 1 -- a level the law's
    % renormalisation removes -- and the bins' stop test read max|alpha2 - 1|:
    % IGG3 on CLEAN data never converged (measured 30 passes, bins 1.72 of
    % their bound; Robust "none" 14). The test reads the change of shape now.
    % Clean data: the robust run must stop as the plain one does, on the same
    % law.
    data = colourFixture(0.2, 0.5);
    base = {'Estimate', "noise", 'Pool', "stream", 'MaxPasses', 30, 'Verbose', false};
    Rp = lsqvce.fit(data, lsqvce.Options(base{:}));
    Rr = lsqvce.fit(data, lsqvce.Options(base{:}, 'Robust', "igg3"));
    tc.verifyTrue(Rp.converged, 'the plain loop converges');
    tc.verifyTrue(Rr.converged, sprintf(['the robust loop converges on clean data (%d passes, ' ...
        'bins %.2f of their bound at the end)'], Rr.history.nIter, Rr.history.binRatio(end)));
    for st = {'x1|a', 'x1|b'}
        Cp = Rp.noise.Cgroups(st{1});  Cr = Rr.noise.Cgroups(st{1});
        % tol: 0.03 on a lag-1 correlation -- the two runs differ by the ~0.3 %
        % of epochs IGG3 rejects on clean data, far inside the 0.08 bound of the
        % planted colour (testPoolStreamGivesEachStreamItsOwnLaw)
        tc.verifyEqual(Cr(2) / Cr(1), Cp(2) / Cp(1), 'AbsTol', 0.03, sprintf('%s: the same law', st{1}));
    end
end

function testHeldArcKeepsItsFactorUnderRobust(tc)
% categories: robustness, contract
    % WP-152 (review 12.84, N2). An arc with as many local parameters as rows
    % has no redundancy: arcVarianceFactors HOLDS it (alpha2 = 1). Under a
    % robust influence the joint loop divided every arc's factor by tau --
    % the held arc's too, which took no ratio and carries no truncation -- so
    % it grew by 1 / tau = 2.7 % per pass and the loop never settled. Worse,
    % rounding made its zero denominator +1e-17 under IGG3: not held at all,
    % clamped to 1e-8, 1e8 times the weight (measured: 9.8e-9, growing
    % 2.7 % per pass, 12 passes, not converged). The held arc's factor must
    % stay 1 on every pass, with and without the influence.
    muteWarn(tc, 'lsqvce:arcVarianceFactors:held');      % the held arc, on purpose
    data = heldArcFixture();
    for rb = ["none", "igg3"]
        trail = containers.Map('KeyType', 'double', 'ValueType', 'any');
        R = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=true, ...
                MaxPasses=12, Robust=rb, Verbose=false, OnPass=@(ctx) keepSigma2(trail, ctx)));
        tc.verifyTrue(R.converged, sprintf('Robust %s converges', rb));
        k = find(strcmp(R.noise.arcKeys, 's1_a_d5_a1'), 1);
        f = cellfun(@(s) s(k), values(trail));
        % tol: 4 eps -- the clamp path multiplies and divides by gammaBar^2
        tc.verifyEqual(f, ones(size(f)), 'RelTol', 4*eps, ...
            sprintf('Robust %s: the held arc''s factor stays 1, pass by pass', rb));
    end
end

function testExcludedRowsAreDeletedFromTheSweep(tc)
% categories: reference, precision
    % WP-152 (review 12.84, N4). An EXCLUDED row (weight 0, zero in A and l) is
    % not an observation. Under a 'diag' law the loader leaves it a zero row
    % without a nugget (independent epochs need none for the SOLVE), and the
    % joint sweep counted it with a zero residual and full redundancy: 10 % of
    % stream a excluded gave its level 0.926 against 0.977 with the rows
    % removed. THE ORACLE: the same stream with those rows REMOVED -- for
    % independent epochs exclusion and removal are the same estimator, so the
    % levels, x and dof agree to rounding.
    R = cell(1, 2);
    for mode = 1:2
        R{mode} = lsqvce.fit(excludedFixture(mode == 1), lsqvce.Options(Estimate="noise", ...
                    Pool="stream", ArcFactors=false, MaxPasses=20, Seed=1, Verbose=false));
    end
    Tx = R{1}.summary();  Tr = R{2}.summary();
    % tol: 1e-10 -- the same estimator, the deleted rows' nugget leaves 1e-12
    % of their weight (measured 5.7e-14)
    tc.verifyEqual(Tx.sigma, Tr.sigma, 'RelTol', 1e-10, 'excluded = removed: the levels');
    % tol: 1e-10 of x's unit -- as above (measured 3.8e-15)
    tc.verifyEqual(R{1}.x, R{2}.x, 'AbsTol', 1e-10, 'excluded = removed: x');
    tc.verifyEqual(R{1}.sol.dof, R{2}.sol.dof, 'excluded = removed: the dof');
end

function keep = keepSigma2(trail, ctx)
%KEEPSIGMA2  OnPass hook: record ctx.sigma2 under the pass number.
    trail(ctx.pass) = ctx.sigma2(:).';  keep = true;
end

function data = heldArcFixture()
%HELDARCFIXTURE  Two 'diag' "scale" streams with arc factors, 5 days x 1 arc x
%   400 epochs; stream a has an arc bias and drift (localParams 'arc'), and
%   its day-5 arc has 2 rows: as many locals as rows, no redundancy (WP-152).
    m = 400;  t = (0:m-1).' * 5;  A = [sin(2*pi*t/3000), cos(2*pi*t/3000), ones(m, 1)];
    rs = RandStream('mt19937ar', 'Seed', 52);
    la = randn(rs, m, 1);  lb = randn(rs, m, 1);
    data = struct();
    data.commonParams = {'common','c1';'common','c2';'common','c3'};
    data.days = 1:5;
    spec = struct('type', 'diag', 'sigma', 1, 'estimate', "scale", 'arcFactors', true);
    data.systems = struct('buildDay', @(d, ot) heldArcDay(ot, A, la, lb, t, d), ...
        'covariance', struct('a', spec, 'b', spec), 'arcs', 1, 'id', 1, 'name', 'held', ...
        'localParams', {{'arc','bias';'arc','drift'}});
end

function D = heldArcDay(ot, A, la, lb, t, d)
%HELDARCDAY  One (day, stream) of heldArcFixture.
    l = la;  if strcmp(ot, 'b'), l = lb; end
    l = A * [1; 2; 3] + circshift(l, 37 * d);          % a different noise draw per day
    B = zeros(numel(t), 10);                            % nl = 2 per arc x 5 days
    if strcmp(ot, 'a'), B(:, 2*d-1:2*d) = [ones(size(t)), t / 1e3]; end
    if d == 5                                           % the short arc: 2 rows, 2 locals
        A = A(1:2, :);  l = l(1:2);  t = t(1:2);  B = B(1:2, :);
    end
    D = struct('A', {{A}}, 'B', {{B}}, 'l', {{l}}, 't', {{t + 86400 * d}});
end

function data = excludedFixture(exclude)
%EXCLUDEDFIXTURE  Two 'diag' "scale" streams, 4 days x 1 arc x 512 epochs;
%   every 10th row of stream a EXCLUDED (exclude = true) or REMOVED (WP-152).
    m = 512;  t = (0:m-1).' * 5;  A = [sin(2*pi*t/3000), cos(2*pi*t/3000), ones(m, 1)];
    rs = RandStream('mt19937ar', 'Seed', 44);
    la = randn(rs, m, 1);  lb = randn(rs, m, 1);
    ex = false(m, 1);  ex(1:10:end) = true;
    data = struct();
    data.commonParams = {'common','c1';'common','c2';'common','c3'};
    data.days = 1:4;
    spec = struct('type', 'diag', 'sigma', 1, 'estimate', "scale", 'arcFactors', false);
    data.systems = struct('buildDay', @(d, ot) excludedDay(ot, A, la, lb, t, ex, d, exclude), ...
        'covariance', struct('a', spec, 'b', spec), 'arcs', 1, 'id', 1, 'name', 'ex');
end

function D = excludedDay(ot, A, la, lb, t, ex, d, exclude)
%EXCLUDEDDAY  One (day, stream) of excludedFixture.
    l = la;  if strcmp(ot, 'b'), l = lb; end
    l = A * [1; 2; 3] + circshift(l, 37 * d);
    tt = t + 86400 * d;
    if strcmp(ot, 'a') && exclude
        A(ex, :) = 0;  l(ex) = 0;                       % an excluded row carries weight 0
        D = struct('A', {{A}}, 'B', {{zeros(size(A, 1), 0)}}, 'l', {{l}}, 't', {{tt}}, 'excluded', {{ex}});
        return
    end
    if strcmp(ot, 'a')
        A = A(~ex, :);  l = l(~ex);  tt = tt(~ex);
    end
    D = struct('A', {{A}}, 'B', {{zeros(size(A, 1), 0)}}, 'l', {{l}}, 't', {{tt}});
end

function testNoiseSmoothDefaultsToOne(tc)
% categories: contract, precision
    % WP-154 (review 12.84, M1). Noise.Smooth defaults to 1. MEASURED in the
    % review's closed loop (3 days x 4 arcs x 256 epochs, AR(1) 0.5 .. 0.99,
    % 10 seeds): x's variance 1.02-1.11 x the oracle's (Smooth 0: 1.10-1.23),
    % the formal sigma 1.00-1.06 x (0: 1.16-1.27), the levels as unbiased, half
    % the passes; 2048-epoch arcs unchanged. The grid tied to the arc length,
    % proposed by the review, biased the level low at high phi (16 bins: -12 %
    % at 0.98) and is not the default. Pinned here: the default reaches the
    % engine, and on the colour fixture the default run settles in fewer passes
    % than Smooth 0 on the same law.
    o = lsqvce.Options();
    tc.verifyEqual(o.Noise.Smooth, 1);
    tc.verifyEqual(o.Noise.Grid, 60, 'the grid stays 60');
    a = o.engineArgs('adjustVCEjoint');
    k = find(strcmp(a.nv(1:2:end), 'Smooth'));
    tc.verifyEqual(a.nv{2*k}, 1, 'the default reaches the joint estimator');
    data = colourFixture(0.2, 0.5);
    base = {'Estimate', "noise", 'Pool', "stream", 'MaxPasses', 30, 'Verbose', false};
    R1 = lsqvce.fit(data, lsqvce.Options(base{:}));
    R0 = lsqvce.fit(data, lsqvce.Options(base{:}, 'Noise', struct('Smooth', 0)));
    tc.verifyTrue(R1.converged && R0.converged);
    tc.verifyLessThan(R1.history.nIter, R0.history.nIter, sprintf( ...
        'Smooth 1 settles sooner (%d passes, Smooth 0: %d)', R1.history.nIter, R0.history.nIter));
    for st = {'x1|a', 'x1|b'}
        C1 = R1.noise.Cgroups(st{1});  C0 = R0.noise.Cgroups(st{1});
        % tol: 0.03 on a lag-1 correlation -- the planted 0.2 / 0.5 are read to
        % 0.08 (testPoolStreamGivesEachStreamItsOwnLaw); the two smoothings agree
        % far inside that
        tc.verifyEqual(C1(2) / C1(1), C0(2) / C0(1), 'AbsTol', 0.03, sprintf('%s: the same law', st{1}));
    end
end

function testPoolDefaultsToStream(tc)
% categories: contract
    % WP-157 (user decision 2026-10-18, measured first). Pool defaults to
    % "stream": every stream its own law. MEASURED in the review's closed loop
    % (10 seeds): pooling identical spectra gained nothing in x (eff 1.00-1.08
    % pooled, 0.90-1.08 per stream; levels unbiased either way), pooling
    % different spectra failed (eff 36.7 for two streams, 132 for four, levels
    % 2-3 x off) -- and "all", the default until then, put streams of different
    % physics in one law, which the Manual says never to do.
    o = lsqvce.Options();
    tc.verifyEqual(o.Pool, "stream");
    a = o.engineArgs('adjustVCEjoint');
    k = find(strcmp(a.nv(1:2:end), 'Pool'));
    tc.verifyEqual(string(a.nv{2*k}), "stream", 'the default reaches the joint estimator');
    R = lsqvce.fit(colourFixture(0.2, 0.5), lsqvce.Options(Estimate="noise", MaxPasses=30, Verbose=false));
    tc.verifyEqual(sort(keys(R.noise.Cgroups)), {'x1|a', 'x1|b'}, 'one law per stream by default');
end

function testPoolGroupOfDifferentSpectraWarns(tc)
% categories: robustness, contract
    % WP-154 (review 12.84, M2). One pool group shares one law: two streams of
    % different spectra in it bias both levels and the formal errors (measured
    % eff 36, levels 1.96 / 0.78, formal errors 4.85 x too large) -- silently.
    % A member not white under the group's law (Ljung-Box p < 1e-3, the last
    % pass; measured p ~ 1e-60 against >= 0.88 for one shared spectrum) warns.
    R = verifyWarning(tc, @() lsqvce.fit(colourFixture(0, 0.9), lsqvce.Options(Estimate="noise", ...
            Pool="all", MaxPasses=30, Verbose=false)), 'lsqvce:adjustVCEjoint:poolNotWhite');
    tc.verifyLessThan(min(R.noise.checks.whitePval), 1e-3, 'the checks carry the evidence');
    R = lsqvce.fit(colourFixture(0.5, 0.5), lsqvce.Options(Estimate="noise", Pool="all", ...
            MaxPasses=30, Verbose=false));
    tc.verifyFalse(any(strcmp(R.meta.warnings, 'lsqvce:adjustVCEjoint:poolNotWhite')), ...
        'one spectrum shared by the group does not warn');
end

function testDriftingLevelWarnsWithoutArcFactors(tc)
% categories: robustness, contract
    % WP-154 (review 12.84, M3). A stream without arc factors carries one level
    % for all its arcs; when the level changes from arc to arc the estimate is
    % a compromise (measured: levels 0.5 / 1 / 2 by day, eff 4.5, the level
    % 1.33 x the truth) -- silently. The arcs' Helmert ratios then scatter far
    % beyond their sampling variance (dispersion ~120 against the bound 4; a
    % homogeneous stream reaches 1.4): it warns and names arc factors.
    data = arcFactorFixture(3);                        % arc 2 of every day 3 x in sigma
    data.systems.covariance.a.arcFactors = false;
    R = verifyWarning(tc, @() lsqvce.fit(data, lsqvce.Options(Estimate="noise", MaxPasses=30, ...
            Verbose=false)), 'lsqvce:adjustVCEjoint:levelVaries');
    tc.verifyGreaterThan(R.noise.checks.dispersion(1), R.noise.checks.dispersionBound(1));
    data = arcFactorFixture(1);
    data.systems.covariance.a.arcFactors = false;
    R = lsqvce.fit(data, lsqvce.Options(Estimate="noise", MaxPasses=30, Verbose=false));
    tc.verifyFalse(any(strcmp(R.meta.warnings, 'lsqvce:adjustVCEjoint:levelVaries')), ...
        'one level for homogeneous arcs does not warn');
    tc.verifyLessThan(R.noise.checks.dispersion(1), R.noise.checks.dispersionBound(1));
end

function testLongMemoryAgainstShortArcsWarns(tc)
% categories: robustness, contract
    % WP-154 (review 12.84, M4). The per-arc local parameters absorb the low
    % frequencies of a law whose memory is long against the arcs: measured on
    % 256-epoch arcs, phi 0.995 (correlation time 200 epochs) read the level
    % 0.88 x the truth, phi 0.999 0.76 -- silently. A correlation time above
    % half the median arc warns.
    R = verifyWarning(tc, @() lsqvce.fit(colourFixture(0.998, 0.5), lsqvce.Options(Estimate="noise", ...
            Pool="stream", MaxPasses=30, Verbose=false)), 'lsqvce:adjustVCEjoint:memoryVsArc');
    c = R.noise.checks;
    k = strcmp(c.lawGroup, 'x1|a');
    tc.verifyGreaterThan(c.corrTime(k), c.arcMedian(k) / 2, 'the checks carry the evidence');
    k = strcmp(c.lawGroup, 'x1|b');
    tc.verifyLessThan(c.corrTime(k), c.arcMedian(k) / 2, 'the short-memory stream is not flagged');
end

function data = levelFixture(phi, seed)
%LEVELFIXTURE  Two CLEAN AR(1) streams, unit innovations, the TRUE law declared
%   ('arma', estimate "scale": only the level is free, and it is 1), 3 days x
%   3 arcs x 256 epochs -- dev/validation/wp115_measure.m's part-A fixture.
    cfg = struct('nc', 4, 'days', 1:3, 'arcs', 3, 'mArc', 256, 'dt', 5, ...
                 'xTrue', [0.8 -0.45 0.3 0.15]);
    data = struct();
    data.commonParams = [repmat({'common'}, cfg.nc, 1), ...
                         arrayfun(@(j) sprintf('c%d', j), (1:cfg.nc).', 'UniformOutput', false)];
    data.days = cfg.days;
    spec = @(p) struct('type', 'arma', 'ar', [1 -p], 'sigma', 1, 'estimate', "scale", 'arcFactors', false);
    data.systems = struct('buildDay', @(day, ot) levelDay(cfg, phi.(ot), seed, day, ot), ...
                          'covariance', struct('a', spec(phi.a), 'b', spec(phi.b)), ...
                          'arcs', cfg.arcs, 'id', 1, 'name', 'levelFixture');
end

function D = levelDay(cfg, phi, seed, day, obsType)
%LEVELDAY  One (day, stream) of levelFixture: stationary AR(1) noise.
    m = cfg.arcs * cfg.mArc;
    t = (day - 1)*86400 + (0:m-1).' * cfg.dt;
    w = 2*pi / (0.5 * 86400);
    D.A = [sin(w*t), cos(w*t), sin(2*w*t), cos(2*w*t)];
    D.B = zeros(m, 0);
    rs = RandStream('mt19937ar', 'Seed', 1e5*seed + 10*day + double(obsType(1)));
    e = filter(1, [1 -phi], randn(rs, m + 2000, 1));
    D.l = D.A * cfg.xTrue(:) + e(2001:end);
    D.t = t;
    D.arcIdx = arrayfun(@(a) (a-1)*cfg.mArc + (1:cfg.mArc), 1:cfg.arcs, 'UniformOutput', false);
    D = arcCells(D);                                 % v6: the cells are the arcs
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

%% ------------------------------------- WP-134: the prior is in the stop rule
function testCommonPriorIsInTheStopRule(tc)
% categories: contract, robustness | tags: diagnostics
    % WP-134 (review 3 finding 1): with every stream at estimate "none" the
    % law and the levels have nothing to move, and the loop declared itself
    % converged on pass 1 -- before the common prior group with estimate
    % "scale" was ever updated (its first Koch-Kusche step is on the pass after
    % it enters the solve): the sigma came back as declared, R.converged true,
    % no warning. The prior is in the stop rule now: it is updated, and a
    % converged run means its last change was below the Tolerance too.
    data = tc.TestData.data;
    for k = 1:numel(data.systems)
        for f = fieldnames(data.systems(k).covariance).'
            data.systems(k).covariance.(f{1}).estimate = "none";
        end
    end
    nc = numel(tc.TestData.Rv.x);
    data.priors = struct('L', eye(nc), 'x0', zeros(nc, 1), 'sigma', 0.01, 'name', 'damping', 'estimate', "scale");
    R = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, MaxPasses=6, ...
                     Noise=struct('Grid', 16), Seed=3, Verbose=false));
    tc.verifyGreaterThanOrEqual(R.history.nIter, 2, 'no convergence before the prior was updated');
    tc.verifyNotEmpty(R.noise.priorHistory, 'the prior was updated');
    tc.verifyNotEqual(R.noise.priors(1).sigma, 0.01, 'the declared sigma moved');
    tc.verifyEqual(numel(R.history.priorDev), R.history.nIter, 'one prior change per pass');
    if R.converged
        % tol: the Tolerance default (1e-2) -- the stop rule's own bound
        tc.verifyLessThanOrEqual(R.history.priorDev(end), 1e-2, 'converged: the prior settled too');
    end
end
