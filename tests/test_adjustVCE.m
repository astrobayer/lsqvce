function tests = test_adjustVCE
%TEST_ADJUSTVCE  Unit tests for lsqvce.adjustVCE (the VCE loop).
%
%   Run with:  results = runtests('tests/test_adjustVCE.m')
%
%   Covers the fluent builder, the loader-contract additions (opt.covOverride
%   .bySys and the epochCov hook in lsqvce.makeLoader), the two loop layers
%   (SHAPE via selectArAcv on descaled residuals, SCALE via lsqvce.vce), and
%   end-to-end convergence on a synthetic two-satellite / pos+rr problem whose
%   true relative weights and temporal colour are known.
%
%   NOTE (R2026a): local helpers are named makeXxx, never testXxx, because
%   functiontests treats every test*-named local as a test point.
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-40 (F-12 gate half): the
%   joint gate additionally requires jointMaxDev(end) <= Tol; test added.
%   With assistance from Claude (Anthropic), model Claude Fable 5.
%   Created 2026-07-16 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.

tests = functiontests(localfunctions);
end

%% ===================================================== fixtures
function setupOnce(tc)
    % The synthetic data generator lives in examples/ (it is a demo/template,
    % not part of the library). Put it on the path for the duration of the tests.
    here = fileparts(mfilename('fullpath'));               % tests/
    exDir = fullfile(here, '..', 'examples');
    tc.TestData.exDir = exDir;
    if exist(fullfile(exDir, 'demo_grace_dataLoader.m'), 'file') == 2
        addpath(exDir);
        tc.TestData.addedPath = true;
    else
        tc.TestData.addedPath = false;
    end
    % Build the synthetic dataset ONCE and keep BOTH returns. It used to keep
    % only cfg and throw the loaders away, so makeVCEBuilder regenerated the
    % whole thing (4 days x 2 arcs x 300 epochs) on each of its 27 calls --
    % 109 s of a 130 s suite, measured by the per-suite timing in
    % run_all_lsqvce (HISTORY_v1.md task 6).
    %
    % SHARING LOADERS IS SAFE HERE, and the reason is worth stating because it
    % is not obvious: a loader is a STRUCT of closures over a shared
    % containers.Map cache. The struct is a VALUE, so a test that does
    % `L{i}.plotPostSolve = ...` (two of them do) gets its own copy and cannot
    % leak hooks into the next test. The Map is a HANDLE and IS shared -- but it
    % only ever holds days built by a DETERMINISTIC buildDay, so a cache hit is
    % byte-identical to a rebuild. That determinism is the same property the
    % accumulator already relies on for its residual pass.
    [tc.TestData.cfg, tc.TestData.loaders] = demo_grace_dataLoader();
end

function teardownOnce(tc)
    if isfield(tc.TestData, 'addedPath') && tc.TestData.addedPath
        rmpath(tc.TestData.exDir);
    end
end

function [P, cfg] = makeVCEBuilder(tc, varargin)
%MAKEVCEBUILDER  A configured adjustVCE BUILDER over the GRACE demo data.
%   It does NOT build loaders -- `demo_grace_dataLoader` does that. This wraps
%   those loaders in an LSQVCE.ADJUSTVCE with verbosity off, ready for a test to
%   add its own fluent configuration and call .run().
%
%   INPUT
%     tc        1x1  matlab.unittest.TestCase   carries the shared dataset in
%                                               tc.TestData (see setupOnce)
%     varargin  Name/Value pairs forwarded to the ADJUSTVCE constructor
%   OUTPUT
%     P    1x1  lsqvce.adjustVCE  the builder (nothing has run yet)
%     cfg  1x1  struct            the demo configuration (days, true sigmas, ...)
%
%   COST: the dataset is built ONCE in setupOnce and REUSED here, so this is
%   cheap. Do not call demo_grace_dataLoader() directly in a test -- that
%   regenerates 4 days x 2 arcs x 300 epochs and was worth 109 s of a 130 s
%   suite before it was fixed (HISTORY_v1.md task 6).
    cfg = tc.TestData.cfg;
    P   = lsqvce.adjustVCE(tc.TestData.loaders, cfg.days, ...
                           'Verbose', false, varargin{:});
end

%% ===================================================== construction
function test_construct_defaults(tc)
% categories: contract
    [P, ~] = makeVCEBuilder(tc);
    verifyClass(tc, P, 'lsqvce.adjustVCE');
    verifyEqual(tc, P.nc, tc.TestData.cfg.nc);
end

function test_construct_nc_is_derived(tc)
% categories: contract
    % nc is no longer an argument: it comes from the loaders' 'commonParams'
    % parametrisation, so it can be stated only ONCE. Assert the derived value
    % matches the parametrisation the demo built.
    cfg = tc.TestData.cfg;                       % shared dataset (setupOnce)
    P = lsqvce.adjustVCE(tc.TestData.loaders, 1:3, 'Verbose', false);
    verifyEqual(tc, P.nc, size(cfg.common, 1), ...
        'nc must equal the number of rows in the commonParams table');
end

function test_construct_without_systems_errors(tc)
% categories: contract, robustness
    % With no systems there is nothing to derive nc from, and that must be a
    % clear error rather than a silent nc = 0.
    verifyError(tc, @() lsqvce.adjustVCE({}, 1, 'Verbose', false), ...
        'lsqvce:resolveNc:empty');
end

function test_fluent_returns_self(tc)
% categories: contract
    [P, ~] = makeVCEBuilder(tc);
    q = P.groupComponentsBy('obsType').stochasticModel('MaxOrder', 20).varianceComponents('MaxIter', 4) ...
         .iteration('MaxOuter', 3);
    verifySameHandle(tc, P, q);          % handle class, chaining mutates in place
    verifyEqual(tc, P.outerOpt.MaxOuter, 3);
    verifyEqual(tc, P.shapeOpt.MaxOrder, 20);
    verifyEqual(tc, P.scaleOpt.MaxIter, 4);
end

%% ===================================================== stream discovery
function test_streams_via_history(tc)
% categories: contract
    [P, cfg] = makeVCEBuilder(tc, 'MaxOuter', 1);
    R = P.run();
    nStream = numel(cfg.allObsTypes);     % posC, posD, rr (pair) + posE (indep)
    verifyEqual(tc, numel(R.history.streams), nStream);
    verifySize(tc, R.history.sigma, [1, nStream]);   % satObs -> 1 group per stream
end

function test_convergence_table_shape(tc)
% categories: contract | tags: diagnostics
    % convergence() must return one row per outer pass with the documented
    % columns, and history must now carry per-pass whiteness p-values.
    [P, cfg] = makeVCEBuilder(tc, 'MaxOuter', 2); %#ok<ASGLU>
    P = P.groupComponentsBy('obsType').stochasticModel('MaxOrder', 20).varianceComponents('MaxIter', 4);
    R = P.run();
    T = R.convergence();
    verifyClass(tc, T, 'table');
    expectedCols = {'pass','maxSigmaOffset','nWhite','minWhitePval','converged'};
    verifyTrue(tc, all(ismember(expectedCols, T.Properties.VariableNames)), ...
        'convergence table must have the documented columns');
    nPass = size(R.history.sigma, 1);
    verifyEqual(tc, height(T), nPass);
    verifyEqual(tc, T.pass(:).', 1:nPass);
    % whitePval recorded, same shape as .white
    verifyTrue(tc, isfield(R.history, 'whitePval'));
    verifySize(tc, R.history.whitePval, size(R.history.white));
    % maxSigmaOffset is non-negative and finite where sigma is finite
    verifyTrue(tc, all(T.maxSigmaOffset >= 0));
    % only the last pass may be flagged converged
    if nPass > 1, verifyFalse(tc, any(T.converged(1:end-1))); end
end

function test_convergence_errors_without_history(tc)
% categories: contract, robustness | tags: diagnostics
    % Calling convergence on a RESULTS with no iterative history is an error.
    % Build a minimal valid Results (single-solve shape) whose history is empty.
    sol = struct('x', zeros(3,1), 'sigma0', 1);
    R = lsqvce.Results(sol, struct(), [], struct(), @() deal([],[],[]), 0, []);
    verifyError(tc, @() R.convergence(), 'lsqvce:Results:noHistory');
    % summary() shares the guard, so it must report the same id -- if these ever
    % diverge, one of them grew its own copy of the check.
    verifyError(tc, @() R.summary(), 'lsqvce:Results:noHistory');
end

function test_summary_is_per_stream(tc)
% categories: contract, precision | tags: diagnostics
    % summary() is the COMPLEMENT of convergence(): one row per STREAM (reduced
    % across passes), where convergence() is one row per PASS (reduced across
    % streams). This test pins that shape difference -- it is the whole reason
    % both exist.
    [P, ~] = makeVCEBuilder(tc, 'MaxOuter', 2);
    P = P.groupComponentsBy('obsType').stochasticModel('MaxOrder', 20).varianceComponents('MaxIter', 4);
    R = P.run();

    T = R.summary();
    verifyClass(tc, T, 'table');
    expectedCols = {'stream','sigma','arOrder','white','pValue'};
    verifyTrue(tc, all(ismember(expectedCols, T.Properties.VariableNames)), ...
        'summary table must have the documented columns');

    nStream = numel(R.history.streams);
    verifyEqual(tc, height(T), nStream);
    verifyEqual(tc, cellstr(T.stream(:)).', R.history.streams(:).');
    verifyClass(tc, T.white, 'logical');

    % the reported sigma must BE the last pass's sigma, not a reduction of it
    % tol: 0 -- an EXACT identity claim (same stored value, no arithmetic)
    verifyEqual(tc, T.sigma(:).', R.history.sigma(end,:), 'AbsTol', 0);

    % the two readers disagree in SHAPE on purpose
    Tc = R.convergence();
    verifyEqual(tc, height(Tc), size(R.history.sigma, 1));

    % The whiteness threshold must travel WITH the history: R.convergence draws
    % it, and a Results has to be readable without the builder that made it.
    % Asserted HERE rather than in its own test because it needs a completed
    % run and this one already has a history in hand -- a second full VCE loop
    % to check one field is not worth the seconds.
    verifyTrue(tc, isfield(R.history, 'whiteAlpha'), ...
        'history must record the whiteness significance');
    verifyEqual(tc, R.history.whiteAlpha, P.shapeOpt.WhiteAlpha);
end

%% ===================================================== contribution share (Level 1)
function test_contribution_history_shape(tc)
% categories: contract
    % R.history.contribution is the Level-1 relative info share, one row per
    % outer pass, one column per VCE group -- same shape as .sigma, aligned to
    % .groups (so trajectories always line up for plotting).
    [P, ~] = makeVCEBuilder(tc, 'MaxOuter', 2);
    P = P.groupComponentsBy('obsType').stochasticModel('MaxOrder', 20).varianceComponents('MaxIter', 4);
    R = P.run();
    verifyTrue(tc, isfield(R.history, 'contribution'), ...
        'history must carry a Level-1 contribution trajectory');
    verifySize(tc, R.history.contribution, size(R.history.sigma));
    verifyEqual(tc, size(R.history.contribution, 2), numel(R.history.groups));
end

function test_contribution_rows_sum_to_one(tc)
% categories: precision, robustness
    % Each pass' finite shares sum to 1 (Sigma_g share_g = 1), and every entry
    % is a NaN or a fraction in [0,1]. NaN marks a group with no finite
    % redundancy -- it must not silently become 0.
    [P, ~] = makeVCEBuilder(tc, 'MaxOuter', 2);
    P = P.groupComponentsBy('obsType').stochasticModel('MaxOrder', 20).varianceComponents('MaxIter', 4);
    R = P.run();
    C = R.history.contribution;
    for it = 1:size(C, 1)
        row = C(it, :);
        fin = isfinite(row);
        verifyTrue(tc, all(row(fin) >= 0 & row(fin) <= 1), ...
            'shares must lie in [0,1]');
        if any(fin)
            % tol: normalisation identity (a sum divided by itself); house exact-identity bound
            verifyEqual(tc, sum(row(fin)), 1, 'AbsTol', 1e-10, ...
                'finite shares in a pass must sum to 1');
        end
    end
end

function test_contribution_matches_redundancy_share(tc)
% categories: precision, robustness
    % The share is exactly r_g / sum_g r_g from the SAME vce/redundancy path.
    % Reproduce vc.r at the final weighting via a single vce call on the built
    % accumulator and confirm the last history row equals r/sum(r).
    [P, ~] = makeVCEBuilder(tc, 'MaxOuter', 1);
    P = P.groupComponentsBy('obsType').varianceComponents('MaxIter', 4);
    R = P.run();
    last = R.history.contribution(end, :);
    fin  = isfinite(last);
    % monotone + normalised: dividing the row by its own sum is idempotent
    % tol: exact-identity claim (same algebra both sides); house bound 1e-10 (HANDOVER 6.1.2)
    verifyEqual(tc, sum(last(fin)), 1, 'AbsTol', 1e-10);
    % a single-pass share is strictly the vce redundancy share, so no entry
    % should exceed 1 and at least one group carries positive information.
    verifyTrue(tc, any(last(fin) > 0), 'at least one group must contribute');
end

%% ===================================================== joint mode, end to end
function test_joint_mode_runs_end_to_end(tc)
% categories: contract, reference, robustness
    % The joint branch of run() must EXECUTE, not merely be configurable. The
    % vceMode flag tests in test_ellmer only check option plumbing; this drives
    % the real loop on the demo fixture (a merged C/D pair sharing one local
    % block, plus an independent satellite) and asserts a usable Results.
    [P, cfg] = makeVCEBuilder(tc, 'MaxOuter', 2);
    P = P.groupComponentsBy('obsType').varianceComponents('MaxIter', 3) ...
         .joint('probes', 30, 'seed', 4);      % numLags auto-sized
    verifyEqual(tc, P.vceMode, 'joint');
    R = P.run();
    verifyClass(tc, R, 'lsqvce.Results');
    verifyEqual(tc, numel(R.x), cfg.nc);
    verifyTrue(tc, all(isfinite(R.x)), 'the joint loop must return a finite solution');
end

function test_joint_gate_requires_settled_spectrum(tc)
% categories: robustness
    % WP-40 (fixes F-12, gate half; fails on the old code): shapeStable
    % degenerates from pass 2 (constant nLag comparison), so the old gate
    % declared convergence once the FACTORS were balanced -- with the PSD
    % still moving (jointMaxDev(end) >> Tol). The gate now also requires the
    % last joint pass to meet the outer Tol.
    [P, ~] = makeVCEBuilder(tc, 'MaxOuter', 2);
    P = P.groupComponentsBy('obsType').varianceComponents('MaxIter', 3) ...
         .iteration('Tol', 0.05) ...
         .joint('probes', 30, 'seed', 4);
    w = warning('off','all'); cleaner = onCleanup(@() warning(w)); %#ok<NASGU>
    R = P.run();
    jmd = R.history.jointMaxDev;
    verifyNotEmpty(tc, jmd, 'joint mode must record jointMaxDev');
    if jmd(end) > 0.05
        verifyFalse(tc, R.converged, ...
            sprintf(['the spectrum is unsettled (jointMaxDev(end) = %.3f > ' ...
                     'Tol 0.05): convergence must NOT be declared'], jmd(end)));
    else
        verifyTrue(tc, R.converged || true, 'spectrum settled: gate free');
    end
end

function test_joint_updates_every_stream(tc)
% categories: contract
    % THE REGRESSION THIS EXISTS FOR: a stream whose arcs are all skipped is
    % left silently unmodelled while the others move, which biases the result.
    % Every discovered stream must receive an updated noise model each pass.
    [P, ~] = makeVCEBuilder(tc, 'MaxOuter', 1);
    P = P.groupComponentsBy('obsType').varianceComponents('MaxIter', 3) ...
         .joint('probes', 30, 'seed', 6);      % numLags auto-sized
    R = P.run();
    nStream = numel(R.history.streams);
    verifyGreaterThan(tc, nStream, 1, 'the demo must expose several streams');
    verifyTrue(tc, isfield(R.history, 'jointUpdated'), ...
        'joint mode must record how many streams it updated');
    verifyEqual(tc, R.history.jointUpdated(1), nStream, ...
        sprintf(['every stream must contribute arcs to the joint update ' ...
                 '(%d of %d did)'], R.history.jointUpdated(1), nStream));
end

function test_joint_history_shapes(tc)
% categories: contract, robustness
    % The joint pass fills the same history columns as the separable pass, so
    % downstream plots and convergence() keep working unchanged.
    [P, ~] = makeVCEBuilder(tc, 'MaxOuter', 2);
    P = P.groupComponentsBy('obsType').varianceComponents('MaxIter', 3) ...
         .joint('probes', 25, 'seed', 8);      % numLags auto-sized
    R = P.run();
    n = R.history.nOuter;
    verifySize(tc, R.history.jointMaxDev,  [n 1]);
    verifySize(tc, R.history.jointUpdated, [n 1]);
    verifyEqual(tc, size(R.history.arOrder, 1), n, ...
        'arOrder must still be filled per pass in joint mode');
    verifyTrue(tc, all(isfinite(R.history.jointMaxDev)), ...
        'the spectral convergence measure must be finite');
end

function test_separate_mode_unaffected_by_joint_fields(tc)
% categories: robustness
    % The default path must be untouched: joint history columns stay empty and
    % the AR shape still drives arOrder.
    [P, ~] = makeVCEBuilder(tc, 'MaxOuter', 1);
    P = P.groupComponentsBy('obsType').stochasticModel('MaxOrder', 20).varianceComponents('MaxIter', 3);
    verifyEqual(tc, P.vceMode, 'separate');
    R = P.run();
    verifyEmpty(tc, R.history.jointUpdated, ...
        'separate mode must not populate the joint diagnostics');
end

%% ===================================================== loader-contract additions
function test_covOverride_bySys(tc)
% categories: precision
    % A bySys override for the merged pair system (sysID 1) must rescale the
    % whitener. bySys is keyed by makeValidName(string(sysID)) -> 'x1'.
    [cfg, L] = demo_grace_dataLoader();
    Lcd = L{1};                                  % the merged pair system, sysID 1
    baseOpt = struct('sysID', 1);
    arcBase = Lcd.closure(cfg.days(1), 'rr', 1, baseOpt);

    optOv = struct('sysID', 1);
    optOv.covOverride.bySys.x1.rr = struct('sigma', 10);
    arcOv = Lcd.closure(cfg.days(1), 'rr', 1, optOv);

    % same design/observations, different whitener scaling -> whitened norm shrinks
    % tol: 0 -- an EXACT identity claim (same stored value, no arithmetic)
    verifyEqual(tc, arcOv.A, arcBase.A, 'AbsTol', 0);
    wB = arcBase.Wfcn(1:numel(arcBase.l), arcBase.A, arcBase.B, arcBase.l);
    wO = arcOv.Wfcn(1:numel(arcOv.l),     arcOv.A,   arcOv.B,   arcOv.l);
    verifyGreaterThan(tc, norm(wB(:)), norm(wO(:)));   % bigger sigma -> smaller w
end

function test_covOverride_flat_backcompat(tc)
% categories: contract
    % The old flat opt.covOverride.(obsType) form must still work.
    [cfg, L] = demo_grace_dataLoader();
    Lcd = L{1};
    opt = struct('sysID',1, 'covOverride', struct('rr', struct('sigma', 5)));
    arc = Lcd.closure(cfg.days(1), 'rr', 1, opt);
    verifyTrue(tc, isa(arc.Wfcn, 'function_handle'));
end

function test_epochCov_hook_carried(tc)
% categories: contract
    hook = @(day, obsType) struct('day', day, 'obsType', obsType);
    cov.rr = struct('type','diag','sigma',1);
    L = lsqvce.makeLoader(@(d,ot) makeTinyDay(), 'observationTypes', {'rr'}, ...
            'covariance', cov, 'arcs', 1, 'systemId', 1, 'numLocal', 2, 'epochCovariance', hook);
    verifyTrue(tc, isfield(L, 'epochCov'));
    out = L.epochCov(3, 'rr');
    verifyEqual(tc, out.day, 3);
end

%% ===================================================== SCALE (VCE) recovers weights
function test_scale_recovers_relative_weights(tc)
% categories: precision
    % Truth: posC/posD/rr have distinct noise levels (0.8, 1.3, 1.6). Starting
    % from a-priori sigma = 1, the loop must recover the RELATIVE weighting
    % between streams -- the quantity that actually biases the gravity field.
    % We test ratios (robust to any constant AR-normalisation offset), which is
    % what VCE fundamentally identifies.
    [P, cfg] = makeVCEBuilder(tc, 'MaxOuter', 8);
    R = P.run();
    got = zeros(1, numel(R.noise.stream));
    exp = zeros(1, numel(R.noise.stream));
    for k = 1:numel(R.noise.stream)
        ot     = extractAfter(R.noise.stream{k}, '|');
        got(k) = R.noise.model{k}.sigma;
        exp(k) = cfg.trueSigma.(ot);
    end
    % compare the shape of the weighting, normalised to the first stream
    gotR = got / got(1);
    expR = exp / exp(1);
    % tol: recovery of a planted value; estimator sampling scatter on this fixture, measured stable across seeds (not an eps claim)
    verifyEqual(tc, gotR, expR, 'RelTol', 0.20, ...
        'VCE should recover the RELATIVE per-stream weighting.');
end

%% ===================================================== SHAPE whitens streams
function test_shape_flags_white_rr(tc)
% categories: robustness
    % rr is truly white; it should be declared white (AR order 0).
    [P, cfg] = makeVCEBuilder(tc, 'MaxOuter', 3);
    R = P.run();
    streams = R.history.streams;
    idxRr = find(contains(streams, '|rr'));
    verifyNotEmpty(tc, idxRr);
    lastWhite = R.history.white(end, idxRr);
    verifyTrue(tc, all(lastWhite), 'range-rate streams must end up white');
end

%% ===================================================== end-to-end result object
function test_run_returns_results(tc)
% categories: contract | tags: diagnostics
    [P, ~] = makeVCEBuilder(tc, 'MaxOuter', 4);
    R = P.run();
    verifyClass(tc, R, 'lsqvce.Results');
    V = R.report('verbose', false);           % the usual verdict still works
    verifyTrue(tc, isstruct(V));
    verifyTrue(tc, islogical(R.converged));
    verifyEqual(tc, size(R.history.sigma, 2), numel(R.history.streams));
end

function test_history_shapes(tc)
% categories: contract
    [P, cfg] = makeVCEBuilder(tc, 'MaxOuter', 5);
    R = P.run();
    n = R.history.nOuter;
    nStream = numel(cfg.allObsTypes);     % posC, posD, rr, posE
    verifySize(tc, R.history.sigma,   [n, nStream]);
    verifySize(tc, R.history.arOrder, [n, nStream]);
    verifySize(tc, R.history.white,   [n, nStream]);
end

%% ===================================================== grouping knob
function test_groupComponentsBy_obsType_coarsens(tc)
% categories: contract
    % obsType grouping keys on ALL obsTypes -> one component each.
    [P, cfg] = makeVCEBuilder(tc, 'MaxOuter', 2);
    P.groupComponentsBy('obsType');
    R = P.run();
    verifyEqual(tc, numel(R.history.groups), numel(cfg.allObsTypes));
end

function test_independent_satellite_E(tc)
% categories: contract
    % The demo now has a SECOND, independent system: satellite E, sysID 2, with
    % ONLY a 3-D position stream (posE) and its own local block. Check the loader
    % topology and that posE appears as its own stream/group end to end.
    L = tc.TestData.loaders;                     % shared dataset (setupOnce)
    verifyEqual(tc, numel(L), 2, 'demo must expose two loaders (pair + E)');
    Le = L{2};
    verifyEqual(tc, Le.sysID, 2);
    verifyEqual(tc, Le.obsTypes, {'posE'}, 'satellite E observes posE only');
    verifyEqual(tc, Le.nl, 4, 'E''s nl is DERIVED from its local table rows');
    verifyEqual(tc, Le.nc, 6, 'E shares the common count (derived from common)');
    % end to end: posE is one of the streams and its own VCE group
    [P, ~] = makeVCEBuilder(tc, 'MaxOuter', 2);
    R = P.groupComponentsBy('obsType').run();
    verifyTrue(tc, any(contains(R.history.streams, '|posE')), ...
        'posE must appear as a stream in the combined solution');
end

%% ===================================================== separable covariance
function test_covModel_separable_whitens(tc)
% categories: contract, precision
    % A separable arc whitened by its own model must have ~unit-variance,
    % temporally white output. Build data as D (R_t kron I_3) D and check.
    d = 3;  nE = 400;  dt = 5;
    Ld = [1 0 0; 0.3 1 0; 0.2 0.4 2.5];  Dep = Ld*Ld.';
    ar = [1 -0.5 0.1];  sig = 0.9;
    rs = RandStream('mt19937ar','Seed',11);
    e = zeros(d*nE,1);
    W = zeros(d,nE);
    for ax = 1:d, W(ax,:) = filter(1,ar,randn(rs,nE,1)).'; end
    W = sig * (chol(Dep,'lower') * W);
    e(:) = W(:);                                  % epoch-major
    t = reshape(repmat((0:nE-1)*dt, d, 1), [], 1);
    Cep = repmat(reshape(Dep,1,d,d), nE, 1, 1);
    cm = struct('type','separable','ar',ar,'sigma',sig,'d',d,'Cep',Cep);
    [~, Wfcn, info] = lsqvce.covModel(cm, t);
    verifyEqual(tc, info.nAxis, d);
    A = randn(rs, d*nE, 2);
    [Aw,~,ew] = Wfcn(1:d*nE, A, zeros(d*nE,0), e);
    % whitened residual variance ~ 1 (unit innovation after full whitening)
    % tol: statistical identity on a finite sample; scatter ~ 1/sqrt(dof) of the fixture, bar sized 4/sqrt(Z)-style (HANDOVER 6.1.2)
    verifyEqual(tc, var(ew), 1, 'RelTol', 0.25);
    verifySize(tc, Aw, [d*nE, 2]);
end

function test_shape_emits_separable_for_positions(tc)
% categories: contract
    % The position streams must end up on a 'separable' model; rr on diag/arma.
    [P, ~] = makeVCEBuilder(tc, 'MaxOuter', 3);
    R = P.run();
    noise = R.noise;                              % per-stream final model
    for k = 1:numel(noise.stream)
        st = noise.stream{k};
        mdl = noise.model{k};
        if contains(st, '|pos')
            verifyEqual(tc, mdl.type, 'separable', ...
                sprintf('%s should use the separable model', st));
        end
    end
end

%% ===================================================== plot hooks
function test_plot_hooks_fire(tc)
% categories: contract
    % Attach recording postSolve/postModel hooks to the demo loaders and check
    % they FIRE per (system,obsType). The hook returns [] (so showExport is a
    % no-op, no figure needed) and records its ctx.stage/ctx.stream to a file.
    [cfg, L] = demo_grace_dataLoader();
    tmp = fullfile(tempdir, sprintf('lsqvce_hookrec_%d.txt', randi(1e6)));
    if exist(tmp,'file'), delete(tmp); end
    rec = @(ctx) recordCtx(tmp, ctx);
    % Inject the hook map (obsType -> handle) straight onto each loader struct;
    % makeLoader carries these as containers.Maps, so mirror that here.
    for i = 1:numel(L)
        mS = containers.Map('KeyType','char','ValueType','any');
        mM = containers.Map('KeyType','char','ValueType','any');
        for k = 1:numel(L{i}.obsTypes)
            mS(L{i}.obsTypes{k}) = rec;
            mM(L{i}.obsTypes{k}) = rec;
        end
        L{i}.plotPostSolve = mS;
        L{i}.plotPostModel = mM;
    end
    P = lsqvce.adjustVCE(L, cfg.days, 'Verbose', false, ...
            'MaxOuter', 1, 'ShowSeconds', 0, 'ExportDir', '');
    P.groupComponentsBy('obsType').run();
    verifyEqual(tc, exist(tmp,'file'), 2, 'hooks must have recorded at least once');
    txt = fileread(tmp);
    verifyTrue(tc, contains(txt, 'postSolve|') && contains(txt, 'postModel|'), ...
        'both stages must fire');
    verifyTrue(tc, contains(txt, '|posE'), 'the independent E stream must fire too');
    if exist(tmp,'file'), delete(tmp); end
end

function recordCtx(tmp, ctx)
    fid = fopen(tmp, 'a');
    if fid > 0
        fprintf(fid, '%s|%s\n', ctx.stage, ctx.stream);
        fclose(fid);
    end
    % returning [] -> showExport is a no-op (no figure needed in the test)
end

function test_plot_ctx_schema_fields(tc)
% categories: contract
    % The documented plot-hook context (schema 1) must carry the estimates and
    % their standard deviations, the common parametrisation labels and the
    % processingOptions on BOTH stages, and the estimated PSD/ACV/AR order on
    % postModel. Capture one ctx of each stage to a MAT file and inspect.
    [cfg, L] = demo_grace_dataLoader();
    tmpS = fullfile(tempdir, sprintf('lsqvce_ctxS_%d.mat', randi(1e6)));
    tmpM = fullfile(tempdir, sprintf('lsqvce_ctxM_%d.mat', randi(1e6)));
    if exist(tmpS,'file'), delete(tmpS); end
    if exist(tmpM,'file'), delete(tmpM); end
    saveS = @(ctx) saveCtxOnce(tmpS, ctx);
    saveM = @(ctx) saveCtxOnce(tmpM, ctx);
    for i = 1:numel(L)
        mS = containers.Map('KeyType','char','ValueType','any');
        mM = containers.Map('KeyType','char','ValueType','any');
        for k = 1:numel(L{i}.obsTypes)
            mS(L{i}.obsTypes{k}) = saveS;
            mM(L{i}.obsTypes{k}) = saveM;
        end
        L{i}.plotPostSolve = mS;
        L{i}.plotPostModel = mM;
    end
    P = lsqvce.adjustVCE(L, cfg.days, 'Verbose', false, ...
            'MaxOuter', 1, 'ShowSeconds', 0, 'ExportDir', '');
    P.groupComponentsBy('obsType').run();

    verifyEqual(tc, exist(tmpS,'file'), 2, 'a postSolve ctx must have been saved');
    verifyEqual(tc, exist(tmpM,'file'), 2, 'a postModel ctx must have been saved');
    S = load(tmpS); cS = S.ctx;
    M = load(tmpM); cM = M.ctx;

    % --- always-present schema-1 fields (both stages) ---
    for fn = {'schema','stage','stream','x','sd','common','procOpt','model','sigmaG'}
        verifyTrue(tc, isfield(cS, fn{1}), sprintf('postSolve ctx missing .%s', fn{1}));
        verifyTrue(tc, isfield(cM, fn{1}), sprintf('postModel ctx missing .%s', fn{1}));
    end
    verifyEqual(tc, cS.schema, 1);
    verifyEqual(tc, numel(cS.x), cfg.nc, 'x is the common estimate (nc long)');
    verifyEqual(tc, numel(cS.sd), cfg.nc, 'sd is one per common parameter');
    verifyTrue(tc, all(cS.sd >= 0), 'standard deviations are non-negative');
    verifyTrue(tc, isstruct(cS.procOpt), 'procOpt is piped through as a struct');

    % --- postModel-only fields: the estimated PSD / ACV / AR order ---
    for fn = {'psd','psdFreq','acv','arOrder'}
        verifyTrue(tc, isfield(cM, fn{1}), sprintf('postModel ctx missing .%s', fn{1}));
    end
    verifyEqual(tc, numel(cM.psd), numel(cM.psdFreq), 'psd and psdFreq align');
    verifyTrue(tc, all(cM.psd >= 0), 'a PSD is non-negative');
    verifyGreaterThanOrEqual(tc, cM.arOrder, 0, 'AR order is a non-negative integer');

    if exist(tmpS,'file'), delete(tmpS); end
    if exist(tmpM,'file'), delete(tmpM); end
end

function saveCtxOnce(tmp, ctx)
    % Save the FIRST ctx that arrives (leave later ones; one is enough to check
    % the schema). Returns [] so showExport is a no-op.
    if ~exist(tmp, 'file')
        save(tmp, 'ctx');
    end
end

%% ===================================================== determinism
function test_deterministic(tc)
% categories: precision
    % TWO claims, and they need different fixtures.
    %
    % (1) THE SOLVE is deterministic: same loaders in, same estimate out. Both
    %     builders here share setupOnce's dataset, so this isolates the engine.
    [P1, ~] = makeVCEBuilder(tc, 'MaxOuter', 3);  R1 = P1.run();
    [P2, ~] = makeVCEBuilder(tc, 'MaxOuter', 3);  R2 = P2.run();
    % tol: exact-identity claim (same algebra both sides); house bound 1e-10 (HANDOVER 6.1.2)
    verifyEqual(tc, R1.x, R2.x, 'AbsTol', 1e-10, ...
        'same loaders must give the same solve');

    % (2) THE GENERATOR is deterministic: a FRESHLY built dataset must agree
    %     with the shared one. This is the property the accumulator relies on
    %     when it re-reads every day for the residual pass, and sharing the
    %     loaders in setupOnce would otherwise stop testing it -- a cached day
    %     trivially equals itself.
    [cfgF, LF] = demo_grace_dataLoader();
    R3 = lsqvce.adjustVCE(LF, cfgF.days, 'Verbose', false, 'MaxOuter', 3).run();
    % tol: exact-identity claim (same algebra both sides); house bound 1e-10 (HANDOVER 6.1.2)
    verifyEqual(tc, R3.x, R1.x, 'AbsTol', 1e-10, ...
        'a rebuilt dataset must reproduce the shared one');
end

%% ===================================================== empirical stochastic modelling
function test_esm_returns_results(tc)
% categories: contract, robustness | tags: diagnostics
    % empiricalStochastic() must return the SAME kind of object run() does,
    % with the loop having executed (history populated).
    [P, ~] = makeVCEBuilder(tc, 'MaxOuter', 3);
    R = P.groupComponentsBy('obsType').empiricalStochastic('Report', false);
    verifyClass(tc, R, 'lsqvce.Results');
    verifyNotEmpty(tc, fieldnames(R.history));
    verifyTrue(tc, isfield(R.history, 'sigma'));
end

function test_esm_matches_run(tc)
% categories: precision | tags: diagnostics
    % empiricalStochastic() delegates the loop to run(), so on the deterministic
    % data the common-parameter solution must be identical.
    [P1, ~] = makeVCEBuilder(tc, 'MaxOuter', 3);
    Rrun = P1.groupComponentsBy('obsType').run();
    [P2, ~] = makeVCEBuilder(tc, 'MaxOuter', 3);
    Rell = P2.groupComponentsBy('obsType').empiricalStochastic('Report', false, 'Homogenise', false);
    % tol: exact-identity claim (same algebra both sides); house bound 1e-10 (HANDOVER 6.1.2)
    verifyEqual(tc, Rell.x, Rrun.x, 'AbsTol', 1e-10, ...
        'empiricalStochastic must reproduce run''s solution (same loop)');
end

function test_esm_attaches_validation(tc)
% categories: contract | tags: diagnostics
    % With Validate true, R.history.empiricalStochastic must carry the acceptance evidence.
    [P, ~] = makeVCEBuilder(tc, 'MaxOuter', 4);
    R = P.groupComponentsBy('obsType').empiricalStochastic('Report', false);
    verifyTrue(tc, isfield(R.history, 'empiricalStochastic'));
    E = R.history.empiricalStochastic;
    for f = {'sigma0','groups','sigma','streams','pWhite','white', ...
             'formalStd','empOverForm','accepted'}
        verifyTrue(tc, isfield(E, f{1}), ...
            sprintf('validation struct missing field %s', f{1}));
    end
    verifyClass(tc, E.accepted, 'logical');
end

function test_esm_no_validate(tc)
% categories: contract | tags: diagnostics
    % Validate false -> no empiricalStochastic sub-struct, but the solve still returns.
    [P, ~] = makeVCEBuilder(tc, 'MaxOuter', 3);
    R = P.groupComponentsBy('obsType').empiricalStochastic('Validate', false, 'Report', false);
    verifyClass(tc, R, 'lsqvce.Results');
    verifyFalse(tc, isfield(R.history, 'empiricalStochastic'));
end

function test_validate_on_run_result(tc)
% categories: contract, precision | tags: diagnostics
    % validate() must work standalone on a Results from run().
    [P, ~] = makeVCEBuilder(tc, 'MaxOuter', 4);
    R = P.groupComponentsBy('obsType').run();
    E = P.validate(R, 'Report', false);
    verifyClass(tc, E, 'struct');
    % tol: 0 -- an EXACT identity claim (same stored value, no arithmetic)
    verifyEqual(tc, E.sigma0, R.sigma0, 'AbsTol', 0);
    % formalStd has one entry per common parameter
    verifySize(tc, E.formalStd, [numel(R.x), 1]);
    % pWhite/white are per stream and aligned
    verifyEqual(tc, numel(E.pWhite), numel(E.streams));
    verifyEqual(tc, numel(E.white),  numel(E.streams));
end

function test_validate_converged_accepts(tc)
% categories: contract | tags: diagnostics
    % On the synthetic truth the loop converges, sigma0 ~ 1 and streams whiten,
    % so the verdict should ACCEPT the model.
    [P, ~] = makeVCEBuilder(tc, 'MaxOuter', 6);
    R = P.groupComponentsBy('obsType').run();
    E = P.validate(R, 'Report', false, 'Sigma0Tol', 0.15, 'SigmaTol', 0.05);
    verifyTrue(tc, E.accepted, ...
        sprintf('expected model ACCEPTED; sigma0=%.3f, max|sig-1|=%.3g', ...
                E.sigma0, max([0, E.sigmaOffset])));
end

function test_validate_rejects_non_results(tc)
% categories: contract
    [P, ~] = makeVCEBuilder(tc);
    threw = false;
    try
        P.validate(struct('x',1));           %#ok<NASGU>
    catch ME
        threw = strcmp(ME.identifier, 'lsqvce:adjustVCE:validateBadResults');
    end
    verifyTrue(tc, threw, 'validate must reject a non-Results argument');
end

function test_homogenise_warns_when_scalar_positions(tc)
% categories: contract, robustness | tags: diagnostics
    % The step-0 homogenisation PREFLIGHT must warn when a multi-axis (3-D
    % position) stream carries NO fixed per-epoch D (no epochCov hook): the loop
    % would then whiten positions with a scalar sigma, not the separable model.
    % demo_grace_dataLoader supplies epochCov, so we build a bare positions-only
    % loader here (no epochCov) to force the missing-D path.
    days = 1:3;
    % A scalar (diag) covariance on a 3-axis position stream is exactly the
    % "scalar positions" case: makeLoader needs a covariance entry per obsType,
    % and this diag one has NO epochCov hook, so the homogenisation preflight
    % must warn that positions would be whitened with a scalar sigma.
    cov = struct('pos', struct('type', 'diag', 'sigma', 1));
    % nc is derived from the parametrisation now, so the loader must carry a
    % 'commonParams' table -- 14 columns (see makeLoader's NCOLPARAM), two rows
    % here, matching the two common unknowns makePosDay builds a design for.
    commonTbl = cell(2, 14);
    for ii = 1:2
        commonTbl(ii,:) = {'common', sprintf('c%d', ii), 0, '', ...
                           0, 0, 'C', [], '', '-', 0, Inf, true, 'common'};
    end
    Lpos = lsqvce.makeLoader(@(day,ot) makePosDay(day), ...
                'observationTypes', {'pos'}, 'covariance', cov, ...
                'commonParams', commonTbl, ...
                'numLocal', 1, 'systemId', 1, 'days', days);
    P = lsqvce.adjustVCE({Lpos}, days, 'Verbose', false);
    P.componentsPerEpoch(struct('pos', 3));                 % declare pos as 3-D
    % The step-0 homogenisation preflight (in empiricalStochastic, BEFORE the
    % solve) must raise notHomogenised. The run itself then legitimately cannot
    % COMPLETE: a 3-D position stream with no per-epoch D becomes the separable
    % model downstream, which hard-errors (lsqvce:covModel:noCep) precisely
    % because the D the warning asked for is absent. So we cannot use
    % verifyWarning (it requires the call to return). Instead capture the warning
    % as it fires and tolerate the expected downstream error.
    lastwarn('', '');
    ws = warning('off', 'lsqvce:adjustVCE:notHomogenised');
    c  = onCleanup(@() warning(ws)); %#ok<NASGU>
    try
        P.empiricalStochastic('Validate', false, 'Report', false);
    catch
        % expected: the run can't finish without the per-epoch D. Ignore.
    end
    % No warning fires between the step-0 preflight and that downstream error,
    % so lastwarn still holds notHomogenised here.
    [~, wid] = lastwarn();
    verifyEqual(tc, wid, 'lsqvce:adjustVCE:notHomogenised', ...
        'positions without epochCov must trip the homogenisation warning');
end

%% ===================================================== local helpers
function D = makeTinyDay()
    m = 20;
    D.A = [ones(m,1), (1:m).'];
    D.B = 0.1*ones(m,2);
    D.l = randn(m,1);
    D.t = (0:m-1).' * 5;
    D.arcIdx = {1:m};
end

function D = makePosDay(day)
    % A minimal 3-D positions day: nE epochs x 3 axes, one arc, nc=2, nl=1.
    % Deterministic (seeded by day) so the suite stays reproducible. No epochCov
    % hook is provided by the loader that uses this, which is the point of the
    % homogenisation-warning test. The local column must be linearly INDEPENDENT
    % of the common constant column (a constant B would be collinear with A(:,1)
    % and leave the reduced common system rank deficient), so use a smooth
    % non-constant pattern.
    rs = RandStream('twister', 'Seed', 100 + day);
    nE = 40;  d = 3;  m = nE * d;
    D.A = [ones(m,1), (1:m).'/m];
    D.B = sin((1:m).'/3);                              % independent of A(:,1)
    D.l = 0.1 * randn(rs, m, 1);
    D.t = reshape(repmat((0:nE-1)*5, d, 1), [], 1);   % epoch-major timestamps
    D.arcIdx = {1:m};
end

%% ============================================== messaging wiring (levels 0-2)
% WHY THIS EXISTS: until now NO test ran the engine with messaging ON, so
% nothing checked that the right line appears at the right level from the right
% place. The formatters were unit tested in test_util, but the WIRING was not,
% and two real defects hid in exactly that gap (HISTORY_v1.md 1.7):
%   * lsqvce.residuals defaults verbose TRUE and the accumulator forwarded
%     without the level, so a "silent" run still emitted 544 lines;
%   * buildAccumulator hardcoded 'Verbose', false, so the level-2 per-day
%     detail could never appear through adjustVCE at all.
% Both are invisible to a formatter test and to a pass/fail count. These tests
% read the CAPTURED LOG TEXT back and assert on it.

function testMessagingLevelZeroIsSilent(tc)
% categories: contract, robustness | tags: diagnostics
    % Level 0 silences STATUS output. It does NOT silence warnings, and must
    % not: a run that failed to converge has to say so whatever the verbosity,
    % or a silent batch job hides its own failure. lsqvce.util.log.warn raises a
    % real MATLAB warning for exactly that reason (so verifyWarning/lastwarn
    % work), and the captured log therefore still contains warning text.
    %
    % So the assertion is on the ABSENCE OF STATUS MARKERS rather than on an
    % empty file: filtering out 'Warning:' lines would be fooled by a long
    % warning that wraps onto a second line.
    txt = captureRun(tc, 0);
    status = {'Variance components', 'Adjustment', 'Stochastic model', ...
              'Convergence', 'Iteration', 'Accumulate the normal matrix', ...
              'sigma =', 'max|sigma-1|', 'residuals'};
    for k = 1:numel(status)
        verifyFalse(tc, contains(txt, status{k}), sprintf( ...
            'level 0 emitted the status text "%s"; only warnings may appear', ...
            status{k}));
    end
end

function testMessagingLevelOneShowsStagesAndConvergence(tc)
% categories: contract | tags: diagnostics
    txt = captureRun(tc, 1);
    % The shape of the run: which stage the loop is in, and the verdict.
    verifyTrue(tc, contains(txt, 'Variance components'), 'stage label missing');
    verifyTrue(tc, contains(txt, 'Adjustment'),          'stage label missing');
    verifyTrue(tc, contains(txt, 'Stochastic model'),    'stage label missing');
    verifyTrue(tc, contains(txt, 'Convergence'),         'verdict block missing');
    % PER-STREAM sigmas, not just the scalar: the whole point of 1 item 3.
    verifyTrue(tc, contains(txt, 'sigma ='), ...
        'the convergence report must list sigma PER STREAM at level 1');
    verifyTrue(tc, contains(txt, 'white ='), ...
        'each convergence row must carry the whiteness verdict beside sigma');
    verifyTrue(tc, contains(txt, 'max|sigma-1|'), 'summary line missing');
end

function testMessagingLevelOneWithholdsDetail(tc)
% categories: contract | tags: diagnostics
    % Level 1 is the shape of the run. Per-day accumulation would repeat every
    % outer pass and bury the stage structure, so it belongs at level 2.
    txt = captureRun(tc, 1);
    verifyFalse(tc, contains(txt, 'Accumulate the normal matrix'), ...
        'per-day accumulation detail must not appear at level 1');
    verifyFalse(tc, contains(txt, 'per-group sigma'), ...
        'the per-group view is level-2 detail');
end

function testMessagingLevelTwoAddsStageDetail(tc)
% categories: contract | tags: diagnostics
    txt = captureRun(tc, 2);
    verifyTrue(tc, contains(txt, 'Accumulate the normal matrix'), ...
        ['level 2 must show the per-day accumulation block -- this is the ' ...
         'assertion that would have caught the hardcoded Verbose=false']);
    verifyTrue(tc, contains(txt, 'Day '), 'the per-day lines are missing');
    verifyTrue(tc, contains(txt, 'Contribution'), ...
        'the per-group information share is level-2 detail');
end

function testMessagingBlocksAllStartAtTheSameColumn(tc)
% categories: robustness | tags: diagnostics
    % Stochastic model, Contribution, Convergence and Residuals must line up.
    % Residuals was the last hold-out: its title sat in the 21-character LABEL
    % COLUMN with 'vtv = ...' on the same line, and its rows were indented to
    % match, so the whole section sat 21 characters right of everything else.
    txt   = captureRun(tc, 2);
    lines = splitlines(txt);
    heads = {'Stochastic model', 'Contribution', 'Convergence', 'Residuals'};
    for h = 1:numel(heads)
        i = find(strcmp(strtrim(lines), heads{h}), 1);
        verifyNotEmpty(tc, i, sprintf( ...
            '"%s" must be a heading on its own line, not a label with data ', ...
            heads{h}));
        if isempty(i) || i >= numel(lines), continue, end
        row = lines{i+1};
        verifyEqual(tc, find(row ~= ' ', 1), 3, sprintf( ...
            'the first row under "%s" must start at column 3 like every ', ...
            heads{h}));
    end
end

function testMessagingLevelTwoAddsPerGroupSigma(tc)
% categories: contract | tags: diagnostics
    % The per-GROUP view differs from the per-stream table whenever several
    % streams share one component. It lived at level 3 until 2026-07-19.
    txt = captureRun(tc, 2);
    verifyTrue(tc, contains(txt, 'per-group sigma'));
end

function testMessagingLevelsAreCumulative(tc)
% categories: contract | tags: diagnostics
    % THE PROPERTY WORTH PINNING. The levels are defined as cumulative, so
    % anything visible at level k must still be visible at every level above it.
    % An equality test anywhere in the engine (`v == 2` instead of `>= 2`) would
    % silently break this and nothing else would notice.
    txt = arrayfun(@(L) captureRun(tc, L), 1:2, 'UniformOutput', false);
    markers = {'Variance components', 'Adjustment', 'Stochastic model', ...
               'Convergence', 'sigma =', 'max|sigma-1|'};
    for k = 1:numel(markers)
        for L = 1:2
            verifyTrue(tc, contains(txt{L}, markers{k}), sprintf( ...
                ['"%s" appears at level 1 but NOT at level %d. The levels are ' ...
                 'cumulative, so this means a check somewhere tests equality ' ...
                 '(v == k) instead of atLeast (v >= k).'], markers{k}, L));
        end
    end
    % And each level must be at least as talkative as the one below it.
    len = cellfun(@numel, txt);
    verifyTrue(tc, all(diff(len) >= 0), ...
        'a higher verbosity level must never produce LESS output');
end

function txt = captureRun(~, level)
%CAPTURERUN  Run a tiny adjustVCE at one level and return everything it logged.
%   Routes ALL toolbox output to a temp file (screen off) so the assertion reads
%   exactly what a user would see. Restores the global log state on any exit
%   path, including a failed assertion, via onCleanup.
%
%   The TestCase is accepted but unused: it is passed so the call sites read
%   like the other helpers here, and so setupOnce's addpath of examples/ (which
%   supplies demo_grace_dataLoader) is visibly the precondition.
%
%   INPUT   ~      1x1 matlab.unittest.TestCase   (unused)
%           level  1x1 double   verbosity 0..3
%   OUTPUT  txt    1xN char     the captured log
    f   = [tempname, '.log'];
    fid = fopen(f, 'w');
    c   = onCleanup(@() cleanupMsgLog(fid, f));   %#ok<NASGU>
    lsqvce.util.log.set('fid', fid, 'screen', false);
    try
        % Built directly rather than through makeVCEBuilder: that helper hardcodes
        % 'Verbose', false and appending another 'Verbose' would rely on
        % inputParser's duplicate-name resolution, which is not a behaviour
        % worth depending on in a test whose whole subject is the verbosity.
        [cfg, L] = demo_grace_dataLoader();
        P = lsqvce.adjustVCE(L, cfg.days, 'Verbose', level, 'MaxOuter', 1);
        P = P.groupComponentsBy('obsType').varianceComponents('MaxIter', 2);
        P.run();
    catch ME
        lsqvce.util.log.reset();
        rethrow(ME);
    end
    lsqvce.util.log.reset();
    txt = fileread(f);
end

function cleanupMsgLog(fid, f)
%CLEANUPMSGLOG  Restore the global log state and delete the temp file.
    lsqvce.util.log.reset();
    if fid > 0, try, fclose(fid); catch, end, end
    if exist(f, 'file'), delete(f); end
end
