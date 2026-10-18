function tests = test_fit
%TEST_FIT  Tests for lsqvce.fit, the one user verb (R4 WP-56, decision R4-D4).
%
%     lsqvce.fit   data struct (the API reference) + lsqvce.Options -> Results;
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
%   Updated 2026-10-16 (CEST, UTC+02:00), WP-145: two mutes that hid nothing
%   removed (demo_fit mutes its own two-pass run, demo_priors converges); the
%   demo capture echoes its warnings. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-16 (CEST, UTC+02:00), WP-142 (code review of WP-141): two
%   tests of prepareDay -- the hand-off (no, five and varargin arguments; the
%   hook sees what buildDay sees, procOpt without 'prepared'; a varargin
%   buildDay takes the hook) and the one-day slot -- and the refusals widened
%   (an lsqvce error passes, handles that do not resolve, the messages name
%   the data fields). With assistance from Claude (Anthropic), model Claude
%   Opus 5.5.
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-141: three tests of prepareDay --
%   one call per resident day for every stream (CacheDays 1, Cache off, Inf),
%   bitwise the same fit with the hook, the refusals and the wrapped throw.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-138 (code review of WP-135): a
%   non-integer day list is lsqvce:fit:badDays, ids of two kinds
%   lsqvce:fit:mixedIdTypes; the declaration test named for what it tests; the
%   report refuses a loop label as Iteration.loop. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-137 (code review of WP-135): the
%   Whitener's case on the noise path, an unsettled Gauss-Newton run's advice;
%   the throwing-buildDay test's name. With assistance from Claude (Anthropic),
%   model Claude Opus 5.5.
%   Updated 2026-10-14 (CEST, UTC+02:00), WP-135 (third review findings 5,
%   11, 12, 13 and the loader lows): duplicate days, one id kind, 'i' / 'j',
%   validity entries, local prior from / names refused; a throwing buildDay
%   named (the bomb a named function); layout without localParams; the
%   calendar of a system without days; a local datum, a local prior row never
%   rejected; systematic over the data rows. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-134 (third review): the robust joint sweep leaves the prior pseudo-arcs out.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-133: zero weights are rejected
%   rows, robust rejections leave the dof, Gauss-Newton residuals and a
%   non-finite trial point (review 3 findings 2, 4, 14, 17). With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
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
%
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-131: a declared physical psd under
%   Estimate "noise"; a stream named prior_* refused; snoop inherits R.diag's
%   test. With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-130: Gauss-Newton sigma0 / vPv /
%   D.x at the returned point; the datum test under a change of parameter
%   units; snoop without the re-solve and its critical value; the mad scale
%   without datum rows. With assistance from Claude (Anthropic), model Claude
%   Opus 5.5.
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-124: references to the data
%   contract point at the API reference. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
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
    % "pooled" = Pool "all" without arc factors (stated: the default is
    % "stream" since WP-157); the fixture pools an ARMA position with a white
    % KBR stream, which WP-154's check reports
    muteWarn(tc, 'lsqvce:adjustVCEjoint:poolNotWhite');
    Rp = lsqvce.fit(data, lsqvce.Options(base{:}, Estimate="noise", Pool="all", ArcFactors=false));
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
    d = data;  d.days = [1 2.5];                                  % WP-138: not noDays
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:badDays');
    d = data;  d.days = {1, 2};                                   % a cell is no day list
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:badDays');
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

function testDeclarationIsRefusedBeforeAnyRead(tc)
% categories: contract
    % WP-135 (review 3 findings 11, 12 and lane 4): refused before any read --
    % the builder THROWS if reached, and the old code reached it (or the
    % accumulator) on every one of these declarations. (WP-138: named for what
    % it tests, not for the review that found it.)
    [data, ~] = declaration();
    bomb = @(varargin) error('test:reached', 'buildDay must not be called during validation');
    data.systems(1).buildDay = bomb;
    % finding 11: a repeated day was accumulated twice
    d = data;  d.days = [1 2 2 3];
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:duplicateDay');
    d = data;  d.systems(1).days = [1 1 2];
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:duplicateDay');
    d = data;  d.systems(1).days = [1 2.5];
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:badDays');
    % finding 12: a character id next to an index default (a number); WP-138:
    % its own id, the accumulator's reason name
    d = data;  d.systems(1).id = 'C';  d.systems(2) = d.systems(1);  d.systems(2).id = [];
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:mixedIdTypes');
    d = data;  d.systems(1).id = NaN;
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:badId');
    % a common validity entry of another type, or a start after its end
    d = data;  d.paramColumns = struct('start', 3, 'end', 4);
    d.commonParams(:, 3:4) = {[]};  d.systems(1).localParams(:, 3:4) = {[]};
    d.commonParams{1, 3} = 'tomorrow';
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:badValidity');
    d.commonParams{1, 3} = 10;  d.commonParams{1, 4} = 5;
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:badValidity');
    % local prior groups: from a pass number, and one stream per group
    d = data;  d.systems(1).priors = struct('L', [1 0 0 0], 'x0', 0, 'sigma', 1, 'from', NaN);
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:badPrior');
    d = data;  d.systems(1).priors = struct('L', {[1 0 0 0], [0 1 0 0]}, 'x0', 0, 'sigma', 1, 'name', 'a');
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:badPrior');
end

function testCharacterIdsAndAThrowingBuildDay(tc)
% categories: contract
    % WP-135: ids 'i' / 'j' parsed as the complex 1i (both: an opaque
    % duplicateSystem); a buildDay that throws is lsqvce:makeLoader:buildDay,
    % with the day and the stream, the user's error the cause (documented,
    % never raised). Past validation, the bomb is the first read.
    % The bomb is a named function: an anonymous @(varargin) error(...) asked
    % for D throws MATLAB:TooManyOutputs, not the user's error (CI, 6.22.0).
    [data, ~] = declaration();
    data.systems(1).buildDay = @throwingBuildDay;
    data.systems(1).id = 'i';  data.systems(2) = data.systems(1);  data.systems(2).id = 'j';
    try
        lsqvce.fit(data, Verbose=false);
        tc.verifyFail('the builder throws');
    catch e
        tc.verifyEqual(e.identifier, 'lsqvce:makeLoader:buildDay');
        tc.verifyEqual(e.cause{1}.identifier, 'test:reached');
        tc.verifyTrue(contains(e.message, 'day') && contains(e.message, 'the user''s error'));
    end
end

function D = throwingBuildDay(varargin) %#ok<STOUT> -- it throws, by design
%THROWINGBUILDDAY  A user's buildDay that fails (testCharacterIdsAndAThrowingBuildDay).
    error('test:reached', 'the user''s error');
end

function testPrepareDayRunsOncePerDayForEveryStream(tc)
% categories: contract
    % WP-141: the per-day hook runs ONCE per build of a day, whatever the
    % number of streams that read the day, and its result reaches buildDay
    % as procOpt.prepared. Read through the loader's own day accessor
    % (L.day): both streams of day 1 cost one call; with CacheDays = 1 the
    % two streams of day 2 evict day 1, so day 1 again is a second call;
    % Cache = false keeps the current day's result (one call for two reads);
    % CacheDays = Inf keeps the built calendar (day 1 stays one call).
    [data, cfg] = declaration();
    s = data.systems;
    base = {'ObservationTypes', fieldnames(s.covariance).', 'Covariance', s.covariance, ...
            'NumLocal', cfg.nl, 'NumCommon', cfg.nc, 'Arcs', cfg.arcsPerDay, 'Days', data.days};
    L0 = lsqvce.engine.makeLoader(@(d, ot) buildDayImpl(cfg, d, ot), base{:});
    for setup = {{true, 1, 2}, {false, 1, 2}, {true, Inf, 1}}
        [useCache, cacheDays, want] = setup{1}{:};
        cnt = containers.Map('KeyType', 'double', 'ValueType', 'double');
        L = lsqvce.engine.makeLoader(@(d, ot, c, l, po) preparedStream(po, ot), base{:}, ...
                'Cache', useCache, 'CacheDays', cacheDays, ...
                'PrepareDay', @(d) prepareCounted(cnt, cfg, d));
        tc.verifyEqual(L.day(1, 'pos'), L0.day(1, 'pos'), 'buildDay reads procOpt.prepared');
        tc.verifyEqual(L.day(1, 'kbr'), L0.day(1, 'kbr'));
        L.day(1, 'pos');
        tc.verifyEqual(cnt(1), 1, 'one call serves every stream and re-read of the day');
        L.day(2, 'pos');  L.day(2, 'kbr');  L.day(1, 'pos');
        tc.verifyEqual(cnt(1), want, sprintf('Cache %d, CacheDays %g: day 1 after day 2', ...
                       useCache, cacheDays));
        tc.verifyEqual(cnt(2), 1);
    end
end

function testPrepareDayKeepsOnlyTheCurrentDay(tc)
% categories: contract
    % WP-142: the hook's result is the CURRENT day's only -- CacheDays keeps
    % the built streams, not the result (it was kept CacheDays days, dead
    % memory once the day's streams were built). With CacheDays = Inf a
    % stream of day 1 that was never built, read after day 2, costs a second
    % call; a built stream costs none; and day 2 again is a second call too.
    [data, cfg] = declaration();
    s = data.systems;
    base = {'ObservationTypes', fieldnames(s.covariance).', 'Covariance', s.covariance, ...
            'NumLocal', cfg.nl, 'NumCommon', cfg.nc, 'Arcs', cfg.arcsPerDay, 'Days', data.days};
    cnt = containers.Map('KeyType', 'double', 'ValueType', 'double');
    L = lsqvce.engine.makeLoader(@(d, ot, c, l, po) preparedStream(po, ot), base{:}, ...
            'CacheDays', Inf, 'PrepareDay', @(d) prepareCounted(cnt, cfg, d));
    L.day(1, 'pos');  L.day(2, 'pos');  L.day(1, 'kbr');
    tc.verifyEqual([cnt(1), cnt(2)], [2, 1], 'day 2 replaced day 1''s result');
    L.day(1, 'pos');
    tc.verifyEqual(cnt(1), 2, 'a built stream is read from the day cache');
    L.day(2, 'kbr');
    tc.verifyEqual(cnt(2), 2, 'day 1 replaced day 2''s result');
end

function testPrepareDayHandOff(tc)
% categories: contract
    % WP-142: the hook gets as many of (day, common, local, procOpt, info) as
    % it declares -- none (WP-141 passed one to @() ...), all five to a
    % varargin hook -- and, through fit, the same common, local and info as
    % buildDay, with procOpt the plain processingOptions (no 'prepared'); a
    % varargin buildDay sees procOpt, so it takes the hook.
    [data, cfg] = declaration();
    s = data.systems;
    base = {'ObservationTypes', fieldnames(s.covariance).', 'Covariance', s.covariance, ...
            'NumLocal', cfg.nl, 'NumCommon', cfg.nc, 'Arcs', cfg.arcsPerDay, 'Days', data.days};
    L0 = lsqvce.engine.makeLoader(@(d, ot) buildDayImpl(cfg, d, ot), base{:});
    P1 = struct('pos', buildDayImpl(cfg, 1, 'pos'), 'kbr', buildDayImpl(cfg, 1, 'kbr'));
    picker = @(d, ot, c, l, po) preparedStream(po, ot);
    L = lsqvce.engine.makeLoader(picker, base{:}, 'PrepareDay', @() P1);
    tc.verifyEqual(L.day(1, 'pos'), L0.day(1, 'pos'), 'a hook without arguments is called with none');
    rec = containers.Map('KeyType', 'char', 'ValueType', 'any');
    L = lsqvce.engine.makeLoader(picker, base{:}, ...
            'PrepareDay', @(varargin) recordedResult(rec, 'args', P1, varargin{:}));
    L.day(1, 'kbr');
    tc.verifyNumElements(rec('args'), 5, 'a varargin hook gets all five');
    L = lsqvce.engine.makeLoader(@(varargin) preparedStream(varargin{5}, varargin{2}), base{:}, ...
            'PrepareDay', @() P1);
    tc.verifyEqual(L.day(1, 'kbr'), L0.day(1, 'kbr'), 'a varargin buildDay takes the hook');
    rec = containers.Map('KeyType', 'char', 'ValueType', 'any');
    d = data;
    d.systems.processingOptions = struct('tag', 'wp142');
    d.systems.prepareDay = @(day, c, l, po, info) recordedResult(rec, sprintf('hook%d', day), ...
        struct('pos', buildDayImpl(cfg, day, 'pos'), 'kbr', buildDayImpl(cfg, day, 'kbr')), c, l, po, info);
    d.systems.buildDay = @(day, ot, c, l, po, info) recordedPick(rec, day, ot, c, l, po, info);
    lsqvce.fit(d, Verbose=false);
    for day = data.days
        hook = rec(sprintf('hook%d', day));  build = rec(sprintf('build%d', day));
        tc.verifyEqual(hook{3}, d.systems.processingOptions, 'the hook sees the plain processingOptions');
        tc.verifyEqual(rmfield(build{3}, 'prepared'), hook{3});
        tc.verifyEqual(hook([1 2 4]), build([1 2 4]), 'common, local and info: what buildDay gets');
    end
end

function P = recordedResult(rec, key, P, varargin)
%RECORDEDRESULT  A prepareDay that returns P and records the arguments it got
%   under KEY in the handle map REC -- the first call per key only, the one
%   made in the same build as the stream it serves (testPrepareDayHandOff).
    if ~isKey(rec, key), rec(key) = varargin; end
end

function D = recordedPick(rec, day, obsType, common, local, procOpt, info)
%RECORDEDPICK  A six-argument buildDay that picks its stream out of
%   procOpt.prepared and records (common, local, procOpt, info) of the
%   day's FIRST build -- the build in which the hook ran (testPrepareDayHandOff).
    key = sprintf('build%d', day);
    if ~isKey(rec, key), rec(key) = {common, local, procOpt, info}; end
    D = procOpt.prepared.(obsType);
end

function testPrepareDayGivesTheSameFitOnceADay(tc)
% categories: metamorphic, contract
    % WP-141: moving the shared work of a day into prepareDay moves work, not
    % numbers -- the same declaration with the day built in the hook and
    % picked per stream gives bitwise the same x and sigma0 -- and with
    % CacheDays = Inf the hook runs exactly once per day of the run.
    [data, cfg] = declaration();
    R0 = lsqvce.fit(data, Verbose=false);
    cnt = containers.Map('KeyType', 'double', 'ValueType', 'double');
    d = data;
    d.systems.buildDay   = @(day, ot, c, l, po) preparedStream(po, ot);
    d.systems.prepareDay = @(day) prepareCounted(cnt, cfg, day);
    R1 = lsqvce.fit(d, Verbose=false);
    tc.verifyEqual(R1.x, R0.x, 'bitwise: the same rows in the same order');
    tc.verifyEqual(R1.sigma0, R0.sigma0);
    remove(cnt, keys(cnt));
    lsqvce.fit(d, lsqvce.Options(Verbose=false, CacheDays=Inf));
    tc.verifyEqual(cell2mat(values(cnt, num2cell(data.days))), ones(1, numel(data.days)), ...
                   'CacheDays Inf: one call per day per run');
end

function testPrepareDayRefusals(tc)
% categories: contract
    % WP-141: a hook that is not a function handle, a buildDay that never sees
    % procOpt (the result would be dropped), a processingOptions that already
    % has 'prepared' (never overwritten) -- refused while building the loader,
    % relabelled under lsqvce:fit by fit; a hook that throws is
    % lsqvce:makeLoader:prepareDay with the day, the user's error the cause.
    [data, cfg] = declaration();
    s = data.systems;
    base = {'ObservationTypes', fieldnames(s.covariance).', 'Covariance', s.covariance, ...
            'NumLocal', cfg.nl, 'NumCommon', cfg.nc, 'Arcs', cfg.arcsPerDay, 'Days', data.days};
    picker = @(d, ot, c, l, po) preparedStream(po, ot);
    tc.verifyError(@() lsqvce.engine.makeLoader(picker, base{:}, 'PrepareDay', 5), ...
                   'lsqvce:makeLoader:prepareDay');
    tc.verifyError(@() lsqvce.engine.makeLoader(@(d, ot) buildDayImpl(cfg, d, ot), base{:}, ...
                   'PrepareDay', @(d) d), 'lsqvce:makeLoader:prepareUnused');
    tc.verifyError(@() lsqvce.engine.makeLoader(picker, base{:}, 'PrepareDay', @(d) d, ...
                   'ProcessingOptions', struct('prepared', 1)), 'lsqvce:makeLoader:preparedField');
    L = lsqvce.engine.makeLoader(picker, base{:}, 'PrepareDay', @throwingPrepareDay);
    try
        L.day(2, 'pos');
        tc.verifyFail('the hook throws');
    catch e
        tc.verifyEqual(e.identifier, 'lsqvce:makeLoader:prepareDay');
        tc.verifyEqual(e.cause{1}.identifier, 'test:reached');
        tc.verifyTrue(contains(e.message, 'day 2'));
    end
    % WP-142: an lsqvce error passes unchanged; a handle that does not resolve
    % is the documented id with the day (nargin threw outside the wrapper),
    % at construction when the hook's check needs buildDay's argument count
    L = lsqvce.engine.makeLoader(picker, base{:}, 'PrepareDay', @lsqvceThrowingPrepareDay);
    tc.verifyError(@() L.day(1, 'pos'), 'lsqvce:test:passThrough');
    L = lsqvce.engine.makeLoader(picker, base{:}, 'PrepareDay', str2func('lsqvce_wp142_noSuchHook'));
    try
        L.day(1, 'pos');
        tc.verifyFail('the hook does not resolve');
    catch e
        tc.verifyEqual(e.identifier, 'lsqvce:makeLoader:prepareDay');
        tc.verifyTrue(contains(e.message, 'day 1'));
    end
    noBuilder = str2func('lsqvce_wp142_noSuchBuilder');
    tc.verifyError(@() lsqvce.engine.makeLoader(noBuilder, base{:}, 'PrepareDay', @() 1), ...
                   'lsqvce:makeLoader:buildDay');
    L = lsqvce.engine.makeLoader(noBuilder, base{:});
    tc.verifyError(@() L.day(1, 'pos'), 'lsqvce:makeLoader:buildDay');
    d = data;  d.systems.prepareDay = 'notAHandle';
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:builderKind');
    d = data;  d.systems.prepareDay = @(day) day;     % the fixture's buildDay is @(d, ot)
    try
        lsqvce.fit(d);
        tc.verifyFail('a buildDay without procOpt must refuse the hook');
    catch e
        tc.verifyEqual(e.identifier, 'lsqvce:fit:prepareUnused');
        tc.verifyEqual(e.cause{1}.identifier, 'lsqvce:makeLoader:prepareUnused');
        tc.verifyTrue(startsWith(e.message, 'data.systems(1): prepareDay is declared'));
    end
    % WP-142: the relabelled message names the field the user wrote
    d = data;  d.systems.buildDay = picker;  d.systems.prepareDay = @(day) day;
    d.systems.processingOptions = struct('prepared', 1);
    try
        lsqvce.fit(d);
        tc.verifyFail('processingOptions.prepared is the user''s and never overwritten');
    catch e
        tc.verifyEqual(e.identifier, 'lsqvce:fit:preparedField');
        tc.verifyEqual(e.cause{1}.identifier, 'lsqvce:makeLoader:preparedField');
        tc.verifyTrue(startsWith(e.message, 'data.systems(1): processingOptions already has'));
    end
end

function D = preparedStream(po, obsType)
%PREPAREDSTREAM  A buildDay that picks its stream out of procOpt.prepared (WP-141).
    D = po.prepared.(obsType);
end

function P = prepareCounted(cnt, cfg, day)
%PREPARECOUNTED  A prepareDay that builds both streams of the day at once and
%   counts its calls per day in the handle map CNT (WP-141 tests).
    if isKey(cnt, day), cnt(day) = cnt(day) + 1; else, cnt(day) = 1; end
    P = struct('pos', buildDayImpl(cfg, day, 'pos'), 'kbr', buildDayImpl(cfg, day, 'kbr'));
end

function P = throwingPrepareDay(varargin) %#ok<STOUT> -- it throws, by design
%THROWINGPREPAREDAY  A prepareDay that fails (testPrepareDayRefusals).
    error('test:reached', 'the user''s error');
end

function P = lsqvceThrowingPrepareDay(varargin) %#ok<STOUT> -- it throws, by design
%LSQVCETHROWINGPREPAREDAY  A prepareDay that fails with an lsqvce id, which
%   passes the wrapper unchanged (testPrepareDayRefusals, WP-142).
    error('lsqvce:test:passThrough', 'a toolbox id names its own cause');
end

function testUnsettledGaussNewtonRunAdvisesItsOwnRecord(tc)
% categories: contract | tags: diagnostics
    % WP-137 (code review of WP-135): the per-loop advice of NOISE_UNSETTLED had
    % its IRLS case tested only. An unsettled Gauss-Newton run (maxPasses, not
    % diverged) is pointed at R.history, its pass record -- R.convergence()
    % refuses a run without the noise loop's history. The fixture of
    % testGaussNewtonDampsAndKeepsTheBestPoint: ten passes from a far start.
    t = linspace(0, 4, 30).';  y = 3*exp(-1.3*t) + 0.05*sin(7*t);
    bf = @(k, x) deal([exp(x(2)*t), x(1)*t.*exp(x(2)*t)], y - x(1)*exp(x(2)*t), []);
    muteWarn(tc, 'lsqvce:adjustNonlinear:maxIter');
    R = lsqvce.fit(struct('batchFcn', bf, 'numBatches', 1, 'x0', [0.1; 1]), ...
                   lsqvce.Options(MaxPasses=10, Tolerance=1e-10, Verbose=false));
    tc.assumeEqual(R.stopReason, "maxPasses", 'ten passes do not settle from this start');
    V = R.report('Verbose', false);
    codes = cellfun(@(f) f.code, V.findings, 'UniformOutput', false);
    k = find(strcmp(codes, 'NOISE_UNSETTLED'), 1);
    tc.assertNotEmpty(k);
    tc.verifyTrue(contains(V.findings{k}.msg, 'Gauss-Newton iteration'), 'it names the loop');
    tc.verifyTrue(contains(V.findings{k}.msg, 'R.history lists every step'));
    tc.verifyFalse(contains(V.findings{k}.msg, 'R.convergence()'));
    % WP-138: Results hands the loop as a code (Iteration.loop) and the report
    % owns the wording -- it branched on Results' display label; a label is no
    % code, and an unknown loop is refused
    it = struct('converged', false, 'stopReason', "maxPasses", 'nIter', 10, ...
                'loop', 'Gauss-Newton iteration');
    tc.verifyError(@() lsqvce.report(R.diag, 'Verbose', false, 'Iteration', it), ...
                   'lsqvce:report:badIteration');
    it = rmfield(it, 'loop');  it.what = 'Gauss-Newton iteration';   % the retired field
    tc.verifyError(@() lsqvce.report(R.diag, 'Verbose', false, 'Iteration', it), ...
                   'lsqvce:report:badIteration');
end

function testLayoutOfASystemWithoutLocals(tc)
% categories: contract, robustness
    % WP-135 (review 3 finding 13): a system without localParams has no plan;
    % lsqvce.layout stopped on plan.nl (the Manual's chapter 6 example).
    [data, ~] = declaration();
    data.systems = rmfield(data.systems, 'localParams');
    lay = lsqvce.layout(data);
    tc.verifyEqual(lay.systems(1).nl, 0);
    tc.verifyEmpty(lay.systems(1).local);
    tc.verifyEqual(lay.systems(1).col(1, 1, []), 0, 'no instance: column 0');
    tc.verifyEqual(lay.colCommon(1, []), 1);
end

function testSystemWithoutDaysClaimsOnlyTheGlobalCalendar(tc)
% categories: contract, robustness
    % WP-135 (review 3, lane 4): a system without its own days reads
    % data.days only, but claimed every calendar day for the common 'day'
    % blocks -- day 3 below got a block nobody observes (singular solve).
    [data, ~] = declaration();
    bomb = @(varargin) error('test:reached', 'layout reads no data');
    data.commonParams = {'common', 1; 'day', 2};
    data.days = 1:2;
    data.systems(1).buildDay = bomb;
    data.systems(1).id = 1;  data.systems(1).days = [];   % no own days: data.days
    s2 = data.systems(1);  s2.id = 2;  s2.days = 1:3;
    s2.arcs = @(day, obsType) 1:(day < 3);           % no arc on day 3
    data.systems(2) = s2;
    lay = lsqvce.layout(data);
    tc.verifyEqual(lay.nc, 3, 'one common row and a day row on days 1 and 2 only');
    tc.verifyEqual(lay.colCommon(2, 3), 0, 'no day-3 block');
end

function testLocalPriorRowsFollowTheCommonPriorsRules(tc)
% categories: contract, robustness | tags: diagnostics
    % WP-135 (review 3 finding 5): a LOCAL prior group got neither the datum
    % rule (WP-127) nor the no-rejection rule (WP-128). (a) Locals 3 and 4
    % enter only as their sum; "y3 - y4 = 0" is a datum the data cannot see:
    % r = 0, and it was an uncontrolled observation (REL_UNCONTROLLED, verdict
    % SUSPECT). (b) A soft local prior 10 sigma off the data's estimate was a
    % snoop candidate and rejected.
    [data, cfg] = declaration();
    data.systems(1).buildDay = @(d, ot) twinLocals(buildDayImpl(cfg, d, ot));
    data.systems(1).priors = struct('L', [0 0 1 -1], 'x0', 0, 'sigma', 1, 'name', 'twin');
    R = lsqvce.fit(data, Verbose=false);
    D = R.diag;
    tc.verifyNumElements(D.priorRows, 1, 'the one prior row');
    tc.verifyEqual(D.datumRows, D.priorRows, 'the blind local prior row is a datum row');
    tc.verifyEqual(D.nUncontrolled, 0, 'no uncontrolled observation');
    tc.verifyEqual(D.nTests, D.n - 1, 'the datum row is not tested');
    % (b) a soft local prior the data contradict
    [data, ~] = declaration();
    R0 = lsqvce.fit(data, Verbose=false);
    y0 = R0.residuals('Verbose', false).sys(1).y;
    data.systems(1).priors = struct('L', [0 1 0 0], 'x0', y0(2) + 10, 'sigma', 1, 'name', 'off');
    R1 = lsqvce.fit(data, Verbose=false);
    pr = R1.diag.priorRows;
    tc.verifyNumElements(pr, 1);
    tc.verifyTrue(ismember(pr, R1.diag.flagIdx), 'the inconsistent prior row is flagged (tested)');
    [~, REJ] = R1.snoop('Verbose', false);
    tc.verifyFalse(any(ismember(REJ.idx, pr)), 'a prior row is a declaration, never rejected');
    % systematic leaves the prior row out: a covariate over the DATA rows works
    t = (1:(R1.diag.n - 1)).';
    S = R1.systematic('Covariates', struct('t', t), 'Verbose', false);
    tc.verifyEqual(S.n, R1.diag.n);
end

function D = twinLocals(D)
% locals 3 and 4 enter every arc identically: only their sum is observed
    for a = 1:numel(D.B), D.B{a}(:, 4) = D.B{a}(:, 3); end
end

function testNotYetPathsAreRefusedExplicitly(tc)
% categories: contract
    % R5 M2 2b: Robust IRLS RUNS on the streaming kind (a policy on the block
    % path). R6 WP-105: its composition with Estimate "noise" RUNS too -- the
    % robust pass inside the joint loop -- so the front door no longer refuses
    % it, and the run carries the factors. What stays refused, explicitly and
    % never as a silent fallback, is Robust on a NONLINEAR batch problem.
    [data, ~] = declaration();
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');   % MaxPasses 2 stops early by design (a contract point)
    R = lsqvce.fit(data, lsqvce.Options(Robust="huber", Estimate="noise", Pool="stream", ...
                                        ArcFactors=false, MaxPasses=2, Verbose=false));
    tc.verifyEqual(R.noise.robust.influence, 'huber');
    tc.verifyEqual(numel(R.noise.robust.factor), 2 * 3 * 2 * 300, 'one factor per observation');
    tc.verifyEqual(size(R.noise.robust.keys, 2), 4);
    tc.verifyEqual(numel(R.history.nDown), R.history.nIter, 'one robust record per pass');
    nl = struct('batchFcn', @(k, x) deal(ones(3, 1), zeros(3, 1), ones(3, 1)), 'numBatches', 2, 'x0', 0);
    tc.verifyEqual(idOf(@() lsqvce.fit(nl, lsqvce.Options(Robust="huber"))), 'lsqvce:fit:notYet');
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

% --- issue #110 / WP-108: robust IRLS against local prior groups ----------
% Three points, not one: each names a single claim, so the runner's failure
% COUNT says which of them fired. Written as one point first, the red run
% could only say "1 failed" -- true, and useless for deciding which
% consequence of the finding is real.
%
% The shared fixture is declaration(): kind A, one system, common (nc 6) and
% local (nl 4) parameters over 2 arcs x 3 days, an ARMA law on pos and a
% diag law on kbr -- already what issue #110's test declaration asks for.

function testRobustStopRuleIsUnitInvariant(tc)
% categories: metamorphic, robustness
    % R7-17: the IRLS stop rule had an absolute floor of 1 on max|x|, so the
    % same problem in SH units (x ~ 1e-9) stopped at pass 2 "converged" about
    % one formal sigma from the fixed point. Rescaling the unknowns (A -> A/k,
    % x -> k x) must not change the pass count, and the estimate must scale.
    rs = RandStream('mt19937ar', 'Seed', 17);
    m = 400;  A = [ones(m, 1), randn(rs, m, 3)];  x = [1; -2; 0.5; 3];
    l = A*x + randn(rs, m, 1);  b = randperm(rs, m, 40);  l(b) = l(b) + 25;   % 10 % blunders
    opts = lsqvce.Options(Robust="igg3", Verbose=false);
    R1 = lsqvce.fit(struct('A', A, 'l', l), opts);
    k  = 1e-9;
    Rk = lsqvce.fit(struct('A', A / k, 'l', l), opts);
    tc.verifyEqual(numel(Rk.history), numel(R1.history), 'the same pass count in any units');
    tc.verifyGreaterThan(numel(R1.history), 2, 'a 10 %-blunder IRLS does not settle in 2 passes');
    % tol: 1e-8 relative -- the same arithmetic on rescaled columns; the stop
    % rule's own Tol (1e-8) bounds how far two runs of equal length can differ
    tc.verifyEqual(Rk.x / k, R1.x, 'RelTol', 1e-8);
end

function testRobustReportSaysTheFitIsRobust(tc)
% categories: contract | tags: diagnostics
    % R7-18: adjustRobust's help promised "report() says so"; report said
    % nothing, and judged the reweighted sigma0 and covariance as least squares.
    rs = RandStream('mt19937ar', 'Seed', 18);
    A = [ones(200, 1), randn(rs, 200, 2)];  l = A*[1; 2; 3] + randn(rs, 200, 1);
    Rr = lsqvce.fit(struct('A', A, 'l', l), lsqvce.Options(Robust="huber", Verbose=false));
    V = Rr.report('Verbose', false);
    codes = cellfun(@(f) f.code, V.findings, 'UniformOutput', false);
    tc.verifyTrue(any(strcmp(codes, 'ROBUST_FIT')), 'a robust fit says so in its verdict');
    R0 = lsqvce.fit(struct('A', A, 'l', l), Verbose=false);
    V0 = R0.report('Verbose', false);
    codes0 = cellfun(@(f) f.code, V0.findings, 'UniformOutput', false);
    tc.verifyFalse(any(strcmp(codes0, 'ROBUST_FIT')), 'a least-squares fit does not');
end

function testRobustReportCountsItsOwnPasses(tc)
% categories: contract | tags: diagnostics
    % WP-125: the settled/unsettled line read history.nIter, which only the
    % noise loop records -- an IRLS run reported "1 pass(es)" and called its
    % loop "the stochastic model". The pass count is the number of records,
    % and the line names the loop that iterated.
    rs = RandStream('mt19937ar', 'Seed', 18);
    A = [ones(200, 1), randn(rs, 200, 2)];  l = A*[1; 2; 3] + randn(rs, 200, 1);
    l(7) = l(7) + 30;                                      % something to re-weight
    Rr = lsqvce.fit(struct('A', A, 'l', l), lsqvce.Options(Robust="huber", Verbose=false));
    tc.assumeGreaterThan(numel(Rr.history), 1, 'IRLS took more than one pass');
    V = Rr.report('Verbose', false);
    codes = cellfun(@(f) f.code, V.findings, 'UniformOutput', false);
    k = find(startsWith(codes, 'NOISE_'), 1);
    tc.assertNotEmpty(k, 'an iterative run says whether its loop settled');
    msg = V.findings{k}.msg;
    tc.verifyTrue(contains(msg, 'robust re-weighting'), 'it names the loop that iterated');
    tc.verifyTrue(contains(msg, sprintf('%d pass', numel(Rr.history))), ...
                  'and counts its passes from the pass record');
end

function testUnsettledRobustRunAdvisesItsOwnOptions(tc)
% categories: contract | tags: diagnostics
    % WP-135 (review 3, lane 4): an unsettled IRLS run was told to raise
    % MaxPasses / Tolerance (they do not govern it) and pointed at
    % R.convergence(), which refused it with "a single solve does not iterate".
    rs = RandStream('mt19937ar', 'Seed', 18);
    A = [ones(200, 1), randn(rs, 200, 2)];  l = A*[1; 2; 3] + randn(rs, 200, 1);
    l(7) = l(7) + 30;
    muteWarn(tc, 'lsqvce:adjustRobust:notConverged');
    Rr = lsqvce.fit(struct('A', A, 'l', l), lsqvce.Options(Robust="huber", ...
                    RobustOpts=struct('MaxIter', 1), Verbose=false));
    tc.assumeFalse(Rr.converged, 'one IRLS pass does not settle');
    V = Rr.report('Verbose', false);
    codes = cellfun(@(f) f.code, V.findings, 'UniformOutput', false);
    k = find(strcmp(codes, 'NOISE_UNSETTLED'), 1);
    tc.assertNotEmpty(k);
    tc.verifyTrue(contains(V.findings{k}.msg, 'RobustOpts.MaxIter'));
    tc.verifyFalse(contains(V.findings{k}.msg, 'R.convergence()'));
    try
        Rr.convergence();
        tc.verifyFail('an IRLS run has no per-stream history');
    catch e
        tc.verifyEqual(e.identifier, 'lsqvce:Results:noHistory');
        tc.verifyTrue(contains(e.message, 'R.history'));
        tc.verifyFalse(contains(e.message, 'a single solve does not iterate'));
    end
end

function testRobustWarnsAtTheBreakdownPoint(tc)
% categories: robustness | tags: diagnostics
    % User decision 2026-10-10 (HANDOVER 12.34 option 2, R7 WP-115): the
    % collapse refusal fires only on a model error or near-total
    % contamination, so a pass that ZEROES more than 45 % of the rows says it
    % is at the breakdown point -- once, changing nothing. 48 % blunders of
    % 25 sigma under IGG3: zeroed, warned, and the fit still returns.
    rs = RandStream('mt19937ar', 'Seed', 45);
    m = 500;  A = [ones(m, 1), randn(rs, m, 2)];  l = A*[1; 2; 3] + randn(rs, m, 1);
    b = randperm(rs, m, 240);  l(b) = l(b) + 25 * sign(randn(rs, 240, 1));
    muteWarn(tc, 'lsqvce:robust:atBreakdown');   % provoked on purpose
    R = lsqvce.fit(struct('A', A, 'l', l), lsqvce.Options(Robust="igg3", Verbose=false));
    tc.verifyTrue(any(strcmp(R.meta.warnings, 'lsqvce:robust:atBreakdown')), ...
                  'the breakdown warning is in the run''s warning inventory');
    tc.verifyTrue(all(isfinite(R.x)), 'diagnostic only: the fit still returns');
    % 10 % blunders: far from breakdown, no warning
    l2 = A*[1; 2; 3] + randn(rs, m, 1);  l2(b(1:50)) = l2(b(1:50)) + 25;
    R2 = lsqvce.fit(struct('A', A, 'l', l2), lsqvce.Options(Robust="igg3", Verbose=false));
    tc.verifyFalse(any(strcmp(R2.meta.warnings, 'lsqvce:robust:atBreakdown')));
end

function testRobustReportsNoFactorForAPriorRow(tc)
% categories: contract
    % DISCRIMINATOR. A local prior group rides as its own stream
    % 'prior_<name>' (makeLoader>priorArc). adjustRobust's
    % lawWhitenedResiduals collected every stream, so the pass reported a
    % factor for pseudo-rows that priorArc never applies one to.
    [d, e, cj] = priorFixture();
    d.systems(1).priors = struct('L', e(cj), 'x0', 0.123, 'sigma', 1e-6, ...
                                 'from', 1, 'name', 'pin', 'estimate', "none");
    R = lsqvce.fit(d, lsqvce.Options(Robust="huber", Verbose=false));
    obsTypes = string(R.history(end).keys(:, 2));
    tc.verifyFalse(any(startsWith(obsTypes, "prior_")), ...
        'the robust pass must report no factor for a prior pseudo-row');
end

function testRobustScaleIsFormedOverDataRowsAlone(tc)
% categories: contract, precision
    % DISCRIMINATOR. A LOOSE group (sigma 1e6) leaves the solution, and so
    % every data residual, where it was; with pseudo-rows excluded the MAD
    % is formed over an identical sample and must reproduce exactly. Loose
    % on purpose: a tight group moves the solution and could not separate
    % the two behaviours.
    [d, e, cj] = priorFixture();
    ropt = lsqvce.Options(Robust="huber", Verbose=false);
    dLoose = d;
    dLoose.systems(1).priors = struct('L', e(cj), 'x0', 0, 'sigma', 1e6, ...
                                      'from', 1, 'name', 'loose', 'estimate', "none");
    % tol: 1e-12 relative -- the same data residuals through the same MAD;
    % only floating-point re-association separates them.
    tc.verifyEqual(lsqvce.fit(dLoose, ropt).history(end).scale, ...
                   lsqvce.fit(d, ropt).history(end).scale, 'RelTol', 1e-12, ...
        'a prior group must not enter the robust scale of the data');
end

function testRobustDoesNotMoveAPinnedLocal(tc)
% categories: contract, robustness
    % REGRESSION GUARD, not a discriminator, and it is expected to pass on
    % the unfixed code too -- precisely because the factor computed for the
    % prior row was never applied to it. It guards the fix.
    [d, e, cj] = priorFixture();
    d.systems(1).priors = struct('L', e(cj), 'x0', 0.123, 'sigma', 1e-6, ...
                                 'from', 1, 'name', 'pin', 'estimate', "none");
    yIrls = lsqvce.fit(d, lsqvce.Options(Robust="huber", Verbose=false)) ...
                .residuals('Verbose', false).sys(1).y;
    yGls  = lsqvce.fit(d, Verbose=false).residuals('Verbose', false).sys(1).y;
    % tol: 1e-6 absolute -- the group's own sigma. On clean data IRLS may
    % move the pinned local only within the pin that holds it.
    tc.verifyEqual(yIrls(2), yGls(2), 'AbsTol', 1e-6, ...
        'IRLS must not move a pinned local relative to the GLS solve');
end

function testRobustFactorVectorCountsDataRowsOnly(tc)
% categories: contract
    % DISCRIMINATOR for consequence (b) of issue #110, the one the other two
    % points do not reach. nAlive, nDown and nZero are all formed over the
    % factor vector, and the collapse check compares nAlive against nc + 1.
    % A prior pseudo-row in that vector pads every one of those counts, so
    % the check could pass on data already past the breakdown point -- a
    % guard that fails OPEN. The vector must be the same length with a prior
    % group as without one. A loose group is used so the solution, and so
    % the data rows themselves, are identical in both runs.
    [d, e, cj] = priorFixture();
    ropt = lsqvce.Options(Robust="huber", Verbose=false);
    nWithout = numel(lsqvce.fit(d, ropt).history(end).factor);
    dP = d;
    dP.systems(1).priors = struct('L', e(cj), 'x0', 0, 'sigma', 1e6, ...
                                  'from', 1, 'name', 'loose', 'estimate', "none");
    nWith = numel(lsqvce.fit(dP, ropt).history(end).factor);
    tc.verifyEqual(nWith, nWithout, ...
        'the factor vector, and so nAlive, must count data rows only');
end

function testRobustWithPriorsAndALawReproducesGlsOnCleanData(tc)
% categories: metamorphic
    % Issue #110 test 1, and the gap WP-107 recorded: nothing ran robust +
    % priors + a law + local parameters in ONE call. declaration() supplies
    % the ARMA law on pos, the diag law on kbr, nc = 6 common and nl = 4
    % local parameters over 2 arcs x 3 days; a LOOSE local group adds the
    % priors without moving the solution. On clean data the robust solve
    % must then reproduce the GLS one.
    %
    % The COMMON prior half of WP-107's sentence is deliberately not here.
    % data.priors reaches adjustRobust through opt.Prior -> acc.setPrior,
    % which is not part of the sweep at all, so it cannot interact with this
    % finding; and no test in this file exercises that schema, so writing
    % one blind would assert a spelling rather than a behaviour.
    [d, e, cj] = priorFixture();
    d.systems(1).priors = struct('L', e(cj), 'x0', 0, 'sigma', 1e6, ...
                                 'from', 1, 'name', 'loose', 'estimate', "none");
    R  = lsqvce.fit(d, lsqvce.Options(Robust="huber", Verbose=false));
    R0 = lsqvce.fit(d, Verbose=false);
    tc.verifyEqual(R.stopReason, "converged");
    % tol: 3e-3 absolute -- the bound testRobustRunsOnTheStreamingKind states,
    % for the same reason: Huber leaves ~5 % of Gaussian rows below factor 1,
    % so the two O(1) estimates differ by a fraction of those rows' share.
    tc.verifyEqual(R.x, R0.x, 'AbsTol', 3e-3);
end

function [data, e, cj] = priorFixture()
%PRIORFIXTURE  declaration(), plus the column selector its prior groups use.
    data = declaration();
    lay  = lsqvce.layout(data);
    nl   = lay.systems(1).nl;
    e    = @(c) double((1:nl) == c);
    cj   = lay.systems(1).col(2, [], []);
end

% --- issue #110 tests 3-6 / WP-109: the rest of the robust declaration ----
% WP-108 answered test 1 and the finding behind test 2. These are the four
% the issue left: they need API surface no point in this file exercised, and
% each states its claim as a RATIO of quantities measured in the same run,
% so no tolerance here is a number guessed at tier B.

function testIgg3RecoversPlantedBlundersWhereGlsDoesNot(tc)
% categories: robustness, reference
    verifyBlundersAreRejected(tc, "igg3");
end

function testBisquareRecoversPlantedBlundersWhereGlsDoesNot(tc)
% categories: robustness, reference
    verifyBlundersAreRejected(tc, "bisquare");
end

function verifyBlundersAreRejected(tc, influence)
%VERIFYBLUNDERSAREREJECTED  Issue #110 test 3, for one redescending influence.
%   ONE claim, two halves that cannot be separated: the robust estimate
%   recovers the planted truth to the accuracy the SAME fit reaches on clean
%   data, and the plain GLS estimate does not. A recovery claim without its
%   control says nothing -- an estimator that ignored the data entirely would
%   pass the first half on its own.
%
%   The margin is measured, not assumed: eClean is this fixture's own GLS
%   error on undisturbed data, so the two bounds below are ratios and carry
%   no absolute scale. 60 blunders of 200 (about 68 sigma on the ARMA 'pos'
%   stream, 1000 on the white 'kbr' one) over 3600 rows: 1.7 %, far inside
%   any breakdown point, and gross enough that the GLS bias is unambiguous.
    [clean, cfg] = declaration();
    dirty = clean;
    dirty.systems(1).buildDay = @(day, ot) plantBlunders(buildDayImpl(cfg, day, ot), 5, 200);
    xTrue  = cfg.x_c;
    eClean = norm(lsqvce.fit(clean, Verbose=false).x - xTrue);
    eGls   = norm(lsqvce.fit(dirty, Verbose=false).x - xTrue);
    R      = lsqvce.fit(dirty, lsqvce.Options(Robust=influence, Verbose=false));
    eRob   = norm(R.x - xTrue);
    tc.verifyEqual(R.stopReason, "converged");
    tc.verifyLessThan(eRob, eGls/10, ...
        sprintf('%s must beat the GLS estimate on planted blunders by an order of magnitude', influence));
    tc.verifyLessThan(eRob, 3*eClean, ...
        sprintf('%s must reach the accuracy this fixture has on clean data', influence));
end

function D = plantBlunders(D, nPerArc, amp)
%PLANTBLUNDERS  Add a gross error to nPerArc rows of every arc of a built day.
%   Deterministic (no RNG, so the planted set is the same on every runner)
%   and spread over the arc's interior, away from the boundaries where the
%   ARMA whitener's start-up would blur the spike into its neighbours.
    for k = 1:numel(D.l)
        m = numel(D.l{k});
        if m < nPerArc + 2, continue, end
        r = round(linspace(1, m, nPerArc + 2));
        D.l{k}(r(2:end-1)) = D.l{k}(r(2:end-1)) + amp;
    end
end

function testRobustFixedPointOnTheBatchKind(tc)
% categories: reference, contract
    % Issue #110 test 4, on kind C (batchFcn with an explicit p) -- the kind
    % the issue names and the one no robust point covered; the kind-B pin
    % testMatrixKindRobustReproducesRobust is its counterpart.
    %
    % The issue asks for 1e-10. Until WP-132 that bound was not reachable: the
    % factors returned were the NEXT pass's (derived from the last estimate),
    % one Lipschitz step from the weights x was solved with. Since WP-132
    % (review low 24) history(end).factor IS the returned solve's weights, so
    % x and the closed form on them agree to rounding.
    [A, l, p] = matrixFixture();
    l(7) = l(7) + 25;  l(22) = l(22) - 18;           % two gross errors
    idx  = {1:15, 16:30};
    bf   = @(k) deal(A(idx{k}, :), l(idx{k}), p(idx{k}));
    data = struct('batchFcn', bf, 'numUnknowns', size(A, 2), 'numBatches', 2);
    R = lsqvce.fit(data, lsqvce.Options(Robust="huber", ...
                   RobustOpts=struct('MaxIter', 40, 'Tol', 1e-12), Verbose=false));
    fac = R.history(end).factor;
    tc.verifyEqual(R.meta.options.kind, 'batch');
    tc.verifyEqual(numel(fac), numel(l), 'one factor per observation, batches in order');
    tc.verifyLessThan(fac([7 22]), 0.2, 'both planted errors are down-weighted');
    Pw = diag(p(:) .* fac(:));  xd = (A.'*Pw*A) \ (A.'*Pw*l);
    % tol: 1e-10 relative -- the same weights on both sides; what is left is
    % the QR solve against the normal-equation closed form on a 30 x 3 design
    % with cond ~ 3, i.e. rounding
    tc.verifyEqual(R.x, xd, 'RelTol', 1e-10);
end

function testArmaAndPsdLawsAgreeNoWorseUnderIrls(tc)
% categories: metamorphic, precision
    % Issue #110 test 5. The same process declared two ways -- the ARMA
    % itself, and the PSD evaluated from it -- is NOT the same model: 'psd'
    % is the SAMPLED spectrum, piecewise linear between 4096 ordinates (and
    % whitened by Noise.Whitener: exact Toeplitz by default since WP-125, an
    % AR(30) fit before), so the two disagree at that level whatever the
    % estimator. The issue asks for agreement
    % "within the AR-order tolerance", which is a number this session cannot
    % know; so the point MEASURES that tolerance as the GLS gap between the
    % two declarations and asks the only question IRLS can answer: does
    % re-weighting AMPLIFY it? It must not, because the factors are a
    % Lipschitz function of residuals that differ only by that gap.
    [dArma, cfg] = declaration();
    dPsd = withPsdLaw(dArma, cfg);
    xA = lsqvce.fit(dArma, Verbose=false).x;
    xP = lsqvce.fit(dPsd,  Verbose=false).x;
    ropt = lsqvce.Options(Robust="huber", Verbose=false);
    rA = lsqvce.fit(dArma, ropt).x;
    rP = lsqvce.fit(dPsd,  ropt).x;
    gapGls  = norm(xA - xP);
    gapIrls = norm(rA - rP);
    % the premise: the two declarations ARE the same process to the AR order,
    % so the GLS gap is small against the estimate itself. Without this the
    % ratio below could pass on two unrelated fits.
    tc.verifyLessThan(gapGls, 1e-2 * norm(xA), 'the psd declaration must model the same process');
    tc.verifyLessThan(gapIrls, 10 * gapGls, ...
        'IRLS must not amplify the difference between two spellings of one law');
end

function testWhitenerReachesADeclaredPsdLaw(tc)
% categories: metamorphic, contract
    % WP-126 (finding F3 of the WP-125 review): Noise.Whitener governs a
    % DECLARED psd law end to end -- Options -> adjustBlock / adjustRobust
    % 'PsdMode' -> the loaders' UserOpt -> covModel. test_options only checks
    % the name-value pair; here the route itself: the same declaration fitted
    % with "toeplitz" and with "ar" must give DIFFERENT estimates (a broken
    % route makes them bitwise equal), both of the same process (the gap is
    % the AR(30) fit, small against the estimate), with plain GLS and IRLS.
    [dArma, cfg] = declaration();
    dPsd = withPsdLaw(dArma, cfg);
    xT = lsqvce.fit(dPsd, lsqvce.Options(Noise=struct('Whitener', "toeplitz"), Verbose=false)).x;
    xR = lsqvce.fit(dPsd, lsqvce.Options(Noise=struct('Whitener', "ar"), Verbose=false)).x;
    tc.verifyNotEqual(xT, xR, 'the Whitener must reach the declared law');
    % tol: the same premise as testArmaAndPsdLawsAgreeNoWorseUnderIrls -- two
    % whiteners of one process differ by the AR(30) fit, well under 1 %
    tc.verifyLessThan(norm(xT - xR), 1e-2 * norm(xT), 'both whiten the same process');
    rT = lsqvce.fit(dPsd, lsqvce.Options(Robust="huber", Noise=struct('Whitener', "toeplitz"), Verbose=false)).x;
    rR = lsqvce.fit(dPsd, lsqvce.Options(Robust="huber", Noise=struct('Whitener', "ar"), Verbose=false)).x;
    tc.verifyNotEqual(rT, rR, 'the Whitener must reach the declared law under IRLS too');
end

function testWhitenerIsCaseFreeOnTheNoisePath(tc)
% categories: metamorphic, contract
    % WP-137 (code review of WP-135, its lane-4 low "Noise.Whitener unvalidated
    % on the noise path"): Options checks the Whitener in any case and its
    % engineArgs lower it, but fit's joint path passed char(Noise.Whitener)
    % raw to jointModelOpt, whose 'PsdMode' takes lowercase only -- "AR" ran
    % under Estimate "none" and died on pass 1 under "noise". The two spellings
    % are now one run, bitwise.
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');   % two passes on purpose (gate rule 7)
    [dArma, cfg] = declaration();
    dPsd = withPsdLaw(dArma, cfg);
    xLo = lsqvce.fit(dPsd, lsqvce.Options(Estimate="noise", MaxPasses=2, Seed=3, ...
                     Noise=struct('Grid', 16, 'Whitener', "ar"), Verbose=false)).x;
    xUp = lsqvce.fit(dPsd, lsqvce.Options(Estimate="noise", MaxPasses=2, Seed=3, ...
                     Noise=struct('Grid', 16, 'Whitener', "AR"), Verbose=false)).x;
    tc.verifyEqual(xUp, xLo, 'the Whitener''s case does not matter on the noise path');
end

function testNoiseEstimatesADeclaredPhysicalPsdLaw(tc)
% categories: contract, robustness
    % WP-131: under Estimate="noise" a declared 'psd' law with a physical
    % frequency axis f (covModel's own EXAMPLE form) takes policy "law"; the
    % estimator's override brings fNorm and a Grid-bin S, but the merge kept
    % the declared f, and covModel takes f whenever it is present: the 4096
    % ordinates of f met the 16 of S and pass 1 died in an unnamed numel
    % error (at equal lengths: the law on the wrong axis, r 0.2x at dt = 5).
    % The override's axis now replaces the declared one (applyCovOverride).
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');   % three passes on purpose (gate rule 7)
    [dArma, cfg] = declaration();
    dPsd = withPsdLaw(dArma, cfg);
    o = lsqvce.Options(Estimate="noise", MaxPasses=3, Noise=struct('Grid', 16), Verbose=false);
    R = lsqvce.fit(dPsd, o);                          % errored on pass 1 before
    tc.verifyEqual(R.stopReason == "converged" || R.stopReason == "maxPasses", true);
    xA = lsqvce.fit(dArma, Verbose=false).x;
    % tol: 5 % -- the premise of testArmaAndPsdLawsAgreeNoWorseUnderIrls (1 %
    % between two fixed spellings of one law), widened because here the law
    % is RE-ESTIMATED on a 16-bin grid over three passes; a law on the wrong
    % axis (the equal-length failure) whitens a different process
    tc.verifyLessThan(norm(R.x - xA), 5e-2 * norm(xA), 'the estimated law models the same process');
end

function testStreamNamedLikeAPriorIsRefused(tc)
% categories: contract
    % WP-131 (user decision): 'prior_*' names the prior pseudo-streams. A data
    % stream so named was taken for one: no arcs were read for it, yet it was
    % listed in R.noise.stream -- data dropped without a word. It is refused by
    % makeLoader (relabelled under lsqvce:fit, the original as the cause).
    d = declaration();
    d.systems(1).covariance.prior_pos = d.systems(1).covariance.pos;
    try
        lsqvce.fit(d, Verbose=false);
        tc.verifyFail('a stream named prior_* must be refused');
    catch e
        tc.verifyEqual(e.identifier, 'lsqvce:fit:priorType');
        tc.verifyEqual(e.cause{1}.identifier, 'lsqvce:makeLoader:priorType');
    end
end

function d = withPsdLaw(d, cfg)
%WITHPSDLAW  The declaration with 'pos' declared as a PSD instead of the ARMA.
%   Only the DECLARED law changes: buildDay still closes over the ARMA cfg,
%   so both runs see the same observations to the last bit. The conversion is
%   the one lsqvce.engine.covModel's own EXAMPLE gives.
    a  = cfg.cov.pos.ar;   dt = cfg.dt;
    f  = linspace(0, 0.5/dt, 4096).';
    H  = 1 ./ (exp(-1i*2*pi*f*dt*(0:numel(a)-1)) * a(:));
    d.systems(1).covariance.pos = struct('type', 'psd', 'f', f, ...
                                         'S', cfg.cov.pos.sigma^2 * dt * 2 * abs(H).^2);
end

function testAnAllRejectingInfluenceRaisesTheCollapseError(tc)
% categories: contract, robustness
    % Issue #110 test 6, the half that is reachable -- and HOW it is reached
    % is the result. adjustRobust refuses a pass when nAlive < nc + 1, but
    % the scale is the normalised MAD, so at least half the sample survives
    % any redescending influence at its default tuning
    % (test_metamorphic>testTheMadScaleLeavesHalfTheSampleAlive proves the
    % bound). On 3600 observations against nc = 6 the guard therefore CANNOT
    % fire at default tuning, however discrepant the data. Driving the
    % tuning to reject essentially everything is the only route to it on a
    % real-sized problem, and that is what this point does.
    %
    % THE TUNING IS ARITHMETIC, NOT A ROUND NUMBER, and the first attempt got
    % it wrong in a way worth recording. igg3 keeps a row alive iff
    % |w| <= k1, and for standardised residuals P(|w| <= k1) ~ 2*k1*phi(0),
    % so the survivors number about 3600 * 2*k1*0.3989:
    %
    %     k1 = 2e-3   ->  E[nAlive] = 5.7   against the threshold 7  -- A COIN
    %                     FLIP, and CI flipped it: 99.84 % rejected was not
    %                     enough, which is this WP's finding biting its own test
    %     k1 = 2e-8   ->  E[nAlive] = 6e-5  -- the guard fires, deterministically
    d = declaration();
    opts = lsqvce.Options(Robust="igg3", RobustOpts=struct('Tuning', [1e-8 2e-8]), ...
                          Verbose=false);
    tc.verifyEqual(idOf(@() lsqvce.fit(d, opts)), 'lsqvce:adjustRobust:collapsed');
end

% --- issue #110 test 6 / WP-110: the degenerate arc, carried by its prior ---
% WP-109 could not write this: declaration()'s locals are SYSTEM-scoped, one
% set for all twelve arcs, so rejecting one arc's data leaves nothing to
% carry. This is the arc-scoped fixture it needed. scopedCfg() had the arc
% scope already but declared no priors, and 'arc' rows are only populated
% under a declared SCHEDULE (scopedBuildDay: info.arcId) -- which is why
% s.arcs is the id vector [1 2] here and not the count 2.

function testADegenerateArcRejectsAWholeArcWorthOfRows(tc)
% categories: robustness, contract
    % The PRECONDITION of the point below, made a claim of its own so the
    % failure count says which half broke. Day 1's arc 1 is wrecked with an
    % alternating +-500 that no column of its design can absorb, so a
    % redescending influence must reject the whole arc: 120 rows of the 480.
    % 25 % contamination, comfortably below any breakdown point.
    [d, cfg, colOf, e] = arcScopedPriorFixture();
    d.systems(1).buildDay = @(day, ot, common, local, po, info) ...
        wreckDayOne(scopedBuildDay(cfg, day, ot, info), day, 500);
    d.systems(1).priors = struct('L', e(colOf(1, 1)), 'x0', 0.25, 'sigma', 0.1, ...
                                 'from', 1, 'name', 'arcPin', 'estimate', "none");
    R = lsqvce.fit(d, lsqvce.Options(Robust="bisquare", Verbose=false));
    tc.verifyEqual(R.stopReason, "converged");
    tc.verifyGreaterThanOrEqual(R.history(end).nZero, cfg.mArc, ...
        'a whole arc of unfittable rows must be zeroed, not merely down-weighted');
end

function testADegenerateArcIsCarriedByItsLocalPrior(tc)
% categories: robustness, reference
    % ISSUE #110 TEST 6. With every data row of one arc rejected, that arc's
    % local parameter has exactly one equation left -- its prior -- and must
    % sit on it. The control is the SAME wrecked declaration solved by plain
    % GLS, which does not reject anything and is dragged far off.
    %
    % The prior is deliberately LOOSE (sigma 0.1 against a data weight some
    % five times larger), the lesson WP-108 paid for: a tight prior would pin
    % the local whether or not the data were rejected, and the point would
    % claim nothing. On CLEAN data this prior leaves the local near 0, not
    % near its x0.
    [d, cfg, colOf, e] = arcScopedPriorFixture();
    c = colOf(1, 1);                                   % day 1, arc 1: the starved column
    d.systems(1).buildDay = @(day, ot, common, local, po, info) ...
        wreckDayOne(scopedBuildDay(cfg, day, ot, info), day, 500);
    d.systems(1).priors = struct('L', e(c), 'x0', 0.25, 'sigma', 0.1, ...
                                 'from', 1, 'name', 'arcPin', 'estimate', "none");
    yRob = lsqvce.fit(d, lsqvce.Options(Robust="bisquare", Verbose=false)) ...
               .residuals('Verbose', false).sys(1).y;
    yGls = lsqvce.fit(d, Verbose=false).residuals('Verbose', false).sys(1).y;
    % tol: 1e-6 absolute -- with its data rejected the column's ONLY equation
    % is the prior, so the solve returns x0 to rounding, not to sigma
    tc.verifyEqual(yRob(c), 0.25, 'AbsTol', 1e-6, ...
        'the arc is carried by its local prior alone');
    tc.verifyGreaterThan(abs(yGls(c) - 0.25), 3 * 0.1, ...
        'and the control: GLS keeps the wrecked rows and is dragged off the prior');
end

function [data, cfg, colOf, e] = arcScopedPriorFixture()
%ARCSCOPEDPRIORFIXTURE  scopedCfg() as a fit declaration with ARC-scoped
%   locals, a declared schedule (without one the 'arc' columns are never
%   populated) and the column accessor a prior group's L needs.
%   2 days x 2 arcs x 120 rows; nc = 2 common, nl = 2*2 arc rows + 1 system.
    cfg  = scopedCfg();
    data = struct();
    data.commonParams = [{'common', 'g1'}; {'common', 'g2'}];
    data.days = 1:2;
    s = struct();
    s.buildDay    = @(d, ot, common, local, po, info) scopedBuildDay(cfg, d, ot, info);
    s.covariance  = cfg.cov;
    s.localParams = [{'arc', 'ic'}; {'system', 'scale'}];
    s.arcs        = [1 2];                       % the SCHEDULE, not the count
    data.systems  = s;
    lay   = lsqvce.layout(data);
    nl    = lay.systems(1).nl;
    colOf = @(day, arc) lay.systems(1).col(1, day, arc);
    e     = @(c) double((1:nl) == c);
end

function D = wreckDayOne(D, day, amp)
%WRECKDAYONE  Make every row of day 1's FIRST arc unfittable, and no other.
%   An alternating +-amp: no single design column can absorb it, so the rows
%   stay discrepant however the local parameter moves. Only day 1, because
%   day 2's arc 1 carries no prior -- starving it too would leave its column
%   with no equation at all, which is a rank defect, not a robustness test.
    if day ~= 1 || isempty(D.l), return, end
    m = numel(D.l{1});
    D.l{1} = D.l{1} + amp * (-1).^(1:m).';
end

function testSignatureRejectsNonStructAndNonOptions(tc)
% categories: contract
    seedFor(tc, 584);   % R7-45: drew from the global stream unseeded (suite-order dependent)
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
    % WP-154: the direct call takes the smoothing fit resolves -- the Options
    % default is 1, the engine's own default 0
    out = lsqvce.engine.adjustVCEjoint(accFcn, 'NumLags', nLag, 'MaxIter', 2, 'Probes', 20, ...
                                'Seed', 7, 'Verbose', false, 'Smooth', opts.Noise.Smooth);
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

function testDiagnosticSourceDoesNotStoreTheRows(tc)
% categories: contract, robustness | tags: diagnostics
    % WP-86 part 1 (R7-24): diagnoseAccum stored every whitened batch at full
    % width -- n x (nc + sum nl) doubles kept alive in R.batchFcn, ~228 GB for
    % one GRAVFIRE-scale month. The source now keeps only each batch's ADDRESS
    % and rebuilds it on demand. Pinned two ways: what the closure holds (minus
    % the accumulator's stores and the one-arc cache) does not grow with the
    % number of rows, and a batch rebuilt after a cache switch is bitwise the
    % one first served (the loader contract the rebuild relies on).
    small = fitWithArcLength(300);
    large = fitWithArcLength(3000);                  % 10x the rows, same arcs
    tc.assertEqual(large.nBatch, small.nBatch, 'same arcs, one batch each (fixture premise)');
    tc.verifyGreaterThan(large.diag.n, 9 * small.diag.n, 'the large fit really has ~10x the rows');
    tc.verifyEqual(heldBytes(large), heldBytes(small), ...
        'the diagnostic source holds addresses, not rows: its size does not follow n');
    [A1, l1] = large.batchFcn(1);
    [~, ~] = large.batchFcn(large.nBatch - 1);       % another arc: the cache switches
    [A2, l2] = large.batchFcn(1);                    % rebuilt from the loader
    tc.verifyEqual(A2, A1, 'a rebuilt batch is the batch first served (exact)');
    tc.verifyEqual(l2, l1);
end

function R = fitWithArcLength(mArc)
%FITWITHARCLENGTH  test_fit's canonical declaration with a chosen arc length.
    cfg = buildCfg();  cfg.mArc = mArc;
    data = struct();
    data.commonParams = [repmat({'common'}, cfg.nc, 1), num2cell((1:cfg.nc).')];
    data.days = 1:3;
    s = struct();
    s.buildDay = @(d, ot) buildDayImpl(cfg, d, ot);
    s.covariance = cfg.cov;
    s.localParams = [repmat({'system'}, cfg.nl, 1), num2cell((1:cfg.nl).')];
    s.arcs = cfg.arcsPerDay;
    data.systems = s;
    R = lsqvce.fit(data, Verbose=false);
end

function b = heldBytes(R)
%HELDBYTES  Bytes the diagnostic source holds beyond the accumulator's own
%   stores and its one-arc cache (both independent of what is being tested).
    ws = functions(R.batchFcn).workspace{1};
    src = rmfield(ws.src, {'stores', 'cache'}); %#ok<NASGU> measured by whos
    w = whos('src');
    b = w.bytes;
end

function testResultsCarryTheVerdictMethods(tc)
% categories: contract | tags: diagnostics
    % WP-86 part 4: everything the verdict tells the user to do is a method on
    % R (lsqvce.engine.verdict behind it) -- systematic, snoop, reliability,
    % ellipse. reliability's quadrants are the report's D.quad counts.
    [data, ~] = declaration();
    R = lsqvce.fit(data, Verbose=false);
    T = R.reliability();
    tc.verifyEqual(height(T), R.diag.n, 'one row per observation');
    tc.verifyEqual(T.r, R.diag.rNum(:));
    tc.verifyEqual(T.w, R.diag.w(:));
    tc.verifyTrue(all(ismember(T.quadrant, ["cleanOutlier", "ambiguous", "verified", "untested"])));
    % the quadrant counts are the report's (D.quad)
    tc.verifyEqual(nnz(T.quadrant == "verified"), R.diag.quad.verified);
    tc.verifyEqual(nnz(T.quadrant == "untested"), R.diag.quad.untested);
    S = R.systematic('Verbose', false);
    tc.verifyTrue(isstruct(S));
    E = R.ellipse(1, 2);
    tc.verifyGreaterThanOrEqual(E.semiMajor, E.semiMinor);
    [D2, REJ] = R.snoop('Verbose', false);
    tc.verifyTrue(isstruct(D2) && isfield(REJ, 'idx'));
end

function testReliabilityNeedsThePerObservationDiagnostics(tc)
% categories: contract | tags: diagnostics
    [data, ~] = declaration();
    R = lsqvce.fit(data, lsqvce.Options(Verbose=false, Diagnostics=struct('Store', "none")));
    tc.verifyError(@() R.reliability(), 'lsqvce:Results:noPerObservation');
end

function testCompareCountsTheEliminatedLocals(tc)
% categories: reference, precision
    % R7-19: compare built df2 = n - numel(R.x) -- the COMMONS only -- while the
    % eliminated locals are parameters of the fit too (DEVELOPMENT, "the
    % redundancy ... counting the eliminated parameters"). s0^2 of the full
    % model was then too small and F too large: 7.4 % false positives at
    % alpha 0.05 (R7 mirror m5_irls_stats.py, 2000 trials). Reference: the
    % closed-form F with r = n - nc - nl, and r is the fit's own redundancy.
    [data, cfg] = declaration();
    Rfull = lsqvce.fit(data, Verbose=false);
    red = data;
    red.commonParams = data.commonParams(1:end-1, :);
    red.systems(1).buildDay = @(d, ot) dropLastCommon(buildDayImpl(cfg, d, ot));
    Rred = lsqvce.fit(red, Verbose=false);
    F = Rfull.compare(Rred);
    n = Rfull.diag.n;
    tc.verifyEqual(F.df2, n - cfg.nc - cfg.nl, 'df2 counts the eliminated locals');
    tc.verifyEqual(F.df2, Rfull.diag.r, 'and is the fit''s own redundancy');
    tc.verifyEqual(F.df1, 1);
    Fref = (Rred.sol.vtv - Rfull.sol.vtv) / (Rfull.sol.vtv / (n - cfg.nc - cfg.nl));
    % tol: 1e-10 relative -- the same two residual sums in a closed form
    tc.verifyEqual(F.F, Fref, 'RelTol', 1e-10);
end

function D = dropLastCommon(D)
%DROPLASTCOMMON  The nested model: the last common column removed from every arc.
    for a = 1:numel(D.A), D.A{a} = D.A{a}(:, 1:end-1); end
end

function testUnseededJointRunRecordsAReplayableSeed(tc)
% categories: contract, reference
    % R7-37: an unseeded joint run drew a fresh private seed every pass and
    % recorded none (R.meta.seeds.jointSeed was []), so it could not be
    % replayed -- against the "seeds in provenance" contract. The driver now
    % pins one seed up front and records it; feeding it back replays the run.
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');   % MaxPasses 2 stops early by design
    [data, cfg] = declaration();
    nLag = cfg.mArc;
    base = {'Estimate', "noise", 'Noise', struct('Grid', nLag), 'MaxPasses', 2, ...
            'Probes', 20, 'Verbose', false};
    R1 = lsqvce.fit(data, lsqvce.Options(base{:}));
    s  = R1.meta.seeds.jointSeed;
    tc.assertTrue(isscalar(s) && isfinite(s), 'an unseeded run records the seed it used');
    R2 = lsqvce.fit(data, lsqvce.Options(base{:}, 'Seed', s));
    % tol: exact -- the same seed fixes the same probes on the same code path
    tc.verifyEqual(R2.x, R1.x, 'replaying the recorded seed reproduces the run');
    tc.verifyEqual(R2.noise.C, R1.noise.C);
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
    %
    % R7-29 (2026-10-10): the test passed Parallel= TWICE. set.Parallel merges
    % onto the DEFAULTS, so the second struct reset Mode to "serial", the
    % accumulator never took the distributed path, and the fused solve was
    % compared with itself. One struct now, and the recorded Mode is
    % asserted. The distributed path needs nDays >= Threshold, and so does
    % mapreduce's parfor, so with the PCT present this opens a pool; a pool
    % opened HERE is closed again in teardown (validation Part 10: a test
    % does not change the machine it measures). Without the PCT mapreduce
    % runs the same distributed tasks serially and says so.
    if exist('gcp', 'file') == 2 && license('test', 'Distrib_Computing_Toolbox')
        if isempty(gcp('nocreate'))
            tc.addTeardown(@() delete(gcp('nocreate')));
        end
    else
        muteWarn(tc, 'lsqvce:mapreduce:serialFallback');   % no PCT: the fallback is the expected path
    end
    [data, ~] = declaration();
    R1 = lsqvce.fit(data, lsqvce.Options(Verbose=false));
    R2 = lsqvce.fit(data, lsqvce.Options(Verbose=false, Parallel=struct('Mode', "parfor", ...
                                         'DistributeDays', true, 'Threshold', 2)));
    % tol: 1e-11 -- serial-vs-distributed reduction order (test_parallel
    % testDistributeDaysEqualsFusedSerial's measured bound); bitwise NOT promised
    tc.verifyEqual(R2.x,      R1.x,      'RelTol', 1e-11);
    tc.verifyEqual(R2.sigma0, R1.sigma0, 'RelTol', 1e-11);
    % adjustBlock records its RESOLVED option struct (lowerCamel, WP-44)
    % R5 M5 5b: every kind records the USER's resolved Options as provenance
    tc.verifyEqual(string(R2.meta.options.Parallel.Mode), "parfor", 'one Parallel struct, Mode not reset (R7-29)');
    tc.verifyEqual(R2.meta.options.Parallel.DistributeDays, true, 'the policy is recorded in the provenance');
    tc.verifyEqual(R2.meta.options.Parallel.Threshold, 2);
    tc.verifyEqual(R1.meta.options.Parallel.DistributeDays, false);
    % the same knobs are accepted by the VCE engine's constructor (signature)
end

%% ================================================== WP-67: the user-path demo
function testDemoFitRuns(tc)
% categories: contract, robustness | tags: diagnostics
    % R5 M4 4b: the demo's two-pass vce run is the joint estimator now and
    % says honestly that two passes did not settle. WP-145: the demo mutes
    % that warning itself (and restores it), so it cannot reach this test --
    % the mute that stood here hid nothing, and muteWarn is not for that
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
    % the estimate: x = argmin sum fac_i p_i v_i^2 with fac the weights the
    % returned solve used (WP-132: history(end).factor)
    Pw = diag(p(:) .* fac(:));  xd = (A.'*Pw*A) \ (A.'*Pw*l);
    % tol: 1e-10 relative -- the same weights on both sides; QR against the
    % normal-equation closed form on a cond ~ 3 design, i.e. rounding
    tc.verifyEqual(R.x, xd, 'RelTol', 1e-10);
    % and the factors are Huber's on the standardised whitened residuals
    vw = sqrt(p(:)) .* (A*R.x - l);
    facDef = lsqvce.kernel.robustWeights(vw, 'huber', []);
    % tol: 1e-6 absolute -- the factors were derived from the PREVIOUS
    % estimate; at Tol 1e-12 on x the two differ by the influence's Lipschitz
    % constant times the last step, well below 1e-6
    tc.verifyEqual(fac, facDef, 'AbsTol', 1e-6);
    tc.verifyEqual(R.meta.options.kind, 'matrix');
    tc.verifyNotEmpty(R.acc, 'kind B robust carries the accumulator too');
end
function testRobustFactorsAreTheReturnedSolvesWeights(tc)
% categories: contract, metamorphic
    % WP-132 (review low 24): R.history(end).factor was the NEXT pass's factor
    % (derived from the returned solve's residuals), so at a maxPasses exit the
    % returned x was NOT the weighted solution on the returned factors -- one
    % IRLS step apart. Now they are the weights the returned solve used: one
    % pass returns the unweighted solve with all factors 1, three passes the
    % weighted solution on the factors pass 2 derived.
    [A, l, p] = matrixFixture();
    l(7) = l(7) + 25;                                   % gross error
    data = struct('A', A, 'l', l, 'p', p);
    muteWarn(tc, 'lsqvce:adjustRobust:notConverged');   % both runs stop early by design
    for maxIter = [1 3]
        R = lsqvce.fit(data, lsqvce.Options(Robust="huber", ...
                       RobustOpts=struct('MaxIter', maxIter), Verbose=false));
        fac = R.history(end).factor;
        Pw = diag(p(:) .* fac(:));  xd = (A.'*Pw*A) \ (A.'*Pw*l);
        % tol: 1e-10 relative -- the same weights on both sides; QR against the
        % normal-equation closed form on a cond ~ 3 design, i.e. rounding
        tc.verifyEqual(R.x, xd, 'RelTol', 1e-10, ...
            sprintf('MaxIter %d: x must be the weighted solution on the returned factors', maxIter));
        if maxIter == 1
            tc.verifyEqual(R.stopReason, "maxPasses");
            tc.verifyEqual(fac, ones(numel(l), 1), 'pass 1 solves unweighted');
        else
            tc.verifyLessThan(fac(7), 0.5, 'by pass 3 the planted error is down-weighted');
        end
    end
end
function testZeroWeightRowsAreRejectedNotObservations(tc)
% categories: reference, contract
    % WP-133 (review 3 finding 17): a zero weight deletes its row's
    % information, but the row stayed in the count -- sigma0 0.753 against
    % 0.949 with 10 of 30 rows at weight 0 -- and in the verdict as a
    % "verified" observation. It keeps its index now and is "rejected": out of
    % the dof, sum(r), the tests and Bonferroni's n. (A negative weight is
    % refused by name: test_errorDoctrine.)
    [A, l, p] = matrixFixture();
    z = [4; 11; 23];  p(z) = 0;  l(z) = l(z) + 100;     % weight 0: their values must not matter
    R = lsqvce.fit(struct('A', A, 'l', l, 'p', p), lsqvce.Options(Verbose=false));
    k = true(numel(l), 1);  k(z) = false;
    P = diag(p(k));  xd = (A(k,:).'*P*A(k,:)) \ (A(k,:).'*P*l(k));  vd = A(k,:)*xd - l(k);
    dof = nnz(k) - size(A, 2);
    % tol: 1e-12 relative -- as testMatrixKindReproducesAdjust (block QR vs
    % dense normal equations, cond ~ 1e1)
    tc.verifyEqual(R.x, xd, 'RelTol', 1e-12, 'the zero-weight rows carry no information');
    tc.verifyEqual(R.sol.dof, dof, 'the zero-weight rows leave the dof');
    tc.verifyEqual(R.sigma0, sqrt((vd.'*P*vd) / dof), 'RelTol', 1e-12);   % tol: as above
    D = R.diag;
    tc.verifyEqual(D.r, dof, 'the verdict counts the same redundancy');
    tc.verifyEqual(D.excludedIdx(:), z, 'their indices are kept');
    tc.verifyEqual(D.quad.rejected, numel(z));
    tc.verifyEqual(D.nTests, nnz(k), 'not tested, not in the Bonferroni n');
    tc.verifyEqual(D.rNum(z), zeros(numel(z), 1));
    % tol: 1e-10 -- sum(r) against n - excluded - u, 30 leverages on cond ~ 1e1
    tc.verifyLessThan(D.sumRRelErr, 1e-10, 'sum(r) counts the live rows');
    T = R.reliability();
    tc.verifyEqual(T.quadrant(z), repmat("rejected", numel(z), 1));
    tc.verifyFalse(any(T.quadrant(k) == "rejected"));
    res = tc.verifyWarningFree(@() R.residuals('Verbose', false));
    tc.verifyEqual(res.dof, dof, 'the residual pass counts the same dof');
end

function testRobustRejectionsLeaveTheDofAndAreLabelled(tc)
% categories: robustness, metamorphic
    % WP-133 (review 3 findings 2 and 14). (2) adjustRobust handed the
    % verdict the SAME containers.Map it then overwrote with the next pass's
    % factors: at a maxPasses exit the diagnostics replayed weights the
    % returned solve never used (vPv against sol.vtv off). (14) A rejected row
    % (factor 0) stayed in the dof and the verdict called it "verified". Now
    % the replay is the returned solve's (vPv equal), the rejected rows are
    % "rejected", out of the dof and sum(r), and the residual pass agrees.
    [A, l, p] = matrixFixture();
    l(7) = l(7) + 25;                                   % one gross error: bisquare zeroes it
    data = struct('A', A, 'l', l, 'p', p);
    muteWarn(tc, 'lsqvce:adjustRobust:notConverged');   % MaxIter 2 stops early by design
    for maxIter = [2 30]
        R = lsqvce.fit(data, lsqvce.Options(Robust="bisquare", ...
                       RobustOpts=struct('MaxIter', maxIter), Verbose=false));
        fac = R.history(end).factor;  z = find(fac == 0);
        tc.verifyTrue(ismember(7, z), sprintf('MaxIter %d: the blunder is rejected', maxIter));
        D = R.diag;
        % tol: 1e-9 relative -- the replay drops the rejected rows (p = 0) the
        % solve kept at weight 1e-12: (1e-6 * 25)^2 against v'v ~ 30
        tc.verifyLessThan(D.vPvRelDiff, 1e-9, sprintf('MaxIter %d: the replay is the returned solve''s', maxIter));
        tc.verifyEqual(R.sol.dof, numel(l) - numel(z) - size(A, 2), 'the rejected rows leave the dof');
        tc.verifyEqual(D.excludedIdx(:), z(:), 'and are named');
        tc.verifyEqual(D.quad.rejected, numel(z));
        % tol: the report's own integrity line; the 1e-12-weight rows move sum(r) by ~1e-12
        tc.verifyLessThan(D.sumRRelErr, 1e-9);
        res = tc.verifyWarningFree(@() R.residuals('Verbose', false));
        % tol: 1e-8 relative -- two routes to the same weighted v'Pv over the same dof
        tc.verifyEqual(res.sigma0, R.sigma0, 'RelTol', 1e-8);
    end
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

function testGaussNewtonDampsAndKeepsTheBestPoint(tc)
% categories: robustness, contract
    % R7 WP-113 (review R7-14): the damping works on the COMMON block, the
    % descent check runs from pass 1, and an unconverged run returns the best
    % ACCEPTED point. Exponential fit y = a exp(b t) (+ a non-zero residual),
    % from two far starts where the shipped loop re-solved a rejected step
    % unchanged until MaxIter and reported x + that step (F 36.8 / 17.5
    % against the best accepted 4.9 / 5.5; mirror mirrors/wp113).
    t = linspace(0, 4, 30).';  y = 3*exp(-1.3*t) + 0.05*sin(7*t);
    bf = @(k, x) deal([exp(x(2)*t), x(1)*t.*exp(x(2)*t)], y - x(1)*exp(x(2)*t), []);
    F  = @(x) sum((y - x(1)*exp(x(2)*t)).^2);
    ref = lsqvce.fit(struct('batchFcn', bf, 'numBatches', 1, 'x0', [1; 0.5]), ...
                     lsqvce.Options(MaxPasses=30, Tolerance=1e-10, Verbose=false));
    for x0 = {[10; -5], [0.1; 1]}
        R = lsqvce.fit(struct('batchFcn', bf, 'numBatches', 1, 'x0', x0{1}), ...
                       lsqvce.Options(MaxPasses=30, Tolerance=1e-10, Verbose=false));
        tc.verifyEqual(R.stopReason, "converged", sprintf('start [%g %g] converges', x0{1}));
        % tol: 1e-9 relative -- both runs stop at TolX 1e-10 formal sd on the same optimum
        tc.verifyEqual(F(R.x), F(ref.x), 'RelTol', 1e-9, 'the far start reaches the optimum');
        ok = [R.history.accepted];  Fs = [R.history(ok).vPv];
        tc.verifyTrue(all(diff([F(x0{1}), Fs]) <= 0), 'every accepted step goes downhill, the first included');
    end
    % a budget too small for the far start: maxPasses, and R.x is the best accepted point
    muteWarn(tc, 'lsqvce:adjustNonlinear:maxIter');
    R = lsqvce.fit(struct('batchFcn', bf, 'numBatches', 1, 'x0', [0.1; 1]), ...
                   lsqvce.Options(MaxPasses=10, Tolerance=1e-10, Verbose=false));
    tc.verifyEqual(R.stopReason, "maxPasses");
    ok = [R.history.accepted];
    % tol: 1e-12 relative -- R.x IS the last accepted iterate, F re-evaluated
    tc.verifyEqual(F(R.x), min([R.history(ok).vPv]), 'RelTol', 1e-12, ...
        'an unconverged run returns the best accepted point, not a rejected step');
end


function testGaussNewtonNeverReturnsAnUnevaluatedPoint(tc)
% categories: robustness, metamorphic
    % WP-127: TolX measures the step in formal-error units; on a plateau the
    % formal errors are huge, so a tiny DAMPED step "converges", and the final
    % UNDAMPED step was applied without evaluating F. Mirror of this loop
    % (the review's, reproduced): a logistic fit from [1 -2 0.1] returned
    % F = 24.7 against an accepted 21.9, from [1 -5 0.1] 49.6 against 28.8,
    % from [0.5 15 0.3] 69.0 against 68.4. Starts whose final linearisation is
    % singular (pivots ~1e-28 deep on the plateau) are refused by name
    % (lsqvce:qr:singularR) and are not used here. Whatever the start, the
    % returned point's F is no worse than the best accepted one (default
    % Tolerance, the setting that exposed it).
    t = linspace(0, 10, 40).';  y = 2 ./ (1 + exp(-(t - 5)/0.7)) + 0.02*sin(3*t);
    g  = @(x) 1 ./ (1 + exp(-(t - x(2))/x(3)));
    bf = @(k, x) deal([g(x), -x(1)*g(x).^2.*exp(-(t - x(2))/x(3))/x(3), ...
                       -x(1)*g(x).^2.*exp(-(t - x(2))/x(3)).*(t - x(2))/x(3)^2], y - x(1)*g(x), []);
    F  = @(x) sum((y - x(1)*g(x)).^2);
    muteWarn(tc, 'lsqvce:adjustNonlinear:maxIter');
    for x0 = {[1; -2; 0.1], [1; -5; 0.1], [0.5; 15; 0.3]}
        R = lsqvce.fit(struct('batchFcn', bf, 'numBatches', 1, 'x0', x0{1}), ...
                       lsqvce.Options(MaxPasses=30, Verbose=false));
        ok = [R.history.accepted];
        Fbest = min([F(x0{1}), R.history(ok).vPv]);
        % tol: 1e-10 relative -- the rounding allowance of the final-step check
        tc.verifyLessThanOrEqual(F(R.x), Fbest * (1 + 1e-10), ...
            sprintf('start [%g %g %g]: R.x is no worse than the best evaluated point', x0{1}));
        % WP-130: sigma0 and the diagnostics are those of R.x, not of the linear
        % prediction at x + dx (sigma0 0.19 against sqrt(F/r) = 0.77 from the
        % first start), and D.x is the estimate, not the final increment.
        % tol: 1e-9 relative -- the same sum of squares, accumulated vs closed form
        tc.verifyEqual(R.sigma0^2 * R.diag.r, F(R.x), 'RelTol', 1e-9, ...
            sprintf('start [%g %g %g]: sigma0^2 r = F(R.x)', x0{1}));
        tc.verifyEqual(R.diag.vPv, F(R.x), 'RelTol', 1e-9, 'the diagnostics replay R.x');
        tc.verifyEqual(R.diag.x(1:3), R.x, 'D.x is the estimate (report, ellipse centre)');
    end
end
function testGaussNewtonResidualsAreThoseOfTheEstimate(tc)
% categories: metamorphic, contract
    % WP-133 (review 3 finding 4): the accumulator a Gauss-Newton fit returns
    % holds the rows linearised AT the estimate, i.e. increments about it, and
    % R.residuals() applied the absolute x to them as if it were the step --
    % sigma0 623x too large and a sigma0Mismatch warning on every nonlinear
    % fit. sol.xLin names the linearisation point; the residual pass now
    % measures the zero step and agrees with R.sigma0.
    rs = RandStream('mt19937ar', 'Seed', 5);
    t  = linspace(0, 1, 40).';
    lobs = exp(0.7*t) - 0.2 + 1e-3*randn(rs, 40, 1);
    bf = @(k, x) deal([t.*exp(x(1)*t), ones(40, 1)], lobs - (exp(x(1)*t) + x(2)), []);
    R = lsqvce.fit(struct('batchFcn', bf, 'numBatches', 1, 'x0', [0.5; 0]), ...
                   lsqvce.Options(MaxPasses=15, Tolerance=1e-9, Verbose=false));
    tc.verifyEqual(R.sol.xLin, R.x, 'the linearisation point is the estimate');
    res = tc.verifyWarningFree(@() R.residuals('Verbose', false));
    % tol: 1e-9 relative -- the same whitened rows at the zero step, summed twice
    tc.verifyEqual(res.sigma0, R.sigma0, 'RelTol', 1e-9);
    w = lobs - (exp(R.x(1)*t) + R.x(2));
    % tol: as above -- the residual pass's mass is F(R.x) itself
    tc.verifyEqual(res.vtv, w.'*w, 'RelTol', 1e-9);
end

function testGaussNewtonRejectsANonFiniteTrialPoint(tc)
% categories: robustness, contract
    % WP-133 (review 3): a trial point where the model is not defined (an
    % exp that overflows on a far step, here x^3 past |x| = 10) stopped the
    % fit with lsqvce:accum:arcNotFinite -- 2 of 40 starts of an exponential
    % decay. It is a REJECTED step now (F = Inf, the damping grows), which is
    % what the damping is for; only a starting point the model cannot
    % evaluate is refused (lsqvce:adjustNonlinear:nonFinite).
    y  = 1 + [0.01; -0.02; 0.015; -0.005; 0.01];        % five looks at x^3, truth 1
    f  = @(x) x.^3 ./ (abs(x) <= 10);                   % Inf past |x| = 10
    bf = @(k, x) deal(3*x^2 ./ (abs(x) <= 10) * ones(5, 1), y - f(x), []);
    R = lsqvce.fit(struct('batchFcn', bf, 'numBatches', 1, 'x0', 0.1), ...
                   lsqvce.Options(MaxPasses=60, Tolerance=1e-10, Verbose=false));
    h = R.history;
    tc.verifyTrue(any(~[h.accepted] & isinf([h.vPv])), ...
        'the first undamped step (to x ~ 33) is rejected with F = Inf, not an abort');
    tc.verifyEqual(R.stopReason, "converged");
    % tol: 1e-8 relative -- the least-squares optimum of x^3 against y is
    % mean(y)^(1/3); TolX 1e-10 formal sd and the final undamped step
    tc.verifyEqual(R.x, mean(y)^(1/3), 'RelTol', 1e-8);
    tc.verifyError(@() lsqvce.fit(struct('batchFcn', bf, 'numBatches', 1, 'x0', 20), ...
                   lsqvce.Options(Verbose=false)), 'lsqvce:adjustNonlinear:nonFinite');
end

function testGaussNewtonPriorIsOnXNotOnTheStep(tc)
% categories: precision, reference
    % R7 WP-113 (review R7-04): a prior K x ~ x0 on a NONLINEAR problem enters
    % every step as K dx ~ x0 - K x. On the increment (the shipped form) the
    % fixed point was not the declared optimum: b reported +0.60 against the
    % regularised +0.30 (sigma 1e-6). Oracle: the stationarity of the DECLARED
    % objective |w|^2_P + |(K x - x0)/sigma|^2, and the tight prior holding.
    rs = RandStream('mt19937ar', 'Seed', 5);
    t  = linspace(0, 1, 40).';
    lobs = exp(0.7*t) - 0.2 + 1e-3*randn(rs, 40, 1);
    p  = 1e6 * ones(40, 1);                          % weights of 1e-3 noise
    bf = @(k, x) deal([t.*exp(x(1)*t), ones(40, 1)], lobs - (exp(x(1)*t) + x(2)), p);
    for c = [0.30 1e-3; 0.30 1e-6; 0 1e-3].'
        data = struct('batchFcn', bf, 'numBatches', 1, 'x0', [0.5; 0], ...
                      'priors', struct('L', [0 1], 'x0', c(1), 'sigma', c(2), 'name', 'b'));
        R = lsqvce.fit(data, lsqvce.Options(MaxPasses=30, Tolerance=1e-10, Verbose=false));
        tc.verifyEqual(R.stopReason, "converged");
        J  = [t.*exp(R.x(1)*t), ones(40, 1)];  w = lobs - (exp(R.x(1)*t) + R.x(2));
        gD = J.' * (p .* w);  gP = [0; 1] * (c(1) - R.x(2)) / c(2)^2;
        % tol: 1e-6 relative to the gradient's parts -- the data and prior pulls
        % cancel at the declared optimum; TolX 1e-10 leaves far less than that
        tc.verifyLessThan(norm(gD + gP), 1e-6 * (norm(gD) + norm(gP) + 1), ...
            sprintf('b ~ %g (sigma %g): the declared objective is stationary at R.x', c(1), c(2)));
        if c(2) < 1e-5
            % tol: 1e-4 absolute -- a 1e-6 prior holds b to its target up to the data's pull
            tc.verifyEqual(R.x(2), c(1), 'AbsTol', 1e-4, 'the tight prior holds');
        end
    end
end

function testDatumCovarianceIsTheMinimumConstraintOne(tc)
% categories: precision, contract
    % R7 WP-114 (R7-05): a datum fixes x for ANY sigma, and now the covariance
    % too -- the minimum-constraint cofactor pinv(N), not s0^2 inv(N + BB'/s^2),
    % whose se grew 0.005 -> 2.48 as the datum sigma went 1 -> 1e3 (review). The
    % parameter cofactor of the diagnostics (D.sigX) follows. And R7-15: with
    % the prior counted, residuals() and redundancy() agree with the solve --
    % no sigma0Mismatch, no identityViolated.
    rs = RandStream('mt19937ar', 'Seed', 5);
    u = 4;  pairs = [1 2; 2 3; 3 4; 4 1; 1 3];
    A = zeros(size(pairs, 1), u);
    for k = 1:size(pairs, 1), A(k, pairs(k, 1)) = -1;  A(k, pairs(k, 2)) = 1; end
    l = A*[0; 1.5; -0.7; 2.2] + 0.01*randn(rs, size(A, 1), 1);
    Np = pinv(A.'*A);
    for sg = [1 1e3 1e6]
        data = struct('A', A, 'l', l, 'priors', struct('L', ones(1, u), 'sigma', sg, 'name', 'datum'));
        R = lsqvce.fit(data, Verbose=false);
        se = sqrt(diag(R.sol.Cxx));
        % tol: 1e-8 relative -- the minimum-constraint cofactor from the re-weighted
        % factor vs the dense pseudo-inverse (mirror: 3e-16 at every sigma)
        tc.verifyEqual(se, sqrt(R.sol.s0sq * diag(Np)), 'RelTol', 1e-8, sprintf('datum sigma %g', sg));
        % tol: as above -- the diagnostics' parameter cofactor is the same one
        tc.verifyEqual(R.diag.sigX(:), sqrt(diag(Np)), 'RelTol', 1e-8, sprintf('D.sigX at datum sigma %g', sg));
        res = tc.verifyWarningFree(@() R.residuals('Verbose', false));
        % tol: 1e-8 relative -- two routes to the same v'Pv over the same dof
        tc.verifyEqual(res.sigma0, R.sigma0, 'RelTol', 1e-8, 'residuals() counts the prior');
        red = tc.verifyWarningFree(@() R.acc.redundancy(R.sol, 'Verbose', false));
        tc.verifyLessThan(red.check, 1e-8, 'sum(r_i) + r_prior = m + nk - u');
    end
    % a REGULARISING prior is not a datum: its covariance stays the solve's
    rs = RandStream('mt19937ar', 'Seed', 6);
    A2 = randn(rs, 40, 6);  l2 = A2 * randn(rs, 6, 1) + randn(rs, 40, 1);
    R = lsqvce.fit(struct('A', A2, 'l', l2, 'priors', struct('L', eye(6), 'sigma', 1, 'name', 'ridge')), Verbose=false);
    % tol: 1e-8 relative -- s0^2 inv(A'A + I) from the factor vs the dense inverse
    tc.verifyEqual(diag(R.sol.Cxx), R.sol.s0sq * diag(inv(A2.'*A2 + eye(6))), 'RelTol', 1e-8);
    res = tc.verifyWarningFree(@() R.residuals('Verbose', false));
    % tol: as above
    tc.verifyEqual(res.sigma0, R.sigma0, 'RelTol', 1e-8, 'the ridge prior is counted');
    red = tc.verifyWarningFree(@() R.acc.redundancy(R.sol, 'Verbose', false));
    tc.verifyLessThan(red.check, 1e-8);
    tc.verifyGreaterThan(red.priorRedundancy, 0);
end

function testDatumGroupBesideARegulariser(tc)
% categories: precision, contract, robustness
    % R7 WP-119 (the residue §12.48 recorded): the datum test runs per declared
    % GROUP. Minimum constraints on four levelled heights beside a ridge on two
    % extra parameters the datum does not touch: before, the stacked prior was
    % not a datum and se(h) grew 0.09 -> 41.6 -> 4.2e4 with the datum sigma
    % (mirror wp119/datum_split.py); now it is s0^2 pinv(N + Ki'Ki) at every
    % sigma, on both backends. A ridge on EVERY parameter, or datum rows inside
    % one partial group, is information and keeps the solve's covariance.
    rs = RandStream('mt19937ar', 'Seed', 7);
    pairs = [1 2; 2 3; 3 4; 4 1; 1 3];
    A = zeros(size(pairs, 1) + 3, 6);
    for k = 1:size(pairs, 1)
        A(k, pairs(k, 1)) = -1;  A(k, pairs(k, 2)) = 1;  A(k, 5) = 0.3*k;  A(k, 6) = (-1)^k;
    end
    A(size(pairs, 1) + 1:end, 5:6) = randn(rs, 3, 2);
    l = A*[0; 1.5; -0.7; 2.2; 0.2; -0.1] + 0.01*randn(rs, size(A, 1), 1);
    datum = [1 1 1 1 0 0];  ridge = [0 0 0 0 1 0; 0 0 0 0 0 1];  sr = 0.5;
    Qref = pinv(A.'*A + ridge.'*ridge / sr^2);
    for bk = ["qr" "normal"]
        for sg = [1 1e3 1e6]
            P = struct('L', {datum, ridge}, 'sigma', {sg, sr}, 'name', {'datum', 'ridge'});
            R = lsqvce.fit(struct('A', A, 'l', l, 'priors', P), Verbose=false, Backend=bk);
            msg = sprintf('%s, datum sigma %g', bk, sg);
            % tol: 1e-8 relative -- the minimum-constraint cofactor of the data plus
            % the ridge from the re-weighted factor vs the dense pseudo-inverse
            % (mirror: 1e-15 at every sigma, both backends)
            tc.verifyEqual(sqrt(diag(R.sol.Cxx)), sqrt(R.sol.s0sq * diag(Qref)), 'RelTol', 1e-8, msg);
            % tol: as above -- the diagnostics' parameter cofactor is the same one
            tc.verifyEqual(R.diag.sigX(:), sqrt(diag(Qref)), 'RelTol', 1e-8, msg);
            % The identities read the SOLVE's cofactor (leverage keeps it, as the
            % datum row must stay absorbed). On the normal backend at datum sigma
            % 1e6 that cofactor comes from N + K'K/s^2 at cond 4.6e12 -- its
            % sum(r) + r_prior is off by 5.6e-6 in the mirror (QR: 0), a property
            % of normal equations, not of the datum part. Checked where it holds.
            if bk == "normal" && sg > 1e3, continue, end
            res = tc.verifyWarningFree(@() R.residuals('Verbose', false));
            % tol: 1e-8 relative -- two routes to the same v'Pv over the same dof
            tc.verifyEqual(res.sigma0, R.sigma0, 'RelTol', 1e-8, msg);
            red = tc.verifyWarningFree(@() R.acc.redundancy(R.sol, 'Verbose', false));
            tc.verifyLessThan(red.check, 1e-8, msg);
        end
    end
    % information, not a datum: the solve's own covariance s0^2 inv(N + K'K)
    info = {struct('L', {datum, eye(6)}, 'sigma', {1, 1}), ...       % a ridge on everything overlaps the datum
            struct('L', [datum; ridge], 'sigma', 1)};                 % ONE group, only partly in the null space
    for k = 1:numel(info)
        K = vertcat(info{k}.L);
        R = lsqvce.fit(struct('A', A, 'l', l, 'priors', info{k}), Verbose=false);
        % tol: 1e-8 relative -- the solve's cofactor from the factor vs the dense inverse
        tc.verifyEqual(diag(R.sol.Cxx), R.sol.s0sq * diag(inv(A.'*A + K.'*K)), 'RelTol', 1e-8, ...
            sprintf('information case %d keeps the solve''s covariance', k));
    end
end

function testDatumSplitOverGroupsIsOneDatum(tc)
% categories: precision, contract, robustness
    % WP-127. (1) A free 2-D distance network whose datum is declared as TWO
    % groups -- translation and rotation. Off-centre, Kt Kr' ~= 0, and WP-119's
    % rule (each group orthogonal to every OTHER group) refused both halves:
    % the covariance fell back to inv(N + K'K/s^2), off pinv(N) by 0.9 at
    % sigma 1 and 9e5 at 1e3 (dev/validation/wp127_datum_split_mirror.py).
    % Blind groups only need orthogonality to the groups the data sees.
    % (2) The datum rows (r = 0 by construction) are constraints, not
    % observations: not uncontrolled, no MDB = Inf, not in the Bonferroni n --
    % every free network's report said SUSPECT.
    rs = RandStream('mt19937ar', 'Seed', 27);
    P  = [20 5; 28 6; 30 13; 24 16; 19 11];                  % not centred; r_i 0.23..0.38
    pairs = nchoosek(1:5, 2);  m = size(pairs, 1);
    A = zeros(m, 10);
    for k = 1:m
        e = P(pairs(k, 2), :) - P(pairs(k, 1), :);  e = e / norm(e);
        A(k, 2*pairs(k, 1) - [1 0]) = -e;  A(k, 2*pairs(k, 2) - [1 0]) = e;
    end
    l  = A * randn(rs, 10, 1) + 0.01 * randn(rs, m, 1);
    Kt = kron(ones(1, 5), eye(2));                           % x and y translation
    Kr = reshape([-P(:, 2).'; P(:, 1).'], 1, []);            % the rotation
    Qref = pinv(A.' * A);
    for bk = ["qr" "normal"]
        for sg = [1 1e3]
            pri = struct('L', {Kt, Kr}, 'sigma', {sg, sg}, 'name', {'trans', 'rot'});
            R = lsqvce.fit(struct('A', A, 'l', l, 'priors', pri), Verbose=false, Backend=bk);
            msg = sprintf('%s, datum sigma %g', bk, sg);
            % tol: 1e-8 relative -- the minimum-constraint cofactor from the
            % re-weighted factor vs the dense pseudo-inverse (mirror: 8e-11 abs)
            tc.verifyEqual(R.diag.sigX(:), sqrt(diag(Qref)), 'RelTol', 1e-8, msg);
            % tol: as above -- the solve's covariance is the same one
            tc.verifyEqual(sqrt(diag(R.sol.Cxx)), sqrt(R.sol.s0sq * diag(Qref)), 'RelTol', 1e-8, msg);
        end
    end
    D = R.diag;
    tc.verifyEqual(D.datumRows(:), (m + 1 : m + 3).', 'the three prior rows are the datum');
    tc.verifyEqual([D.nTests, D.nUncontrolled], [m, 0], 'datum rows are neither tests nor uncontrolled');
    V = R.report('Verbose', false);
    codes = cellfun(@(f) f.code, V.findings, 'UniformOutput', false);
    tc.verifyFalse(any(ismember(codes, {'REL_UNCONTROLLED', 'REL_MDBINF'})), ...
        'a datum is not an uncontrolled observation');
    T = R.reliability();
    tc.verifyEqual(T.quadrant(m + 1 : end), repmat("datum", 3, 1));
end

function testDatumTestIsInvariantToParameterUnits(tc)
% categories: metamorphic, robustness
    % WP-130: the datum test compared |S Kg'| with the norm of the WHOLE
    % system, so a parameter in small units passed as "blind": with a prior
    % group of its own it got the minimum-constraint covariance, sd 0
    % (dev/validation/wp130_datum_scale_mirror.py). Now the product is held
    % against its own rounding bound |S||Kg'|. Metamorphic: rescaling
    % parameter 6 by c (column A(:,6) c, its prior sigma / c) rescales its sd
    % by 1/c and leaves every other sd unchanged; the datum on 1..4 stays one.
    rs = RandStream('mt19937ar', 'Seed', 130);
    pairs = [1 2; 2 3; 3 4; 4 1; 1 3];
    A = zeros(size(pairs, 1) + 3, 6);
    for k = 1:size(pairs, 1)
        A(k, pairs(k, 1)) = -1;  A(k, pairs(k, 2)) = 1;  A(k, 5) = 0.3*k;  A(k, 6) = (-1)^k;
    end
    A(size(pairs, 1) + 1:end, 5:6) = randn(rs, 3, 2);
    l = A*[0; 1.5; -0.7; 2.2; 0.2; -0.1] + 0.01*randn(rs, size(A, 1), 1);
    c = 1e-9;  As = A;  As(:, 6) = c * A(:, 6);
    pri = @(s6) struct('L', {[1 1 1 1 0 0], [0 0 0 0 0 1]}, 'sigma', {1, s6}, 'name', {'datum', 'p6'});
    R1 = lsqvce.fit(struct('A', A,  'l', l, 'priors', pri(0.5)),     Verbose=false);
    R2 = lsqvce.fit(struct('A', As, 'l', l, 'priors', pri(0.5 / c)), Verbose=false);
    sd1 = sqrt(diag(R1.sol.Cxx));  sd2 = sqrt(diag(R2.sol.Cxx));
    % tol: 1e-6 relative -- one QR of a design whose column 6 is 1e-9 of the
    % others, against one of the unscaled design (eps * 1/c is not involved:
    % the scaling is exact in the factor's column)
    tc.verifyEqual(sd2(6) * c, sd1(6), 'RelTol', 1e-6, 'sd(x6) rescales by 1/c (it was 0)');
    % tol: as above
    tc.verifyEqual(sd2(1:5), sd1(1:5), 'RelTol', 1e-6, 'the other parameters do not see the units');
    tc.verifyGreaterThan(sd2(6), 0);
    % tol: as above -- the diagnostics' cofactor from the same factor
    tc.verifyEqual(R2.diag.sigX(6) * c, R1.diag.sigX(6), 'RelTol', 1e-6, 'the diagnostics'' cofactor too');
end

function testSnoopKeepsTheDatumAndTheOriginalIndex(tc)
% categories: precision, contract, robustness
    % WP-128. (1) snoop's clean re-solve re-accumulated a prior as rows of
    % DATA (the setPrior branch never ran): with a datum the cofactor became
    % inv(N + K'K/s^2), sigX 250 against the minimum-constraint 0.164 at
    % sigma 1e3 on this net (mirror). (2) The re-diagnosis numbered the
    % SURVIVING rows while D.rejectedIdx counts the original ones -- after
    % one rejection "observation 399" was 400. (3) A prior row is never a
    % rejection candidate.
    rs = RandStream('mt19937ar', 'Seed', 5);
    pairs = repmat([1 2; 2 3; 3 4; 4 1; 1 3; 2 4], 7, 1);  m = size(pairs, 1);
    A = zeros(m, 4);  A(sub2ind(size(A), (1:m).', pairs(:, 1))) = -1;
    A(sub2ind(size(A), (1:m).', pairs(:, 2))) = 1;
    l = A*[0; 1.5; -0.7; 2.2] + 0.01*randn(rs, m, 1);  l(3) = l(3) + 0.2;   % a 20-sigma blunder
    Ar = A([1:2, 4:m], :);
    R = lsqvce.fit(struct('A', A, 'l', l, 'priors', struct('L', ones(1, 4), 'sigma', 1e3, 'name', 'datum')), ...
                   Verbose=false);
    [D, REJ] = R.snoop('Verbose', false);
    tc.verifyEqual(REJ.idx, 3, 'the blunder, by its original index');
    tc.verifyEqual(D.rejectedIdx, 3);
    % tol: 1e-8 relative -- the minimum-constraint cofactor of the surviving
    % rows from the re-weighted factor vs the dense pseudo-inverse
    tc.verifyEqual(D.sigX(:), sqrt(diag(pinv(Ar.' * Ar))), 'RelTol', 1e-8, ...
        'the clean re-solve keeps the datum (minimum-constraint cofactor)');
    tc.verifyEqual(D.obsIdx(:), [1:2, 4:m, m + 1].', 'survivors, then the prior row, by original index');
    tc.verifyEqual(D.datumRows, m + 1, 'the datum row keeps its original index');
    tc.verifyFalse(ismember(3, D.flagIdx), 'the rejected observation is not flagged again');
end

function testSnoopAfterRejectionWithoutTheReSolve(tc)
% categories: precision, contract, metamorphic
    % WP-130. (1) With Reaccumulate=false snoop returned the INPUT D (the
    % rejected observation still flagged, no .rejectedIdx) and sol.s0sq / n /
    % r / vPv from before the rejection (the blunder still in s0): the 1/r_i
    % identity was computed and dropped. Now D is re-diagnosed on the
    % downdated cofactor over the survivors and sol is after rejection; both
    % agree with the clean re-solve. (2) snoop's Bonferroni n counted the
    % datum row, diagnose's did not: the first-iteration critical value is
    % now R.diag.kCrit exactly.
    rs = RandStream('mt19937ar', 'Seed', 5);
    pairs = repmat([1 2; 2 3; 3 4; 4 1; 1 3; 2 4], 7, 1);  m = size(pairs, 1);
    A = zeros(m, 4);  A(sub2ind(size(A), (1:m).', pairs(:, 1))) = -1;
    A(sub2ind(size(A), (1:m).', pairs(:, 2))) = 1;
    l = A*[0; 1.5; -0.7; 2.2] + 0.01*randn(rs, m, 1);  l(3) = l(3) + 0.2;   % a 20-sigma blunder
    R = lsqvce.fit(struct('A', A, 'l', l, 'priors', struct('L', ones(1, 4), 'sigma', 1e3, 'name', 'datum')), ...
                   Verbose=false);
    [Dt, REJt, solT] = R.snoop('Verbose', false);
    [Df, REJf, solF] = R.snoop('Verbose', false, 'Reaccumulate', false);
    tc.assertEqual([REJt.idx, REJf.idx], [3 3], 'the blunder, both ways');
    % tol: 1e-8 relative -- the 1/r_i identity on the downdated cofactor vs the clean re-solve
    tc.verifyEqual(solF.s0sq, solT.s0sq, 'RelTol', 1e-8, 's0^2 after rejection');
    tc.verifyEqual([solF.n, solF.r], [solT.n, solT.r], 'n and r after rejection');
    tc.verifyEqual(Df.rejectedIdx, 3);
    tc.verifyFalse(ismember(3, Df.flagIdx), 'the rejected observation is not flagged again');
    tc.verifyEqual(Df.obsIdx(:), Dt.obsIdx(:), 'one index space either way');
    % tol: as above -- the downdated cofactor vs the re-solved one
    tc.verifyEqual(Df.sigX, Dt.sigX, 'RelTol', 1e-8);
    % tol: as above
    tc.verifyEqual(Df.rNum(:), Dt.rNum(:), 'AbsTol', 1e-8, 'the redundancies of the survivors');
    % (2) tol: 1e-12 relative -- the same tau of the same alpha / nTests and r
    tc.verifyEqual(REJt.history(1, 4), R.diag.kCrit, 'RelTol', 1e-12, ...
        'snoop tests at diagnose''s critical value (the datum row is no test)');
end

function testSnoopInheritsTheDiagnosticsTest(tc)
% categories: contract
    % WP-131 (user decision): R.snoop tested at its own defaults (Alpha 1e-3,
    % Pope) whatever R.diag had used -- an observation R.diag flagged at
    % Alpha 0.01 could pass R.snoop "clean". Alpha, Bonferroni and Test now
    % default to R.diag's; an explicit pair still wins.
    rs = RandStream('mt19937ar', 'Seed', 5);
    pairs = repmat([1 2; 2 3; 3 4; 4 1; 1 3; 2 4], 7, 1);  m = size(pairs, 1);
    A = zeros(m, 4);  A(sub2ind(size(A), (1:m).', pairs(:, 1))) = -1;
    A(sub2ind(size(A), (1:m).', pairs(:, 2))) = 1;
    l = A*[0; 1.5; -0.7; 2.2] + 0.01*randn(rs, m, 1);  l(3) = l(3) + 0.2;
    data = struct('A', A, 'l', l, 'priors', struct('L', ones(1, 4), 'sigma', 1e3, 'name', 'datum'));
    R = lsqvce.fit(data, lsqvce.Options(Diagnostics=struct('Alpha', 0.01), Verbose=false));
    [~, REJ] = R.snoop('Verbose', false);
    % tol: 1e-12 relative -- the same tau of the same alpha / nTests and r
    tc.verifyEqual(REJ.history(1, 4), R.diag.kCrit, 'RelTol', 1e-12, 'snoop tests at R.diag''s Alpha');
    [~, REJ2] = R.snoop('Verbose', false, 'Alpha', 1e-3);
    tc.verifyGreaterThan(REJ2.history(1, 4), REJ.history(1, 4), 'an explicit Alpha still wins');
end

function testSnoopReDiagnosisKeepsTheCallersStorage(tc)
% categories: contract | tags: diagnostics
    % WP-135 (review 3, lane 1): snoop's re-diagnosis forwarded only alpha /
    % Bonferroni / test / thresholds: Store "flagged" came back "all" (five
    % n-sized vectors at the scale "flagged" exists for). And under Test "mad"
    % the report labels the MDB as the Baarda value (user decision): it is the
    % normal delta0, and the bias that hides from the t(r/2) threshold is larger.
    rs = RandStream('mt19937ar', 'Seed', 5);
    pairs = repmat([1 2; 2 3; 3 4; 4 1; 1 3; 2 4], 7, 1);  m = size(pairs, 1);
    A = zeros(m, 4);  A(sub2ind(size(A), (1:m).', pairs(:, 1))) = -1;
    A(sub2ind(size(A), (1:m).', pairs(:, 2))) = 1;
    l = A*[0; 1.5; -0.7; 2.2] + 0.01*randn(rs, m, 1);  l(3) = l(3) + 0.2;
    data = struct('A', A, 'l', l, 'priors', struct('L', ones(1, 4), 'sigma', 1e3, 'name', 'datum'));
    R = lsqvce.fit(data, lsqvce.Options(Diagnostics=struct('Store', "flagged"), Verbose=false));
    tc.verifyEqual(R.diag.store, 'flagged');
    [D2, REJ] = R.snoop('Verbose', false);
    tc.verifyEqual(REJ.idx, 3, 'the blunder is rejected, so snoop re-diagnoses');
    tc.verifyEqual(D2.store, 'flagged', 'the re-diagnosis keeps R.diag''s Store');
    Rm = lsqvce.fit(data, lsqvce.Options(Diagnostics=struct('Test', "mad"), Verbose=false));
    V = Rm.report('Verbose', false);
    codes = cellfun(@(f) f.code, V.findings, 'UniformOutput', false);
    k = find(strcmp(codes, 'REL_HEAD'), 1);
    tc.assertNotEmpty(k);
    tc.verifyTrue(contains(V.findings{k}.msg, 'Baarda'), 'the MDB under the MAD test is labelled');
end

function testMadScaleLeavesTheDatumRowsOut(tc)
% categories: metamorphic, contract
    % WP-130: Test "mad" formed its scale in a separate sweep that kept the
    % datum rows (r = 0, but qv ~ +-eps by rounding, so a near-zero raw score
    % about half the time) -- the WP-127 exclusion sat in the fused path only;
    % the scale came out low (median 0.85 of the fused one, 0.37 at worst,
    % in the reviewer's mirror). The two sweeps form the SAME raw scores, so
    % the scales are bitwise equal; snoop's 'mad' scale is that one too
    % (its first |w| is diagnose's). Over two backends and two sigmas, twelve
    % datum rows: a rounding pattern that hides the old defect has p ~ 2^-12.
    rs = RandStream('mt19937ar', 'Seed', 27);
    P  = [20 5; 28 6; 30 13; 24 16; 19 11];
    pairs = nchoosek(1:5, 2);  m = size(pairs, 1);
    A = zeros(m, 10);
    for k = 1:m
        e = P(pairs(k, 2), :) - P(pairs(k, 1), :);  e = e / norm(e);
        A(k, 2*pairs(k, 1) - [1 0]) = -e;  A(k, 2*pairs(k, 2) - [1 0]) = e;
    end
    l  = A * randn(rs, 10, 1) + 0.01 * randn(rs, m, 1);
    Kt = kron(ones(1, 5), eye(2));  Kr = reshape([-P(:, 2).'; P(:, 1).'], 1, []);
    for bk = ["qr" "normal"]
        for sg = [1 1e3]
            pri = struct('L', {Kt, Kr}, 'sigma', {sg, sg}, 'name', {'trans', 'rot'});
            R = lsqvce.fit(struct('A', A, 'l', l, 'priors', pri), Verbose=false, Backend=bk);
            Dm = lsqvce.engine.verdict.diagnose(R.source, R.batchFcn, 'nBatch', R.nBatch, 'Test', "mad");
            msg = sprintf('%s, datum sigma %g', bk, sg);
            tc.verifyEqual(Dm.robustScale, R.diag.robustScale, ['the mad sweep is the fused one: ' msg]);
            tc.verifyEqual(Dm.fhat, R.diag.fhat, msg);
        end
    end
    % snoop 'mad' on the levelling net with a datum: its first |w| is diagnose's
    rs = RandStream('mt19937ar', 'Seed', 5);
    pairs = repmat([1 2; 2 3; 3 4; 4 1; 1 3; 2 4], 7, 1);  m = size(pairs, 1);
    A = zeros(m, 4);  A(sub2ind(size(A), (1:m).', pairs(:, 1))) = -1;
    A(sub2ind(size(A), (1:m).', pairs(:, 2))) = 1;
    l = A*[0; 1.5; -0.7; 2.2] + 0.01*randn(rs, m, 1);  l(3) = l(3) + 0.2;
    R = lsqvce.fit(struct('A', A, 'l', l, 'priors', struct('L', ones(1, 4), 'sigma', 1e3, 'name', 'datum')), ...
                   Verbose=false);
    Dm = lsqvce.engine.verdict.diagnose(R.source, R.batchFcn, 'nBatch', R.nBatch, 'Test', "mad");
    [~, REJ] = R.snoop('Verbose', false, 'Test', "mad");
    tc.assertEqual(REJ.idx(1), 3, 'the blunder');
    % tol: 1e-10 relative -- the same raw score over the same MAD, two sweeps
    tc.verifyEqual(REJ.history(1, 2), abs(Dm.w(3)), 'RelTol', 1e-10, ...
        'snoop''s mad scale is diagnose''s (no datum score)');
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
    muteWarn(tc, 'lsqvce:fit:priorNotActive');   % provoked on purpose (gate rule 7)
    R2 = lsqvce.fit(data, Verbose=false);
    tc.verifyEqual(R2.diag.n, R0.diag.n);
    % WP-128: and SAYS so -- only data.priors were counted, a late local
    % group vanished without a word
    tc.verifyTrue(any(strcmp(R2.meta.warnings, 'lsqvce:fit:priorNotActive')), ...
        'a late local group is reported, as a late common one is');
    % WP-127: and through the joint loop. Pass 1 lists the group's pseudo-arc
    % (the arc identity is pinned) but it carries no rows: an INACTIVE arc,
    % not a 'gap' skip; from pass 2 it is one more stream in the per-pass
    % record -- the record's width changed and the run died before
    data.systems(1).priors = struct('L', e(cj), 'x0', target, 'sigma', 1e-6, 'from', 2, ...
                                    'name', 'late', 'estimate', "none");
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');
    Rj = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, ...
                    MaxPasses=3, Noise=struct('Grid', 16), Seed=3, Verbose=false));
    h = Rj.history;
    tc.verifyTrue(any(contains(h.streams, 'prior_late')), 'the late group is a stream of the record');
    tc.verifyEqual(size(h.sigma, 2), numel(h.streams), 'one column per pinned stream');
    tc.verifyFalse(any(contains(string(h.skips(:, 2)), 'prior_late')), 'an inactive prior arc is not a skip');
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

function testLawOnFixedStructureIsRefused(tc)
% categories: contract
    % R7-28: a USER-tier identifier no test had fired. A per-row sigma is a
    % fixed structure: it can carry a level ("scale"), not a spectral shape
    % ("law"); fit refuses the request by name before anything is solved.
    [data, ~] = declaration();
    data.systems(1).covariance.kbr.sigma = [0.2; 0.2];
    data.systems(1).covariance.kbr.estimate = "law";
    tc.verifyEqual(idOf(@() lsqvce.fit(data, lsqvce.Options(Estimate="noise", Verbose=false))), ...
                   'lsqvce:fit:lawOnFixedStructure');
end

function testVarLawUnderNoiseIsRefused(tc)
% categories: contract
    % R7-07: the joint sweep treated a 'var' law's d-interleaved block-Toeplitz
    % Cfcn as a per-axis temporal law over epochs (a wrong Sigma, silently) and
    % applied the estimated level as a std to a covariance. Refused by name.
    [data, ~] = declaration();
    data.systems(1).covariance.kbr = struct('type', 'var', 'ar', 0.3, 'sigma', 0.04);
    tc.verifyEqual(idOf(@() lsqvce.fit(data, lsqvce.Options(Estimate="noise", Verbose=false))), ...
                   'lsqvce:fit:varUnderNoise');
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
    % WP-145: no mute -- demo_priors converges (measured: it raises nothing)
    out = evalc('[R, s] = demo_priors();');
    echoWarnings(out);   % WP-145: its warnings reach the gate

    tc.verifyClass(R, 'lsqvce.Results');
    tc.verifyTrue(contains(out, 'lsqvce.layout says'), 'it prints the layout census');
    tc.verifyEqual([s.nc, s.nl], [7, 13], 'nc and nl are derived from the two tables');
    tc.verifyEqual(s.nullDim, 3, 'one null direction per day');
    % R7 WP-112 (R7-16): the refusal is the rank decision on the pivots,
    % singularR -- the normal-form fold that dropped the null directions (and
    % made it rankDeficient) is gone from the QR accumulator path
    tc.verifyEqual(s.unregularised, 'lsqvce:qr:singularR', ...
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

function testRobustJointSweepLeavesThePriorArcsOut(tc)
% categories: contract, robustness
    % WP-134 (review 3 low): the robust joint sweep counted a local prior
    % group's pseudo-arc as live epochs -- in the collapse guard, the
    % breakdown fraction and the hard solution's tau (the WP-108 defect that
    % adjustRobust had fixed). The sweep leaves them out: their factor stays 1
    % (no key), and R.noise.robust lists the data epochs only.
    [d, e, cj] = priorFixture();
    d.systems(1).priors = struct('L', e(cj), 'x0', 0, 'sigma', 1e6, ...
                                 'from', 1, 'name', 'loose', 'estimate', "none");
    muteWarn(tc, 'lsqvce:adjustVCEjoint:notConverged');   % MaxPasses 2: a contract point
    R = lsqvce.fit(d, lsqvce.Options(Robust="huber", Estimate="noise", Pool="stream", ...
                                     ArcFactors=false, MaxPasses=2, Verbose=false));
    obsT = string(R.noise.robust.keys(:, 2));
    tc.verifyFalse(any(startsWith(obsT, "prior_")), 'no prior pseudo-arc in the robust sweep');
    tc.verifyEqual(numel(R.noise.robust.factor), 2 * 3 * 2 * 300, 'one factor per data epoch');
end

%% ================================================ WP-150: the front-door verdicts
function testRobustStopRuleIsPerParameterSigma(tc)
% categories: contract, precision
    % WP-150 (F2): the stop rule max|dx| < Tol * max(max|x|, max sigma_x)
    % let one parameter's size set the threshold for all -- a coordinate of
    % 6.4e6 known to 1e-4 stopped IRLS "converged" after 2 passes at 22.6
    % formal sd from its fixed point (reproduced). Every parameter's step is
    % now measured in its own formal sigma.
    rng(1);  n = 1e4;  t = linspace(-1, 1, n).';  A = [ones(n, 1), t];
    l = A*[6.4e6; 2] + 1e-3*randn(n, 1);
    k = randperm(n, 1000);  l(k) = l(k) + 0.05 + 0.05*rand(1000, 1);
    R1 = lsqvce.fit(struct('A', A, 'l', l), Robust="huber", Verbose=false);
    % the reference asks for 1e-12 sigma, which this problem cannot resolve:
    % observations of 6.4e6 against noise of 1e-3 leave ~10 digits in the
    % residuals, so the Huber weights jitter at ~1e-5 sigma -- its budget runs out
    muteWarn(tc, 'lsqvce:adjustRobust:notConverged');
    R2 = lsqvce.fit(struct('A', A, 'l', l), Robust="huber", ...
                    RobustOpts=struct('Tol', 1e-12, 'MaxIter', 200), Verbose=false);
    tc.verifyTrue(R1.converged);
    % tol: the default Tol is 1e-8 of a formal sigma per step; the distance
    % to a fixed point ten thousand times tighter stays far below 0.01 sd
    % (measured 1.8e-7; the old rule left 22.6)
    tc.verifyLessThan(max(abs(R1.x - R2.x) ./ sqrt(diag(R2.sol.Cxx))), 0.01);
    tc.verifyLessThan(R1.history(end).dx, 1e-8, 'history.dx is the stop quantity, in formal sigma');
end

function testExcludedRowsStayOutOfTheRobustScale(tc)
% categories: contract, robustness
    % WP-150 (F4): rows of weight 0 have a whitened residual of 0 by
    % construction, and they entered the MAD -- with 40 % zero weights the
    % robust scale came out 0.33 of the true one and IGG3 rejected a third of
    % the good rows (reproduced). The scale is the live rows' alone.
    rng(2);  n = 200;  A = [ones(n, 1), (1:n).'/n];  l = A*[1; 2] + randn(n, 1);
    p = ones(n, 1);  p(1:80) = 0;
    R = lsqvce.fit(struct('A', A, 'l', l, 'p', p), Robust="igg3", Verbose=false);
    Rl = lsqvce.fit(struct('A', A(81:end, :), 'l', l(81:end)), Robust="igg3", Verbose=false);
    % tol: the same 120 live rows -- the scale and the estimate must not move
    % because 80 zero-weight rows ride along (measured 0.956 both ways)
    tc.verifyEqual(R.history(1).scale, Rl.history(1).scale, 'RelTol', 1e-10);
    tc.verifyEqual(R.x, Rl.x, 'RelTol', 1e-10);
end

function testHeldScalePriorWarns(tc)
% categories: contract
    % WP-150 (F3): a prior group with estimate "scale" in a run that cannot
    % estimate a prior's sigma was applied with its declared sigma, silently
    % (a group at 1e-3 shrank x to ~1e-4 against a truth of 1). It warns.
    rng(4);  A = randn(100, 5);  l = A*ones(5, 1) + randn(100, 1);
    d = struct('A', A, 'l', l, 'priors', struct('L', eye(5), 'sigma', 1e-3, 'estimate', "scale"));
    tc.verifyWarning(@() lsqvce.fit(d, Verbose=false), 'lsqvce:fit:priorScaleHeld');
    d.priors.estimate = "none";
    tc.verifyWarningFree(@() lsqvce.fit(d, Verbose=false), 'a declared "none" is no surprise');
end

function testSnoopOnGaussNewtonReturnsEstimates(tc)
% categories: contract
    % WP-150 (F5): a Gauss-Newton fit's rows are increments about R.x, and
    % R.snoop() returned the increment as the estimate (reproduced: [0.021
    % 0.028] for a fit at [2.985 0.676]). Both snoop paths -- the closed-form
    % downdate and the clean re-solve -- now return R.x + increment.
    rs = RandStream('mt19937ar', 'Seed', 5);
    t = (0:.1:5).';  y = 3*exp(-.7*t) + 1e-2*randn(rs, size(t));  y(20) = y(20) + 1;
    bf = @(k, x) deal([exp(-x(2)*t), -x(1)*t.*exp(-x(2)*t)], y - x(1)*exp(-x(2)*t));
    R = lsqvce.fit(struct('batchFcn', bf, 'numBatches', 1, 'x0', [1; .5]), Verbose=false);
    for re = [true false]
        [D, ~, sol] = R.snoop(Reaccumulate=re, Verbose=false);
        % tol: one blunder out of 51 rows moves the estimate by a few of its
        % formal errors (~1e-2), not back to zero
        tc.verifyEqual(sol.x, [3; 0.7], 'AbsTol', 0.05, sprintf('Reaccumulate %d: an estimate', re));
        tc.verifyEqual(D.x(1:2), sol.x, 'AbsTol', 1e-12, 'D.x agrees');
    end
end

function testRefusedVtPvIsDiagnosedWithItsOwnSigma0(tc)
% categories: contract, robustness
    % WP-150 (F6): after the normal backend refused a cancelled v'Pv,
    % diagnose studentised with the fallback sigma0 and flagged 28 rows where
    % the QR backend flags the one blunder (reproduced; R.sigma0 3.27 on
    % both). The verdict re-diagnoses with its own pass's s0^2.
    m = 300;  tt = (0:m-1).' * 5;  A = [sin(2*pi*tt/3000), cos(2*pi*tt/3000), ones(m, 1)];
    rs = RandStream('mt19937ar', 'Seed', 6);
    l = A * [1e5; -2e5; 3e5] + 3 * randn(rs, m, 1);  l(50) = l(50) + 25;
    muteWarn(tc, 'lsqvce:normalEq:cancellation');            % refused on purpose
    Rn = lsqvce.fit(struct('A', A, 'l', l), lsqvce.Options(Backend="normal", Verbose=0));
    Rq = lsqvce.fit(struct('A', A, 'l', l), lsqvce.Options(Backend="qr", Verbose=0));
    tc.verifyFalse(Rn.diag.s0Fallback, 'no fallback to the a-priori sigma0');
    tc.verifyEqual(Rn.diag.flagIdx, Rq.diag.flagIdx, 'the two backends flag the same rows');
    tc.verifyEqual(Rn.diag.flagIdx, 50);
end
