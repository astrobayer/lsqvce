function tests = test_adjustBlock
%TEST_ADJUSTBLOCK  Tests for the recover/Results/makeLoader access layer.
%
%     lsqvce.engine.adjustBlock      one-call block-arrow adjustment -> Results
%     lsqvce.Results      solution + diagnostics + provenance object
%     lsqvce.engine.makeLoader   day-builder -> accumulator-ready Loader
%
%   Both the legacy closure path (via synthLoader) and the new
%   makeLoader path must recover the synthetic truth and return sigma0 ~ 1.
%
%   Run: runtests('tests/test_adjustBlock')
    tests = functiontests(localfunctions);
end

function testAdjustBlockResultHasNoIterationHistory(tc)
% categories: contract, robustness | tags: diagnostics
    % adjustBlock is a SINGLE SOLVE: it returns a Results but never calls
    % setIterative, so .history stays empty. The two iteration readers must
    % therefore refuse CLEANLY -- naming the engine that does populate a history
    % -- rather than erroring on a missing field or returning an empty table
    % that looks like a converged run with zero passes.
    %
    % This is the case a real user hits (the sibling test in test_adjustVCE uses
    % a hand-built Results); adjustBlock is where an empty history occurs
    % naturally.
    cfg = synthCfg();
    loader = @(d,ot,a,uo) lsqvce.engine.synthLoader(cfg, 1, d, ot, a, uo);
    sys = struct('sysID',1, 'nl',cfg.nl, 'nc',cfg.nc, 'loader',loader);
    R = lsqvce.engine.adjustBlock(sys, 1:3);

    tc.verifyClass(R, 'lsqvce.Results');
    tc.verifyEmpty(R.history);
    tc.verifyTrue(R.converged);      % a single solve IS done, by definition
    tc.verifyError(@() R.convergence(), 'lsqvce:Results:noHistory');
    tc.verifyError(@() R.summary(),     'lsqvce:Results:noHistory');
end

%% ===================================================== R4 G-07 (Level L2)
function testSystemWithoutLocalUnknownsIsValid(tc)
% categories: contract, precision, robustness | tags: diagnostics
    % R4 G-07 (user decision 2026-08-31): nl = 0 is a valid system. The
    % streaming path used to refuse it at addSystem ('positive'); now the
    % local block is empty, pre-elimination is the identity and every local
    % consumer returns the consistent empty shape. Reference: the SAME day
    % data (same A, l, noise) run with ONE dummy local column carrying zero
    % signal -- the two GLS estimates differ by that column's projection
    % only, i.e. within the formal error.
    cfg1 = buildCfg();  cfg1.nl = 1;  cfg1.y1 = 0;              % dummy carries nothing
    buildOne  = @(d, ot) buildDayImpl(cfg1, d, ot);
    buildZero = @(d, ot) setfield(buildDayImpl(cfg1, d, ot), 'B', ...
                    zeros(cfg1.arcsPerDay*cfg1.mArc, 0));        %#ok<SFLD>  same A, l
    L0 = lsqvce.engine.makeLoader(buildZero, ...
            'observationTypes',{'pos','kbr'}, 'covariance',cfg1.cov, 'arcs',cfg1.arcsPerDay, ...
            'systemId',1, 'numLocal',0, 'numCommon',cfg1.nc, 'days',1:3);
    R0 = lsqvce.engine.adjustBlock({L0}, []);
    tc.verifyClass(R0, 'lsqvce.Results');
    tc.verifyEqual(numel(R0.x), cfg1.nc);
    tc.verifyLessThan(max(abs(R0.x - cfg1.x_c)), 0.15);       % recovers the truth
    tc.verifyLessThan(abs(R0.sigma0 - 1), 0.15);
    y = R0.acc.backSubstitute(1, R0.x);
    tc.verifyEqual(size(y), [0 1], 'no locals -> a 0x1 back-substitution');
    Cyy = R0.acc.localCovariance(1, R0.sol.Cxx);
    tc.verifyEqual(size(Cyy), [0 0], 'no locals -> a 0x0 local covariance');
    st = R0.acc.systemStore(1);
    tc.verifyEqual(st.nl, 0);
    V = evalc('R0.report();');                                 % the verdict still prints
    tc.verifyTrue(~isempty(V));

    L1 = lsqvce.engine.makeLoader(buildOne, ...
            'observationTypes',{'pos','kbr'}, 'covariance',cfg1.cov, 'arcs',cfg1.arcsPerDay, ...
            'systemId',1, 'numLocal',1, 'numCommon',cfg1.nc, 'days',1:3);
    R1 = lsqvce.engine.adjustBlock({L1}, []);
    % tol: 3 sigma -- same data; the estimates differ by the dummy column's
    % projection, whose size is the formal error of x (Cxx), not a rounding
    % bound; elementwise |dx| <= 3*sqrt(diag(Cxx)).
    se = sqrt(diag(R0.sol.Cxx));
    tc.verifyTrue(all(abs(R0.x - R1.x) <= 3*se), ...
        'nl = 0 and one zero-signal dummy local must agree within 3 sigma');
end

%% ===================================================== R4 G-06 (Level L1)
function testRawClosureHonoursPerSatelliteOverride(tc)
% categories: contract, robustness | tags: diagnostics
    % R4 finding G-06: adjustVCE emits its reweighting ONLY as
    % opt.covOverride.bySys.(satKey).(obsType) (adjustVCE.m covOverrideOpt),
    % while the raw-closure template synthLoader read only the flat
    % opt.covOverride.(obsType) -- a template-derived loader therefore
    % ignored every variance component under adjustVCE, and no VCE suite
    % drove a raw closure, so nothing caught it. Both shapes now resolve
    % through the ONE package-private applyCovOverride; this test drives the
    % template with each shape and pins that they whiten identically.
    cfg = synthCfg();
    loader = @(d,ot,a,uo) lsqvce.engine.synthLoader(cfg, 1, d, ot, a, uo);
    uo0 = struct('sysID', 1);
    b0  = loader(1, 'kbr', 1, uo0);                    % base: diag sigma 0.2
    idx = 1:numel(b0.l);
    [~, ~, w0] = b0.Wfcn(idx, b0.A, b0.B, b0.l, uo0);

    satKey = matlab.lang.makeValidName(char(string(1)));
    uoS = uo0;  uoS.covOverride.bySys.(satKey).kbr = struct('sigma', 0.6);
    bS  = loader(1, 'kbr', 1, uoS);
    [~, ~, wS] = bS.Wfcn(idx, bS.A, bS.B, bS.l, uoS);

    uoF = uo0;  uoF.covOverride.kbr = struct('sigma', 0.6);
    bF  = loader(1, 'kbr', 1, uoF);
    [~, ~, wF] = bF.Wfcn(idx, bF.A, bF.B, bF.l, uoF);

    % tol: 1e-14 -- l/0.6 versus (l/0.2)/3 differ by one extra rounding of a
    % scalar division (2 ulp); a dropped override would give a factor 3.
    tc.verifyEqual(wS, w0/3, 'RelTol', 1e-14, ...
        'the per-satellite (bySys) override must reach a raw closure');
    % tol: exact -- both shapes run the same merge and the same whitener
    tc.verifyEqual(wF, wS, ...
        'flat and per-satellite overrides must whiten identically');

    % A bySys entry for ANOTHER satellite must leave this loader untouched.
    uoX = uo0;  uoX.covOverride.bySys.(matlab.lang.makeValidName('2')).kbr = struct('sigma', 0.6);
    bX  = loader(1, 'kbr', 1, uoX);
    [~, ~, wX] = bX.Wfcn(idx, bX.A, bX.B, bX.l, uoX);
    % tol: exact -- no override applies, same code path as the base call
    tc.verifyEqual(wX, w0, 'an override keyed on another satellite must not apply');
end

%% ============================================================ recover paths
function testRecoverSynthClosure(tc)
% categories: contract, precision
    cfg = synthCfg();
    loader = @(d,ot,a,uo) lsqvce.engine.synthLoader(cfg, 1, d, ot, a, uo);
    sys = struct('sysID',1, 'nl',cfg.nl, 'nc',cfg.nc, 'loader',loader);
    R = lsqvce.engine.adjustBlock(sys, 1:3);
    tc.verifyClass(R, 'lsqvce.Results');
    tc.verifyClass(R.acc, 'lsqvce.engine.BlockAccumulator');
    tc.verifyLessThan(max(abs(R.x - cfg.x_c)), 0.15);   % recovers the truth
    tc.verifyLessThan(abs(R.sigma0 - 1), 0.15);         % noise model is correct
    tc.verifyTrue(R.converged);
end

function testRecoverMakeLoaderCell(tc)
% categories: precision
    cfg = buildCfg();
    L = lsqvce.engine.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos','kbr'}, 'covariance',cfg.cov, 'arcs',cfg.arcsPerDay, ...
            'systemId',1, 'numLocal',cfg.nl, 'numCommon',cfg.nc, 'days',1:3);
    R = lsqvce.engine.adjustBlock({L}, []);                % days come from the loader
    tc.verifyLessThan(max(abs(R.x - cfg.x_c)), 0.15);
    tc.verifyLessThan(abs(R.sigma0 - 1), 0.15);
end

function testResultsAccessors(tc)
% categories: robustness | tags: diagnostics
    cfg = buildCfg();
    L = lsqvce.engine.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos','kbr'}, 'covariance',cfg.cov, 'arcs',cfg.arcsPerDay, ...
            'systemId',1, 'numLocal',cfg.nl, 'numCommon',cfg.nc);
    R = lsqvce.engine.adjustBlock({L}, 1:2);
    tc.verifyEqual(numel(R.x), cfg.nc);
    tc.verifyGreaterThan(R.redundancy(), 0);
    V = R.report('verbose', false);                     % silent; returns the verdict
    tc.verifyNotEmpty(V);
    tc.verifyTrue(isfield(V, 'verdict'));               % the level-0..4 struct
    res = R.residuals('Verbose', false);                % re-reads via acc
    tc.verifyTrue(isstruct(res));
end

%% ============================================================ makeLoader unit
function testMakeLoaderMeta(tc)
% categories: contract | tags: provenance
    cfg = buildCfg();
    L = lsqvce.engine.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos','kbr'}, 'covariance',cfg.cov, 'arcs',2, 'numLocal',4, 'numCommon',6);
    m = L.meta();
    tc.verifyEqual(m.obsTypes, {'pos','kbr'});
    tc.verifyEqual(m.arcs{1}, 1:2);
    m2 = L.closure('meta', 1, [], struct());            % accumulator protocol
    tc.verifyEqual(m2.obsTypes, {'pos','kbr'});
    b = L.closure(1, 'pos', 1, struct());               % one arc batch
    tc.verifyTrue(all(isfield(b, {'A','B','l','t','Cfcn','Wfcn','filterMemory'})));
    tc.verifyEqual(size(b.A,1), cfg.mArc);
end

function testMakeLoaderValidation(tc)
% categories: contract
    cfg = buildCfg();
    tc.verifyError(@() lsqvce.engine.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
        'covariance',cfg.cov, 'arcs',2), 'lsqvce:makeLoader:noObsTypes');
    tc.verifyError(@() lsqvce.engine.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
        'observationTypes',{'nope'}, 'covariance',cfg.cov), 'lsqvce:makeLoader:noCov');
    bad = @(d,ot) struct('A',randn(10,6));              % missing B,l,t
    Lb = lsqvce.engine.makeLoader(bad, 'observationTypes',{'pos'}, 'covariance',cfg.cov, 'numLocal',4, 'numCommon',6);
    tc.verifyError(@() Lb.closure(1,'pos',1,struct()), 'lsqvce:makeLoader:badDay');
end

%% ==================================================== option names
function testMakeLoaderCanonicalNames(tc)
% categories: contract
    % Every option has ONE name. Building with them must warn about nothing.
    cfg = buildCfg();
    ws = warning('off', 'all'); c = onCleanup(@() warning(ws)); %#ok<NASGU>
    lastwarn('', '');
    L = lsqvce.engine.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos'}, 'covariance',cfg.cov, ...
            'systemId',3, 'numLocal',4, 'numCommon',6);
    tc.verifyEqual(L.sysID, 3);          % carried onto the internal field
    tc.verifyEqual(L.nl, 4);
    tc.verifyEqual(L.nc, 6);
    tc.verifyEqual(L.obsTypes, {'pos'});
end

function testMakeLoaderShortAliasesRemoved(tc)
% categories: contract
    % The terse aliases (sysID/nl/nc/cov/obsTypes/...) were REMOVED 2026-07-18.
    % inputParser has PartialMatching off, so an old name is now an unmatched
    % parameter and must ERROR rather than silently doing nothing.
    cfg = buildCfg();
    % WP-39 c1: the unknown-name rejection moved from inputParser into the
    % typed arguments signature -- the id is now MATLAB:TooManyInputs (the
    % arguments-block unknown-NV id). The protective contract (a removed
    % terse alias must ERROR, never silently no-op) is unchanged and
    % release-noted with the b1-r3 rule.
    tc.verifyError(@() lsqvce.engine.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos'}, 'covariance',cfg.cov, ...
            'sysID',3, 'numLocal',4, 'numCommon',6), ...
        'MATLAB:TooManyInputs');
end

%% ============================================================ recover errors
function testRecoverErrors(tc)
% categories: contract
    tc.verifyError(@() lsqvce.engine.adjustBlock({1}, 1:3), 'lsqvce:adjustBlock:badCell');
    tc.verifyError(@() lsqvce.engine.adjustBlock(struct('sysID',1,'nl',4,'nc',6), 1:3), ...
                   'lsqvce:adjustBlock:badStruct');
end

function testShapeIsValidatedBeforeNcIsDerived(tc)
% categories: contract | tags: diagnostics
    % REGRESSION GUARD (2026-07-18). adjustBlock derives nc with lsqvce.engine.resolveNc
    % and validates the shape of `systems` with normalizeSystems. The ORDER is
    % load-bearing: deriving nc first makes every malformed `systems` report the
    % generic lsqvce:resolveNc:missing ("no system declares a common-parameter
    % count") instead of the specific id that names what is actually wrong.
    %
    % This shipped as a real bug -- the required order was written in the
    % handover but never applied to the code -- so it is pinned here rather than
    % left to prose. The distinguishing feature of these inputs is that they are
    % malformed AND declare no nc, so whichever step runs first decides the id.
    tc.verifyError(@() lsqvce.engine.adjustBlock({1}, 1:3), ...
        'lsqvce:adjustBlock:badCell', ...
        'a non-Loader cell entry must report badCell, not resolveNc:missing');

    tc.verifyError(@() lsqvce.engine.adjustBlock(struct('sysID',1,'nl',4), 1:3), ...
        'lsqvce:adjustBlock:badStruct', ...
        'a spec missing .loader must report badStruct, not resolveNc:missing');

    % The complement: a WELL-FORMED systems argument that genuinely declares no
    % nc must still reach resolveNc and report the nc problem.
    tc.verifyError(@() lsqvce.engine.adjustBlock(struct('sysID',1,'nl',4, ...
                        'loader',@(varargin) []), 1:3), ...
        'lsqvce:resolveNc:missing', ...
        'a well-formed spec without nc must still report resolveNc:missing');
end

%% ====================================================== parametrisation pipe
function testParamBackCompatTwoArg(tc)
% categories: contract
    % An existing 2-arg buildDay must keep working unchanged when common/local
    % are supplied to makeLoader (callBuildDay passes only what it declares).
    cfg = buildCfg();
    L = lsqvce.engine.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos'}, 'covariance',cfg.cov, 'arcs',2, ...
            'commonParams', paramRows(6), 'localParams', paramRows(4));
    b = L.closure(1, 'pos', 1, struct());
    tc.verifyEqual(size(b.A,1), cfg.mArc);        % builds fine, ignores extras
end

function testParamPipedToFiveArgBuilder(tc)
% categories: contract
    % A 5-arg builder must RECEIVE common/local/procOpt exactly as passed.
    cfg = buildCfg();
    common = paramRows(6);  local = paramRows(4);
    po = struct('flag', 42);
    seen = containers.Map('KeyType','char','ValueType','any');   % v6: a day carries
    grab = @(d,ot,c,l,p) captureBuild(cfg,d,ot,c,l,p,seen);      % no extra fields
    L = lsqvce.engine.makeLoader(grab, 'observationTypes',{'pos'}, 'covariance',cfg.cov, 'arcs',2, ...
            'commonParams', common, 'localParams', local, 'processingOptions', po);
    L.day(1, 'pos');
    tc.verifyEqual(seen('common'), common);
    tc.verifyEqual(seen('local'),  local);
    tc.verifyEqual(seen('po').flag, 42);
end

function testParamDerivesNcNl(tc)
% categories: contract
    % nc = rows(common), nl = rows(local), carried on L.
    cfg = buildCfg();
    L = lsqvce.engine.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos'}, 'covariance',cfg.cov, ...
            'commonParams', paramRows(6), 'localParams', paramRows(4));
    tc.verifyEqual(L.nc, 6);
    tc.verifyEqual(L.nl, 4);
end

function testParamNlConflictErrors(tc)
% categories: contract
    % Passing both 'nl' and 'local' with disagreeing counts must error.
    cfg = buildCfg();
    tc.verifyError(@() lsqvce.engine.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
        'observationTypes',{'pos'}, 'covariance',cfg.cov, 'numLocal',5, 'localParams', paramRows(4)), ...
        'lsqvce:makeLoader:nlConflict');
end

function testParamColumnCountValidated(tc)
% categories: contract
    % v3 (L10 WP-72): the fixed 14-column shape is RETIRED -- any RxC cell is
    % accepted (nc = rows, C the caller's convention), while a non-cell is
    % still rejected with paramType. Both halves of the relaxation pinned.
    cfg = buildCfg();
    L = lsqvce.engine.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
        'observationTypes',{'pos'}, 'covariance',cfg.cov, 'commonParams', cell(6,10));
    tc.verifyEqual(L.nc, 6, 'nc must be the row count, whatever C is');
    % a non-cell is stopped by the SIGNATURE (MATLAB validation id, b1-r3);
    % the body's paramType guards what the signature cannot: a cell that is
    % not a 2-D matrix.
    tc.verifyError(@() lsqvce.engine.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
        'observationTypes',{'pos'}, 'covariance',cfg.cov, 'commonParams', cell(2,2,2)), ...
        'lsqvce:makeLoader:paramType');
end

function testScopePlanOffsetsAreGapSafe(tc)
% categories: contract, robustness
    % L11 WP-76: the SCOPED plan (engine tier, 'Scoped') pins the binding
    % layout against a hand enumeration on a RAGGED schedule: day 2 loses arc
    % 2, day 3 loses everything. Per day [arc blocks, sorted ids | day
    % block], trailing system block; a day without arcs gets NO blocks.
    cfg = buildCfg();
    tbl = [{'arc','ic1'}; {'arc','ic2'}; {'day','bias'}; {'system','scale'}];
    L = lsqvce.engine.makeLoader(@(d,ot) buildDayImpl(cfg, d, ot), ...
        'observationTypes', {'pos'}, 'covariance', cfg.cov, ...
        'LocalParams', tbl, 'Scoped', true, ...
        'Arcs', @(day, obsType) scheduleOf(day), 'Days', 1:3);
    p = L.localPlan;
    % hand enumeration: d1 arcs [1 2] -> blocks at 0,2 + day at 4;
    % d2 arc [1] -> block at 5 + day at 7; d3 nothing; system at 8; nl = 9
    tc.verifyEqual(p.nl, 9);
    tc.verifyEqual(L.nl, 9, 'nl is DERIVED from the plan');
    tc.verifyEqual(p.arcId{1},  [1 2]);   tc.verifyEqual(p.offset{1}, [0 2]);
    tc.verifyEqual(p.dayOffset(1), 4);    tc.verifyEqual(p.dayCount(1), 1);
    tc.verifyEqual(p.arcId{2},  1);       tc.verifyEqual(p.offset{2}, 5);
    tc.verifyEqual(p.dayOffset(2), 7);    tc.verifyEqual(p.dayCount(2), 1);
    tc.verifyEqual(p.arcId{3},  zeros(1,0), 'day 3 has no arcs');
    tc.verifyEqual(p.dayCount(3), 0, 'a day without arcs gets no day block');
    tc.verifyEqual(p.sysOffset, 8);       tc.verifyEqual(p.sysCount, 1);
    tc.verifyError(@() lsqvce.engine.makeLoader(@(d,ot) buildDayImpl(cfg, d, ot), ...
        'observationTypes', {'pos'}, 'covariance', cfg.cov, ...
        'LocalParams', [{'weekly','x'}], 'Scoped', true, 'Days', 1:3), ...
        'lsqvce:makeLoader:badScope');
end

function ids = scheduleOf(day)
    switch day
        case 1, ids = [1 2];
        case 2, ids = 1;
        otherwise, ids = zeros(1, 0);
    end
end

%% ====================================================== plot-hook options
function rec = recorderHandle()
% A tiny MUTABLE recorder usable from inside a plot-hook closure. containers.Map
% is a HANDLE object, so writes through the closure are visible to the caller
% without any nested-function sharing. Keys are 1..N insertion order; .log below
% reads them back as a cell in order.
    m = containers.Map('KeyType','double','ValueType','any');
    rec = struct('map', m);
end

function ctxOut = recordCtx(rec, ctx)
% Append ctx to the recorder and return [] (no figure -> showExport no-ops).
    n = rec.map.Count + 1;
    rec.map(n) = ctx;
    ctxOut = [];
end

function logged = recorderLog(rec)
% Read the recorded ctx structs back in insertion order.
    ks = sort(cell2mat(keys(rec.map)));
    logged = cell(1, numel(ks));
    for i = 1:numel(ks), logged{i} = rec.map(ks(i)); end
end

function cfg = buildCfg()
% cfg for the makeLoader buildDay path.
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
    D.l = D.A*cfg.x_c + D.B*cfg.y1 + e;
    D.arcIdx = {1:cfg.mArc, cfg.mArc + (1:cfg.mArc)};
    D = arcCells(D);                                % v6: the cells are the arcs
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

function C = paramRows(nRow)
% A valid 14-column parametrisation array with nRow unknowns (one per row).
    C = cell(nRow, 14);
    for i = 1:nRow
        C{i,1} = i;              % e.g. an index in column 1; rest left empty
    end
end

function D = captureBuild(cfg, day, obsType, common, local, procOpt, seen)
% A 5-arg builder that both builds a day AND records what it received in the
% handle Map SEEN, so the test can assert the parametrisation reached
% buildDay intact (v6: a day struct carries no field but the contract's).
    D = buildDayImpl(cfg, day, obsType);
    seen('common') = common;  seen('local') = local;  seen('po') = procOpt;   %#ok<NASGU>
end

function cfg = synthCfg()
% cfg in the exact shape synthLoader expects.
    cfg = struct('nc',6, 'nl',4, 'arcsPerDay',2, 'mArc',300, 'dt',5, ...
                 'seed',42, 'collinear',false, 'collinearEps',1e-7);
    cfg.obsTypes = {'pos','kbr'};
    rs = RandStream('mt19937ar','Seed',3);
    cfg.x_c = randn(rs, 6, 1);
    cfg.y   = containers.Map(1, randn(rs, 4, 1));
    cfg.cov.pos = struct('type','arma','ar',[1 -1.6 0.9 -0.25],'ma',1,'sigma',0.7);
    cfg.cov.kbr = struct('type','diag','sigma',0.2);
end

%% ============================ WP-39 stage (a): typed signature (F-07)
%% =========================== WP-39 stage (a): typed-signature enforcement
function test_block_signature_types_and_enum_values(tc)
% categories: contract
    % WP-39 (F-07), premise CORRECTED in round a1-2: MATLAB matches
    % comma-syntax name-value NAMES case-insensitively (and partially)
    % even under arguments blocks, so lowercase 'backend' still WORKS --
    % the true breaks are (1) enum VALUES, now case-sensitive lowercase
    % under mustBeMember where strcmpi accepted 'QR' (the fails-on-old
    % lever), and (2) typed constraints firing at the signature.
    ok = false;
    try
        lsqvce.engine.adjustBlock({}, [], 'Backend', 'QR');   % enum VALUE casing
    catch e
        ok = startsWith(e.identifier, 'MATLAB:');
    end
    verifyTrue(tc, ok, ...
        'Backend value ''QR'' must be rejected (case-sensitive mustBeMember)');
    ok = false;
    try
        lsqvce.engine.adjustBlock({}, [], 'RegLambda', -1);   % typed constraint
    catch e
        ok = startsWith(e.identifier, 'MATLAB:');
    end
    verifyTrue(tc, ok, 'negative RegLambda must be rejected by the signature');
    ok = false;
    try
        lsqvce.engine.adjustBlock({}, [], 'Backend', 'qr', 'nosuchopt', 1);
    catch e
        ok = startsWith(e.identifier, 'MATLAB:');
    end
    verifyTrue(tc, ok, 'an unknown option name must be rejected');
end

function test_signature_typed_arguments(tc)
% categories: contract
    % WP-39 (fails on the old code): the typed `arguments` block rejects
    % what inputParser tolerated -- unknown names, WRONG-CASE names (the
    % old lowercase 'backend' etc.), and wrong types/enum casings.
    sys = struct('sysID', 1, 'nl', 2, 'nc', 3, ...
                 'loader', @(varargin) struct('A',[],'B',[],'l',[]));
    mustErr = { ...
        @() lsqvce.engine.adjustBlock(sys, 1, 'backend', 'qr'), ...     % old casing
        @() lsqvce.engine.adjustBlock(sys, 1, 'Backend', 'QR'), ...     % enum casing
        @() lsqvce.engine.adjustBlock(sys, 1, 'Backend', 'cholesky'),...% bad enum
        @() lsqvce.engine.adjustBlock(sys, 1, 'RegLambda', -1), ...     % sign
        @() lsqvce.engine.adjustBlock(sys, 1, 'NoSuchOption', 1) };
    for k = 1:numel(mustErr)
        ok = false;
        try, mustErr{k}(); catch, ok = true; end
        verifyTrue(tc, ok, sprintf('signature case %d must error', k));
    end
end


%% ============================================================ run provenance
function testResultsMetaContract(tc)
% categories: contract | tags: provenance
% WP-44 (T-05): every Results carries .meta from the ONE composer. Asserts
% the documented field list (SUPERSET lock: fields are add-only, never
% renamed), the field types, and immutability. try/catch identifier pattern
% throughout (the WP-43 r1 lesson -- never the verifyError output form).
    cfg = buildCfg();
    L = lsqvce.engine.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos','kbr'}, 'covariance',cfg.cov, ...
            'arcs',cfg.arcsPerDay, 'systemId',1, ...
            'numLocal',cfg.nl, 'numCommon',cfg.nc);
    R = lsqvce.engine.adjustBlock({L}, 1:2);

    m = R.meta;
    tc.verifyClass(m, 'struct');
    documented = {'toolboxVersion','matlabRelease','createdUTC','seeds','options'};
    tc.verifyEmpty(setdiff(documented, fieldnames(m)), ...
        'meta must carry AT LEAST the documented field list (add-only lock)');
    tc.verifyClass(m.toolboxVersion, 'char');
    tc.verifyClass(m.matlabRelease,  'char');
    tc.verifyClass(m.createdUTC,     'char');
    tc.verifyClass(m.seeds,          'struct');
    tc.verifyClass(m.options,        'struct');
    % ISO 8601 with trailing Z, as documented
    tc.verifyTrue(~isempty(regexp(m.createdUTC, ...
        '^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$', 'once')), ...
        'createdUTC must be ISO 8601 with a trailing Z');
    % the version comes from Contents.m, the ONE source
    txt = fileread(fullfile(fileparts(which('lsqvce.Contents')), 'Contents.m'));
    tok = regexp(txt, '^%\s*Version\s+(\S+)', 'tokens', 'once', 'lineanchors');
    tc.verifyEqual(m.toolboxVersion, tok{1}, ...
        'toolboxVersion must equal the Contents.m version line');
    % adjustBlock resolves no RNG seed: empty seeds, resolved options present
    tc.verifyEmpty(fieldnames(m.seeds));
    tc.verifyTrue(isfield(m.options, 'backend'), ...
        'options must be the RESOLVED adjustBlock option struct');

    % immutability
    caught = MException.empty;
    try
        R.meta = struct();
    catch e
        caught = e;
    end
    tc.verifyNotEmpty(caught, 'meta must be immutable');
    tc.verifyEqual(caught.identifier, 'MATLAB:class:SetProhibited');
end

function testResultsMetaDeterminismAndBareCtor(tc)
% categories: robustness | tags: provenance
% WP-44: (1) two runs with identical inputs/options give identical meta
% EXCEPT .createdUTC; (2) a bare 7-argument construction still carries a
% complete meta with empty seeds/options -- never a missing field.
    cfg = buildCfg();
    mk = @() lsqvce.engine.adjustBlock({lsqvce.engine.makeLoader( ...
            @(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos','kbr'}, 'covariance',cfg.cov, ...
            'arcs',cfg.arcsPerDay, 'systemId',1, ...
            'numLocal',cfg.nl, 'numCommon',cfg.nc)}, 1:2);
    m1 = mk().meta;   m2 = mk().meta;
    m1c = rmfield(m1, 'createdUTC');   m2c = rmfield(m2, 'createdUTC');
    tc.verifyEqual(m1c, m2c, ...
        'identical runs must give identical meta up to createdUTC');

    sol = struct('x', zeros(3,1), 'sigma0', 1);
    Rb = lsqvce.Results(sol, struct(), [], struct(), @() deal([],[],[]), 0, []);
    tc.verifyClass(Rb.meta, 'struct');
    tc.verifyEmpty(setdiff({'toolboxVersion','matlabRelease','createdUTC', ...
        'seeds','options'}, fieldnames(Rb.meta)), ...
        'a bare construction must still carry the full meta');
    tc.verifyEmpty(fieldnames(Rb.meta.seeds));
    tc.verifyEmpty(fieldnames(Rb.meta.options));
end

function testReportPrintsTheProvenanceLine(tc)
% categories: contract | tags: diagnostics, provenance
% WP-44: R.report() forwards meta and the printed level-0 block carries the
% ONE provenance line; the finding is machine-readable (code 'PROVENANCE');
% a bare lsqvce.report(D) prints none of it.
    cfg = buildCfg();
    L = lsqvce.engine.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos','kbr'}, 'covariance',cfg.cov, ...
            'arcs',cfg.arcsPerDay, 'systemId',1, ...
            'numLocal',cfg.nl, 'numCommon',cfg.nc);
    R = lsqvce.engine.adjustBlock({L}, 1:2);

    out = evalc('V = R.report();');
    tc.verifyTrue(contains(out, 'run: lsqvce'), ...
        'the printed report must carry the provenance line');
    tc.verifyTrue(contains(out, R.meta.toolboxVersion));
    codes = cellfun(@(f) f.code, V.findings, 'UniformOutput', false);
    tc.verifyTrue(any(strcmp(codes, 'PROVENANCE')), ...
        'the finding list must carry the PROVENANCE note');

    % the standalone report has no provenance -- no line, no finding
    out0 = evalc('V0 = lsqvce.report(R.diag);');
    tc.verifyFalse(contains(out0, 'run: lsqvce'));
    codes0 = cellfun(@(f) f.code, V0.findings, 'UniformOutput', false);
    tc.verifyFalse(any(strcmp(codes0, 'PROVENANCE')));
end

%% ===================================== WP-52: native ragged arc schedules
function testArcResolverAnswersPerDay(tc)
% categories: contract | tags: diagnostics
% WP-52: the 'meta' probe carries the day, so an @(day,obsType) resolver
% gives each day its OWN arc list -- the mechanism a ragged schedule needs.
% Asserted through the loader's public closure, which is what the engines
% actually drive.
    plan = {[1 2 3], [1 3], [2 3 4]};          % day 1..3, arcs absent BY ID
    L = raggedLoader(tc, plan, []);
    for d = 1:3
        m = L.closure('meta', d, [], struct());
        tc.verifyEqual(m.arcs{1}, plan{d}, ...
            sprintf('day %d must report its own arc list', d));
    end
    % the static form still ignores the day, exactly as before
    Ls = raggedLoader(tc, 4, []);
    tc.verifyEqual(Ls.closure('meta', 2, [], struct()).arcs{1}, 1:4);
end

function testEmptyCellIsADataGapNotAnError(tc)
% categories: contract, robustness
% v6 (R6 WP-102): the cells are the arcs. Under a declared schedule the day
% returns one cell per SCHEDULED id, and an EMPTY cell is that id's data gap
% (empty batch, which every engine skips). An id the schedule does not name
% stays the positional out-of-range error -- a mis-declared static schedule
% is a different mistake from a ragged one -- and a day returning the wrong
% NUMBER of cells is refused by name.
    L = gapLoader(tc);                          % schedule [1 2 3]; day 2 lacks arc 2
    got = L.closure(2, 'pot', 3, struct());     % day 2 HAS arc 3
    tc.verifyNotEmpty(got);
    tc.verifyNotEmpty(got.l);
    gap = L.closure(2, 'pot', 2, struct());     % day 2's cell for arc 2 is empty
    tc.verifyEmpty(gap, 'an empty cell must come back as a data gap');
    m = L.closure('meta', 2, [], struct());
    tc.verifyEqual(m.arcs{1}, [1 2 3], 'the meta probe answers the SCHEDULE, gaps included');

    caught = MException.empty;
    Ls = staticPositionalLoader(tc);
    try
        Ls.closure(1, 'pot', 9, struct());
    catch e
        caught = e;
    end
    tc.verifyNotEmpty(caught, 'the positional guard must survive');
    tc.verifyEqual(caught.identifier, 'lsqvce:makeLoader:badArc');

    % three arcs scheduled, two cells returned: refused by name
    Lc = lsqvce.engine.makeLoader(@(d, ot) simpleDay(d, ot), ...
        'observationTypes', {'pot'}, ...
        'covariance', struct('pot', struct('type','diag','sigma',1e-3)), ...
        'arcs', 3, 'days', 1:2, 'numCommon', 4, 'numLocal', 1);
    tc.verifyError(@() Lc.closure(1, 'pot', 1, struct()), 'lsqvce:makeLoader:arcCount');
end

function testCellsAreTheScheduleWhenNoneIsDeclared(tc)
% categories: contract, metamorphic
% v6 (R6 WP-102): 'Arcs' is optional. Without it the 'meta' probe reads the
% day and answers 1:nArc from its cells, and the loader is BITWISE the same
% loader as the one declaring the count it would have discovered -- the
% declaration states nothing the data does not already state.
    cfg = buildCfg();
    Ld = lsqvce.engine.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos','kbr'}, 'covariance',cfg.cov, 'numLocal',4, 'numCommon',6);
    Ls = lsqvce.engine.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos','kbr'}, 'covariance',cfg.cov, 'arcs',cfg.arcsPerDay, ...
            'numLocal',4, 'numCommon',6);
    m = Ld.closure('meta', 1, [], struct());
    tc.verifyEqual(m.arcs, {1:cfg.arcsPerDay, 1:cfg.arcsPerDay}, 'discovered from the cells');
    for a = 1:cfg.arcsPerDay
        bd = Ld.closure(1, 'pos', a, struct());  bs = Ls.closure(1, 'pos', a, struct());
        tc.verifyEqual(bd.A, bs.A);  tc.verifyEqual(bd.l, bs.l);  tc.verifyEqual(bd.t, bs.t);
    end
    Rd = lsqvce.engine.adjustBlock({Ld}, 1:2);
    Rs = lsqvce.engine.adjustBlock({Ls}, 1:2);
    tc.verifyEqual(Rd.x, Rs.x, 'discovered == declared, bitwise');
    tc.verifyEqual(Rd.sigma0, Rs.sigma0);
    % the string form is gone: refused at the door, by name
    tc.verifyError(@() lsqvce.engine.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos'}, 'covariance',cfg.cov, 'arcs','fromDay', ...
            'numLocal',4, 'numCommon',6), 'lsqvce:makeLoader:badArcs');
end

function testDayContractRefusesTheV5Shape(tc)
% categories: contract, robustness
% v6: A, B, l, t must be 1 x nArc cells agreeing in number, every arc's four
% entries sharing a row count; the v5 fields arcIdx / arcId and any other
% stray field are refused BY NAME, at the day, with day and stream in the
% message; Cep must be a matching cell.
    cov = struct('pot', struct('type','diag','sigma',1));
    mk  = @(bd) lsqvce.engine.makeLoader(bd, 'observationTypes', {'pot'}, 'covariance', cov, ...
                                         'numCommon', 2, 'numLocal', 0);
    m = 10;  A = ones(m, 2);  B = zeros(m, 0);  l = ones(m, 1);  t = (1:m).';
    ok = struct('A', {{A}}, 'B', {{B}}, 'l', {{l}}, 't', {{t}});
    L = mk(@(d, ot) ok);
    tc.verifyEqual(size(L.closure(1, 'pot', 1, struct()).A), [m 2], 'the single-arc day is {A}');

    v5 = struct('A', A, 'B', B, 'l', l, 't', t);                   % matrices, not cells
    L = mk(@(d, ot) v5);
    tc.verifyError(@() L.closure(1, 'pot', 1, struct()), 'lsqvce:makeLoader:badDay');
    withIdx = ok;  withIdx.arcIdx = {1:m};                            % the v5 field
    L = mk(@(d, ot) withIdx);
    try
        L.closure(1, 'pot', 1, struct());
        tc.verifyFail('arcIdx must be refused');
    catch e
        tc.verifyEqual(e.identifier, 'lsqvce:makeLoader:badDay');
        tc.verifySubstring(e.message, 'arcIdx');
        tc.verifySubstring(e.message, 'retired');
    end
    uneven = ok;  uneven.B = {B, B};                                  % 1 vs 2 cells
    L = mk(@(d, ot) uneven);
    tc.verifyError(@() L.closure(1, 'pot', 1, struct()), 'lsqvce:makeLoader:dayShape');
    short = ok;  short.l = {l(1:m-1)};                                % rows disagree inside arc 1
    L = mk(@(d, ot) short);
    tc.verifyError(@() L.closure(1, 'pot', 1, struct()), 'lsqvce:makeLoader:dayShape');
    cepMat = ok;  cepMat.Cep = ones(m, 1);                            % Cep not a cell
    L = mk(@(d, ot) cepMat);
    tc.verifyError(@() L.closure(1, 'pot', 1, struct()), 'lsqvce:makeLoader:badDayCep');
    cepShort = ok;  cepShort.Cep = {ones(m - 1, 1)};                  % Cep does not cover the arc
    L = mk(@(d, ot) cepShort);
    tc.verifyError(@() L.closure(1, 'pot', 1, struct()), 'lsqvce:makeLoader:badDayCep');
    % an empty day (no cells at all) is valid: the meta probe reports no arcs
    empty = struct('A', {{}}, 'B', {{}}, 'l', {{}}, 't', {{}});
    L = mk(@(d, ot) empty);
    tc.verifyEqual(L.closure('meta', 1, [], struct()).arcs{1}, zeros(1, 0));
end

function testLocalsPerArcDerivesNlAndOffsets(tc)
% categories: contract | tags: provenance
% WP-52: makeLoader owns the per-arc local bookkeeping. nl is DERIVED from
% the schedule (so it cannot disagree), the offsets are contiguous in
% (day, arc) order, and buildDay receives its day's slice as arcInfo.
    plan = {[1 2 3], [1 3], [2 3 4]};           % 8 arcs over 3 days
    seen = containers.Map('KeyType','double','ValueType','any');
    L = raggedLoader(tc, plan, 2, seen);        % 2 locals per arc
    tc.verifyEqual(L.nl, 16, 'nl must be derived: 8 arcs x 2 locals');
    tc.verifyEqual(L.localPlan.nArcTot, 8);

    L.day(2, 'pot');                            % triggers buildDay for day 2
    info = seen(2);
    tc.verifyEqual(info.arcId, [1 3]);
    tc.verifyEqual(info.count, [2 2]);
    tc.verifyEqual(info.offset, [6 8], ...
        'day 2 follows day 1''s three arcs: offsets 6 and 8');
    tc.verifyEqual(info.ordinal, [4 5]);
    tc.verifyEqual(info.nl, 16);
    tc.verifyEqual(info.offsetOf(3), 8, ...
        'offsetOf must map an arc id to its column offset');

    % a per-arc count: locals may differ from arc to arc
    Lv = raggedLoader(tc, plan, @(d, a) 1 + mod(a, 2));
    tc.verifyEqual(Lv.nl, sum(cellfun(@(ids) sum(1 + mod(ids, 2)), plan)));
end

function testLocalsPerArcRejectsUnknowableSchedules(tc)
% categories: contract, robustness
% nl must be known at construction, so the plan refuses the two cases where
% it is not: no day list, and no declared schedule (v6: the cells would then
% be the schedule, and they exist only once the day is read).
    caught = MException.empty;
    try
        lsqvce.engine.makeLoader(@(d, ot) simpleDay(d, ot), ...
            'observationTypes', {'pot'}, ...
            'covariance', struct('pot', struct('type','diag','sigma',1)), ...
            'arcs', 2, 'numCommon', 2, 'localsPerArc', 2);   % no 'days'
    catch e
        caught = e;
    end
    tc.verifyNotEmpty(caught);
    tc.verifyEqual(caught.identifier, 'lsqvce:makeLoader:localPlan');

    caught = MException.empty;
    try
        lsqvce.engine.makeLoader(@(d, ot) simpleDay(d, ot), ...
            'observationTypes', {'pot'}, ...
            'covariance', struct('pot', struct('type','diag','sigma',1)), ...
            'days', 1:2, 'numCommon', 2, 'localsPerArc', 2);   % no 'arcs': would have to READ
    catch e
        caught = e;
    end
    tc.verifyEqual(caught.identifier, 'lsqvce:makeLoader:localPlan');
end

function testRaggedScheduleSolvesAndSecondPassAgrees(tc)
% categories: precision | tags: diagnostics
% The end-to-end promise: a ragged schedule (arcs absent by id, every arc a
% different length, locals owned by makeLoader) recovers the planted truth,
% AND the residual pass -- which re-enumerates the arcs independently --
% agrees with the accumulation. That agreement is the sharp edge: a
% resolver that answered differently on the second sweep would show up here.
    plan = {[1 2 3], [1 3], [2 3 4]};
    L = raggedLoader(tc, plan, 2);
    R = lsqvce.engine.adjustBlock({L}, 1:3);
    % The recovery is asserted against the FORMAL errors, not a hand-picked
    % band: a fixed RelTol here would only encode this fixture's noise draw
    % (the first attempt guessed 1e-6 and CI measured 1.5e-2 on the worst
    % coefficient -- the harmonics are only partially separated over three
    % short days). |x - xTrue| < 4 sigma IS the claim worth making.
    sig = sqrt(diag(R.sol.Cxx)) * R.sigma0;
    z = (R.x - (1:numel(R.x)).' / 10) ./ sig;
    % tol: 4-sigma bright line on 4 coefficients (P(false alarm) < 3e-4)
    tc.verifyLessThan(max(abs(z)), 4, ...
        'the ragged schedule must recover the planted truth');
    % tol: sigma0 from a correct model over ~1000 observations, 1/sqrt(2*dof)
    % scatter ~ 2e-2; band 0.15 is ~7 sd
    tc.verifyEqual(R.sigma0, 1, 'AbsTol', 0.15, ...
        'a correct ragged model must give sigma0 ~ 1');
    res = R.acc.residuals(R.sol, 'verbose', false);
    tc.verifyFalse(res.supersedesSol, ...
        'the two passes must enumerate the ragged schedule identically');
end

%% -------------------------------------------------- WP-52 fixture helpers
function L = raggedLoader(tc, arcPlan, localsPerArc, seen)   %#ok<INUSD>
%RAGGEDLOADER  A loader over a ragged schedule. arcPlan is a cell (one arc-id
%   vector per day) resolved by an @(day,obsType) function, or a static form.
%   `seen` (optional Map) records the arcInfo each buildDay call received.
    if nargin < 4, seen = containers.Map('KeyType','double','ValueType','any'); end
    nc = 4;
    if iscell(arcPlan)
        arcsSpec = @(d, ot) arcPlan{d};
        days = 1:numel(arcPlan);
    else
        arcsSpec = arcPlan;
        days = 1:3;
    end
    args = {'observationTypes', {'pot'}, ...
            'covariance', struct('pot', struct('type','diag','sigma',1e-3)), ...
            'arcs', arcsSpec, 'days', days, 'numCommon', nc, 'systemId', 1};
    if isempty(localsPerArc)
        args = [args, {'numLocal', 1}];
        builder = @(d, ot) raggedDay(d, ot, arcPlan, nc, [], seen);
    else
        args = [args, {'localsPerArc', localsPerArc}];
        builder = @(d, ot, cm, lo, po, ai) raggedDay(d, ot, arcPlan, nc, ai, seen);
    end
    L = lsqvce.engine.makeLoader(builder, args{:});
end

function L = gapLoader(tc)   %#ok<INUSD>
%GAPLOADER  Schedule [1 2 3] on every day; day 2 carries arcs 1 and 3 only and
%   returns an EMPTY cell for arc 2 (the v6 data gap).
    plan = {[1 2 3], [1 3], [1 2 3]};
    L = lsqvce.engine.makeLoader(@(d, ot) raggedDay(d, ot, plan, 4, [], [], [1 2 3]), ...
        'observationTypes', {'pot'}, ...
        'covariance', struct('pot', struct('type','diag','sigma',1e-3)), ...
        'arcs', [1 2 3], 'days', 1:3, 'numCommon', 4, 'numLocal', 1);
end

function L = staticPositionalLoader(tc)   %#ok<INUSD>
%STATICPOSITIONALLOADER  A plain STATIC loader: two arcs a day, addressed by
%   the declared count (the pre-WP-52 behaviour and its guard).
    L = lsqvce.engine.makeLoader(@(d, ot) simpleDay(d, ot), ...
        'observationTypes', {'pot'}, ...
        'covariance', struct('pot', struct('type','diag','sigma',1e-3)), ...
        'arcs', 2, 'days', 1:2, 'numCommon', 4, 'numLocal', 1);
end

function D = raggedDay(day, obsType, arcPlan, nc, arcInfo, seen, schedule)   %#ok<INUSL>
%RAGGEDDAY  One day of the ragged fixture: only the arcs this day carries,
%   every arc a DIFFERENT length, locals placed at the offsets makeLoader
%   handed over (when it owns the plan). v6: the day's cells follow SCHEDULE
%   (default: the ids it carries); a scheduled id it lacks is an empty cell.
    if nargin < 7, schedule = []; end
    if iscell(arcPlan), ids = arcPlan{day}; else, ids = 1:arcPlan; end
    if ~isempty(arcInfo), seen(day) = arcInfo; end                  %#ok<NASGU>
    lens = 40 + 7 * ids + 3 * day;              % every arc a different length
    mTot = sum(lens);
    if isempty(arcInfo), nl = 1; else, nl = arcInfo.nl; end
    A = zeros(mTot, nc);  B = zeros(mTot, nl);  l = zeros(mTot, 1);
    xTrue = (1:nc).' / 10;
    row0 = 0;
    for k = 1:numel(ids)
        rows = row0 + (1:lens(k));
        tau = (0:lens(k)-1).' / lens(k);        % arc-local, for the LOCALS
        % The commons are harmonics of a GLOBAL clock, so they are not
        % absorbed by the per-arc bias/trend -- without this the system is
        % rank deficient by construction (the locals span [1, tau] exactly).
        tg = (row0 + (0:lens(k)-1).' + 300 * day) / 250;
        A(rows, :) = [cos(2*pi*tg), sin(2*pi*tg), cos(4*pi*tg), sin(4*pi*tg)];
        if ~isempty(arcInfo)
            off = arcInfo.offsetOf(ids(k));
            cnt = arcInfo.count(arcInfo.arcId == ids(k));
            B(rows, off + 1) = 1;
            if cnt >= 2, B(rows, off + 2) = tau; end
        else
            B(rows, 1) = 1;
        end
        l(rows) = A(rows,:) * xTrue;
        row0 = row0 + lens(k);
    end
    rs = RandStream('mt19937ar', 'Seed', 5200 + day);
    l = l + 1e-3 * randn(rs, mTot, 1);
    D = struct('A', A, 'B', B, 'l', l, 't', (1:mTot).', ...
               'arcIdx', {arcSlices(lens)}, 'arcId', ids);
    D = arcCells(D, schedule);                       % v6: the cells are the arcs
end

function c = arcSlices(lens)
%ARCSLICES  Row-index vectors for consecutive arcs of the given lengths.
    c = cell(1, numel(lens));  row0 = 0;
    for k = 1:numel(lens)
        c{k} = row0 + (1:lens(k));  row0 = row0 + lens(k);
    end
end

function D = simpleDay(day, obsType)   %#ok<INUSD>
%SIMPLEDAY  Two equal arcs: the plain static case (two cells).
    m = 60;  nc = 4;
    tau = (0:m-1).' / m;
    A = [ones(m,1), tau, cos(2*pi*tau), sin(2*pi*tau)];
    D = struct('A', A, 'B', ones(m,1), 'l', A * ((1:nc).'/10), ...
               't', (1:m).', 'arcIdx', {{1:m/2, m/2+1:m}});
    D = arcCells(D);                                 % v6: the cells are the arcs
end
