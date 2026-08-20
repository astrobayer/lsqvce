function tests = test_adjustBlock
%TEST_ADJUSTBLOCK  Tests for the recover/Results/makeLoader access layer.
%
%     lsqvce.adjustBlock      one-call block-arrow adjustment -> Results
%     lsqvce.Results      solution + diagnostics + provenance object
%     lsqvce.makeLoader   day-builder -> accumulator-ready Loader
%
%   Both the legacy closure path (via lsqvce.synthLoader) and the new
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
    loader = @(d,ot,a,uo) lsqvce.synthLoader(cfg, 1, d, ot, a, uo);
    sys = struct('sysID',1, 'nl',cfg.nl, 'nc',cfg.nc, 'loader',loader);
    R = lsqvce.adjustBlock(sys, 1:3);

    tc.verifyClass(R, 'lsqvce.Results');
    tc.verifyEmpty(R.history);
    tc.verifyTrue(R.converged);      % a single solve IS done, by definition
    tc.verifyError(@() R.convergence(), 'lsqvce:Results:noHistory');
    tc.verifyError(@() R.summary(),     'lsqvce:Results:noHistory');
end

%% ============================================================ recover paths
function testRecoverSynthClosure(tc)
% categories: contract, precision
    cfg = synthCfg();
    loader = @(d,ot,a,uo) lsqvce.synthLoader(cfg, 1, d, ot, a, uo);
    sys = struct('sysID',1, 'nl',cfg.nl, 'nc',cfg.nc, 'loader',loader);
    R = lsqvce.adjustBlock(sys, 1:3);
    tc.verifyClass(R, 'lsqvce.Results');
    tc.verifyClass(R.acc, 'lsqvce.BlockAccumulator');
    tc.verifyLessThan(max(abs(R.x - cfg.x_c)), 0.15);   % recovers the truth
    tc.verifyLessThan(abs(R.sigma0 - 1), 0.15);         % noise model is correct
    tc.verifyTrue(R.converged);
end

function testRecoverMakeLoaderCell(tc)
% categories: precision
    cfg = buildCfg();
    L = lsqvce.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos','kbr'}, 'covariance',cfg.cov, 'arcs',cfg.arcsPerDay, ...
            'systemId',1, 'numLocal',cfg.nl, 'numCommon',cfg.nc, 'days',1:3);
    R = lsqvce.adjustBlock({L}, []);                % days come from the loader
    tc.verifyLessThan(max(abs(R.x - cfg.x_c)), 0.15);
    tc.verifyLessThan(abs(R.sigma0 - 1), 0.15);
end

function testResultsAccessors(tc)
% categories: robustness | tags: diagnostics
    cfg = buildCfg();
    L = lsqvce.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos','kbr'}, 'covariance',cfg.cov, 'arcs',cfg.arcsPerDay, ...
            'systemId',1, 'numLocal',cfg.nl, 'numCommon',cfg.nc);
    R = lsqvce.adjustBlock({L}, 1:2);
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
    L = lsqvce.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
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
    tc.verifyError(@() lsqvce.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
        'covariance',cfg.cov, 'arcs',2), 'lsqvce:makeLoader:noObsTypes');
    tc.verifyError(@() lsqvce.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
        'observationTypes',{'nope'}, 'covariance',cfg.cov), 'lsqvce:makeLoader:noCov');
    bad = @(d,ot) struct('A',randn(10,6));              % missing B,l,t
    Lb = lsqvce.makeLoader(bad, 'observationTypes',{'pos'}, 'covariance',cfg.cov, 'numLocal',4, 'numCommon',6);
    tc.verifyError(@() Lb.closure(1,'pos',1,struct()), 'lsqvce:makeLoader:badDay');
end

%% ==================================================== option names
function testMakeLoaderCanonicalNames(tc)
% categories: contract
    % Every option has ONE name. Building with them must warn about nothing.
    cfg = buildCfg();
    ws = warning('off', 'all'); c = onCleanup(@() warning(ws)); %#ok<NASGU>
    lastwarn('', '');
    L = lsqvce.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
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
    tc.verifyError(@() lsqvce.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos'}, 'covariance',cfg.cov, ...
            'sysID',3, 'numLocal',4, 'numCommon',6), ...
        'MATLAB:TooManyInputs');
end

%% ============================================================ recover errors
function testRecoverErrors(tc)
% categories: contract
    tc.verifyError(@() lsqvce.adjustBlock({1}, 1:3), 'lsqvce:adjustBlock:badCell');
    tc.verifyError(@() lsqvce.adjustBlock(struct('sysID',1,'nl',4,'nc',6), 1:3), ...
                   'lsqvce:adjustBlock:badStruct');
end

function testShapeIsValidatedBeforeNcIsDerived(tc)
% categories: contract | tags: diagnostics
    % REGRESSION GUARD (2026-07-18). adjustBlock derives nc with lsqvce.resolveNc
    % and validates the shape of `systems` with normalizeSystems. The ORDER is
    % load-bearing: deriving nc first makes every malformed `systems` report the
    % generic lsqvce:resolveNc:missing ("no system declares a common-parameter
    % count") instead of the specific id that names what is actually wrong.
    %
    % This shipped as a real bug -- the required order was written in the
    % handover but never applied to the code -- so it is pinned here rather than
    % left to prose. The distinguishing feature of these inputs is that they are
    % malformed AND declare no nc, so whichever step runs first decides the id.
    tc.verifyError(@() lsqvce.adjustBlock({1}, 1:3), ...
        'lsqvce:adjustBlock:badCell', ...
        'a non-Loader cell entry must report badCell, not resolveNc:missing');

    tc.verifyError(@() lsqvce.adjustBlock(struct('sysID',1,'nl',4), 1:3), ...
        'lsqvce:adjustBlock:badStruct', ...
        'a spec missing .loader must report badStruct, not resolveNc:missing');

    % The complement: a WELL-FORMED systems argument that genuinely declares no
    % nc must still reach resolveNc and report the nc problem.
    tc.verifyError(@() lsqvce.adjustBlock(struct('sysID',1,'nl',4, ...
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
    L = lsqvce.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
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
    grab = @(d,ot,c,l,p) captureBuild(cfg,d,ot,c,l,p);
    L = lsqvce.makeLoader(grab, 'observationTypes',{'pos'}, 'covariance',cfg.cov, 'arcs',2, ...
            'commonParams', common, 'localParams', local, 'processingOptions', po);
    D = L.day(1, 'pos');                          % whole day carries the echoes
    tc.verifyEqual(D.gotCommon, common);
    tc.verifyEqual(D.gotLocal,  local);
    tc.verifyEqual(D.gotPO.flag, 42);
end

function testParamDerivesNcNl(tc)
% categories: contract
    % nc = rows(common), nl = rows(local), carried on L.
    cfg = buildCfg();
    L = lsqvce.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos'}, 'covariance',cfg.cov, ...
            'commonParams', paramRows(6), 'localParams', paramRows(4));
    tc.verifyEqual(L.nc, 6);
    tc.verifyEqual(L.nl, 4);
end

function testParamNlConflictErrors(tc)
% categories: contract
    % Passing both 'nl' and 'local' with disagreeing counts must error.
    cfg = buildCfg();
    tc.verifyError(@() lsqvce.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
        'observationTypes',{'pos'}, 'covariance',cfg.cov, 'numLocal',5, 'localParams', paramRows(4)), ...
        'lsqvce:makeLoader:nlConflict');
end

function testParamColumnCountValidated(tc)
% categories: contract
    % A parametrisation array without 14 columns must be rejected.
    cfg = buildCfg();
    tc.verifyError(@() lsqvce.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
        'observationTypes',{'pos'}, 'covariance',cfg.cov, 'commonParams', cell(6,10)), ...
        'lsqvce:makeLoader:paramCols');
end

%% ====================================================== plot-hook options
function testPlotHookSingleHandleBroadcast(tc)
% categories: contract, robustness
    % A single handle is broadcast to every obsType of the loader.
    cfg = buildCfg();
    h = @(ctx) [];
    L = lsqvce.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos','kbr'}, 'covariance',cfg.cov, 'numLocal',4, 'numCommon',6, ...
            'plotPostSolve', h);
    tc.verifyClass(L.plotPostSolve, 'containers.Map');
    tc.verifyTrue(isKey(L.plotPostSolve, 'pos') && isKey(L.plotPostSolve, 'kbr'));
    tc.verifyEqual(L.plotPostModel.Count, uint64(0));   % none given -> empty Map
end

function testPlotHookPerObsTypeMap(tc)
% categories: contract
    % A per-obsType map is passed through, keeping only the named obsTypes.
    cfg = buildCfg();
    m = containers.Map({'pos'}, {@(ctx) []});
    L = lsqvce.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos','kbr'}, 'covariance',cfg.cov, 'numLocal',4, 'numCommon',6, ...
            'plotPostModel', m);
    tc.verifyTrue(isKey(L.plotPostModel, 'pos'));
    tc.verifyFalse(isKey(L.plotPostModel, 'kbr'));       % not provided
end

function testPlotHookBadHandleErrors(tc)
% categories: contract
    % A non-handle value in a map must be rejected.
    muteWarn(tc, 'lsqvce:adjustBlock:plotHookFailed');   % provoked on purpose here (WP-47, gate rule 7)
    cfg = buildCfg();
    m = containers.Map({'pos'}, {42});
    tc.verifyError(@() lsqvce.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
        'observationTypes',{'pos'}, 'covariance',cfg.cov, 'numLocal',4, 'numCommon',6, 'plotPostSolve', m), ...
        'lsqvce:makeLoader:badPlotHandle');
end

%% ======================================== recover postSolve plot hooks
function testRecoverFiresPostSolveHook(tc)
% categories: contract
    % A postSolve hook on the loader must fire once per obsType after recover's
    % single solve, with a well-formed ctx. The hook records into a handle
    % object and returns [] (no figure), so it needs no display.
    muteWarn(tc, 'lsqvce:adjustBlock:plotHookFailed');   % provoked on purpose here (WP-47, gate rule 7)
    cfg = buildCfg();
    rec = recorderHandle();
    h   = @(ctx) recordCtx(rec, ctx);
    L = lsqvce.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos','kbr'}, 'covariance',cfg.cov, 'numLocal',4, 'numCommon',6, ...
            'plotPostSolve', h);
    R = lsqvce.adjustBlock({L}, 1:2);
    tc.verifyClass(R, 'lsqvce.Results');
    % Both obsTypes fired exactly once.
    logged = recorderLog(rec);
    tc.verifyEqual(numel(logged), 2);
    stages = cellfun(@(c) c.stage, logged, 'uni', 0);
    tc.verifyTrue(all(strcmp(stages, 'postSolve')));
    obsSeen = sort(cellfun(@(c) c.obsType, logged, 'uni', 0));
    tc.verifyEqual(obsSeen, {'kbr','pos'});
    % ctx is well-formed: pass==1, stream key 'satKey|obsType', carries R.
    c1 = logged{1};
    tc.verifyEqual(c1.pass, 1);
    tc.verifyTrue(contains(c1.stream, '|'));
    tc.verifyEqual(c1.stream, sprintf('%s|%s', lsqvce.satKeyOf(c1.sysID), c1.obsType));
    tc.verifyEqual(c1.x, R.x);
end

function testRecoverRawClosureIgnoresPlotOpts(tc)
% categories: contract, precision
    % Loader contract: a RAW 4-arg closure carries none of the optional Loader
    % fields, so the postSolve plot seam must be a silent no-op even when
    % ShowSeconds/ExportDir are set -- it must never look for .plotPostSolve on a
    % bare handle and error. (Guards the "raw closures are first-class" claim.)
    cfg = synthCfg();
    loader = @(d,ot,a,uo) lsqvce.synthLoader(cfg, 1, d, ot, a, uo);
    sys = struct('sysID',1, 'nl',cfg.nl, 'nc',cfg.nc, 'loader',loader);
    R = lsqvce.adjustBlock(sys, 1:3, 'ShowSeconds', 0, 'ExportDir', '');
    tc.verifyClass(R, 'lsqvce.Results');
    tc.verifyLessThan(max(abs(R.x - cfg.x_c)), 0.15);
end

function testRecoverNoHookIsNoOp(tc)
% categories: contract
    % A loader with no plot hook (the common case) must solve exactly as before
    % and never touch showExport. Just assert it runs and returns Results.
    cfg = buildCfg();
    L = lsqvce.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos'}, 'covariance',cfg.cov, 'numLocal',4, 'numCommon',6);
    R = lsqvce.adjustBlock({L}, 1:2);
    tc.verifyClass(R, 'lsqvce.Results');
    tc.verifyFalse(isKey(L.plotPostSolve, 'pos'));   % none registered
end

function testRecoverHookFailureNeverBreaksSolve(tc)
% categories: contract, precision, robustness
    % A throwing hook must be logged and skipped -- the solve still returns a
    % valid Results. Mirrors the adjustVCE guarantee.
    cfg = buildCfg();
    h = @(ctx) error('boom:intentional', 'intentional');
    L = lsqvce.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos'}, 'covariance',cfg.cov, 'numLocal',4, 'numCommon',6, ...
            'plotPostSolve', h);
    R = lsqvce.adjustBlock({L}, 1:2);   % must not throw
    tc.verifyClass(R, 'lsqvce.Results');
    tc.verifyLessThan(max(abs(R.x - cfg.x_c)), 0.2);
end

%% =================================================================== helpers
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

function D = captureBuild(cfg, day, obsType, common, local, procOpt)
% A 5-arg builder that both builds a day AND echoes what it received, so the
% test can assert the parametrisation reached buildDay intact.
    D = buildDayImpl(cfg, day, obsType);
    D.gotCommon = common;
    D.gotLocal  = local;
    D.gotPO     = procOpt;
end

function cfg = synthCfg()
% cfg in the exact shape lsqvce.synthLoader expects.
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
        lsqvce.adjustBlock({}, [], 'Backend', 'QR');   % enum VALUE casing
    catch e
        ok = startsWith(e.identifier, 'MATLAB:');
    end
    verifyTrue(tc, ok, ...
        'Backend value ''QR'' must be rejected (case-sensitive mustBeMember)');
    ok = false;
    try
        lsqvce.adjustBlock({}, [], 'RegLambda', -1);   % typed constraint
    catch e
        ok = startsWith(e.identifier, 'MATLAB:');
    end
    verifyTrue(tc, ok, 'negative RegLambda must be rejected by the signature');
    ok = false;
    try
        lsqvce.adjustBlock({}, [], 'Backend', 'qr', 'nosuchopt', 1);
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
        @() lsqvce.adjustBlock(sys, 1, 'backend', 'qr'), ...     % old casing
        @() lsqvce.adjustBlock(sys, 1, 'Backend', 'QR'), ...     % enum casing
        @() lsqvce.adjustBlock(sys, 1, 'Backend', 'cholesky'),...% bad enum
        @() lsqvce.adjustBlock(sys, 1, 'RegLambda', -1), ...     % sign
        @() lsqvce.adjustBlock(sys, 1, 'NoSuchOption', 1) };
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
    L = lsqvce.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos','kbr'}, 'covariance',cfg.cov, ...
            'arcs',cfg.arcsPerDay, 'systemId',1, ...
            'numLocal',cfg.nl, 'numCommon',cfg.nc);
    R = lsqvce.adjustBlock({L}, 1:2);

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
    mk = @() lsqvce.adjustBlock({lsqvce.makeLoader( ...
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
    L = lsqvce.makeLoader(@(d,ot) buildDayImpl(cfg,d,ot), ...
            'observationTypes',{'pos','kbr'}, 'covariance',cfg.cov, ...
            'arcs',cfg.arcsPerDay, 'systemId',1, ...
            'numLocal',cfg.nl, 'numCommon',cfg.nc);
    R = lsqvce.adjustBlock({L}, 1:2);

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
