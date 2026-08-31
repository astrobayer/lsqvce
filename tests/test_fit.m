function tests = test_fit
%TEST_FIT  Tests for lsqvce.fit, the one user verb (R4 WP-56, decision R4-D4).
%
%     lsqvce.fit   data struct (DATA_CONTRACT.md) + lsqvce.Options -> Results;
%                  builds the loaders through makeLoader and dispatches on
%                  opts.Model to adjustBlock / adjustVCE
%
%   Claims under test: (1) EXACTNESS of the dispatch -- fit reproduces the
%   direct driver call to solver precision; (2) every rule of the data
%   contract fires its documented lsqvce:fit id BEFORE any data is read;
%   (3) the declaration is order-free (fields, systems); (4) nl = 0 and raw
%   loaders pass through.
%
%   Run: runtests('tests/test_fit')
%
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
    Rd = lsqvce.adjustBlock({L}, data.days, a.nv{:});
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

function testFitVceReproducesTheFluentChain(tc)
% categories: reference, precision | tags: diagnostics
    % WP-58 parity: fit(data, Options(Model="vce", ...)) versus the fluent
    % configurator chain it replaces. Exact trace mode (vce default) makes
    % the separate loop deterministic; two outer passes keep the CI cost low.
    [data, cfg] = declaration();
    opts = lsqvce.Options(Model="vce", MaxOuter=2, VceMaxIter=3, Verbose=false);
    Rf = lsqvce.fit(data, opts);
    L  = loaderOf(cfg, data.systems(1), data);
    Rd = lsqvce.adjustVCE({L}, data.days, 'MaxOuter', 2, 'Verbose', false) ...
            .varianceComponents('MaxIter', 3) ...
            .run();
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
    d = data;  d.common.nc = 6;
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:unknownField');
    d = rmfield(data, 'systems');
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:noSystems');
    d = data;  d.systems = struct([]);
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:noSystems');
    d = data;  d.systems = rmfield(d.systems, 'covariance');
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:missingField');
    d = data;  d.systems(1).systemId = 'CD';
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:missingField');
    d = data;  d.systems(1).localsPerArc = 2;                    % now TWO ways
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:localCount');
    d = data;  d.systems(1).numLocal = [];                       % now NO way
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:localCount');
    d = data;  d.common.commonParams = cell(6, 14);              % both ways
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:commonCount');
    d = rmfield(data, 'common');
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:commonCount');
    d = rmfield(data, 'days');
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:noDays');
    d = data;  d.days = [1 2.5];
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:noDays');
    d = data;  d.systems(2) = d.systems(1);                      % same systemId
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:duplicateSystem');
    d = data;  d.systems(1).loader = bomb;                       % both kinds
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:builderKind');
    d = data;  d.systems(1).buildDay = 'notAHandle';
    tc.verifyEqual(idOf(@() lsqvce.fit(d)), 'lsqvce:fit:builderKind');
    d = data;  d.systems(2) = d.systems(1);  d.systems(2).systemId = 2;
    d.systems(2).buildDay = [];  d.systems(2).loader = bomb;    % mixed kinds
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
    % Model "ellmer" and the batch-path Robust are Options the engine space
    % has but fit does not route yet (WP-56b): an explicit refusal, never a
    % silent fallback to another model.
    [data, ~] = declaration();
    tc.verifyEqual(idOf(@() lsqvce.fit(data, lsqvce.Options(Model="ellmer"))), 'lsqvce:fit:notYet');
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
    % G-07 via the front door: numLocal = 0 with B = zeros(m, 0).
    [data, cfg] = declaration();
    cfg.nl = 0;
    m = cfg.arcsPerDay * cfg.mArc;
    data.systems(1).buildDay = @(d, ot) setfield(buildDayImpl(cfg, d, ot), 'B', zeros(m, 0)); %#ok<SFLD>
    data.systems(1).numLocal = 0;
    R = lsqvce.fit(data);
    tc.verifyClass(R, 'lsqvce.Results');
    tc.verifyEqual(numel(R.x), cfg.nc);
    tc.verifyLessThan(max(abs(R.x - cfg.x_c)), 0.15);
end

function testRawLoaderEscapeHatch(tc)
% categories: contract, reference
    % ENGINE tier: .loader (raw 4-arg closure) instead of .buildDay maps onto
    % adjustBlock's struct-array form. Same call, same result.
    cfg = synthCfg();
    raw = @(d, ot, a, uo) lsqvce.synthLoader(cfg, 1, d, ot, a, uo);
    data.common.numCommon = cfg.nc;
    data.days = 1:3;
    s.loader = raw;  s.observationTypes = cfg.obsTypes;  s.covariance = cfg.cov;
    s.systemId = 1;  s.numLocal = cfg.nl;
    data.systems = s;
    Rf = lsqvce.fit(data);
    Rd = lsqvce.adjustBlock(struct('sysID', 1, 'nl', cfg.nl, 'nc', cfg.nc, 'loader', raw), 1:3);
    % tol: 1e-12 -- identical code path (see the reference tests)
    tc.verifyEqual(Rf.x, Rd.x, 'RelTol', 1e-12);
    tc.verifyEqual(Rf.sigma0, Rd.sigma0, 'RelTol', 1e-12);
end

%% ================================================================ helpers
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
% The canonical single-system declaration over the makeLoader-path fixture.
    cfg = buildCfg();
    data = struct();
    data.common.numCommon = cfg.nc;
    data.days = 1:3;
    s = struct();
    s.buildDay = @(d, ot) buildDayImpl(cfg, d, ot);
    s.observationTypes = {'pos', 'kbr'};
    s.covariance = cfg.cov;
    s.systemId = 1;
    s.numLocal = cfg.nl;
    s.arcs = cfg.arcsPerDay;
    data.systems = s;
end

function L = loaderOf(~, s, data)
% What fit builds for one system -- written out, so the tests compare two
% independent spellings of the same declaration.
    L = lsqvce.makeLoader(s.buildDay, 'ObservationTypes', s.observationTypes, ...
            'Covariance', s.covariance, 'SystemId', s.systemId, 'NumLocal', s.numLocal, ...
            'Arcs', s.arcs, 'Days', data.days, 'NumCommon', data.common.numCommon);
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

function cfg = synthCfg()
% cfg in the exact shape lsqvce.synthLoader expects (mirrors test_adjustBlock).
    cfg = struct('nc',6, 'nl',4, 'arcsPerDay',2, 'mArc',300, 'dt',5, ...
                 'seed',42, 'collinear',false, 'collinearEps',1e-7);
    cfg.obsTypes = {'pos','kbr'};
    rs = RandStream('mt19937ar','Seed',3);
    cfg.x_c = randn(rs, 6, 1);
    cfg.y   = containers.Map(1, randn(rs, 4, 1));
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
