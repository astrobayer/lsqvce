function tests = test_streamUtils
%TEST_STREAMUTILS  Direct tests for the functions moved into +lsqvce this session:
%
%     lsqvce.engine.selectArAcv         automatic AR-order selection by ACV fidelity
%     lsqvce.kernel.whitenessCorr          whiteness statistic
%     lsqvce.engine.collectResiduals    stack a residual field over arcs
%     lsqvce.kernel.residualAcv temporal ACV of D^-1 v (pos + rr paths)
%
%   Reference properties/tolerances established in validation/v9verify.py (14).
%
%   Run:
%     results = runtests('tests/test_streamUtils');  disp(table(results))
%
%   See also LSQVCE.ENGINE.SELECTARACV, LSQVCE.KERNEL.RESIDUALACV, TEST_V8
    tests = functiontests(localfunctions);
end

function setupOnce(~)
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(root);  addpath(fullfile(root, 'examples'));   % fixture_grace (moved tests, R5 M4 4c)
end

%% ===================================================== selectArAcv
function testSelectArAcvRecoversOrder(tc)
% categories: precision | tags: diagnostics
    a0 = [1 -1.6 0.9 -0.25];  s0 = 0.7;
    r  = s0^2 * lsqvce.kernel.armaACV(a0, 1, 200);
    [a, s2, p, sel] = lsqvce.engine.selectArAcv(r, 5000);
    tc.verifyEqual(p, 3);
    % tol: a(1)=1 exact by construction (1e-12 house); coefficient/sigma
    % recovery at 1e-6 = v9verify.py measured level on the exact ACV
    tc.verifyEqual(a(1), 1, 'AbsTol', 1e-12);
    tc.verifyEqual(a, a0, 'AbsTol', 1e-6);
    tc.verifyEqual(sqrt(s2), s0, 'AbsTol', 1e-6);         % tol: same source
    tc.verifyEqual(sel.verdict, 'good');
    tc.verifyLessThan(sel.acvRelErr, 1e-3);
    tc.verifyEqual(numel(sel.aic), 121);            % ladder 0..pMax(=120)
end

function testSelectArAcvWhiteReturnsZero(tc)
% categories: precision | tags: diagnostics
    r = [1; zeros(49,1)];                           % delta ACV = white
    [a, s2, p, sel] = lsqvce.engine.selectArAcv(r, 5000);
    tc.verifyEqual(p, 0);
    tc.verifyEqual(a, 1);
    % tol: delta ACV -> s2 = r0 = 1 exactly; house round-trip bound 1e-12
    tc.verifyEqual(s2, 1, 'AbsTol', 1e-12);
    tc.verifyEqual(sel.verdict, 'white');
end

function testSelectArAcvDefaultsAndValidation(tc)
% categories: contract
    r = 0.5^2 * lsqvce.kernel.armaACV([1 -0.6], 1, 64);
    [~, ~, p] = lsqvce.engine.selectArAcv(r, []);          % nData [] -> default 1e4
    tc.verifyEqual(p, 1);
    tc.verifyError(@() lsqvce.engine.selectArAcv([], 100), 'MATLAB:selectArAcv:expectedVector');
end

%% ===================================================== whitenessCorr
function testMaxLagCorrWhiteVsCorrelated(tc)
% categories: contract
    rs = RandStream('mt19937ar','Seed',9);
    w  = randn(rs, 4000, 1);
    tc.verifyLessThan(lsqvce.kernel.whitenessCorr(w, 5), 3/sqrt(4000));
    ar1 = filter(1, [1 -0.8], randn(rs, 4000, 1));
    mc  = lsqvce.kernel.whitenessCorr(ar1, 5);
    tc.verifyGreaterThan(mc, 3/sqrt(4000));
    tc.verifyGreaterThanOrEqual(mc, 0);
    tc.verifyLessThanOrEqual(mc, 1);
end

function testMaxLagCorrShortSeries(tc)
% categories: contract
    tc.verifyEqual(lsqvce.kernel.whitenessCorr(7, 5), 0);     % 1 element -> 0
end

%% ===================================================== collectResiduals
function testCollectResidualsFilters(tc)
% categories: robustness
    res = makeRes();
    wC = lsqvce.engine.collectResiduals(res, 'C', 'pos', 'w');
    tc.verifyEqual(numel(wC), 18);                  % 2 arcs x 3 epochs x 3 axes (C/pos)
    vAll = lsqvce.engine.collectResiduals(res, [], 'pos', 'v');
    tc.verifyEqual(numel(vAll), 36);                % C+D pooled (4 arcs x 9)
    empty = lsqvce.engine.collectResiduals(res, 'C', 'nope');
    tc.verifyEmpty(empty);
    wDef = lsqvce.engine.collectResiduals(res, 'D', 'pos'); % default field 'w'
    tc.verifyEqual(numel(wDef), 18);
end

function testCollectResidualsRejectsMissingField(tc)
% categories: contract
% WP-43: the documented rejection lsqvce:collectResiduals:badRes was never
% exercised directly (audit check 16 found the gap). A residual structure
% whose matching arcs lack the requested field must fire that id -- and only
% that id -- naming the offending field.
    res = makeRes();
    % try/catch capture (WP-39 pattern): verifyError's output-argument form
    % errored on the function-based TestCase in CI (R2026a hosted runner) --
    % round-1 red, both new err-capturing tests. Captured lesson: assert the
    % id AND the message through the caught exception, never via the output.
    caught = MException.empty;
    try
        lsqvce.engine.collectResiduals(res, 'C', 'pos', 'zz');
    catch e
        caught = e;
    end
    tc.verifyNotEmpty(caught, 'the absent field must be rejected');
    tc.verifyEqual(caught.identifier, 'lsqvce:collectResiduals:badRes');
    tc.verifySubstring(caught.message, 'zz');   % the message names the field
end

%% ===================================================== residualAcv
function testResidualAcvRejectsBadDayFcnAndEmptyStream(tc)
% categories: contract
% WP-43: both documented rejections of residualAcv, unexercised until the
% check-16 count. dayFcn must be callable; a stream with no matching
% observations must refuse with noObs rather than return an empty ACV.
    rs = RandStream('mt19937ar','Seed',7);
    [res, dayFcn] = makePosRes(rs, eye(3), 20, 1);
    tc.verifyError(@() lsqvce.kernel.residualAcv(res, 'C', 'pos', 3, 5, 42), ...
                   'lsqvce:residualAcv:dayFcn');
    tc.verifyError(@() lsqvce.kernel.residualAcv(res, 'C', 'nope', 3, 5, dayFcn), ...
                   'lsqvce:residualAcv:noObs');
end

function testDescaledAcvPosWhitens(tc)
% categories: precision
% pos: v = L w (coloured only instantaneously by D=LL'), de-scale -> white ACV.
    rs = RandStream('mt19937ar','Seed',3);
    d = 3;  nE = 600;  maxLag = 40;  L = [1 0 0; 0.3 0.8 0; -0.2 0.1 0.6];
    [res, dayFcn] = makePosRes(rs, L, nE, 4);       % 4 arcs
    [acv, nSamp] = lsqvce.kernel.residualAcv(res, 'C', 'pos', d, maxLag, dayFcn);
    % tol: statistical identity on a finite sample; scatter ~ 1/sqrt(dof) of the fixture, bar sized 4/sqrt(Z)-style (HANDOVER 6.1.2)
    tc.verifyEqual(acv(1), 1, 'RelTol', 0.1);       % r0 ~ 1
    tc.verifyLessThan(abs(acv(2))/acv(1), 0.1);     % |r1/r0| small
    tc.verifyEqual(nSamp, 4*d*nE);
end

function testDescaledAcvRrWhitens(tc)
% categories: precision
% rr: v = sigma .* w, de-scale by sigma -> white.
    rs = RandStream('mt19937ar','Seed',5);
    nE = 800;  maxLag = 30;
    [res, dayFcn] = makeRrRes(rs, nE, 3);
    [acv, nSamp] = lsqvce.kernel.residualAcv(res, 'C', 'rr', 1, maxLag, dayFcn);
    % tol: statistical identity on a finite sample; scatter ~ 1/sqrt(dof) of the fixture, bar sized 4/sqrt(Z)-style (HANDOVER 6.1.2)
    tc.verifyEqual(acv(1), 1, 'RelTol', 0.1);
    tc.verifyLessThan(abs(acv(2))/acv(1), 0.1);
    tc.verifyEqual(nSamp, 3*nE);
end

function testSatKeyOf(tc)
% categories: contract
    tc.verifyEqual(lsqvce.engine.satKeyOf('C'), 's_C');
    tc.verifyEqual(lsqvce.engine.satKeyOf("D"), 's_D');
    tc.verifyEqual(lsqvce.engine.satKeyOf(1),   's_1');
    tc.verifyTrue(isvarname(lsqvce.engine.satKeyOf(7)));   % always a legal field name
end

%% =================================================================== helpers
function res = makeRes()
% Minimal res: C has 2 pos arcs (3 epochs each), D has 2 pos arcs.
    res.sys(1).sysID = 'C';
    res.sys(1).obs   = [arc('pos',1,1,3), arc('pos',1,2,3)];
    res.sys(2).sysID = 'D';
    res.sys(2).obs   = [arc('pos',1,1,3), arc('pos',1,2,3)];
end

function a = arc(obsType, day, arcId, nE)
    a = struct('obsType',obsType, 'day',day, 'arc',arcId, ...
               'v', (1:3*nE).', 'w', (1:3*nE).');
end

function [res, dayFcn] = makePosRes(rs, L, nE, nArc)
    d = size(L,1);  Cep = zeros(nE, d, d);
    for i = 1:nE, Cep(i,:,:) = L*L.'; end
    obs = repmat(arc('pos',1,1,nE), 1, nArc);
    for k = 1:nArc
        w = randn(rs, d, nE);  v = L*w;             % colour by D only
        obs(k) = struct('obsType','pos', 'day',1, 'arc',k, ...
                        'v', reshape(v, [], 1), 'w', reshape(v, [], 1));
    end
    res.sys(1).sysID = 'C';  res.sys(1).obs = obs;
    dayStruct.Cep = Cep;
    for k = 1:nArc, dayStruct.arc(k).epochs = 1:nE; dayStruct.arc(k).rows = 1:nE; end
    dayFcn = @(s,dd) dayStruct;
end

function [res, dayFcn] = makeRrRes(rs, nE, nArc)
    sigma = 0.5 + rand(rs, nE, 1);                  % per-epoch std
    obs = repmat(arc('rr',1,1,nE), 1, nArc);
    for k = 1:nArc
        w = randn(rs, nE, 1);  v = sigma .* w;
        obs(k) = struct('obsType','rr', 'day',1, 'arc',k, 'v', v, 'w', v);
    end
    res.sys(1).sysID = 'C';  res.sys(1).obs = obs;
    dayStruct.sigma = sigma;
    for k = 1:nArc, dayStruct.arc(k).rows = 1:nE; dayStruct.arc(k).epochs = 1:nE; end
    dayFcn = @(s,dd) dayStruct;
end

%% ================= moved from test_adjustVCE (R5 M4 4c): loader / covModel
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
    [~, Wfcn, info] = lsqvce.engine.covModel(cm, t);
    verifyEqual(tc, info.nAxis, d);
    A = randn(rs, d*nE, 2);
    [Aw,~,ew] = Wfcn(1:d*nE, A, zeros(d*nE,0), e);
    % whitened residual variance ~ 1 (unit innovation after full whitening)
    % tol: statistical identity on a finite sample; scatter ~ 1/sqrt(dof) of the fixture, bar sized 4/sqrt(Z)-style (HANDOVER 6.1.2)
    verifyEqual(tc, var(ew), 1, 'RelTol', 0.25);
    verifySize(tc, Aw, [d*nE, 2]);
end

function test_covOverride_bySys(tc)
% categories: precision
    % A bySys override for the merged pair system (sysID 1) must rescale the
    % whitener. bySys is keyed by makeValidName(string(sysID)) -> 'x1'.
    [cfg, L] = fixture_grace();
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
    [cfg, L] = fixture_grace();
    Lcd = L{1};
    opt = struct('sysID',1, 'covOverride', struct('rr', struct('sigma', 5)));
    arc = Lcd.closure(cfg.days(1), 'rr', 1, opt);
    verifyTrue(tc, isa(arc.Wfcn, 'function_handle'));
end

function test_epochCov_hook_carried(tc)
% categories: contract
    hook = @(day, obsType) struct('day', day, 'obsType', obsType);
    cov.rr = struct('type','diag','sigma',1);
    L = lsqvce.engine.makeLoader(@(d,ot) makeTinyDay(), 'observationTypes', {'rr'}, ...
            'covariance', cov, 'arcs', 1, 'systemId', 1, 'numLocal', 2, 'epochCovariance', hook);
    verifyTrue(tc, isfield(L, 'epochCov'));
    out = L.epochCov(3, 'rr');
    verifyEqual(tc, out.day, 3);
end

%% ===================================================== SCALE (VCE) recovers weights

function D = makeTinyDay()
    m = 20;
    D.A = [ones(m,1), (1:m).'];
    D.B = 0.1*ones(m,2);
    D.l = randn(m,1);
    D.t = (0:m-1).' * 5;
    D.arcIdx = {1:m};
end

