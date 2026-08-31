function tests = test_streamUtils
%TEST_STREAMUTILS  Direct tests for the functions moved into +lsqvce this session:
%
%     lsqvce.engine.selectArAcv         automatic AR-order selection by ACV fidelity
%     lsqvce.kernel.whitenessCorr          whiteness statistic
%     lsqvce.engine.collectResiduals    stack a residual field over arcs
%     lsqvce.kernel.residualAcv temporal ACV of D^-1 v (pos + rr paths)
%     lsqvce.engine.normalizeSatList    satellite-list normalisation
%     lsqvce.engine.resolveSatSpecs     per-satellite spec resolution
%
%   Reference properties/tolerances established in validation/v9verify.py (14).
%
%   Run:
%     results = runtests('tests/test_streamUtils');  disp(table(results))
%
%   See also LSQVCE.ENGINE.SELECTARACV, LSQVCE.KERNEL.RESIDUALACV, TEST_V8
    tests = functiontests(localfunctions);
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

%% ===================================================== normalizeSatList
function testNormalizeSatList(tc)
% categories: contract, robustness
    tc.verifyEqual(lsqvce.engine.normalizeSatList([1 2]), {1, 2});
    tc.verifyEqual(lsqvce.engine.normalizeSatList('CD'), {'C', 'D'});
    tc.verifyEqual(lsqvce.engine.normalizeSatList({'C','D'}), {'C', 'D'});
    tc.verifyEqual(lsqvce.engine.normalizeSatList("C"), {'C'});
    tc.verifyError(@() lsqvce.engine.normalizeSatList([]), 'lsqvce:normalizeSatList:empty');
end

%% ===================================================== resolveSatSpecs
function testResolveSatSpecsOverrides(tc)
% categories: contract
    defaults = struct('nl', 12, 'maxLag', 128);
    us(1).satId = 'D';  us(1).nl = 9;               % override only D.nl
    [S, specFor] = lsqvce.engine.resolveSatSpecs({'C','D'}, us, defaults);

    tc.verifyEqual(specFor('C').nl, 12);            % default kept
    tc.verifyEqual(specFor('D').nl, 9);             % overridden
    tc.verifyEqual(specFor('C').maxLag, 128);       % default kept
    tc.verifyEqual(S(1).key, 's_C');                % via lsqvce.engine.satKeyOf
    tc.verifyEqual(S(2).key, 's_D');
    tc.verifyError(@() specFor('Z'), 'lsqvce:resolveSatSpecs:noSpec');
end

function testResolveSatSpecsNumericKeyAndErrors(tc)
% categories: contract
    [S, ~] = lsqvce.engine.resolveSatSpecs({1, 2}, [], struct('nl', 5));
    tc.verifyEqual(S(1).key, 's_1');                % numeric id -> valid field
    tc.verifyEqual(S(1).nl, 5);
    bad.notSatId = 1;
    tc.verifyError(@() lsqvce.engine.resolveSatSpecs({1}, bad, struct('nl',5)), ...
                   'lsqvce:resolveSatSpecs:badSatSpec');
end

%% ===================================================== satKeyOf (public)
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
