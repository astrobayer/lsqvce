function tests = test_streamUtils
%TEST_STREAMUTILS  Direct tests for the functions moved into +lsqvce this session:
%
%     lsqvce.selectArAcv         automatic AR-order selection by ACV fidelity
%     lsqvce.whitenessCorr          whiteness statistic
%     lsqvce.collectResiduals    stack a residual field over arcs
%     lsqvce.residualAcv temporal ACV of D^-1 v (pos + rr paths)
%     lsqvce.normalizeSatList    satellite-list normalisation
%     lsqvce.resolveSatSpecs     per-satellite spec resolution
%
%   Reference properties/tolerances established in validation/v9verify.py (14).
%
%   Run:
%     results = runtests('tests/test_streamUtils');  disp(table(results))
%
%   See also LSQVCE.SELECTARACV, LSQVCE.RESIDUALACV, TEST_V8
    tests = functiontests(localfunctions);
end

%% ===================================================== selectArAcv
function testSelectArAcvRecoversOrder(tc)
    a0 = [1 -1.6 0.9 -0.25];  s0 = 0.7;
    r  = s0^2 * lsqvce.armaACV(a0, 1, 200);
    [a, s2, p, sel] = lsqvce.selectArAcv(r, 5000);
    tc.verifyEqual(p, 3);
    tc.verifyEqual(a(1), 1, 'AbsTol', 1e-12);
    tc.verifyEqual(a, a0, 'AbsTol', 1e-6);
    tc.verifyEqual(sqrt(s2), s0, 'AbsTol', 1e-6);
    tc.verifyEqual(sel.verdict, 'good');
    tc.verifyLessThan(sel.acvRelErr, 1e-3);
    tc.verifyEqual(numel(sel.aic), 121);            % ladder 0..pMax(=120)
end

function testSelectArAcvWhiteReturnsZero(tc)
    r = [1; zeros(49,1)];                           % delta ACV = white
    [a, s2, p, sel] = lsqvce.selectArAcv(r, 5000);
    tc.verifyEqual(p, 0);
    tc.verifyEqual(a, 1);
    tc.verifyEqual(s2, 1, 'AbsTol', 1e-12);
    tc.verifyEqual(sel.verdict, 'white');
end

function testSelectArAcvDefaultsAndValidation(tc)
    r = 0.5^2 * lsqvce.armaACV([1 -0.6], 1, 64);
    [~, ~, p] = lsqvce.selectArAcv(r, []);          % nData [] -> default 1e4
    tc.verifyEqual(p, 1);
    tc.verifyError(@() lsqvce.selectArAcv([], 100), 'MATLAB:selectArAcv:expectedVector');
end

%% ===================================================== whitenessCorr
function testMaxLagCorrWhiteVsCorrelated(tc)
    rs = RandStream('mt19937ar','Seed',9);
    w  = randn(rs, 4000, 1);
    tc.verifyLessThan(lsqvce.whitenessCorr(w, 5), 3/sqrt(4000));
    ar1 = filter(1, [1 -0.8], randn(rs, 4000, 1));
    mc  = lsqvce.whitenessCorr(ar1, 5);
    tc.verifyGreaterThan(mc, 3/sqrt(4000));
    tc.verifyGreaterThanOrEqual(mc, 0);
    tc.verifyLessThanOrEqual(mc, 1);
end

function testMaxLagCorrShortSeries(tc)
    tc.verifyEqual(lsqvce.whitenessCorr(7, 5), 0);     % 1 element -> 0
end

%% ===================================================== collectResiduals
function testCollectResidualsFilters(tc)
    res = makeRes();
    wC = lsqvce.collectResiduals(res, 'C', 'pos', 'w');
    tc.verifyEqual(numel(wC), 18);                  % 2 arcs x 3 epochs x 3 axes (C/pos)
    vAll = lsqvce.collectResiduals(res, [], 'pos', 'v');
    tc.verifyEqual(numel(vAll), 36);                % C+D pooled (4 arcs x 9)
    empty = lsqvce.collectResiduals(res, 'C', 'nope');
    tc.verifyEmpty(empty);
    wDef = lsqvce.collectResiduals(res, 'D', 'pos'); % default field 'w'
    tc.verifyEqual(numel(wDef), 18);
end

%% ===================================================== residualAcv
function testDescaledAcvPosWhitens(tc)
% pos: v = L w (coloured only instantaneously by D=LL'), de-scale -> white ACV.
    rs = RandStream('mt19937ar','Seed',3);
    d = 3;  nE = 600;  maxLag = 40;  L = [1 0 0; 0.3 0.8 0; -0.2 0.1 0.6];
    [res, dayFcn] = makePosRes(rs, L, nE, 4);       % 4 arcs
    [acv, nSamp] = lsqvce.residualAcv(res, 'C', 'pos', d, maxLag, dayFcn);
    tc.verifyEqual(acv(1), 1, 'RelTol', 0.1);       % r0 ~ 1
    tc.verifyLessThan(abs(acv(2))/acv(1), 0.1);     % |r1/r0| small
    tc.verifyEqual(nSamp, 4*d*nE);
end

function testDescaledAcvRrWhitens(tc)
% rr: v = sigma .* w, de-scale by sigma -> white.
    rs = RandStream('mt19937ar','Seed',5);
    nE = 800;  maxLag = 30;
    [res, dayFcn] = makeRrRes(rs, nE, 3);
    [acv, nSamp] = lsqvce.residualAcv(res, 'C', 'rr', 1, maxLag, dayFcn);
    tc.verifyEqual(acv(1), 1, 'RelTol', 0.1);
    tc.verifyLessThan(abs(acv(2))/acv(1), 0.1);
    tc.verifyEqual(nSamp, 3*nE);
end

%% ===================================================== normalizeSatList
function testNormalizeSatList(tc)
    tc.verifyEqual(lsqvce.normalizeSatList([1 2]), {1, 2});
    tc.verifyEqual(lsqvce.normalizeSatList('CD'), {'C', 'D'});
    tc.verifyEqual(lsqvce.normalizeSatList({'C','D'}), {'C', 'D'});
    tc.verifyEqual(lsqvce.normalizeSatList("C"), {'C'});
    tc.verifyError(@() lsqvce.normalizeSatList([]), 'lsqvce:normalizeSatList:empty');
end

%% ===================================================== resolveSatSpecs
function testResolveSatSpecsOverrides(tc)
    defaults = struct('nl', 12, 'maxLag', 128);
    us(1).satId = 'D';  us(1).nl = 9;               % override only D.nl
    [S, specFor] = lsqvce.resolveSatSpecs({'C','D'}, us, defaults);

    tc.verifyEqual(specFor('C').nl, 12);            % default kept
    tc.verifyEqual(specFor('D').nl, 9);             % overridden
    tc.verifyEqual(specFor('C').maxLag, 128);       % default kept
    tc.verifyEqual(S(1).key, 's_C');                % via lsqvce.satKeyOf
    tc.verifyEqual(S(2).key, 's_D');
    tc.verifyError(@() specFor('Z'), 'lsqvce:resolveSatSpecs:noSpec');
end

function testResolveSatSpecsNumericKeyAndErrors(tc)
    [S, ~] = lsqvce.resolveSatSpecs({1, 2}, [], struct('nl', 5));
    tc.verifyEqual(S(1).key, 's_1');                % numeric id -> valid field
    tc.verifyEqual(S(1).nl, 5);
    bad.notSatId = 1;
    tc.verifyError(@() lsqvce.resolveSatSpecs({1}, bad, struct('nl',5)), ...
                   'lsqvce:resolveSatSpecs:badSatSpec');
end

%% ===================================================== satKeyOf (public)
function testSatKeyOf(tc)
    tc.verifyEqual(lsqvce.satKeyOf('C'), 's_C');
    tc.verifyEqual(lsqvce.satKeyOf("D"), 's_D');
    tc.verifyEqual(lsqvce.satKeyOf(1),   's_1');
    tc.verifyTrue(isvarname(lsqvce.satKeyOf(7)));   % always a legal field name
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
