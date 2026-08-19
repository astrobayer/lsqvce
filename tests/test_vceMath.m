function tests = test_vceMath
%TEST_VCEMATH  Variance components, regularisation, multi-axis, order selection.
%
%   Covers everything added in v3:
%     * exact per-observation leverage (the core of VCE and outlier detection)
%     * Helmert VCE per arc
%     * data snooping
%     * basis-agnostic regularisation of the COMMON parameters
%     * separable multi-axis whitening
%     * cross-correlation diagnostic
%     * automatic AR order selection
%     * block-diagonal Schur for arc-wise local parameters
%
%   Run:  results = runtests('tests/test_vceMath');  disp(table(results))
%
%   See also TEST_PRECISION, TEST_RELIABILITY, TEST_SPEED
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    tc.TestData.cfg = makeCfg();
end

%% ================================================== leverage (the core)
function testLeverageEqualsBruteForceHatMatrix(tc)
% categories: precision, reference
% h_ii from the block-arrow inverse must equal the diagonal of the explicit
% joint hat matrix. This is the foundation of VCE and outlier detection: if it
% is wrong, both give plausible but wrong answers.
    rs = RandStream('mt19937ar','Seed',3);
    nc = 8;  nl = 4;  nSys = 3;

    A = cell(1,nSys);  B = cell(1,nSys);
    for s = 1:nSys
        m = 200 + 30*s;
        A{s} = randn(rs, m, nc);
        B{s} = randn(rs, m, nl);
    end

    % --- explicit joint hat matrix (the reference) -------------------------
    ntot = nc + nSys*nl;
    Aw = [];
    for s = 1:nSys
        m   = size(A{s},1);
        blk = zeros(m, ntot);
        blk(:, 1:nc) = A{s};
        blk(:, nc + (s-1)*nl + (1:nl)) = B{s};
        Aw = [Aw; blk];        %#ok<AGROW>
    end
    N     = Aw.'*Aw;
    hTrue = diag(Aw * (N \ Aw.'));

    % --- via the toolbox --------------------------------------------------
    Ncc_red = zeros(nc);
    Nll = cell(1,nSys);  Nlc = cell(1,nSys);
    for s = 1:nSys
        Nll{s}  = B{s}.'*B{s};
        Nlc{s}  = B{s}.'*A{s};
        Ncc_red = Ncc_red + A{s}.'*A{s} - Nlc{s}.'*(Nll{s}\Nlc{s});
    end

    hEst = [];
    for s = 1:nSys
        blk  = lsqvce.leverageBlocks(Ncc_red, Nll{s}, Nlc{s});
        hEst = [hEst; lsqvce.leverage(A{s}, B{s}, blk)];      %#ok<AGROW>
    end

    tc.verifyEqual(hEst, hTrue, 'AbsTol', 1e-9);
    tc.verifyEqual(sum(hEst), ntot, 'RelTol', 1e-10);      % sum(h) = u
end

function testRedundancyIdentity(tc)
% categories: precision
% sum(r_i) = m - u exactly. The strongest global check available.
    rs = RandStream('mt19937ar','Seed',5);
    nc = 6;  nl = 3;  nSys = 4;  m = 150;

    Ncc_red = zeros(nc);
    A = cell(1,nSys); B = cell(1,nSys); Nll = cell(1,nSys); Nlc = cell(1,nSys);
    for s = 1:nSys
        A{s} = randn(rs, m, nc);
        B{s} = randn(rs, m, nl);
        Nll{s} = B{s}.'*B{s};
        Nlc{s} = B{s}.'*A{s};
        Ncc_red = Ncc_red + A{s}.'*A{s} - Nlc{s}.'*(Nll{s}\Nlc{s});
    end

    h = [];
    for s = 1:nSys
        blk = lsqvce.leverageBlocks(Ncc_red, Nll{s}, Nlc{s});
        h   = [h; lsqvce.leverage(A{s}, B{s}, blk)];          %#ok<AGROW>
    end

    r = 1 - h;
    tc.verifyEqual(sum(r), nSys*m - (nc + nSys*nl), 'AbsTol', 1e-7);
    tc.verifyGreaterThanOrEqual(min(h), -1e-10);
    tc.verifyLessThanOrEqual(max(h),  1 + 1e-10);
end

%% ================================================== VCE
function testVCERecoversVarianceComponents(tc)
% categories: precision
% The end-to-end claim: Helmert VCE recovers the true per-group sigma from a
% deliberately wrong start.
    rs = RandStream('mt19937ar','Seed',77);
    nc = 25;  nl = 4;  nSys = 2;  nGrp = 4;  mG = 350;
    TRUE = [1.0 2.5 0.6 1.5];

    xt = randn(rs, nc, 1);
    yt = arrayfun(@(~) randn(rs, nl, 1), 1:nSys, 'UniformOutput', false);

    Ad = cell(nSys,nGrp); Bd = cell(nSys,nGrp); ld = cell(nSys,nGrp);
    for s = 1:nSys
        for g = 1:nGrp
            Ad{s,g} = randn(rs, mG, nc);
            Bd{s,g} = randn(rs, mG, nl);
            ld{s,g} = Ad{s,g}*xt + Bd{s,g}*yt{s} + TRUE(g)*randn(rs, mG, 1);
        end
    end

    sig = ones(1, nGrp);
    for it = 1:6
        % --- accumulate with the current sigma ----------------------------
        Ncc_red = zeros(nc);  ncc_red = zeros(nc,1);
        S = struct('Nll',{},'Nlc',{},'nl',{});
        for s = 1:nSys
            Ncc = zeros(nc); Nlc = zeros(nl,nc); Nll = zeros(nl);
            nc_ = zeros(nc,1); nl_ = zeros(nl,1);
            for g = 1:nGrp
                w = 1/sig(g);
                Ai = Ad{s,g}*w;  Bi = Bd{s,g}*w;  li = ld{s,g}*w;
                Ncc = Ncc + Ai.'*Ai;   Nlc = Nlc + Bi.'*Ai;   Nll = Nll + Bi.'*Bi;
                nc_ = nc_ + Ai.'*li;   nl_ = nl_ + Bi.'*li;
            end
            Ncc_red = Ncc_red + Ncc - Nlc.'*(Nll\Nlc);
            ncc_red = ncc_red + nc_ - Nlc.'*(Nll\nl_);
            S(s) = struct('Nll',Nll, 'Nlc',Nlc, 'nl',nl_);
        end
        x = Ncc_red \ ncc_red;

        % --- redundancy and vtv per group ---------------------------------
        rg = zeros(1,nGrp);  vtv = zeros(1,nGrp);
        for s = 1:nSys
            blk = lsqvce.leverageBlocks(Ncc_red, S(s).Nll, S(s).Nlc);
            y   = S(s).Nll \ (S(s).nl - S(s).Nlc*x);
            for g = 1:nGrp
                w = 1/sig(g);
                Ai = Ad{s,g}*w;  Bi = Bd{s,g}*w;  li = ld{s,g}*w;
                h  = lsqvce.leverage(Ai, Bi, blk);
                rg(g)  = rg(g)  + numel(li) - sum(h);
                v      = li - Ai*x - Bi*y;
                vtv(g) = vtv(g) + v.'*v;
            end
        end
        sig = sig .* sqrt(vtv ./ max(rg,1));
    end

    M = nSys*nGrp*mG;  u = nc + nSys*nl;
    tc.verifyEqual(sum(rg), M - u, 'AbsTol', 1e-5);        % identity holds
    tc.verifyEqual(sig ./ TRUE, ones(1,nGrp), 'AbsTol', 0.15);
end

function testOutlierDetection(tc)
% categories: contract
% Data snooping with the correct sqrt(r_i) normalisation.
    rs = RandStream('mt19937ar','Seed',9);
    nc = 18;  nl = 4;  m = 700;

    A = randn(rs, m, nc);  B = randn(rs, m, nl);
    l = A*randn(rs,nc,1) + B*randn(rs,nl,1) + randn(rs, m, 1);
    bad = [23 111 400 650];
    l(bad) = l(bad) + 11;

    Nll = B.'*B;  Nlc = B.'*A;
    Ncc_red = A.'*A - Nlc.'*(Nll\Nlc);
    ncc_red = A.'*l - Nlc.'*(Nll\(B.'*l));
    x = Ncc_red \ ncc_red;
    y = Nll \ (B.'*l - Nlc*x);

    blk = lsqvce.leverageBlocks(Ncc_red, Nll, Nlc);
    h   = lsqvce.leverage(A, B, blk);
    v   = l - A*x - B*y;
    w   = abs(v) ./ sqrt(max(1-h, eps));

    flagged = find(w > 4.0);
    tc.verifyTrue(all(ismember(bad, flagged.')), 'missed an injected outlier');
    tc.verifyLessThanOrEqual(numel(setdiff(flagged.', bad)), 2, 'too many false alarms');
end

function testRawResidualWouldMissHighLeverageOutlier(tc)
% categories: precision, robustness
% WHY the sqrt(r_i) normalisation matters: an observation with high leverage is
% fitted almost exactly by construction, so its RAW residual is small even when
% it is badly wrong. This test shows the normalisation catching what the raw
% residual misses.
    rs = RandStream('mt19937ar','Seed',13);
    nc = 5;  m = 40;

    A = randn(rs, m, nc);
    A(1,:) = 8*A(1,:);                      % row 1 has very high leverage
    B = zeros(m, 0);
    l = A*randn(rs, nc, 1) + 0.1*randn(rs, m, 1);
    l(1) = l(1) + 3;                        % a real blunder on the leverage point

    N = A.'*A;
    x = N \ (A.'*l);
    v = l - A*x;

    blk = lsqvce.leverageBlocks(N, eye(0), zeros(0,nc));
    h   = lsqvce.leverage(A, B, blk);

    tc.verifyGreaterThan(h(1), 0.5, 'row 1 should have high leverage');

    % the raw residual is SMALL despite the blunder ...
    tc.verifyLessThan(abs(v(1)), max(abs(v)), ...
        'the blunder should be masked in the raw residual');

    % ... but the normalised residual exposes it
    w = abs(v) ./ sqrt(max(1-h, eps));
    [~, k] = max(w);
    tc.verifyEqual(k, 1, 'the normalised residual must rank the blunder first');
end

%% ================================================== common-part prior
function testPriorEqualsAugmentedLeastSquares(tc)
% categories: precision
% The prior must be algebraically exact: identical to stacking K*x ~ x0 as extra
% observation rows.
    rs = RandStream('mt19937ar','Seed',17);
    nc = 10;  nl = 4;  m = 300;

    A = randn(rs, m, nc);  B = randn(rs, m, nl);  l = randn(rs, m, 1);
    K   = diag(1./(1:nc).^1.5);
    x0  = 0.1*randn(rs, nc, 1);
    sig = 0.7;

    % --- explicit augmentation (the reference) ----------------------------
    Aa = [ A,      B ;
           K/sig,  zeros(nc, nl) ];
    la = [ l ; x0/sig ];
    xf = Aa \ la;

    % --- via the prior ----------------------------------------------------
    pri = lsqvce.CommonPrior(K, x0, sig);
    [Nk, nk] = pri.normals();

    Nll = B.'*B;  Nlc = B.'*A;
    Ncc_red = A.'*A - Nlc.'*(Nll\Nlc) + Nk;
    ncc_red = A.'*l - Nlc.'*(Nll\(B.'*l)) + nk;
    x = Ncc_red \ ncc_red;

    tc.verifyEqual(x, xf(1:nc), 'RelTol', 1e-8);
end

function testPriorAbsorbedDoFBounds(tc)
% categories: precision
% absorbed = nk - r_K must lie in [0, nk], and must respond correctly to sigma:
% a strong prior absorbs nearly everything, a weak one almost nothing.
    rs = RandStream('mt19937ar','Seed',19);
    nc = 12;  m = 200;
    A = randn(rs, m, nc);
    K = eye(nc);

    absorbed = zeros(1,2);
    sigs = [0.01 100];
    for k = 1:2
        pri = lsqvce.CommonPrior(K, [], sigs(k));
        [Nk, ~] = pri.normals();
        Lc = chol(A.'*A + Nk, 'lower');
        [rK, absorbed(k)] = pri.redundancy(Lc);

        tc.verifyGreaterThanOrEqual(absorbed(k), -1e-8);
        tc.verifyLessThanOrEqual(absorbed(k), pri.nk + 1e-8);
        tc.verifyEqual(rK, pri.nk - absorbed(k), 'AbsTol', 1e-10);
    end

    tc.verifyGreaterThan(absorbed(1), 0.9*nc, 'a strong prior must absorb nearly all DoF');
    tc.verifyLessThan(   absorbed(2), 0.1*nc, 'a weak prior must absorb almost none');
end

function testPriorIsBasisAgnostic(tc)
% categories: precision, robustness
% K is just a matrix. Identity, Kaula, a mascon Laplacian, a subset selector --
% all must work without the toolbox knowing what basis it is.
    rs = RandStream('mt19937ar','Seed',23);
    nc = 16;  m = 120;
    A = randn(rs, m, nc);

    L = zeros(nc-1, nc);                            % mascon smoothness
    for i = 1:nc-1
        L(i,i) = -1;  L(i,i+1) = 1;
    end
    Ksub = [zeros(4, nc-4), eye(4)];                % subset only

    Ks = {eye(nc), diag(1./(1:nc).^2), L, Ksub};
    nm = {'identity','kaula','laplacian','subset'};

    for k = 1:numel(Ks)
        pri = lsqvce.CommonPrior(Ks{k}, [], 1.0);
        [Nk, ~] = pri.normals();
        Lc = chol(A.'*A + Nk, 'lower');
        [rK, absorbed] = pri.redundancy(Lc);

        tc.verifyTrue(isfinite(rK) && isfinite(absorbed), nm{k});
        tc.verifyGreaterThanOrEqual(absorbed, -1e-8, nm{k});
        tc.verifyLessThanOrEqual(absorbed, pri.nk + 1e-8, nm{k});
    end
end

function testPriorRescuesRankDeficientCommon(tc)
% categories: robustness
% A datum defect in the common block: without a prior the normal matrix is
% singular; with one it is solvable.
    rs = RandStream('mt19937ar','Seed',29);
    nc = 8;  m = 60;
    A = randn(rs, m, nc);
    A(:,end) = A(:,end-1);                          % exactly collinear

    [~, fail] = chol(A.'*A);
    tc.verifyGreaterThan(fail, 0, 'the common block should be singular here');

    pri = lsqvce.CommonPrior(eye(nc), [], 1.0);
    [Nk, ~] = pri.normals();
    [~, fail2] = chol(A.'*A + Nk);
    tc.verifyEqual(fail2, 0, 'the prior must make it solvable');
end

function testPriorRejectsBadTarget(tc)
% categories: contract
    K = eye(6);
    tc.verifyError(@() lsqvce.CommonPrior(K, ones(3,1)), 'lsqvce:prior:badTarget');
    tc.verifyError(@() lsqvce.CommonPrior(K, [], -1),    'MATLAB:lsqvce:CommonPrior:expectedPositive');
end

%% ================================================== separable multi-axis
function testSeparableWhitenEqualsDenseCholesky(tc)
% categories: precision
% The separable whitener must reproduce the dense Cholesky of the full
% (nEpoch*d) x (nEpoch*d) covariance.
    rs = RandStream('mt19937ar','Seed',31);
    d = 3;  nE = 60;  m = nE*d;
    a = [1 -1.2 0.4];  s = 0.8;
    r = s^2 * lsqvce.armaACV(a, 1, nE);

    Rt  = toeplitz(r/r(1));                          % temporal correlation
    sig = 0.5 + rand(rs, nE, d);

    C = diag(reshape(sig.', m, 1)) * kron(Rt, eye(d))*r(1) * diag(reshape(sig.', m, 1));

    A = randn(rs, m, 5);  B = randn(rs, m, 2);  l = randn(rs, m, 1);

    W  = lsqvce.decorrFactor(C);
    Ad = W*A;  Bd = W*B;  ld = W*l;

    [Aw, Bw, lw] = lsqvce.separableWhiten(a, s, r, sig, d, A, B, l, true);

    % the factors may differ by a rotation; the NORMALS must not
    tc.verifyEqual(Aw.'*Aw, Ad.'*Ad, 'RelTol', 1e-8);
    tc.verifyEqual(Aw.'*lw, Ad.'*ld, 'RelTol', 1e-8);
    tc.verifyEqual(Bw.'*Bw, Bd.'*Bd, 'RelTol', 1e-8);
end

function testSeparableWithFullEpochCovariance(tc)
% categories: precision
% The realistic case: a full 3x3 per-epoch covariance from kinematic PPP, where
% the vertical is worse than horizontal AND correlated with it.
    rs = RandStream('mt19937ar','Seed',37);
    d = 3;  nE = 40;  m = nE*d;
    a = [1 -0.9];  s = 1.0;
    r = s^2 * lsqvce.armaACV(a, 1, nE);
    Rt = toeplitz(r/r(1));

    Ce = zeros(nE, d, d);
    Db = zeros(m, m);
    for k = 1:nE
        Q = 0.2*randn(rs,d,d) + diag([1 1 2.5]);
        Ck = Q*Q.';
        Ce(k,:,:) = Ck;
        ib = (k-1)*d + (1:d);
        Db(ib, ib) = chol(Ck, 'lower');
    end
    C = Db * (kron(Rt, eye(d))*r(1)) * Db.';

    A = randn(rs, m, 4);  B = randn(rs, m, 2);  l = randn(rs, m, 1);
    W = lsqvce.decorrFactor(C);

    [Aw, ~, ~] = lsqvce.separableWhiten(a, s, r, Ce, d, A, B, l, true);
    tc.verifyEqual(Aw.'*Aw, (W*A).'*(W*A), 'RelTol', 1e-8);
end

function testSeparableRejectsPartialEpoch(tc)
% categories: contract
    a = [1 -0.5];  r = lsqvce.armaACV(a, 1, 4);
    sig = ones(10, 3);
    A = zeros(29, 2);  B = zeros(29, 1);  l = zeros(29, 1);   % 29 is not 3*k
    tc.verifyError(@() lsqvce.separableWhiten(a, 1, r, sig, 3, A, B, l, true), ...
                   'lsqvce:separable:notWholeEpochs');
end

%% ================================================== cross-correlation
function testCrossCorrWhiteIsOK(tc)
% categories: contract | tags: diagnostics
    rs = RandStream('mt19937ar','Seed',41);
    c = lsqvce.crossCorr(randn(rs, 4000*3, 1), 3, 'MaxLag', 15);
    tc.verifyEqual(c.verdict, 'ok');
end

function testCrossCorrDetectsCrossAxis(tc)
% categories: contract | tags: diagnostics
% R(0) with off-diagonals -> the per-epoch covariance is wrong.
    rs = RandStream('mt19937ar','Seed',43);
    nE = 4000;  d = 3;
    M = [1 .6 .3; .6 1 .4; .3 .4 1];
    W = randn(rs, nE, d) * chol(M);
    c = tc.verifyWarning(@() lsqvce.crossCorr(reshape(W.', nE*d, 1), d), ...
                         'lsqvce:crossCorr:cross_axis_zero_lag');
    tc.verifyEqual(c.verdict, 'cross-axis-zero-lag');
end

function testCrossCorrDetectsARTooLow(tc)
% categories: contract | tags: diagnostics
% R(tau>0) with a nonzero DIAGONAL -> leftover temporal colour.
    rs = RandStream('mt19937ar','Seed',47);
    nE = 4000;  d = 3;
    W = filter(1, [1 -0.6], randn(rs, nE, d), [], 1);
    W = W ./ std(W, 0, 1);
    c = tc.verifyWarning(@() lsqvce.crossCorr(reshape(W.', nE*d, 1), d), ...
                         'lsqvce:crossCorr:ar_order_low');
    tc.verifyEqual(c.verdict, 'ar-order-low');
end

function testCrossCorrDetectsNeedsVAR(tc)
% categories: contract | tags: diagnostics
% R(tau>0) with OFF-diagonals -> genuine cross-axis TEMPORAL correlation.
% Separability fails; a VAR is required. This is the case the toolbox cannot
% handle, and it must say so rather than silently mis-weighting.
    rs = RandStream('mt19937ar','Seed',53);
    nE = 6000;  d = 3;
    Z = randn(rs, nE, d);
    W = Z;
    W(2:end,1) = W(2:end,1) + 0.5*Z(1:end-1,2);      % axis1(k) <- axis2(k-1)
    W = W ./ std(W, 0, 1);
    c = tc.verifyWarning(@() lsqvce.crossCorr(reshape(W.', nE*d, 1), d), ...
                         'lsqvce:crossCorr:needs_var');
    tc.verifyEqual(c.verdict, 'needs-var');
end

%% ================================================== AR order selection
function testSelectAROrderFindsSufficientOrder(tc)
% categories: contract | tags: diagnostics
    dt = 5;
    a  = [1 -1.6 0.9 -0.25];
    f  = linspace(1e-9, 0.5/dt, 4096).';
    H  = 1 ./ (exp(-1i*2*pi*f*dt*(0:3)) * a(:));
    cm = struct('type','psd', 'f',f, 'S', 0.7^2*dt*2*abs(H).^2);

    info = lsqvce.selectArOrder(cm, dt, 'AcvTol', 1e-3, 'Verbose', false);
    tc.verifyEqual(info.verdict, 'good');
    tc.verifyGreaterThanOrEqual(info.p, 3);
    tc.verifyLessThan(info.acvRelErr, 1e-3);
end

function testSelectAROrderReportsCurve(tc)
% categories: contract
    dt = 5;
    a  = [1 -1.6 0.9 -0.25];
    f  = linspace(1e-9, 0.5/dt, 2048).';
    H  = 1 ./ (exp(-1i*2*pi*f*dt*(0:3)) * a(:));
    cm = struct('type','psd', 'f',f, 'S', dt*2*abs(H).^2);

    info = lsqvce.selectArOrder(cm, dt, 'OrderGrid', [2 5 10 20 40], 'Verbose', false);
    tc.verifyGreaterThanOrEqual(numel(info.curve), 2);
    e = [info.curve.acvRelErr];
    tc.verifyLessThan(e(end), e(1), 'the ACV error must fall with order');
    tc.verifyTrue(all(isfield(info.curve, {'aic','bic','maxRoot'})));
end

function testSelectAROrderFlagsHopelessSpectrum(tc)
% categories: contract | tags: diagnostics
% A razor-sharp spectral line: AR cannot represent it at low order and must say
% so rather than returning a silently bad filter.
    dt = 1;
    f  = linspace(1e-6, 0.5, 8192).';
    S  = 1e-4*ones(size(f)) + 5*exp(-((f-0.1)/2e-4).^2);
    cm = struct('type','psd', 'f',f, 'S',S);

    info = lsqvce.selectArOrder(cm, dt, 'OrderGrid', [2 5 10 20], ...
                             'AcvTol', 1e-3, 'Verbose', false);
    tc.verifyTrue(any(strcmp(info.verdict, {'marginal','bad'})));
end

%% ================================================== block-diagonal Schur
function testBlockSchurEqualsDenseSchur(tc)
% categories: precision
% Arc-wise local parameters make Nll block diagonal. The block path must give
% exactly the same reduced normals as the dense path.
    rs = RandStream('mt19937ar','Seed',59);
    nc = 8;  nArc = 20;  nlA = 6;
    nl = nArc*nlA;

    Ncc = zeros(nc);  Nlc = zeros(nl,nc);  Nll = zeros(nl);
    nc_ = zeros(nc,1);  nl_ = zeros(nl,1);
    for a = 1:nArc
        mA = 60;
        Aa = randn(rs, mA, nc);
        Ba = randn(rs, mA, nlA);
        la = randn(rs, mA, 1);
        ib = (a-1)*nlA + (1:nlA);

        Ncc = Ncc + Aa.'*Aa;   nc_ = nc_ + Aa.'*la;
        Nll(ib,ib) = Ba.'*Ba;
        Nlc(ib,:)  = Ba.'*Aa;
        nl_(ib)    = Ba.'*la;
    end

    % dense
    Zd = Nll \ [Nlc, nl_];
    Ncc_dense = Ncc - Nlc.'*Zd(:,1:nc);
    ncc_dense = nc_ - Nlc.'*Zd(:,end);

    % block
    Zb = zeros(nl, nc+1);
    for a = 1:nArc
        ib = (a-1)*nlA + (1:nlA);
        Zb(ib,:) = Nll(ib,ib) \ [Nlc(ib,:), nl_(ib)];
    end
    Ncc_blk = Ncc - Nlc.'*Zb(:,1:nc);
    ncc_blk = nc_ - Nlc.'*Zb(:,end);

    tc.verifyEqual(Ncc_blk, Ncc_dense, 'RelTol', 1e-10);
    tc.verifyEqual(ncc_blk, ncc_dense, 'RelTol', 1e-10);
end

function testBlockSchurIsFaster(tc)
% categories: speed
% O(nArc*blk^3) versus O((nArc*blk)^3): the block Schur elimination avoids the
% dense factorisation of the full local matrix. The asymptotic win only shows up
% once the dense cube dominates LAPACK's constant factor, so this uses a LARGE
% local dimension; on tiny problems MATLAB loop overhead hides the advantage and
% the comparison is meaningless. A 1.5x margin keeps it from flaking on a busy
% machine while still catching a real regression (block should be many times
% faster here, not merely 1.5x).
    rs = RandStream('mt19937ar','Seed',61);
    nc = 12;  nArc = 200;  nlA = 8;
    nl = nArc*nlA;                                   % 1600 local params

    Nll = zeros(nl);  Nlc = randn(rs, nl, nc);  nl_ = randn(rs, nl, 1);
    for a = 1:nArc
        ib = (a-1)*nlA + (1:nlA);
        Q = randn(rs, nlA, nlA);
        Nll(ib,ib) = Q*Q.' + nlA*eye(nlA);
    end

    tDense = timeit(@() Nll \ [Nlc, nl_]);
    tBlock = timeit(@() blockSolve(Nll, [Nlc, nl_], nArc, nlA));

    tc.verifyLessThan(tBlock, tDense / 1.5, ...
        sprintf('block %.2f ms vs dense %.2f ms (want block < dense/1.5)', ...
                tBlock*1e3, tDense*1e3));
end

%% ================================================== helpers
function Z = blockSolve(Nll, RHS, nArc, nlA)
    Z = zeros(size(RHS));
    for a = 1:nArc
        ib = (a-1)*nlA + (1:nlA);
        Z(ib,:) = Nll(ib,ib) \ RHS(ib,:);
    end
end

function cfg = makeCfg()
    rs = RandStream('mt19937ar','Seed',42);
    cfg.nc   = 6;
    cfg.nl   = 4;
    cfg.sats = [1 2 3];
    cfg.x_c  = randn(rs, cfg.nc, 1);
    cfg.y    = containers.Map('KeyType','double','ValueType','any');
    for s = cfg.sats, cfg.y(s) = randn(rs, cfg.nl, 1); end
    cfg.obsTypes     = {'pos','kbr'};
    cfg.arcsPerDay   = 2;
    cfg.mArc         = 400;
    cfg.dt           = 5;
    cfg.seed         = 42;
    cfg.collinear    = false;
    cfg.collinearEps = 1e-7;
    cfg.cov.pos = struct('type','arma','ar',[1 -1.6 0.9 -0.25],'ma',1,'sigma',0.7);
    cfg.cov.kbr = struct('type','diag','sigma',0.2);
end
