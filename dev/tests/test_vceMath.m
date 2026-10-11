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
%
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-86 part 3: the leverage
%   points run on lsqvce.engine.BlockCofactor (leverage / leverageBlocks
%   retired). With assistance from Claude (Anthropic), model Claude Opus 5.5.
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

    Q = cofactorOf(Ncc_red, Nll, Nlc);
    hEst = [];
    for s = 1:nSys
        hEst = [hEst; Q.leverage(s, A{s}, B{s})];      %#ok<AGROW>
    end

    % tol: block-computed leverages vs the dense joint hat diagonal -- same
    % algebra through Schur blocks; eps*cond headroom (house bounds)
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

    Q = cofactorOf(Ncc_red, Nll, Nlc);
    h = [];
    for s = 1:nSys
        h = [h; Q.leverage(s, A{s}, B{s})];          %#ok<AGROW>
    end

    r = 1 - h;
    % tol: redundancy partition identity sum(1-h)=n-u; fp accumulation over
    % n leverages, 1e-7 headroom
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
        Q = cofactorOf(Ncc_red, {S.Nll}, {S.Nlc});
        for s = 1:nSys
            y   = S(s).Nll \ (S(s).nl - S(s).Nlc*x);
            for g = 1:nGrp
                w = 1/sig(g);
                Ai = Ad{s,g}*w;  Bi = Bd{s,g}*w;  li = ld{s,g}*w;
                h  = Q.leverage(s, Ai, Bi);
                rg(g)  = rg(g)  + numel(li) - sum(h);
                v      = li - Ai*x - Bi*y;
                vtv(g) = vtv(g) + v.'*v;
            end
        end
        sig = sig .* sqrt(vtv ./ max(rg,1));
    end

    M = nSys*nGrp*mG;  u = nc + nSys*nl;
    % tol: group redundancies summed over the Helmert loop; identity up to
    % fp accumulation across nSys*nGrp groups, 1e-5 headroom
    tc.verifyEqual(sum(rg), M - u, 'AbsTol', 1e-5);        % identity holds
    % tol: recovery of a planted value; estimator sampling scatter on this fixture, measured stable across seeds (not an eps claim)
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

    Q   = cofactorOf(Ncc_red, {Nll}, {Nlc});
    h   = Q.leverage(1, A, B);
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

    Q   = lsqvce.engine.BlockCofactor(chol((N + N.')/2), []);   % no locals
    h   = Q.leverage(0, A, B);

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
    pri = lsqvce.engine.CommonPrior(K, x0, sig);
    [Nk, nk] = pri.normals();

    Nll = B.'*B;  Nlc = B.'*A;
    Ncc_red = A.'*A - Nlc.'*(Nll\Nlc) + Nk;
    ncc_red = A.'*l - Nlc.'*(Nll\(B.'*l)) + nk;
    x = Ncc_red \ ncc_red;

    % tol: Schur-reduced solve vs the dense joint solve of the same system;
    % difference ~ cond*eps, 1e-8 house headroom
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
        pri = lsqvce.engine.CommonPrior(K, [], sigs(k));
        [Nk, ~] = pri.normals();
        Lc = chol(A.'*A + Nk, 'lower');
        [rK, absorbed(k)] = pri.redundancy(Lc);

        tc.verifyGreaterThanOrEqual(absorbed(k), -1e-8);
        tc.verifyLessThanOrEqual(absorbed(k), pri.nk + 1e-8);
        % tol: dense trace/quadratic-form reference of the same algebra; eps*cond headroom (house bounds)
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
        pri = lsqvce.engine.CommonPrior(Ks{k}, [], 1.0);
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

    % Singularity is asserted SCALE-FREE: whether chol's dpotrf reports
    % failure on a to-eps-singular matrix is BLAS/LAPACK-regime dependent
    % (the trailing pivot rounds to +/- tiny; the 2026-08-30 hosted-runner
    % image rounds it positive). The mathematics, not the factoriser's
    % boundary behaviour, is the contract (WP-47 r1e REGIME doctrine).
    % tol: 1e-12 -- 4 orders above eps(1) and >10 orders below the next
    % eigenvalue ratio of this seeded fixture (~0.1); negative min(lam)
    % (rounding) passes the same bound.
    lam = eig(A.'*A);
    tc.verifyLessThan(min(lam)/max(lam), 1e-12, ...
        'the common block should be singular here (eps-level eigenvalue ratio)');

    pri = lsqvce.engine.CommonPrior(eye(nc), [], 1.0);
    [Nk, ~] = pri.normals();
    [~, fail2] = chol(A.'*A + Nk);
    tc.verifyEqual(fail2, 0, 'the prior must make it solvable');
end

function testPriorRejectsBadTarget(tc)
% categories: contract
    K = eye(6);
    tc.verifyError(@() lsqvce.engine.CommonPrior(K, ones(3,1)), 'lsqvce:prior:badTarget');
    tc.verifyError(@() lsqvce.engine.CommonPrior(K, [], -1),    'MATLAB:lsqvce:engine:CommonPrior:expectedPositive');
end

%% ================================================== separable multi-axis
function testSeparableWhitenEqualsDenseCholesky(tc)
% categories: precision
% The separable whitener must reproduce the dense Cholesky of the full
% (nEpoch*d) x (nEpoch*d) covariance.
    rs = RandStream('mt19937ar','Seed',31);
    d = 3;  nE = 60;  m = nE*d;
    a = [1 -1.2 0.4];  s = 0.8;
    r = s^2 * lsqvce.kernel.armaACV(a, 1, nE);

    Rt  = toeplitz(r/r(1));                          % temporal correlation
    sig = 0.5 + rand(rs, nE, d);

    C = diag(reshape(sig.', m, 1)) * kron(Rt, eye(d))*r(1) * diag(reshape(sig.', m, 1));

    A = randn(rs, m, 5);  B = randn(rs, m, 2);  l = randn(rs, m, 1);

    W  = lsqvce.kernel.decorrFactor(C);
    Ad = W*A;  Bd = W*B;  ld = W*l;

    [Aw, Bw, lw] = lsqvce.kernel.separableWhiten(a, s, r, sig, d, A, B, l, true);

    % the factors may differ by a rotation; the NORMALS must not
    % tol: rotation-invariant quadratic forms of the same whitening;
    % eps*cond headroom (house bounds, HANDOVER 6.1.2)
    tc.verifyEqual(Aw.'*Aw, Ad.'*Ad, 'RelTol', 1e-8);
    tc.verifyEqual(Aw.'*lw, Ad.'*ld, 'RelTol', 1e-8);
    tc.verifyEqual(Bw.'*Bw, Bd.'*Bd, 'RelTol', 1e-8);     % tol: same source
end

function testSeparableWithFullEpochCovariance(tc)
% categories: precision
% The realistic case: a full 3x3 per-epoch covariance from kinematic PPP, where
% the vertical is worse than horizontal AND correlated with it.
    rs = RandStream('mt19937ar','Seed',37);
    d = 3;  nE = 40;  m = nE*d;
    a = [1 -0.9];  s = 1.0;
    r = s^2 * lsqvce.kernel.armaACV(a, 1, nE);
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
    W = lsqvce.kernel.decorrFactor(C);

    [Aw, ~, ~] = lsqvce.kernel.separableWhiten(a, s, r, Ce, d, A, B, l, true);
    % tol: separable whitener vs dense decorrelation of the same covariance;
    % quadratic forms at eps*cond, 1e-8 house headroom
    tc.verifyEqual(Aw.'*Aw, (W*A).'*(W*A), 'RelTol', 1e-8);
end

function testSeparableRejectsPartialEpoch(tc)
% categories: contract
    a = [1 -0.5];  r = lsqvce.kernel.armaACV(a, 1, 4);
    sig = ones(10, 3);
    A = zeros(29, 2);  B = zeros(29, 1);  l = zeros(29, 1);   % 29 is not 3*k
    tc.verifyError(@() lsqvce.kernel.separableWhiten(a, 1, r, sig, 3, A, B, l, true), ...
                   'lsqvce:separable:notWholeEpochs');
end

%% ================================================== cross-correlation
function testCrossCorrWhiteIsOK(tc)
% categories: contract | tags: diagnostics
    rs = RandStream('mt19937ar','Seed',41);
    c = lsqvce.kernel.crossCorr(randn(rs, 4000*3, 1), 3, 'MaxLag', 15);
    tc.verifyEqual(c.verdict, 'ok');
end

function testCrossCorrDetectsCrossAxis(tc)
% categories: contract | tags: diagnostics
% R(0) with off-diagonals -> the per-epoch covariance is wrong.
    rs = RandStream('mt19937ar','Seed',43);
    nE = 4000;  d = 3;
    M = [1 .6 .3; .6 1 .4; .3 .4 1];
    W = randn(rs, nE, d) * chol(M);
    c = tc.verifyWarning(@() lsqvce.kernel.crossCorr(reshape(W.', nE*d, 1), d), ...
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
    c = tc.verifyWarning(@() lsqvce.kernel.crossCorr(reshape(W.', nE*d, 1), d), ...
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
    c = tc.verifyWarning(@() lsqvce.kernel.crossCorr(reshape(W.', nE*d, 1), d), ...
                         'lsqvce:crossCorr:needs_var');
    tc.verifyEqual(c.verdict, 'needs-var');
end

%% ================================================== AR order selection
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

    % tol: independent dense reference, same-precision algebra; eps*cond headroom (house bounds, HANDOVER 6.1.2)
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

    % R7-42: two sequential timeit calls (drift lands on the later side, no
    % rotation) replaced by THE timing instrument; same 1.5x bound
    rec = timePair(@() blockSolve(Nll, [Nlc, nl_], nArc, nlA), @() Nll \ [Nlc, nl_], ...
                   'Label', 'block/dense local solve');
    tc.verifyLessThan(rec.ratio, 1/1.5, ...
        sprintf('block/dense = %.3f [%.3f %.3f] (want < 1/1.5)', rec.ratio, rec.band));
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

%% ===================================================================== helpers
function Q = cofactorOf(Ncc_red, Nll, Nlc)
%COFACTOROF  The block cofactor of a block-arrow system given its reduced
%   commons normal Ncc_red and, per system s, Nll{s} / Nlc{s}; the local
%   columns of system s follow the commons in system order (WP-86 part 3:
%   lsqvce.engine.BlockCofactor replaces leverageBlocks + leverage).
    nc = size(Ncc_red, 1);
    blocks = struct('cols', {}, 'Nll', {}, 'Nlc', {});
    col = nc;
    for s = 1:numel(Nll)
        nl = size(Nll{s}, 1);
        blocks(s) = struct('cols', col + (1:nl), 'Nll', Nll{s}, 'Nlc', Nlc{s});
        col = col + nl;
    end
    Q = lsqvce.engine.BlockCofactor(chol((Ncc_red + Ncc_red.') / 2), blocks);
end
