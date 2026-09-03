function tests = test_metamorphic
%TEST_METAMORPHIC  Invariances the estimators must obey (R5 M7, the seventh category).
%   tests = test_metamorphic
%
%   A metamorphic point changes the INPUT in a way whose effect on the
%   OUTPUT is known exactly, and asserts that relation -- no reference
%   solution needed. The relations here are the ones a geodesist checks by
%   hand: scaling the observations scales the estimate and sigma0, scaling
%   the weights leaves the estimate alone, a batch partition is invisible,
%   a whitener is linear, a spectral round trip is the identity, robust
%   factors depend on the residuals' scale only through the MAD.
%
%   See also TEST_INSTRUMENTS, TEST_FIT, TEST_ELLMER
%
%   Created 2026-09-03 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M7, WP-89.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    here = fileparts(mfilename('fullpath'));
    root = fileparts(fileparts(here));
    addpath(root);  addpath(fullfile(root, 'examples'));
    rs = RandStream('mt19937ar', 'Seed', 71);
    A = [ones(60, 1), randn(rs, 60, 4)];  x = [2; -1; 0.5; 3; -2];
    p = 0.5 + rand(rs, 60, 1);
    l = A*x + randn(rs, 60, 1) ./ sqrt(p);
    tc.TestData.A = A;  tc.TestData.l = l;  tc.TestData.p = p;
end

function testScalingTheObservationsScalesTheEstimate(tc)
% categories: metamorphic
    % l -> c*l  =>  x -> c*x, sigma0 -> |c|*sigma0, Cxx -> c^2 Cxx (the block
    % adapter lsqvce.engine.batchSystem, lsqvce.engine.adjustBlock, the
    % diagnostics lsqvce.engine.diagnose / lsqvce.engine.diagnoseAccum)
    A = tc.TestData.A;  l = tc.TestData.l;  p = tc.TestData.p;  c = -3.7;
    R1 = lsqvce.fit(struct('A', A, 'l', l,   'p', p), Verbose=false);
    R2 = lsqvce.fit(struct('A', A, 'l', c*l, 'p', p), Verbose=false);
    % tol: 1e-12 relative -- the same arithmetic on scaled numbers (no
    % cancellation is introduced by a scalar factor)
    tc.verifyEqual(R2.x, c*R1.x, 'RelTol', 1e-12);
    % tol: as above
    tc.verifyEqual(R2.sigma0, abs(c)*R1.sigma0, 'RelTol', 1e-12);
    % tol: as above
    tc.verifyEqual(R2.sol.Cxx, c^2*R1.sol.Cxx, 'RelTol', 1e-12);
    % tol: as above
    tc.verifyEqual(R2.diag.rNum, R1.diag.rNum, 'RelTol', 1e-12, 'redundancy numbers are scale-free');
end

function testScalingTheWeightsLeavesTheEstimate(tc)
% categories: metamorphic
    % p -> k*p leaves x and the redundancies unchanged; sigma0 scales by
    % sqrt(k) (the a-posteriori factor absorbs the common weight)
    A = tc.TestData.A;  l = tc.TestData.l;  p = tc.TestData.p;  k = 25;
    R1 = lsqvce.fit(struct('A', A, 'l', l, 'p', p),   Verbose=false);
    R2 = lsqvce.fit(struct('A', A, 'l', l, 'p', k*p), Verbose=false);
    % tol: 1e-12 relative -- as above
    tc.verifyEqual(R2.x, R1.x, 'RelTol', 1e-12);
    % tol: as above
    tc.verifyEqual(R2.sigma0, sqrt(k)*R1.sigma0, 'RelTol', 1e-12);
    % tol: as above
    tc.verifyEqual(R2.sol.Cxx, R1.sol.Cxx, 'RelTol', 1e-11, 'Cxx = s0sq * Qxx is weight-scale free');
    % tol: as above
    tc.verifyEqual(R2.diag.rNum, R1.diag.rNum, 'RelTol', 1e-12);
end

function testBatchPartitionIsInvisible(tc)
% categories: metamorphic
    % the same rows in 1, 3 or 7 batches (lsqvce.engine.batchSystem: batch k =
    % day k; lsqvce.engine.BlockAccumulator merges days; lsqvce.engine.QRSystem
    % folds) give one solution
    A = tc.TestData.A;  l = tc.TestData.l;  p = tc.TestData.p;  n = numel(l);
    R1 = lsqvce.fit(struct('A', A, 'l', l, 'p', p), Verbose=false);
    for nb = [3 7]
        edges = round(linspace(0, n, nb + 1));
        bf = @(k) deal(A(edges(k)+1:edges(k+1), :), l(edges(k)+1:edges(k+1)), p(edges(k)+1:edges(k+1)));
        Rk = lsqvce.fit(struct('batchFcn', bf, 'numUnknowns', 5, 'numBatches', nb), Verbose=false);
        % tol: 1e-11 relative -- QR folds in a different order; cond(A) ~ 10
        tc.verifyEqual(Rk.x, R1.x, 'RelTol', 1e-11, sprintf('%d batches', nb));
        % tol: as above
        tc.verifyEqual(Rk.sigma0, R1.sigma0, 'RelTol', 1e-11);
    end
end

function testWhitenerIsLinearAndScales(tc)
% categories: metamorphic
    % lsqvce.kernel.armaWhiten is linear in its rows and scales inversely with
    % the innovation std; lsqvce.kernel.armaACV scales with sigma^2
    a = [1 -0.6 0.2];  b = 1;  m = 200;
    rs = RandStream('mt19937ar', 'Seed', 3);
    X = randn(rs, m, 3);  Y = randn(rs, m, 3);  l = randn(rs, m, 1);
    r1 = lsqvce.kernel.armaACV(a, b, 32);
    [WX, ~, Wl]   = lsqvce.kernel.armaWhiten(a, b, 1,   r1, X,     zeros(m, 0), l, true);
    [WY, ~, ~]    = lsqvce.kernel.armaWhiten(a, b, 1,   r1, Y,     zeros(m, 0), l, true);
    [WXY, ~, ~]   = lsqvce.kernel.armaWhiten(a, b, 1,   r1, 2*X+Y, zeros(m, 0), l, true);
    [WX3, ~, Wl3] = lsqvce.kernel.armaWhiten(a, b, 3,   r1, X,     zeros(m, 0), l, true);
    % tol: 1e-12 relative -- linear filters on O(1) numbers
    tc.verifyEqual(WXY, 2*WX + WY, 'RelTol', 1e-12, 'linearity');
    % tol: as above
    tc.verifyEqual(WX3, WX/3, 'RelTol', 1e-12, 'the innovation std divides the rows');
    % tol: as above
    tc.verifyEqual(Wl3, Wl/3, 'RelTol', 1e-12);
    r3 = 9 * lsqvce.kernel.armaACV(a, b, 32);
    % tol: as above
    tc.verifyEqual(r3, 9*r1, 'RelTol', 1e-12, 'the ACV scales with sigma^2');
end

function testSpectralRoundTripIsTheIdentity(tc)
% categories: metamorphic
    % lsqvce.kernel.dctI forward then inverse is the identity (orthonormal),
    % and lsqvce.kernel.psd2acv of a scaled PSD is the scaled ACV
    rs = RandStream('mt19937ar', 'Seed', 5);
    c = exp(-(0:63).'/10) .* (1 + 0.05*randn(rs, 64, 1));
    S = lsqvce.kernel.dctI(c, 'Norm', 'ortho');
    c2 = lsqvce.kernel.dctI(S, 'Norm', 'ortho');
    % tol: 1e-12 absolute -- an orthonormal transform applied twice
    tc.verifyEqual(c2, c, 'AbsTol', 1e-12);
    f = linspace(0, 0.1, 40).';  Sp = 1 + 0.5*cos(2*pi*f/0.1);
    r1 = lsqvce.kernel.psd2acv(f, Sp, 5, 10);
    r4 = lsqvce.kernel.psd2acv(f, 4*Sp, 5, 10);
    % tol: as above
    tc.verifyEqual(r4, 4*r1, 'RelTol', 1e-12, 'psd2acv is linear in the PSD');
end

function testRobustFactorsAreScaleFree(tc)
% categories: metamorphic
    % lsqvce.kernel.robustWeights: residuals scaled by any c > 0 give the SAME
    % factors (the MAD scale absorbs c) and a scale multiplied by c
    rs = RandStream('mt19937ar', 'Seed', 9);
    vw = randn(rs, 300, 1);  vw(7) = 12;  vw(40) = -9;
    for infl = ["huber", "igg3", "bisquare"]
        [f1, s1] = lsqvce.kernel.robustWeights(vw, infl, []);
        [f2, s2] = lsqvce.kernel.robustWeights(0.037*vw, infl, []);
        % tol: 1e-12 -- the standardised residuals are identical to rounding
        tc.verifyEqual(f2, f1, 'AbsTol', 1e-12, char(infl));
        % tol: as above
        tc.verifyEqual(s2, 0.037*s1, 'RelTol', 1e-12);
    end
end

function testPermutingTheRowsPermutesTheDiagnostics(tc)
% categories: metamorphic
    % a row permutation of (A, l, p) permutes lsqvce.engine.diagnose's
    % per-observation outputs and leaves the estimate alone
    A = tc.TestData.A;  l = tc.TestData.l;  p = tc.TestData.p;  n = numel(l);
    rs = RandStream('mt19937ar', 'Seed', 11);  q = randperm(rs, n);
    R1 = lsqvce.fit(struct('A', A, 'l', l, 'p', p), Verbose=false);
    R2 = lsqvce.fit(struct('A', A(q, :), 'l', l(q), 'p', p(q)), Verbose=false);
    % tol: 1e-11 relative -- one QR in a different row order
    tc.verifyEqual(R2.x, R1.x, 'RelTol', 1e-11);
    % tol: as above
    tc.verifyEqual(R2.diag.v, R1.diag.v(q), 'RelTol', 1e-10, 'AbsTol', 1e-14);
    % tol: as above
    tc.verifyEqual(R2.diag.rNum, R1.diag.rNum(q), 'RelTol', 1e-10);
end
