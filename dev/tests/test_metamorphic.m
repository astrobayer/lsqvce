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
    pathBefore = path;   % R7-45: a standalone run (run_lsqvce / runtests) leaves the path as found
    tc.addTeardown(@() path(pathBefore));
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
    % diagnostics lsqvce.engine.verdict.diagnose / lsqvce.engine.verdict.source)
    A = tc.TestData.A;  l = tc.TestData.l;  p = tc.TestData.p;  c = -3.7;
    R1 = lsqvce.engine.fitDeclaration(struct('A', A, 'l', l,   'p', p), Verbose=false);
    R2 = lsqvce.engine.fitDeclaration(struct('A', A, 'l', c*l, 'p', p), Verbose=false);
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
    R1 = lsqvce.engine.fitDeclaration(struct('A', A, 'l', l, 'p', p),   Verbose=false);
    R2 = lsqvce.engine.fitDeclaration(struct('A', A, 'l', l, 'p', k*p), Verbose=false);
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
    R1 = lsqvce.engine.fitDeclaration(struct('A', A, 'l', l, 'p', p), Verbose=false);
    for nb = [3 7]
        edges = round(linspace(0, n, nb + 1));
        bf = @(k) deal(A(edges(k)+1:edges(k+1), :), l(edges(k)+1:edges(k+1)), p(edges(k)+1:edges(k+1)));
        Rk = lsqvce.engine.fitDeclaration(struct('batchFcn', bf, 'numUnknowns', 5, 'numBatches', nb), Verbose=false);
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
    % the ACV belongs to the process: with innovation std 3 it is 9*r1 (the
    % exact start-up rows use it; the filter rows use s directly)
    [WX3, ~, Wl3] = lsqvce.kernel.armaWhiten(a, b, 3, 9*r1, X,     zeros(m, 0), l, true);
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

function testTheMadScaleLeavesHalfTheSampleAlive(tc)
% categories: metamorphic, robustness
    % R6 WP-109 (issue #110 test 6): WHY THE COLLAPSE GUARD IS ALMOST
    % UNREACHABLE, as a property of the scale rather than of any one fit.
    %
    % lsqvce.engine.adjustRobust refuses a pass when nAlive = sum(fac > 0)
    % falls below nc + 1. But robustWeights standardises by the normalised
    % MAD, s = 1.4826 * median(|vw - med|), and by the definition of a median
    % at least half the sample lies within MAD of med. A redescending
    % influence zeroes only beyond k1*s = 4.45*MAD (igg3) or c*s = 6.95*MAD
    % (bisquare), so that half survives whenever |med| <= 3.4*MAD -- and
    %
    %     nAlive >= ceil(n/2)   for ANY sample, however contaminated.
    %
    % WP-110 CORRECTS WHAT WP-109 CONCLUDED FROM THIS BOUND. The bound is
    % real, and it is what the four samples below verify. The CONSEQUENCE
    % 12.33 drew from it -- "the guard can fire only for n < 2*(nc + 1), and
    % cannot fire at any contamination" -- is FALSE, because it ignored the
    % bound's own premise. The premise |med| <= 3.4*MAD fails exactly at the
    % 50% breakdown point (measured: it holds at 49% contamination and
    % fails at 51%), so this bound says no more than "we are below
    % breakdown". Above it nAlive falls with the clean fraction, and the
    % guard does fire. The three points below pin that; 12.34 has the
    % measurement.
    %
    % The four samples are adversarial on purpose: clean, half gross, a
    % bimodal set two orders apart, and a Cauchy one with no variance at all.
    rs = RandStream('mt19937ar', 'Seed', 109);
    n  = 400;
    sets = { randn(rs, n, 1), ...
             [randn(rs, n/2, 1); 80*randn(rs, n/2, 1)], ...
             [zeros(n/2, 1); 1e4*ones(n/2, 1)], ...
             tan(pi*(rand(rs, n, 1) - 0.5)) };
    for k = 1:numel(sets)
        vw = sets{k};
        for infl = ["igg3", "bisquare"]
            [fac, s] = lsqvce.kernel.robustWeights(vw, infl, []);
            mad = s / 1.4826;
            tc.verifyLessThanOrEqual(abs(median(vw)), 3.4 * mad, ...
                'the stated premise: the residuals are not offset far from zero');
            tc.verifyGreaterThanOrEqual(sum(fac > 0), ceil(n/2), ...
                sprintf('sample %d, %s: the MAD scale must leave half the sample alive', k, infl));
        end
    end
end

function testTheHalfSampleBoundsPremiseFailsExactlyAtBreakdown(tc)
% categories: metamorphic, robustness
    % R6 WP-110, and the CORRECTION of WP-109. The bound above is conditional
    % on |med| <= 3.4*MAD, and that premise is not a technicality: it is the
    % breakdown point itself, wearing a different hat. Below 50% the median
    % sits in the clean cluster and the ratio is O(1); above it the median
    % moves into the contaminated cluster and the ratio explodes.
    %
    % Deterministic by construction -- no RandStream, so the numbers here are
    % the same arithmetic on every runner and were predicted exactly before
    % the run (0.50 below, 518.89 above).
    n = 1000;
    [below, madB] = contaminated(n, 490);         % 49% contaminated
    [above, madA] = contaminated(n, 510);         % 51% contaminated
    tc.verifyLessThanOrEqual(abs(median(below)) / madB, 3.4, ...
        'below breakdown the premise of the half-sample bound holds');
    tc.verifyGreaterThan(abs(median(above)) / madA, 3.4, ...
        'above breakdown it fails -- which is why the bound says nothing there');
end

function testTheHalfSampleBoundHoldsOnlyBelowBreakdown(tc)
% categories: metamorphic, robustness
    % The other half of the correction: the bound itself flips with its
    % premise, at the same place. 510 alive at 49% contamination, 490 at
    % 51% -- so nAlive >= ceil(n/2) is a RESTATEMENT of being below
    % breakdown, not an independent limit on what the collapse guard can see.
    n = 1000;
    below = contaminated(n, 490);
    above = contaminated(n, 510);
    tc.verifyGreaterThanOrEqual(sum(lsqvce.kernel.robustWeights(below, "igg3", []) > 0), ...
        ceil(n/2), 'below breakdown at least half the sample survives');
    tc.verifyLessThan(sum(lsqvce.kernel.robustWeights(above, "igg3", []) > 0), ...
        ceil(n/2), 'above breakdown it does not');
end

function testContaminationLeavesExactlyTheCleanRowsAlive(tc)
% categories: metamorphic, robustness
    % What the collapse guard actually sees under contamination: nAlive is
    % the CLEAN count, to the row. The guard compares that against an
    % ABSOLUTE nc + 1, so on n rows it needs the clean count below nc + 1 --
    % 99.8% contamination on test_fit's 3600-row fixture against nc = 6.
    % "More than half the data is discrepant", which the error message
    % names, is nowhere near the trigger: at 60% contamination 400 of 1000
    % rows are still alive.
    %
    % Standardisation is by |vw|, not |vw - med| (the Yang convention the
    % kernel documents), which is why the contaminated MAJORITY is the part
    % that gets zeroed above breakdown, not the clean minority.
    n = 1000;
    for k = [600 900 995]
        vw = contaminated(n, k);
        for infl = ["igg3", "bisquare"]
            tc.verifyEqual(sum(lsqvce.kernel.robustWeights(vw, infl, []) > 0), n - k, ...
                sprintf('%s, %d contaminated: the clean rows, and only those, stay alive', infl, k));
        end
    end
end

function testAUniformResidualOffsetRejectsEveryRow(tc)
% categories: metamorphic, robustness
    % The route to the collapse guard that DOES exist on a full-sized problem,
    % and it is not an outlier problem at all: residuals offset bodily from
    % zero by several MAD, with no contamination whatsoever. Because the
    % influence standardises vw/s rather than (vw - med)/s, a systematic
    % offset -- an unabsorbed bias, a datum defect -- rejects the whole
    % sample. That is exactly the case adjustRobust's message should name,
    % and after WP-110 it does.
    n = 1000;
    vw = 6 + linspace(-1, 1, n).';            % MAD ~ 0.5, so the offset is ~12 MAD
    for infl = ["igg3", "bisquare"]
        tc.verifyEqual(sum(lsqvce.kernel.robustWeights(vw, infl, []) > 0), 0, ...
            sprintf('%s: a bodily offset rejects every row', infl));
    end
end

function [vw, mad] = contaminated(n, k)
%CONTAMINATED  n law-whitened residuals of which k are shifted to +1000.
%   Deterministic: each block is a linspace, so there is no RandStream and
%   the numbers are identical on every runner. The spread keeps the MAD away
%   from zero, which would otherwise take robustWeights' std fallback.
    vw = [linspace(-1, 1, n - k).'; 1000 + linspace(-1, 1, k).'];
    [~, s] = lsqvce.kernel.robustWeights(vw, "igg3", []);
    mad = s / 1.4826;
end

function testPermutingTheRowsPermutesTheDiagnostics(tc)
% categories: metamorphic
    % a row permutation of (A, l, p) permutes lsqvce.engine.verdict.diagnose's
    % per-observation outputs and leaves the estimate alone
    A = tc.TestData.A;  l = tc.TestData.l;  p = tc.TestData.p;  n = numel(l);
    rs = RandStream('mt19937ar', 'Seed', 11);  q = randperm(rs, n);
    R1 = lsqvce.engine.fitDeclaration(struct('A', A, 'l', l, 'p', p), Verbose=false);
    R2 = lsqvce.engine.fitDeclaration(struct('A', A(q, :), 'l', l(q), 'p', p(q)), Verbose=false);
    % tol: 1e-11 relative -- one QR in a different row order
    tc.verifyEqual(R2.x, R1.x, 'RelTol', 1e-11);
    % tol: as above
    tc.verifyEqual(R2.diag.v, R1.diag.v(q), 'RelTol', 1e-10, 'AbsTol', 1e-14);
    % tol: as above
    tc.verifyEqual(R2.diag.rNum, R1.diag.rNum(q), 'RelTol', 1e-10);
end
