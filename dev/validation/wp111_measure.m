function wp111_measure(nSeed, outFile)
%WP111_MEASURE  The tier-A measurements WP-111 waits for (R7-01, R7-02, R7-20).
%
%   WP111_MEASURE(nSeed, outFile)
%
%   A MEASUREMENT, not a test: it asserts nothing and prints tables. It runs
%   on the self-hosted CI runner through .github/workflows/measure.yml; the
%   log and the diary are the record, the numbers go into HANDOVER by hand.
%   Every input is synthetic and seeded.
%
%   PART A  R7-02 / R7-20 -- the shipped joint estimator on ONE clean AR(1)
%     stream, law estimated (Grid 60, MaxPasses 10, the shipped defaults):
%     phi 0.5 (9 arcs), 0.8 (20 arcs), 0.9 (40 arcs) x 256 epochs. Per seed:
%     the stop reason, passes, the level r(0) = sigma^2 sigma0^2 against the
%     truth 1/(1-phi^2), lag-1 correlation, and the coefficients' z-scores
%     (x - xTrue) / sd. The review's mirror: phi 0.8 runs away to 1e7 x the
%     level and stops "converged" with |z| up to 37.
%
%   PART B  the CANDIDATE (user decision 2026-10-11: positive components) on
%     the SAME data, dense, standalone: the law is a positive spectrum on L
%     bins, piecewise linear in frequency; its covariance r = K S is the
%     exact integral (K closed form, no wrap); every bin's component
%     toeplitz(K(:, j)) is the ACV of a hat spectrum, hence PSD, so the
%     Helmert factor alpha2_j = (K' Omega)_j / (K' s)_j is >= 0 by
%     construction. Exact traces (dense), full lags. Mirror: all converge,
%     level 0.97-1.04 x truth, |z| <= 2.03.
%
%   PART C  R7-01 -- the declared PSD law's ACV through covModel 'psd'
%     'toeplitz' on even and odd arcs: r(m-1) / r(0) (a circulant wraps it
%     to r(1) or r(0)), against the closed form, and Toeplitz cond.
%
%   INPUT
%     nSeed    scalar   seeds per phi (default 6)
%     outFile  char     optional diary file
%
%   See also lsqvce.fit, lsqvce.engine.covModel
%
%   Created 2026-10-11 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- R7 WP-111: the tier-A measurement.
    if nargin < 1 || isempty(nSeed), nSeed = 6; end
    if nargin >= 2 && ~isempty(outFile)
        diary(outFile);
        cleanDiary = onCleanup(@() diary('off'));    %#ok<NASGU>
    end
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    oldPath = addpath(root);
    restorePath = onCleanup(@() path(oldPath));       %#ok<NASGU>
    w0 = warning('off', 'all');
    restoreWarn = onCleanup(@() warning(w0));         %#ok<NASGU>
    fprintf('WP-111 measurement, %s, MATLAB %s\n', ...
            char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), version);
    measureWrap();
    cases = [0.5 9; 0.8 20; 0.9 40];
    for c = 1:size(cases, 1)
        measureLoop(cases(c, 1), cases(c, 2), nSeed);
    end
end

%% ---------------------------------------------------------------- PART A + B
function measureLoop(phi, nArc, nSeed)
%MEASURELOOP  The shipped estimator and the candidate on the same data.
    m = 256;  L = 60;  xt = [0.8 -0.45 0.3 0.15];  truth = 1 / (1 - phi^2);
    fprintf('\n=== phi %.1f, %d arcs x %d, truth r0 %.3f, rho1 %.2f\n', phi, nArc, m, truth, phi);
    K = hatK(L, m);
    for s = 1:nSeed
        [A, l] = arcData(phi, nArc, m, s, xt);
        % PART A: the shipped estimator
        try
            tic;
            R = lsqvce.engine.fitDeclaration(fixture(A, l), lsqvce.Options('Estimate', "noise", 'Pool', "stream", ...
                                                         'MaxPasses', 10, 'Verbose', false));
            tA = toc;
            C = R.noise.Cgroups('x1|a');
            lvl = (R.noise.model{1}.sigma * R.sigma0)^2;
            zA = (R.x(:) - xt(:)) ./ sqrt(diag(R.sol.Cxx));
            fprintf(['seed %d SHIPPED  : %-9s passes %2d  r0 %.4g (x%.3g)  rho1 %.3f  ' ...
                     '|z|max %.2f  %.1fs\n'], s, char(R.stopReason), R.history.nIter, lvl, ...
                    lvl / truth, C(2) / C(1), max(abs(zA)), tA);
        catch err
            fprintf('seed %d SHIPPED  : FAILED %s (%s)\n', s, err.message, err.identifier);
        end
        % PART B: the candidate
        tic;  h = candidate(A, l, K, 30, 1e-2);  tB = toc;
        zB = (h.x - xt(:)) ./ sqrt(diag(h.Q));
        fprintf(['seed %d CANDIDATE: %-9s passes %2d  r0 %.4g (x%.3g)  rho1 %.3f  ' ...
                 '|z|max %.2f  min alpha2 %.3f  %.1fs\n'], s, h.stop, h.it, h.r(1), h.r(1) / truth, ...
                h.r(2) / h.r(1), max(abs(zB)), min(h.a2), tB);
    end
end

function h = candidate(A, l, K, maxPass, tol)
%CANDIDATE  Positive-component spectral VCE, dense (equal arc lengths).
    [m, L] = size(K);  nArc = numel(A);  S = ones(L, 1);
    h = struct('stop', 'maxPasses', 'it', maxPass);
    for it = 1:maxPass
        r = K * S;
        Sig = toeplitz(r);  W = chol(Sig);
        N = 0;  n = 0;
        for k = 1:nArc
            Ab = W' \ A{k};  lb = W' \ l{k};
            N = N + Ab' * Ab;  n = n + Ab' * lb;
        end
        x = N \ n;  Q = inv(N);
        Si = inv(Sig);  Si = (Si + Si') / 2;
        q = zeros(L, 1);  t = zeros(L, 1);
        for k = 1:nArc
            et = Sig \ (l{k} - A{k} * x);
            B = Si * A{k};
            Rbar = Si - B * Q * B';
            q = q + K' * lagSums(et * et');
            t = t + K' * lagSums(Rbar);
        end
        a2 = q ./ t;
        S = a2 .* S;
        if max(abs(a2 - 1)) <= tol
            h.stop = 'converged';  h.it = it;
            break
        end
    end
    h.r = K * S;  h.x = x;  h.Q = Q;  h.a2 = a2;
end

function v = lagSums(M)
%LAGSUMS  v(n+1) = sum of M over the +n and -n diagonals (n = 0..m-1).
    m = size(M, 1);  v = zeros(m, 1);
    v(1) = trace(M);
    for n = 1:m-1
        v(n+1) = sum(diag(M, -n)) + sum(diag(M, n));
    end
end

function K = hatK(L, m)
%HATK  r = K S: the exact ACV (sample lags 0..m-1) of the piecewise-linear
%   spectrum through S on f_j = j / (2(L-1)), two-sided, unit sample.
    d = 1 / (2 * (L - 1));  fj = (0:L-1) * d;  tau = (0:m-1)';
    w = 2 * pi * tau;
    env = (1 - cos(w * d)) ./ (d * w.^2);  env(1) = d / 2;
    c = 2 * ones(1, L);  c([1 L]) = 1;
    K = 2 * (env .* cos(w * fj)) .* c;
end

function [A, l] = arcData(phi, nArc, m, seed, xt)
%ARCDATA  nArc stationary AR(1) arcs (unit innovations) with a 4-term signal.
    rs = RandStream('mt19937ar', 'Seed', 1000 * seed + round(100 * phi));
    A = cell(1, nArc);  l = cell(1, nArc);
    w = 2 * pi / (0.5 * 86400);
    for k = 1:nArc
        t = ((k - 1) * m + (0:m-1)') * 5;
        A{k} = [sin(w*t), cos(w*t), sin(2*w*t), cos(2*w*t)];
        e = filter(1, [1 -phi], randn(rs, m + 2000, 1));
        l{k} = A{k} * xt(:) + e(2001:end);
    end
end

function data = fixture(A, l)
%FIXTURE  The same arcs as ONE day of ONE stream, law estimated.
    nArc = numel(A);  m = size(A{1}, 1);
    data = struct();
    data.commonParams = [repmat({'common'}, 4, 1), arrayfun(@(j) sprintf('c%d', j), (1:4).', ...
                                                            'UniformOutput', false)];
    data.days = 1;
    D = struct('A', {A}, 'B', {repmat({zeros(m, 0)}, 1, nArc)}, 'l', {l}, ...
               't', {arrayfun(@(k) ((k - 1) * m + (0:m-1)') * 5, 1:nArc, 'UniformOutput', false)});
    data.systems = struct('buildDay', @(day, ot) D, ...
                          'covariance', struct('a', struct('type', 'diag', 'sigma', 1, ...
                                                           'estimate', "law", 'arcFactors', false)), ...
                          'arcs', nArc, 'id', 1, 'name', 'wp111');
end

%% ---------------------------------------------------------------- PART C
function measureWrap()
%MEASUREWRAP  R7-01: the 'psd' 'toeplitz' ACV on even / odd arcs.
    fprintf('\n=== PART C (R7-01): AR(1) phi 0.5 as a 60-bin PSD law through covModel ''psd'' ''toeplitz''\n');
    L = 60;  phi = 0.5;
    fN = (0:L-1)' / (2 * (L - 1));
    S = 2 ./ abs(1 - phi * exp(-2i * pi * fN)).^2;   % ONE-sided, covModel's convention (r = int_0^fNyq S cos)
    for m = [256 255 2111]
        cm = struct('type', 'psd', 'fNorm', fN, 'S', S);
        [Cfcn, ~, ~] = lsqvce.engine.covModel(cm, (0:m-1)' * 5, struct('psdMode', 'toeplitz', 'BatchSize', 20000));
        r = Cfcn(1:m, struct());  r = r(:, 1);
        Kc = hatK(L, m) * (S / 2);                     % hatK takes the two-sided density
        fprintf(['m %4d | shipped: r(m-1)/r(0) %.4f  r(1)/r(0) %.4f  cond %.2e | closed form: ' ...
                 'r(m-1)/r(0) %.2e  r(1)/r(0) %.4f  cond %.2e | truth r(1)/r(0) %.4f\n'], m, r(m) / r(1), ...
                r(2) / r(1), cond(toeplitz(r)), Kc(m) / Kc(1), Kc(2) / Kc(1), cond(toeplitz(Kc)), phi);
    end
end
