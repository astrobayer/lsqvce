function wp115_measure(nSeed, outFile)
%WP115_MEASURE  The tier-A measurements WP-115 part 3 waits for (R7-12, R7-10).
%
%   WP115_MEASURE(nSeed, outFile)
%
%   A MEASUREMENT, not a test: it asserts nothing and prints tables. It runs
%   on the self-hosted CI runner through .github/workflows/measure.yml (the
%   log and the uploaded text file are the record; the numbers go into
%   HANDOVER by hand, never as a data file into the repository). Every input
%   is synthetic and seeded.
%
%   PART A  R7-12 -- is the robust joint loop's level biased low?
%     Two clean AR(1) streams (phi 0.2 and 0.5, unit innovations, 3 days x
%     3 arcs x 256 epochs, the colour fixture of test_vcePolicy), NO
%     blunders. Per seed four fits through lsqvce.fit, Estimate = "noise":
%       law   / plain   law estimated, Robust "none"
%       law   / igg3    law estimated, Robust "igg3"
%       scale / plain   the TRUE law declared ('arma'), only the level free
%       scale / igg3    the same, robust
%     Reported per stream: the effective level sigma*sigma0 (1 = unbiased in
%     scale mode) and the robust/plain ratio on the same data (law mode).
%     The review's mirror read 0.83-0.86 (law fixed) / ~0.93 (law free); one
%     MATLAB reading was 1.007.
%
%   PART B  R7-10 -- the Koch fold against deletion and the additive nugget
%     Dense GLS on one arc, m = 256, A = [1, trend, sinusoid], AR(1) law R
%     (unit marginal variance), phi in {0.5, 0.9}, k epochs with robust
%     factor f. The three covariances of a down-weighted epoch set:
%       fold     C = D R D, D_ii = 1/sqrt(f_i)       (inflateCep today)
%       nugget   C = R + diag(1/f_i - 1)             (same marginal variance)
%       delete   the epochs removed                  (the f -> 0 reference)
%     True noise R (clean data down-weighted, IGG3's tail case). Reported:
%     formal/true sd (1 = honest), true sd / deletion's sd, and the lag-0
%     level expectation E[v'C^-1 v]/(m-n) under the true noise (1 = the
%     level an unbiased Helmert step would return).
%     Then the COST of the nugget per arc: dense chol of R + diag(c) against
%     a Woodbury update on the AR(1) precision (tridiagonal), which is the
%     low-rank path a nugget arc would need, with their agreement.
%
%   INPUT
%     nSeed    scalar   Part A seeds (default 20)
%     outFile  char     optional diary file (default: none)
%
%   See also lsqvce.fit, lsqvce.Options
%
%   Created 2026-10-11 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- R7 WP-115 part 3: the tier-A
%   measurement.
    if nargin < 1 || isempty(nSeed), nSeed = 20; end
    if nargin >= 2 && ~isempty(outFile)
        diary(outFile);
        cleanDiary = onCleanup(@() diary('off'));    %#ok<NASGU>
    end
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    oldPath = addpath(root);
    restorePath = onCleanup(@() path(oldPath));       %#ok<NASGU>
    w0 = warning('off', 'lsqvce:adjustVCEjoint:notConverged');
    restoreWarn = onCleanup(@() warning(w0));         %#ok<NASGU>
    fprintf('WP-115 part 3 measurement, %s, MATLAB %s\n', ...
            char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), version);
    measureFold();
    measureNuggetCost();
    measureLevel(nSeed);
end

%% ------------------------------------------------------------ PART A
function measureLevel(nSeed)
%MEASURELEVEL  R7-12: the effective level of two clean streams, per seed.
    phi = struct('a', 0.2, 'b', 0.5);
    runs = {'law', "none"; 'law', "igg3"; 'scale', "none"; 'scale', "igg3"};
    lvl = NaN(nSeed, size(runs, 1), 2);  pass = NaN(nSeed, size(runs, 1));
    zeroed = NaN(nSeed, 2);  secs = NaN(nSeed, size(runs, 1));
    fprintf('\n=== PART A (R7-12): effective level sigma*sigma0, clean AR(1) streams a (phi 0.2), b (phi 0.5)\n');
    for s = 1:nSeed
        for r = 1:size(runs, 1)
            data = levelFixture(phi, s, runs{r, 1});
            opt = lsqvce.Options('Estimate', "noise", 'Pool', "stream", 'MaxPasses', 16, ...
                                 'Verbose', false, 'Robust', runs{r, 2});
            try
                tic;  R = lsqvce.engine.fitDeclaration(data, opt);  secs(s, r) = toc;
            catch err
                fprintf('seed %d run %s/%s FAILED: %s (%s)\n', s, runs{r, 1}, runs{r, 2}, ...
                        err.message, err.identifier);
                continue
            end
            for j = 1:2
                k = find(strcmp(R.noise.stream, sprintf('x1|%s', char('a' + j - 1))), 1);
                lvl(s, r, j) = R.noise.model{k}.sigma * R.sigma0;
            end
            pass(s, r) = R.history.nIter;
            if r == 4 && ~isempty(R.noise.robust.factor)
                zeroed(s, :) = [mean(R.noise.robust.factor == 0), mean(R.noise.robust.factor)];
            end
        end
        fprintf('seed %2d  a: %s   b: %s   passes %s\n', s, ...
                mat2str(squeeze(lvl(s, :, 1)), 4), mat2str(squeeze(lvl(s, :, 2)), 4), mat2str(pass(s, :)));
    end
    fprintf('\ncolumns: law/plain  law/igg3  scale/plain  scale/igg3\n');
    for j = 1:2
        L = lvl(:, :, j);
        fprintf('stream %c  mean %s  sd %s  (n = %s)\n', char('a' + j - 1), ...
                mat2str(mean(L, 1, 'omitnan'), 4), mat2str(std(L, 0, 1, 'omitnan'), 3), ...
                mat2str(sum(isfinite(L), 1)));
        q = L(:, 2) ./ L(:, 1);   p = L(:, 4) ./ L(:, 3);
        fprintf('          robust/plain, law mode   mean %.4f  sd %.4f  se %.4f\n', ...
                mean(q, 'omitnan'), std(q, 'omitnan'), std(q, 'omitnan') / sqrt(sum(isfinite(q))));
        fprintf('          robust/plain, scale mode mean %.4f  sd %.4f  se %.4f\n', ...
                mean(p, 'omitnan'), std(p, 'omitnan'), std(p, 'omitnan') / sqrt(sum(isfinite(p))));
    end
    fprintf('scale/igg3: zeroed fraction %.4f, mean factor %.4f; seconds per fit %s\n', ...
            mean(zeroed(:, 1), 'omitnan'), mean(zeroed(:, 2), 'omitnan'), mat2str(mean(secs, 1, 'omitnan'), 3));
end

function data = levelFixture(phi, seed, mode)
%LEVELFIXTURE  Two clean AR(1) streams, no locals, seeded on (seed, day, stream).
    cfg = struct('nc', 4, 'days', 1:3, 'arcs', 3, 'mArc', 256, 'dt', 5, ...
                 'xTrue', [0.8 -0.45 0.3 0.15]);
    data = struct();
    data.commonParams = [repmat({'common'}, cfg.nc, 1), ...
                         arrayfun(@(j) sprintf('c%d', j), (1:cfg.nc).', 'UniformOutput', false)];
    data.days = cfg.days;
    if strcmp(mode, 'law')
        spec = @(p) struct('type', 'diag', 'sigma', 1, 'estimate', "law", 'arcFactors', false);
    else
        spec = @(p) struct('type', 'arma', 'ar', [1 -p], 'sigma', 1, 'estimate', "scale", ...
                           'arcFactors', false);
    end
    data.systems = struct('buildDay', @(day, ot) levelDay(cfg, phi.(ot), seed, day, ot), ...
                          'covariance', struct('a', spec(phi.a), 'b', spec(phi.b)), ...
                          'arcs', cfg.arcs, 'id', 1, 'name', 'wp115level');
end

function D = levelDay(cfg, phi, seed, day, obsType)
%LEVELDAY  One (day, stream): v6 cells, one per arc, stationary AR(1) noise.
    m = cfg.arcs * cfg.mArc;
    t = (day - 1)*86400 + (0:m-1).' * cfg.dt;
    w = 2*pi / (0.5 * 86400);
    A = [sin(w*t), cos(w*t), sin(2*w*t), cos(2*w*t)];
    rs = RandStream('mt19937ar', 'Seed', 1e5*seed + 10*day + double(obsType(1)));
    e = filter(1, [1 -phi], randn(rs, m + 2000, 1));
    l = A * cfg.xTrue(:) + e(2001:end);
    D = struct('A', {cell(1, cfg.arcs)}, 'B', {cell(1, cfg.arcs)}, ...
               'l', {cell(1, cfg.arcs)}, 't', {cell(1, cfg.arcs)});
    for a = 1:cfg.arcs
        r = (a-1)*cfg.mArc + (1:cfg.mArc);
        D.A{a} = A(r, :);  D.B{a} = zeros(cfg.mArc, 0);  D.l{a} = l(r);  D.t{a} = t(r);
    end
end

%% ------------------------------------------------------------ PART B
function measureFold()
%MEASUREFOLD  R7-10: formal/true sd and level expectation, fold / nugget / delete.
    m = 256;  t = (0:m-1).' * 5;
    A = [ones(m, 1), (t - mean(t)) / (max(t) - min(t)), sin(2*pi*t/1500)];
    fprintf('\n=== PART B (R7-10): one arc m = %d, A = [1 trend sinusoid], truth = the declared AR(1) law\n', m);
    fprintf('    sd ratios on the sinusoid coefficient; lvl = E[v''C^-1 v]/(m-n) under the true noise\n');
    rs = RandStream('mt19937ar', 'Seed', 9);
    for phi = [0.5 0.9]
        R = toeplitz(phi .^ (0:m-1));
        for f = [1e-12 0.3]
            for k = [1 5 13 51]
                z = sort(randperm(rs, m - 10, k) + 5);
                keep = setdiff(1:m, z);
                fv = ones(m, 1);  fv(z) = f;
                Cf = R ./ sqrt(fv) ./ sqrt(fv.');            % fold
                Cn = R + diag(1 ./ fv - 1);                   % nugget
                [qf, ef, lf] = gls(A, Cf, R);
                [qn, en, ln] = gls(A, Cn, R);
                [~, ed] = gls(A(keep, :), R(keep, keep), R(keep, keep));
                fprintf(['phi %.1f f %-6.0e k %2d | fold: formal/true %.3f true/del %.3f lvl %.3f' ...
                         ' | nugget: formal/true %.3f true/del %.3f lvl %.3f\n'], phi, f, k, ...
                        sqrt(qf(3)/ef(3)), sqrt(ef(3)/ed(3)), lf, ...
                        sqrt(qn(3)/en(3)), sqrt(en(3)/ed(3)), ln);
            end
        end
    end
end

function [q, e, lv] = gls(A, C, Rt)
%GLS  Formal and true variances of the GLS under C when the noise is Rt, and
%   the level a Helmert step under C expects: E[v'C^-1 v]/(m-n).
    G = (A' / C * A) \ (A' / C);
    Q = inv(A' / C * A);
    E = G * Rt * G';
    P = eye(size(A, 1)) - A * G;
    lv = trace(C \ (P * Rt * P')) / (size(A, 1) - size(A, 2));
    q = diag(Q);  e = diag(E);
end

function measureNuggetCost()
%MEASURENUGGETCOST  The per-arc cost of a nugget arc: dense chol vs Woodbury.
    fprintf('\n=== PART B cost: per-arc normal contribution A''C^-1 A, nc = 4, AR(1) phi 0.9\n');
    fprintf('    toeplitz = the fold''s path today (one chol of R, per arc a descale); dense = chol(R + diag c);\n');
    fprintf('    woodbury = tridiagonal precision of R plus a k x k capacitance; err = max rel diff to dense\n');
    phi = 0.9;
    for m = [256 1024 4096]
        t = (0:m-1).' * 5;
        A = [sin(2*pi*t/4e4), cos(2*pi*t/4e4), sin(4*pi*t/4e4), cos(4*pi*t/4e4)];
        R = toeplitz(phi .^ (0:m-1));
        Pr = arPrecision(phi, m);                       % inv(R), tridiagonal, sparse
        for frac = [0.02 0.2]
            rs = RandStream('mt19937ar', 'Seed', m);
            k = max(1, round(frac * m));
            z = sort(randperm(rs, m, k));
            c = 1 ./ (0.3 + 0.7 * rand(rs, k, 1)) - 1;  c(1:ceil(k/3)) = 1e12;   % partial and zeroed
            Rc = chol(R);
            tT = timeit(@() foldNormal(Rc, A, ones(m, 1)));
            tD = timeit(@() denseNormal(R, A, z, c));
            tW = timeit(@() woodburyNormal(Pr, A, z, c));
            Nd = denseNormal(R, A, z, c);  Nw = woodburyNormal(Pr, A, z, c);
            fprintf('m %5d k %4d | toeplitz %8.2e s  dense %8.2e s  woodbury %8.2e s | err %.1e\n', ...
                    m, k, tT, tD, tW, max(abs(Nd(:) - Nw(:))) / max(abs(Nd(:))));
        end
    end
end

function N = foldNormal(Rc, A, d)
%FOLDNORMAL  The fold path: descale, one triangular solve with chol(R).
    Ab = Rc' \ (A .* d);
    N = Ab' * Ab;
end

function N = denseNormal(R, A, z, c)
%DENSENORMAL  The nugget by brute force.
    C = R;  idx = sub2ind(size(C), z, z);  C(idx) = C(idx) + c(:).';
    Ab = chol(C)' \ A;
    N = Ab' * Ab;
end

function N = woodburyNormal(Pr, A, z, c)
%WOODBURYNORMAL  inv(R + U diag(c) U') = Pr - Pr U (diag(1/c) + U' Pr U)^-1 U' Pr.
    PA = Pr * A;
    K = diag(1 ./ c) + full(Pr(z, z));
    N = A' * PA - PA(z, :)' * (K \ PA(z, :));
end

function P = arPrecision(phi, m)
%ARPRECISION  inv(toeplitz(phi.^(0:m-1))): tridiagonal.
    d = (1 + phi^2) * ones(m, 1);  d([1 m]) = 1;
    e = -phi * ones(m, 1);
    P = spdiags([e d e], -1:1, m, m) / (1 - phi^2);
end
