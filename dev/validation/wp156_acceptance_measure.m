function T = wp156_acceptance_measure(outFile)
%WP156_ACCEPTANCE_MEASURE  WP-156's acceptance bar, measured on the full loop.
%
%   T = WP156_ACCEPTANCE_MEASURE(outFile)
%
%   A MEASUREMENT, not a test: it asserts nothing and prints; it runs on the
%   self-hosted CI runner through .github/workflows/measure.yml and the job
%   log is the record (HANDOVER 12.91 quotes it). The bar (user decision
%   2026-10-18): paired against Robust "none" on clean data over 20 seeds,
%   every stream's variance and sigma0^2 within 0.5 % for phi <= 0.98, IGG3
%   and bisquare; with 2 % blunders of 2 sigma at phi 0.9 the red stream's
%   level within 2 % and >= 95 % of its blunders rejected; convergence in
%   every seed. The fixture is the design study's (lsqvce.fit, Estimate
%   "noise"): 3 days x 4 arcs x 256 epochs at 5 s, 6 commons, one bias per
%   arc shared by the streams; stream a white sigma 0.5, stream b AR(1) phi
%   sigma 1.2; both declared white sigma 1 and estimated ("law").
%     PART A  clean, phi_b 0.5 / 0.9 / 0.98: per influence the paired
%             variance ratio - 1 of each stream and of sigma0^2 (mean +- se),
%             converged runs, runs with a frozen cycle
%     PART B  2 % x 2 sigma blunders on both streams, phi_b 0.9, 12 seeds:
%             the red stream's variance against the clean plain run, its
%             lag-1 correlation, the share of its blunders rejected
%
%   INPUT   outFile  char  optional diary file
%   OUTPUT  T  struct, the per-seed numbers of both parts
%
%   See also wp156_selection_measure
%
%   Created 2026-10-18 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- WP-156 round 2, the acceptance bar.
    if nargin > 0 && ~isempty(outFile)
        diary(outFile);  cleanDiary = onCleanup(@() diary('off')); %#ok<NASGU>
    end
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    p0 = path;  restore = onCleanup(@() path(p0)); %#ok<NASGU>
    addpath(root);
    w0 = warning('off', 'lsqvce:adjustVCEjoint:notConverged');  cw = onCleanup(@() warning(w0)); %#ok<NASGU>
    fprintf('WP-156 acceptance measurement, %s, MATLAB %s\n', char(datetime('now')), version);
    base = {'Estimate', "noise", 'Pool', "stream", 'ArcFactors', false, 'MaxPasses', 30, 'Verbose', 0};
    T = struct('clean', struct([]), 'blunders', struct([]));

    fprintf('\n===== PART A: clean, paired against Robust "none", 20 seeds =====\n');
    fprintf('%-6s %-9s | %-22s | %-22s | %-22s | conv  frozen\n', 'phi_b', 'influence', 'var a - 1', 'var b - 1', 'sigma0^2 - 1');
    for phi = [0.5 0.9 0.98]
        R0 = cell(1, 20);
        for sd = 1:20
            R0{sd} = lsqvce.fit(fixture(phi, sd, [0 0]), lsqvce.Options(base{:}, 'Seed', sd));
        end
        for infl = ["igg3", "bisquare"]
            r = nan(20, 3);  cv = 0;  fz = 0;
            for sd = 1:20
                Rr = lsqvce.fit(fixture(phi, sd, [0 0]), lsqvce.Options(base{:}, 'Seed', sd, 'Robust', infl));
                r(sd, :) = [varOf(Rr, 'a') / varOf(R0{sd}, 'a'), varOf(Rr, 'b') / varOf(R0{sd}, 'b'), ...
                            (Rr.sigma0 / R0{sd}.sigma0)^2] - 1;
                cv = cv + Rr.converged;  fz = fz + (Rr.noise.robust.frozen > 0);
            end
            fprintf('%-6.2f %-9s | %s | %s | %s | %2d/20 %2d\n', phi, infl, pm(r(:, 1)), pm(r(:, 2)), pm(r(:, 3)), cv, fz);
            T.clean(end+1).phi = phi;  T.clean(end).influence = infl;  T.clean(end).ratios = r; %#ok<AGROW>
            T.clean(end).converged = cv;  T.clean(end).frozen = fz;
        end
    end

    fprintf('\n===== PART B: 2 %% x 2 sigma blunders on both streams, phi_b 0.9, 12 seeds =====\n');
    fprintf('%-9s | %-22s | %-12s | %-14s | conv\n', 'influence', 'var b / clean plain - 1', 'lag-1 b', 'blunders b');
    Rc = cell(1, 12);
    for sd = 1:12
        Rc{sd} = lsqvce.fit(fixture(0.9, sd, [0 0]), lsqvce.Options(base{:}, 'Seed', sd));
    end
    for infl = ["none", "igg3"]
        r = nan(12, 1);  ph = r;  hit = 0;  pl = 0;  cv = 0;
        for sd = 1:12
            Rr = lsqvce.fit(fixture(0.9, sd, [0.02 2]), lsqvce.Options(base{:}, 'Seed', sd, 'Robust', infl));
            r(sd) = varOf(Rr, 'b') / varOf(Rc{sd}, 'b') - 1;
            Cb = Rr.noise.Cgroups(lsqvce.util.streamKey(1, 'b'));  ph(sd) = Cb(2) / Cb(1);
            cv = cv + Rr.converged;
            if infl ~= "none"
                [h, n] = caught(Rr, sd);  hit = hit + h;  pl = pl + n;
            end
        end
        fprintf('%-9s | %s | %.3f+-%.3f | %5d of %-5d | %2d/12\n', infl, pm(r), mean(ph), std(ph) / sqrt(12), hit, pl, cv);
        T.blunders(end+1).influence = infl;  T.blunders(end).varRatio = r; %#ok<AGROW>
        T.blunders(end).lag1 = ph;  T.blunders(end).caught = [hit pl];  T.blunders(end).converged = cv;
    end
    fprintf('\nDone, %s\n', char(datetime('now')));
end

function s = pm(x)
    s = sprintf('%+7.3f %% +- %5.3f', 100 * mean(x), 100 * std(x) / sqrt(numel(x)));
end

function v = varOf(R, st)
%VAROF  A stream's variance, sigma x sigma0, squared.
    v = (R.noise.model{find(strcmp(R.noise.stream, lsqvce.util.streamKey(1, st)), 1)}.sigma * R.sigma0)^2;
end

function [hit, n] = caught(R, sd)
%CAUGHT  Stream b's planted blunders at factor 0 (the robust pass's keys).
    f = R.noise.robust.factor;  k = R.noise.robust.keys;  off = 0;  hit = 0;  n = 0;
    for j = 1:size(k, 1)
        if strcmp(k{j, 2}, 'b')
            idx = blunderRows(sd, k{j, 3}, k{j, 4}, 2, 0.02);
            hit = hit + nnz(f(off + idx) == 0);  n = n + numel(idx);
        end
        off = off + 256;
    end
end

function data = fixture(phiB, seed, blunder)
%FIXTURE  The design study's closed loop (scratch harness mcScenario).
    cfg = struct('days', 1:3, 'arcs', 4, 'm', 256, 'dt', 5, 'nc', 6, 'streams', {{'a', 'b'}}, ...
                 'level', [0.5 1.2], 'phi', [0 phiB], 'seed', seed, 'blunder', blunder);
    rsx = RandStream('mt19937ar', 'Seed', 90000 + seed);
    xTrue = randn(rsx, cfg.nc, 1);
    bias = randn(rsx, numel(cfg.days) * cfg.arcs, 1);
    data = struct();
    data.commonParams = [repmat({'common'}, cfg.nc, 1), ...
        arrayfun(@(j) sprintf('c%d', j), (1:cfg.nc).', 'UniformOutput', false)];
    data.days = cfg.days;
    spec = struct('type', 'diag', 'sigma', 1, 'estimate', "law", 'arcFactors', false);
    s = struct();
    s.buildDay = @(day, ot, common, local, po, info) fixtureDay(cfg, xTrue, bias, day, ot, info);
    s.covariance = struct('a', spec, 'b', spec);
    s.arcs = cfg.arcs;  s.id = 1;  s.name = 'wp156';
    s.localParams = {'arc', 'bias'};
    data.systems = s;
end

function D = fixtureDay(cfg, xTrue, bias, day, obsType, info)
    k = find(strcmp(cfg.streams, obsType), 1);
    di = find(cfg.days == day, 1);
    D = struct('A', {cell(1, cfg.arcs)}, 'B', {cell(1, cfg.arcs)}, 'l', {cell(1, cfg.arcs)}, 't', {cell(1, cfg.arcs)});
    w = 2 * pi / (0.5 * 86400);
    for a = 1:cfg.arcs
        t = (day - 1) * 86400 + ((a - 1) * cfg.m + (0:cfg.m - 1)).' * cfg.dt;
        A = [sin(w*t), cos(w*t), sin(2*w*t), cos(2*w*t), sin(3*w*t), cos(3*w*t)];
        phi = cfg.phi(k);
        rs = RandStream('mt19937ar', 'Seed', mod(cfg.seed * 1000003 + day * 7919 + a * 104729 + k * 13, 2^31 - 1));
        nB = ceil(10 / max(1 - abs(phi), 1e-3));
        e = filter(1, [1 -phi], cfg.level(k) * sqrt(1 - phi^2) * randn(rs, cfg.m + nB, 1));
        e = e(nB+1:end);
        if cfg.blunder(1) > 0
            idx = blunderRows(cfg.seed, day, a, k, cfg.blunder(1));
            rb = RandStream('mt19937ar', 'Seed', mod(cfg.seed * 31337 + day * 101 + a * 7 + k, 2^31 - 1));
            e(idx) = e(idx) + cfg.blunder(2) * cfg.level(k) * sign(randn(rb, numel(idx), 1) + eps);
        end
        B = zeros(cfg.m, info.nl);
        B(:, info.offsetOf(a) + 1) = 1;
        g = (di - 1) * cfg.arcs + a;
        D.A{a} = A;  D.B{a} = B;  D.t{a} = t;
        D.l{a} = A * xTrue + bias(g) + e;
    end
end

function idx = blunderRows(seed, day, arc, k, fraction)
%BLUNDERROWS  The planted epochs of one arc of stream k (of 256).
    rs = RandStream('mt19937ar', 'Seed', mod(seed * 7777 + day * 31 + arc * 1009 + k * 17, 2^31 - 1));
    idx = sort(randperm(rs, 256, round(fraction * 256))).';
end
