function bench_m0(opt)
%BENCH_M0  R5 M0 (WP-79): the three measurements that decide the M2 deletions.
%   bench_m0
%   bench_m0('Quick', true)          % smaller sizes, ~2 min, for a smoke run
%   bench_m0('OutFile', path)        % where the record goes
%
%   Tier B (no MATLAB in the container): the user runs this locally and
%   uploads dev/validation/bench_m0_result.txt; the numbers are pinned in
%   dev/DEVELOPMENT.md and cited by the M2 deletion commits (HANDOVER
%   §12.5 WP-79). Nothing here is a pass/fail test: every line is a RATIO
%   with its band from the ONE timing instrument dev/tests/private/timePair
%   (called here through a local copy of its logic, because a private
%   function is not reachable from dev/validation -- the copy is deleted
%   with this script after M2; it is not a second instrument).
%
%   THE THREE QUESTIONS
%     (a) PrefetchLoader: does overlapping the loader with the accumulation
%         beat plain serial accumulation, and does either beat
%         DistributeDays=true over a parfor pool? If prefetch never wins,
%         H-01 class (iii) deletes it (R5-D5).
%     (b) exact leverage sweep vs matrix-free probe sweep at nc in
%         {500, 2000, 4000}: engine.redundancy (exact per-observation r_i)
%         against engine.streamContribution (Hutchinson .r with 100 probes).
%         Decides the v5 sweep's default route per nc (theory_r5 §3.3).
%     (c) block-accumulator overhead on tiny arcs: BlockAccumulator with
%         nl = 0 and 10-row arcs against the flat Accumulator on the same
%         rows. H-02's caution: kinds B/C become adapters into the block
%         accumulator only if this ratio is acceptable (<= 3 predicted).
%
%   INPUT (Name-Value)
%     'Quick'    logical  false   halve every size (smoke run)
%     'OutFile'  char     dev/validation/bench_m0_result.txt
%   OUTPUT  none; the record is appended to OutFile and echoed.
%
%   ERRORS  none of its own; a failing measurement is RECORDED in the file
%   with the error identifier and the bench continues (this is a bench in
%   dev/validation, not package code -- the catch binds and reports).
%
%   See also TIMEPAIR, TEST_SPEED, LSQVCE.ENGINE.PREFETCHLOADER,
%            LSQVCE.ENGINE.REDUNDANCY, LSQVCE.ENGINE.STREAMCONTRIBUTION
%
%   Created 2026-09-01 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M0, WP-79.

arguments
    opt.Quick (1,1) logical = false
    opt.OutFile (1,:) char = fullfile(fileparts(mfilename('fullpath')), 'bench_m0_result.txt')
end
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    pathBefore = path();  restore = onCleanup(@() path(pathBefore)); %#ok<NASGU>
    addpath(root);

    fid = fopen(opt.OutFile, 'a');
    closeIt = onCleanup(@() fclose(fid)); %#ok<NASGU>
    say(fid, '\n===== bench_m0  %s  %s %s  quick=%d =====\n', ...
        datestr(now, 'yyyy-mm-dd HH:MM'), computer(), version('-release'), opt.Quick); %#ok<TNOW1,DATST>
    q = 1 + opt.Quick;                                  % size divisor

    % ---------------------------------------------------------- (a) prefetch
    try
        cfg  = makeCfg(200 / q, 4, 1000 / q, 4);        % nc, arcs/day, mArc, nl
        days = 1:8;
        f    = binder(cfg, 1);
        plain = @() runBlock(cfg, days, f, {});
        pf    = @() runPrefetch(cfg, days, f);
        rec = pairTime(pf, plain, 5, 'prefetch/plain');
        say(fid, '(a) prefetch / plain serial      : ratio %.3f  band [%.3f %.3f]  tA %.2fs tB %.2fs\n', ...
            rec.ratio, rec.band, rec.tA, rec.tB);
        if hasPCT()
            par = @() runBlock(cfg, days, f, {'Parallel','parfor','DistributeDays',true,'ParallelThreshold',1});
            rec2 = pairTime(par, plain, 5, 'parfor-days/plain');
            say(fid, '(a) parfor DistributeDays / plain: ratio %.3f  band [%.3f %.3f]  (pool: %d workers)\n', ...
                rec2.ratio, rec2.band, poolSize());
        else
            say(fid, '(a) parfor DistributeDays / plain: SKIPPED (no PCT)\n');
        end
    catch ME
        say(fid, '(a) FAILED: %s -- %s\n', ME.identifier, strrep(ME.message, newline, ' '));
    end

    % ------------------------------------------------ (b) exact vs probe sweep
    for nc = [500 2000 4000] / q
        try
            cfg  = makeCfg(nc, 4, max(2*nc/4, 1000), 4);   % rows per arc >= nc/2, 4 arcs/day
            days = 1:4;
            f    = binder(cfg, 1);
            acc  = lsqvce.engine.BlockAccumulator(cfg.nc, 'Verbose', false);
            acc.addSystem(1, cfg.nl, f, days);
            sol  = acc.solve();
            exact = @() lsqvce.engine.redundancy(acc, sol, 'GroupBy', 'arc', 'Verbose', false);
            probe = @() lsqvce.engine.streamContribution(acc, sol, 'GroupBy', 'arc', 'Probes', 100, 'Seed', 1);
            rec = pairTime(exact, probe, 3, sprintf('exact/probe nc=%d', nc));
            say(fid, '(b) nc=%5d  exact sweep / probe sweep: ratio %.3f  band [%.3f %.3f]  exact %.2fs probe %.2fs\n', ...
                nc, rec.ratio, rec.band, rec.tA, rec.tB);
        catch ME
            say(fid, '(b) nc=%5d FAILED: %s -- %s\n', nc, ME.identifier, strrep(ME.message, newline, ' '));
        end
    end

    % ------------------------------------------ (c) block vs flat on tiny arcs
    try
        u = 50;  nArc = 2000 / q;  mArc = 10;
        cfg  = makeCfg(u, nArc, mArc, 0);              % nl = 0: a batch problem
        days = 1;
        f    = binder(cfg, 1);
        block = @() runBlock(cfg, days, f, {'Backend','normal'});
        flat  = @() runFlat(cfg, f, days, u);
        rec = pairTime(block, flat, 5, 'block/flat 10-row arcs');
        say(fid, '(c) block / flat, %d arcs x %d rows, u=%d: ratio %.3f  band [%.3f %.3f]  block %.2fs flat %.2fs\n', ...
            nArc, mArc, u, rec.ratio, rec.band, rec.tA, rec.tB);
    catch ME
        say(fid, '(c) FAILED: %s -- %s\n', ME.identifier, strrep(ME.message, newline, ' '));
    end
    say(fid, '===== end bench_m0 =====\n');
end

% ---------------------------------------------------------------- helpers
function say(fid, varargin)
    fprintf(fid, varargin{:});  fprintf(1, varargin{:});
end

function tf = hasPCT()
    tf = license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
end

function n = poolSize()
    p = gcp('nocreate');
    if isempty(p), n = 0; else, n = p.NumWorkers; end
end

function cfg = makeCfg(nc, arcsPerDay, mArc, nl)
    rs = RandStream('mt19937ar', 'Seed', 42);
    cfg.nc = nc;  cfg.nl = nl;  cfg.sats = 1;
    cfg.x_c = randn(rs, nc, 1);
    cfg.y = containers.Map('KeyType', 'double', 'ValueType', 'any');
    cfg.y(1) = randn(rs, nl, 1);                    % nl = 0 -> 0x1 (B is m x 0)
    cfg.obsTypes = {'kbr'};
    cfg.arcsPerDay = arcsPerDay;  cfg.mArc = mArc;  cfg.dt = 5;  cfg.seed = 42;
    cfg.collinear = false;  cfg.collinearEps = 1e-7;
    cfg.cov.kbr = struct('type', 'diag', 'sigma', 0.2);
end

function f = binder(cfg, s)
    f = @(d, ot, a, uo) lsqvce.engine.synthLoader(cfg, s, d, ot, a, uo);
end

function sol = runBlock(cfg, days, f, nv)
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Verbose', false, nv{:});
    acc.addSystem(1, cfg.nl, f, days);
    sol = acc.solve();
end

function sol = runPrefetch(cfg, days, f)
    pf  = lsqvce.engine.PrefetchLoader(f, days, 'Enable', 'auto');
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Verbose', false);
    acc.addSystem(1, cfg.nl, pf.handle(), days);
    sol = acc.solve();
    delete(pf);
end

function sol = runFlat(cfg, f, days, u)
% The flat path on the SAME rows: every (day, arc) becomes one batch.
    acc = lsqvce.engine.Accumulator(u, 'method', 'normal');
    uo  = struct('sysID', 1, 'nc', cfg.nc, 'nl', cfg.nl, 'BatchSize', 20000);
    for d = days
        meta = f('meta', d, [], uo);
        for a = meta.arcs{1}(:).'
            uo.obsType = meta.obsTypes{1};  uo.day = d;  uo.arc = a;
            arc = f(d, uo.obsType, a, uo);
            acc.add(arc.A, arc.l, []);
        end
    end
    sol = acc.solve();
end

function rec = pairTime(fA, fB, R, label)
% Local copy of dev/tests/private/timePair's logic (rotation, per-repeat
% ratio, median, band) -- see the header note; deleted with this script.
    fA();  fB();
    tA = zeros(R, 1);  tB = zeros(R, 1);
    for r = 1:R
        if mod(r, 2) == 1
            t0 = tic; fA(); tA(r) = toc(t0);  t0 = tic; fB(); tB(r) = toc(t0);
        else
            t0 = tic; fB(); tB(r) = toc(t0);  t0 = tic; fA(); tA(r) = toc(t0);
        end
    end
    ratios = tA ./ tB;
    rec = struct('ratio', median(ratios), 'band', [min(ratios) max(ratios)], ...
                 'tA', median(tA), 'tB', median(tB), 'label', label);
end
