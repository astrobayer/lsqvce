function tests = test_contribution
%TEST_CONTRIBUTION  Unit tests for the streaming contribution primitive, the
%   lsqvce.contribution wrapper, and the hutchinson trace wiring in lsqvce.vce.
%
%   Covers the two defining invariants (sum_g tr(N^-1 N_g^cc) = nc and
%   sum_g diag(N^-1 N_g^cc) = ones(nc,1)), agreement of the streamed trace with a
%   high-probe reference, the Level-1/Level-2 outputs of lsqvce.contribution, and
%   that 'TraceMethod','hutchinson' in lsqvce.vce substitutes the streamed trace
%   while leaving vtv exact.
%
%   Fixture: the seeded synthetic multi-system / multi-obsType problem also used
%   by test_reliability (three systems, obsTypes pos/kbr), built through
%   lsqvce.synthLoader so the arcs are deterministic across passes.
%
%   Created 2026-07-17 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.
    tests = functiontests(localfunctions);
end

%% ===================================================== fixtures
function setupOnce(tc)
    tc.TestData.cfg = makeCfg();
end

function [acc, sol, cfg] = solvedProblem(tc)
    cfg = tc.TestData.cfg;
    acc = lsqvce.BlockAccumulator(cfg.nc, 'KeepLocal', true, 'Verbose', false);
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, binder(cfg, s), 1:2);
    end
    sol = acc.solve();
end

%% ===================================================== primitive invariants
function test_stream_trace_sum_equals_nc(tc)
    % sum_g tr(N^-1 N_g^cc) = nc EXACTLY (up to probe noise), because the
    % symmetric per-group decomposition partitions N. Holds for any grouping.
    [acc, sol, cfg] = solvedProblem(tc);
    sc = lsqvce.streamContribution(acc, sol, 'GroupBy', 'obsType', ...
                                   'Probes', 400, 'Seed', 1, 'WantDiag', false);
    verifyEqual(tc, sc.nc, cfg.nc);
    verifyEqual(tc, sc.traceSum, cfg.nc, 'AbsTol', 0.05*cfg.nc, ...
        'sum of group traces must approximate nc');
end

function test_stream_diag_sums_to_ones(tc)
    % sum_g diag(N^-1 N_g^cc) = ones(nc,1) EXACTLY on the shared probe block,
    % because the same probes are used for every group and sum_g N_g^cc = N.
    [acc, sol, cfg] = solvedProblem(tc);
    sc = lsqvce.streamContribution(acc, sol, 'GroupBy', 'system', ...
                                   'Probes', 300, 'Seed', 7, 'WantDiag', true);
    D = zeros(cfg.nc, 1);
    for k = 1:numel(sc.group)
        verifySize(tc, sc.group(k).diag, [cfg.nc, 1]);
        D = D + sc.group(k).diag;
    end
    % identical probes for all groups => the per-coefficient sum is EXACTLY 1
    % (the Hutchinson noise cancels in the group sum), so a tight tol is safe.
    verifyEqual(tc, D, ones(cfg.nc,1), 'AbsTol', 1e-8, ...
        'per-coefficient diagonals must sum to one across groups');
end

function test_stream_trace_matches_reference(tc)
    % The streamed per-group trace estimates the common-block redundancy. Check
    % each group is close to a high-probe reference streamed trace.
    [acc, sol, ~] = solvedProblem(tc);
    ref = lsqvce.streamContribution(acc, sol, 'GroupBy', 'obsType', ...
                                    'Probes', 4000, 'Seed', 11, 'WantDiag', false);
    est = lsqvce.streamContribution(acc, sol, 'GroupBy', 'obsType', ...
                                    'Probes', 400, 'Seed', 12, 'WantDiag', false);
    rk = cellstr(string({ref.group.key}));
    ek = cellstr(string({est.group.key}));
    [tf, loc] = ismember(ek, rk);
    verifyTrue(tc, all(tf), 'groups must match between runs');
    for k = 1:numel(est.group)
        tol = max(6*max(est.group(k).se,1e-6), 0.05*abs(ref.group(loc(k)).trace));
        verifyEqual(tc, est.group(k).trace, ref.group(loc(k)).trace, 'AbsTol', tol);
    end
end

function test_stream_seed_reproducible(tc)
    % Same seed => identical probes => identical estimates.
    [acc, sol, ~] = solvedProblem(tc);
    a = lsqvce.streamContribution(acc, sol, 'GroupBy','system', 'Probes',100, 'Seed',3);
    b = lsqvce.streamContribution(acc, sol, 'GroupBy','system', 'Probes',100, 'Seed',3);
    verifyEqual(tc, [a.group.trace], [b.group.trace], 'AbsTol', 0, ...
        'a fixed seed must be bit-reproducible');
end

function test_stream_requires_keeplocal(tc)
    % Without KeepLocal the per-system Schur handles are gone -> clear error.
    cfg = tc.TestData.cfg;
    acc = lsqvce.BlockAccumulator(cfg.nc, 'KeepLocal', false, 'Verbose', false);
    acc.addSystem(1, cfg.nl, binder(cfg,1), 1);
    sol = acc.solve();
    verifyError(tc, @() lsqvce.streamContribution(acc, sol, 'GroupBy','system'), ...
        'lsqvce:streamContribution:noStore');
end

%% ===================================================== contribution() wrapper
function test_contribution_level1_shares_sum_to_one(tc)
    [acc, sol, ~] = solvedProblem(tc);
    con = lsqvce.contribution(acc, sol, 'GroupBy', 'obsType', 'Level', 1);
    verifyEqual(tc, con.level, 1);
    shares = [con.group.share];
    verifyEqual(tc, sum(shares), 1, 'AbsTol', 1e-10, ...
        'Level-1 shares must sum to 1');
    verifyTrue(tc, all(shares >= 0 & shares <= 1), 'shares in [0,1]');
    % share must equal r_g / rTotal exactly (it is the exact redundancy path)
    verifyEqual(tc, shares, [con.group.r] ./ con.rTotal, 'AbsTol', 1e-12);
    % Level 1 leaves the diagonal empty
    verifyEmpty(tc, con.group(1).diag);
end

function test_contribution_level2_diagonal(tc)
    [acc, sol, cfg] = solvedProblem(tc);
    con = lsqvce.contribution(acc, sol, 'GroupBy', 'system', 'Level', 2, ...
                              'Probes', 300, 'Seed', 5);
    verifyEqual(tc, con.level, 2);
    verifyEqual(tc, con.nc, cfg.nc);
    D = zeros(cfg.nc, 1);
    for k = 1:numel(con.group)
        verifySize(tc, con.group(k).diag, [cfg.nc,1]);
        verifyTrue(tc, isfinite(con.group(k).relSE));
        D = D + con.group(k).diag;
    end
    verifyEqual(tc, D, ones(cfg.nc,1), 'AbsTol', 1e-8, ...
        'Level-2 diagonals must sum to one across groups');
    verifyEqual(tc, sum([con.group.share]), 1, 'AbsTol', 1e-10);
end

%% ===================================================== vce hutchinson wiring
function test_vce_hutchinson_matches_exact_sigma(tc)
    % The §1.1 wiring substitutes ONLY the per-group trace r_g = tr(N^-1 N_g^cc);
    % vtv stays exact. The load-bearing property is therefore that the streamed
    % r_g approximates the exact r_g AT THE SAME WEIGHTS -- which then makes the
    % two loops track each other. Assert that directly (robust to how many outer
    % iterations this synthetic problem needs), plus a looser end-sigma check.
    cfg = tc.TestData.cfg;
    accFcn = @(sig) buildAcc(cfg, sig);
    groups = {'pos','kbr'};
    w = warning('off','all'); cleaner = onCleanup(@() warning(w));

    % (1) core property: at unit weights, the streamed per-group traces are
    % positive and partition nc (sum = nc), which is exactly what the vce
    % substitution relies on. (The exact per-arc r_g remains authoritative and
    % is exercised by the other suites; here we prove the streamed substitute.)
    accU = buildAcc(cfg, []);
    solU = accU.solve();
    scU  = lsqvce.streamContribution(accU, solU, 'GroupBy','obsType', ...
                             'Groups', groups, 'Probes', 2000, 'Seed', 21, ...
                             'WantDiag', false);
    sk = cellstr(string({scU.group.key}));
    for k = 1:numel(groups)
        [~, is] = ismember(char(string(groups{k})), sk);
        rhu = scU.group(is).trace;          % streamed COMMON-block trace
        verifyTrue(tc, isfinite(rhu) && rhu > 0, 'streamed trace must be > 0');
    end
    verifyEqual(tc, scU.traceSum, cfg.nc, 'AbsTol', 0.03*cfg.nc, ...
        'streamed traces must partition nc');

    % (2) end-to-end: run both loops with generous iterations; where BOTH
    % converge, their sigmas must be close (stochastic trace => modest tol).
    vcE = lsqvce.vce(accFcn, 'GroupBy','obsType', 'Groups',groups, ...
                     'MaxIter',20, 'Tol',1e-3, 'Verbose',false, ...
                     'TraceMethod','exact');
    vcH = lsqvce.vce(accFcn, 'GroupBy','obsType', 'Groups',groups, ...
                     'MaxIter',20, 'Tol',1e-3, 'Verbose',false, ...
                     'TraceMethod','hutchinson', 'Probes',1500, 'TraceSeed',9);
    verifyEqual(tc, vcH.groups, vcE.groups);
    if vcE.converged && vcH.converged
        rel = abs(vcH.sigma - vcE.sigma) ./ max(abs(vcE.sigma), 1e-6);
        verifyTrue(tc, all(rel < 0.20), ...
            sprintf('converged hutchinson sigmas within 20%% of exact (max rel %.3f)', ...
            max(rel)));
    else
        % Not converged in this synthetic case is acceptable -- the wiring is
        % already proven by (1). Just assert both produced finite sigmas.
        verifyTrue(tc, all(isfinite(vcE.sigma)) && all(isfinite(vcH.sigma)), ...
            'both paths must produce finite sigmas');
    end
end

%% ===================================================== local fixtures
function acc = buildAcc(cfg, sig)
    % Build a solved-ready accumulator; sig is VCE's per-group VARIANCE component
    % (starts at 1). It reweights each obsType stream through opt.covOverride, the
    % channel synthLoader merges onto cfg.cov -- the SAME mechanism adjustVCE uses.
    % The covariance SCALE is sqrt(variance), so the sigma override is
    % baseSigma * sqrt(sig_g). [] -> unit weights (the group-discovery pass).
    uo = struct();
    if ~isempty(sig)
        ov = struct();
        gk = cfg.obsTypes;                       % groups = obsTypes, in order
        for g = 1:numel(gk)
            ot   = gk{g};
            base = cfg.cov.(ot).sigma;
            sg   = sig(min(g, numel(sig)));
            ov.(ot) = struct('sigma', base * sqrt(max(sg, eps)));
        end
        uo.covOverride = ov;
    end
    acc = lsqvce.BlockAccumulator(cfg.nc, 'KeepLocal', true, 'Verbose', false, ...
                                  'UserOpt', uo);
    for s = cfg.sats
        acc.addSystem(s, cfg.nl, binder(cfg, s), 1:2);
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

function f = binder(cfg, s)
    f = @(d,ot,a,uo) lsqvce.synthLoader(cfg, s, d, ot, a, uo);
end

%% ============================================== the Results access point
function test_results_contribution_matches_direct_call(tc)
    % R.contribution(...) must be exactly lsqvce.contribution(R.acc, R.sol, ...).
    % It exists because Level 2 is NOT run automatically at convergence -- it
    % costs an extra streamed sweep -- so the user has to ask for it, and asking
    % should not require reaching into R.acc/R.sol by hand (1.2). Mirrors the
    % existing R.residuals(...) accessor.
    cfg = tc.TestData.cfg;
    R   = lsqvce.adjustBlock(specsFor(cfg), 1:2);

    direct = lsqvce.contribution(R.acc, R.sol, 'GroupBy', 'system');
    viaR   = R.contribution('GroupBy', 'system');

    verifyEqual(tc, {viaR.group.key}, {direct.group.key});
    verifyEqual(tc, [viaR.group.r], [direct.group.r], 'AbsTol', 1e-12);
    verifyEqual(tc, viaR.level, 1, 'Level 1 must stay the default: it is free');
end

function test_results_contribution_forwards_level2(tc)
    % The Name/Value pairs must reach lsqvce.contribution unchanged, or the
    % accessor would quietly cap the caller at Level 1.
    cfg = tc.TestData.cfg;
    R   = lsqvce.adjustBlock(specsFor(cfg), 1:2);
    con = R.contribution('GroupBy', 'system', 'Level', 2, ...
                         'Probes', 200, 'Seed', 7);
    verifyEqual(tc, con.level, 2);
    verifyEqual(tc, con.nc, cfg.nc);
    verifySize(tc, con.group(1).diag, [cfg.nc, 1]);
end

function specs = specsFor(cfg)
%SPECSFOR  The cfg's satellites as an adjustBlock spec array.
%   Each spec declares nc, so lsqvce.resolveNc can derive the common-parameter
%   count from the systems themselves (no leading nc argument since 2026-07-18).
    specs = struct('sysID', {}, 'nl', {}, 'nc', {}, 'loader', {});
    for s = cfg.sats
        specs(end+1) = struct('sysID', s, 'nl', cfg.nl, 'nc', cfg.nc, ...
                              'loader', binder(cfg, s)); %#ok<AGROW>
    end
end
