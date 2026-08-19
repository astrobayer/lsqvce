function tests = test_contribution
%TEST_CONTRIBUTION  Unit tests for the streaming contribution primitive, the
%   lsqvce.contribution wrapper, and the hutchinson trace wiring in lsqvce.vce.
%
%   Covers the two defining invariants (sum_g tr(N^-1 N_g^cc) = nc and
%   sum_g diag(N^-1 N_g^cc) = ones(nc,1)), agreement of the streamed trace with a
%   high-probe reference, the Level-1/Level-2 outputs of lsqvce.contribution, and
%   that 'TraceMethod','hutchinson' in lsqvce.vce divides by the streamed
%   redundancy r_g = n_g - local_g - tr(N^-1 N_g^cc) (WP-32, fixes F-01)
%   while leaving vtv exact.
%
%   Fixture: the seeded synthetic multi-system / multi-obsType problem also used
%   by test_reliability (three systems, obsTypes pos/kbr), built through
%   lsqvce.synthLoader so the arcs are deterministic across passes.
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-40 (F-10): the freezeR
%   freeze-and-warn test. With assistance from Claude (Anthropic), model
%   Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-36: stream global-RNG
%   invariance + the vce traceSeedUsed replay test (final pass == loop
%   probes). With assistance from Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-18 (CEST, UTC+02:00), WP-32: hutchinson wiring prose
%   corrected + two new tests (streamed-redundancy parity/partition; the
%   F-01 regression that FAILS on the pre-fix code). With assistance from
%   Claude (Anthropic), model Claude Fable 5.
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
    % The §1.1 wiring estimates the REDUNDANCY r_g = n_g - local_g -
    % tr(N^-1 N_g^cc) matrix-free (WP-32; only the cc-trace is stochastic);
    % vtv stays exact. The load-bearing property is that the streamed r_g
    % approximates the exact r_g AT THE SAME WEIGHTS -- which then makes the
    % two loops track each other. Assert that directly (robust to how many
    % outer iterations this synthetic problem needs), plus an end-sigma check.
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

%% ===================================== WP-32: streamed redundancy (F-01 fix)
function test_streamed_redundancy_partitions_and_matches_exact(tc)
    % Spec tests (1)+(2) of WP-32 on the standing block-arrow fixture, with
    % obsType groups that CUT ACROSS every system (the general sub-system
    % split). (a) Parity: per group the streamed r must sit within 3x its own
    % Monte-Carlo se of the exact leverage-sum redundancy. (b) Partitions:
    % the EXACT pieces must partition exactly -- sum_g n_g = m (row count),
    % sum_g localAbsorbed_g = sum_s nl_s (block identity: whole systems),
    % sum_g trace_g ~ nc -- hence sum_g r_g ~ m - u.
    cfg = tc.TestData.cfg;
    w = warning('off','all'); cleaner = onCleanup(@() warning(w)); %#ok<NASGU>
    acc = buildAcc(cfg, []);
    sol = acc.solve();
    groups = {'pos','kbr'};
    sc = lsqvce.streamContribution(acc, sol, 'GroupBy','obsType', ...
             'Groups', groups, 'Probes', 2000, 'Seed', 21, 'WantDiag', false);
    red = lsqvce.redundancy(acc, sol, 'GroupBy','obsType', ...
             'Collect', false, 'Verbose', false);

    sk = cellstr(string({sc.group.key}));
    rk = cellstr(string({red.group.key}));
    for k = 1:numel(groups)
        [~, is] = ismember(char(string(groups{k})), sk);
        [~, ir] = ismember(char(string(groups{k})), rk);
        tol = 3*sc.group(is).rSE + 1e-6*sc.group(is).n;
        verifyEqual(tc, sc.group(is).r, red.group(ir).r, 'AbsTol', tol, ...
            sprintf('group %s: streamed r within 3x MC se of exact r', ...
                    groups{k}));
    end

    nSys = numel(cfg.sats);
    u    = cfg.nc + nSys*cfg.nl;
    m    = sum([sc.group.n]);
    verifyEqual(tc, sum([sc.group.localAbsorbed]), nSys*cfg.nl, ...
        'AbsTol', 1e-6, 'exact local dof must partition to sum_s nl_s');
    rTol = 3*sqrt(sum([sc.group.rSE].^2));
    verifyEqual(tc, sum([sc.group.r]), m - u, 'AbsTol', rTol, ...
        'streamed redundancies must sum to m - u');
end

function test_vce_hutchinson_regression_F01(tc)
    % Spec test (3) + the v10review.py T2 port: this test FAILS on the
    % pre-WP-32 code, where 'hutchinson' divided by tr(N^-1 N_g^cc). There
    % (T2a/b/c, measured): sum_g denominator ~ nc = 6 instead of m-u
    % (hundreds), and sigma inflated 8.4x / diverged to 1e23. On the fixed
    % code both loops estimate the SAME quantity, so they must both converge
    % and agree -- NO escape hatch.
    cfg = tc.TestData.cfg;
    accFcn = @(sig) buildAcc(cfg, sig);
    groups = {'pos','kbr'};
    w = warning('off','all'); cleaner = onCleanup(@() warning(w)); %#ok<NASGU>

    vcE = lsqvce.vce(accFcn, 'GroupBy','obsType', 'Groups',groups, ...
                     'MaxIter',20, 'Tol',1e-3, 'Verbose',false, ...
                     'TraceMethod','exact');
    vcH = lsqvce.vce(accFcn, 'GroupBy','obsType', 'Groups',groups, ...
                     'MaxIter',20, 'Tol',1e-3, 'Verbose',false, ...
                     'TraceMethod','hutchinson', 'Probes',1500, 'TraceSeed',9);

    verifyTrue(tc, vcE.converged, 'exact loop must converge on the fixture');
    verifyTrue(tc, vcH.converged, ...
        'hutchinson loop must converge (F-01: the old denominator diverged)');
    % T2a port: BOTH paths report the redundancy, so vc.r sums agree (~m-u);
    % the old code summed to ~nc = 6 here and fails this by orders of magnitude.
    verifyEqual(tc, sum(vcH.r)/sum(vcE.r), 1, 'AbsTol', 0.20, ...
        'hutchinson vc.r must sum to ~m-u like the exact path');
    % T2b/c port: the sigmas agree; in particular NO >5x inflation.
    rel = abs(vcH.sigma - vcE.sigma) ./ max(abs(vcE.sigma), 1e-6);
    verifyTrue(tc, all(rel < 0.20), ...
        sprintf('hutchinson sigma within 20%% of exact (max rel %.3f)', max(rel)));
    verifyTrue(tc, all(vcH.sigma ./ vcE.sigma < 5), ...
        'the F-01 inflation (>5x) must be gone');
end

%% ================================ WP-36: one probe generator, honest RNG (F-06)
function test_stream_default_seed_global_rng_untouched(tc)
    % streamContribution's half of the F-06 invariance lever (the other
    % consumers are covered in test_traceEstimate): an EMPTY Seed must leave
    % the global RNG state bit-identical. Fails on the pre-fix code (probes
    % came from the global stream while the comment claimed otherwise).
    cfg = tc.TestData.cfg;
    w = warning('off','all'); cleaner = onCleanup(@() warning(w)); %#ok<NASGU>
    acc = buildAcc(cfg, []);
    sol = acc.solve();
    rng(31337, 'twister');
    before = rng;
    lsqvce.streamContribution(acc, sol, 'GroupBy','obsType', ...
        'Groups', {'pos','kbr'}, 'Probes', 64, 'WantDiag', false);
    after = rng;
    verifyEqual(tc, after.State, before.State, ...
        'an empty Seed must never touch the global stream');
end

%% ============================ WP-39 stage (b): typed-signature enforcement
function test_stageB_signatures_contribution_family(tc)
    % WP-39 stage (b) (F-07): enum VALUES are the fails-on-old levers --
    % the a1-2 lesson says assert the SPECIFIC validator id
    % (MATLAB:validators:mustBeMember), which the old inputParser code
    % could not produce (it threw InputParser:ArgumentFailedValidation).
    dummy = @(s) s;
    ok = false;
    try
        lsqvce.vce(dummy, 'TraceMethod', 'Exact');
    catch e
        ok = strcmp(e.identifier, 'MATLAB:validators:mustBeMember');
    end
    verifyTrue(tc, ok, 'vce TraceMethod ''Exact'' must fail mustBeMember');
    ok = false;
    try
        lsqvce.streamContribution([], [], 'Probes', 0);
    catch e
        ok = startsWith(e.identifier, 'MATLAB:validators:');
    end
    verifyTrue(tc, ok, 'streamContribution Probes = 0 must fail the signature');
    ok = false;
    try
        lsqvce.contribution([], [], 'Level', 3);
    catch e
        ok = strcmp(e.identifier, 'MATLAB:validators:mustBeMember');
    end
    verifyTrue(tc, ok, 'contribution Level = 3 must fail mustBeMember');
    ok = false;
    try
        lsqvce.sharedProbe(eye(3), 'Distribution', 'Rademacher');
    catch e
        ok = strcmp(e.identifier, 'MATLAB:validators:mustBeMember');
    end
    verifyTrue(tc, ok, 'sharedProbe Distribution ''Rademacher'' must fail');
end

function test_vce_freezeR_freezes_and_warns(tc)
    % WP-40 (fixes F-10; fails on the old code, which had no 'freezeR' and
    % SILENTLY floored r_g to 1): with freezeR = Inf every group is below
    % the estimability threshold -- sigma must stay at its start BITWISE
    % (frozen, upd = 1) and the low-redundancy warning must fire once.
    cfg = tc.TestData.cfg;
    accFcn = @(sig) buildAcc(cfg, sig);
    w = warning('off','all'); cleaner = onCleanup(@() warning(w)); %#ok<NASGU>
    fired = false;
    ws = warning('query', 'lsqvce:vce:lowRedundancy'); %#ok<NASGU>
    lastwarn('', '');
    vc = lsqvce.vce(accFcn, 'GroupBy','obsType', 'Groups',{'pos','kbr'}, ...
                    'MaxIter', 3, 'Verbose', false, 'FreezeR', Inf);
    [~, wid] = lastwarn;
    fired = strcmp(wid, 'lsqvce:vce:lowRedundancy');
    verifyTrue(tc, fired, 'the low-redundancy warning must fire');
    verifyEqual(tc, vc.sigma, ones(size(vc.sigma)), ...
        'all groups frozen: sigma must stay at the start values bitwise');
end

function test_vce_hutchinson_final_reuses_loop_probes(tc)
    % WP-36 in vce: an EMPTY TraceSeed is pinned once up front, exposed as
    % vc.traceSeedUsed, and a replay with that seed reproduces vc.r and
    % sigma BIT-IDENTICALLY -- the "final pass probes with the loop's Z"
    % guarantee made observable. Fails on the pre-fix code (no
    % traceSeedUsed field; the final pass drew fresh global-stream probes).
    cfg = tc.TestData.cfg;
    accFcn = @(sig) buildAcc(cfg, sig);
    groups = {'pos','kbr'};
    w = warning('off','all'); cleaner = onCleanup(@() warning(w)); %#ok<NASGU>
    vc1 = lsqvce.vce(accFcn, 'GroupBy','obsType', 'Groups',groups, ...
                     'MaxIter',6, 'Tol',1e-3, 'Verbose',false, ...
                     'TraceMethod','hutchinson', 'Probes',300);
    verifyTrue(tc, isscalar(vc1.traceSeedUsed), ...
        'the pinned probe seed must be exposed');
    vc2 = lsqvce.vce(accFcn, 'GroupBy','obsType', 'Groups',groups, ...
                     'MaxIter',6, 'Tol',1e-3, 'Verbose',false, ...
                     'TraceMethod','hutchinson', 'Probes',300, ...
                     'TraceSeed', vc1.traceSeedUsed);
    verifyEqual(tc, vc2.r,     vc1.r,     'replay must reproduce vc.r bitwise');
    verifyEqual(tc, vc2.sigma, vc1.sigma, 'replay must reproduce sigma bitwise');
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
