function tests = test_lsqvce_combine
%TEST_LSQVCE_COMBINE  Unit tests for multi-dataset combination by elimination.
%
%   Run with:  results = runtests('tests/test_lsqvce_combine.m')
%
%   The two tests that matter most are the TRAP tests: they verify that the
%   toolbox does NOT commit the two classic elimination errors, both of which
%   produce plausible numbers that are wrong in the optimistic direction.
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.

tests = functiontests(localfunctions);
end

%% =================================================== fixtures
function S = arrowProblem(seed, uc, ul, nobs) %#ok<*DEFNU>
    % An arrow-structured problem: uc common parameters, ul(i) locals per
    % system, nobs(i) observations per system.
    rs = RandStream('mt19937ar', 'Seed', seed);
    nsys = numel(ul);
    utot = uc + sum(ul);

    S.uc = uc;  S.ul = ul;  S.nsys = nsys;  S.utot = utot;
    S.xt = randn(rs, utot, 1);
    S.n  = sum(nobs);

    Afull = zeros(S.n, utot);
    pfull = zeros(S.n, 1);
    r0 = 0;  off = uc;
    S.sys = cell(nsys,1);
    for i = 1:nsys
        ni = nobs(i);
        AC = randn(rs, ni, uc);
        AL = randn(rs, ni, ul(i));
        p  = 0.5 + rand(rs, ni, 1);
        S.sys{i} = struct('AC', AC, 'AL', AL, 'p', p, 'off', off, ...
                          'ul', ul(i), 'n', ni);
        Afull(r0+(1:ni), 1:uc) = AC;
        Afull(r0+(1:ni), off+(1:ul(i))) = AL;
        pfull(r0+(1:ni)) = p;
        r0 = r0 + ni;  off = off + ul(i);
    end
    l = Afull*S.xt + (1./sqrt(pfull)) .* randn(rs, S.n, 1);
    r0 = 0;
    for i = 1:nsys
        ni = S.sys{i}.n;
        S.sys{i}.l = l(r0+(1:ni));
        r0 = r0 + ni;
    end
    S.Afull = Afull;  S.pfull = pfull;  S.l = l;
end

function cmb = buildCombine(S, method)
    cmb = lsqvce.Combine(S.uc, 'method', method);
    for i = 1:S.nsys
        st = S.sys{i};
        % one batch per system for simplicity; the batching is tested elsewhere
        src = @(k) deal(st.AC, st.AL, st.l, st.p);
        cmb.addSystem(st.ul, src, 1, 'name', sprintf('sys%d', i));
    end
end

function [x, Qxx, v, r, vPv] = denseRef(S)
    P = diag(S.pfull);
    N = S.Afull' * P * S.Afull;
    Qxx = inv(N);
    x = Qxx * (S.Afull' * P * S.l);
    v = S.Afull*x - S.l;
    vPv = v' * P * v;
    Qvv = inv(P) - S.Afull*Qxx*S.Afull';
    r = diag(Qvv * P);
end

%% ================================================ CORRECTNESS
function test_matches_the_full_solution(tc)
% categories: precision | tags: diagnostics
    S = arrowProblem(1, 5, [4 6 3], [60 80 50]);
    [xr, Qxx, vr, rr, vPvr] = denseRef(S);

    for m = {'qr','normal'}
        cmb = buildCombine(S, m{1});
        sol = cmb.solve();
        D   = cmb.diagnose(sol);

        % tol: independent dense reference, same-precision algebra; eps*cond headroom (house bounds, HANDOVER 6.1.2)
        verifyEqual(tc, sol.xC, xr(1:S.uc), 'RelTol', 1e-8, ...
            sprintf('%s: common parameters wrong', m{1}));
        verifyEqual(tc, sol.QCC, Qxx(1:S.uc,1:S.uc), 'RelTol', 1e-7, ...
            sprintf('%s: Qxx_CC must be the reduced inverse', m{1}));
        % tol: independent dense solve of the same system; difference ~ cond*eps, house bounds add headroom (HANDOVER 6.1.2)
        verifyEqual(tc, D.vPv, vPvr, 'RelTol', 1e-8);
    end
end

function test_backsubstituted_locals_are_exact(tc)
% categories: precision
    S = arrowProblem(2, 4, [5 3], [70 60]);
    [xr, ~, ~, ~, ~] = denseRef(S);

    cmb = buildCombine(S, 'qr');
    sol = cmb.solve();

    off = S.uc;
    for i = 1:S.nsys
        ul = S.ul(i);
        % tol: independent dense reference, same-precision algebra; eps*cond headroom (house bounds, HANDOVER 6.1.2)
        verifyEqual(tc, sol.xL{i}, xr(off+(1:ul)), 'RelTol', 1e-8, ...
            sprintf('back-substituted locals of system %d are wrong', i));
        off = off + ul;
    end
end

function test_local_sigmas_include_the_inflation_term(tc)
% categories: precision
    % Q_ii = Nll^-1 + S QCC S'.  The SECOND TERM IS NOT OPTIONAL: Nll^-1 alone
    % is the covariance with the common parameters held FIXED, but xC was
    % ESTIMATED and its uncertainty propagates into the locals. Dropping it
    % makes the local sigmas too optimistic.
    S = arrowProblem(3, 5, [6 4], [80 70]);
    [~, Qxx, ~, ~, ~] = denseRef(S);

    cmb = buildCombine(S, 'qr');
    sol = cmb.solve();

    off = S.uc;
    for i = 1:S.nsys
        ul = S.ul(i);
        sig_ref = sqrt(diag(Qxx(off+(1:ul), off+(1:ul))));
        % tol: independent dense solve of the same system; difference ~ cond*eps, house bounds add headroom (HANDOVER 6.1.2)
        verifyEqual(tc, sol.sigL{i}, sig_ref, 'RelTol', 1e-7, ...
            sprintf(['local sigmas of system %d do not match the full Qxx -- ' ...
                     'the S*QCC*S'''' inflation term is probably missing'], i));
        off = off + ul;
    end
end

function test_residuals_match(tc)
% categories: precision | tags: diagnostics
    S = arrowProblem(4, 4, [5 5 5], [60 60 60]);
    [~, ~, vr, ~, ~] = denseRef(S);

    cmb = buildCombine(S, 'qr');
    sol = cmb.solve();
    D   = cmb.diagnose(sol);

    v = vertcat(D.v{:});
    % tol: independent dense reference, same-precision algebra; eps*cond headroom (house bounds, HANDOVER 6.1.2)
    verifyEqual(tc, v, vr, 'AbsTol', 1e-8);
end

%% ===================================== THE TWO TRAPS
function test_redundancy_numbers_include_the_cross_term(tc)
% categories: precision | tags: diagnostics
    % TRAP 1. a'Qxx a spans BOTH blocks:
    %     a'Qxx a = aC' Q_CC aC + 2 aC' Q_Ci aL + aL' Q_ii aL
    % Dropping the cross term (using only Q_CC, which is what the reduced
    % system hands you) makes a'Qxx a too small, hence r_i too LARGE.
    % Measured: sum(r) inflated by 7.3%, r_i biased +0.066 on average.
    S = arrowProblem(5, 5, [4 6 3], [60 80 50]);
    [~, ~, ~, rr, ~] = denseRef(S);

    cmb = buildCombine(S, 'qr');
    sol = cmb.solve();
    D   = cmb.diagnose(sol);

    r = vertcat(D.rNum{:});
    % tol: independent dense reference, same-precision algebra; eps*cond headroom (house bounds, HANDOVER 6.1.2)
    verifyEqual(tc, r, rr, 'AbsTol', 1e-8, ...
        ['redundancy numbers do not match the full solution -- the cross-term ' ...
         'in a''Qxx a has probably been dropped']);
end

function test_integrity_identity_uses_u_total(tc)
% categories: precision | tags: diagnostics
    % TRAP 2. The redundancy is n - u_TOTAL, counting the ELIMINATED locals.
    % They were removed from the SOLVE, not from the MODEL: they still cost
    % degrees of freedom. Using n - uC deflates s0^2 by ~7% here -- and in the
    % optimistic direction, which is the worst kind of wrong.
    S = arrowProblem(6, 5, [4 6 3], [60 80 50]);

    cmb = buildCombine(S, 'qr');
    sol = cmb.solve();
    D   = cmb.diagnose(sol);

    verifyEqual(tc, sol.r, S.n - S.utot, ...
        'redundancy must be n - u_total, NOT n - u_common');
    verifyNotEqual(tc, sol.r, S.n - S.uc, ...
        'redundancy must NOT be n - u_common');

    % and the identity must hold
    % tol: probe noise CANCELS in the group sum (identical probes); residual is fp accumulation over nc terms, 1e-8 leaves wide headroom
    verifyEqual(tc, D.sumR, D.r, 'RelTol', 1e-8, ...
        'sum(r_i) must equal n - u_total');
    verifyLessThan(tc, D.sumRRelErr, 1e-8);
end

%% ================================================ ERRORS
function test_singular_local_block_errors(tc)
% categories: contract, robustness
    % A system whose own data do not determine its own locals cannot be
    % eliminated. You cannot remove a parameter the data never constrained.
    rs = RandStream('mt19937ar','Seed', 11);
    uc = 3;  ul = 4;  n = 20;
    AC = randn(rs, n, uc);
    AL = randn(rs, n, ul);
    AL(:, 4) = AL(:, 1);              % exactly collinear -> N_ll singular
    l  = randn(rs, n, 1);

    cmb = lsqvce.Combine(uc);
    src = @(k) deal(AC, AL, l, ones(n,1));
    verifyError(tc, @() cmb.addSystem(ul, src, 1), ...
        'lsqvce:combine:singularLocal');
end

function test_common_block_undetermined_errors(tc)
% categories: contract, robustness
    % If NO system constrains a common direction, the reduced system is
    % singular. That is a parameterisation defect, not a numerical accident.
    rs = RandStream('mt19937ar','Seed', 13);
    uc = 4;  ul = 3;  n = 30;
    cmb = lsqvce.Combine(uc);
    for i = 1:2
        AC = randn(rs, n, uc);
        AC(:, 4) = 0;                 % NOBODY sees common parameter 4
        AL = randn(rs, n, ul);
        l  = randn(rs, n, 1);
        cmb.addSystem(ul, @(k) deal(AC, AL, l, ones(n,1)), 1);
    end
    verifyError(tc, @() cmb.solve(), 'lsqvce:combine:singularCommon');
end

function test_wrong_AL_width_errors(tc)
% categories: contract
    rs = RandStream('mt19937ar','Seed', 17);
    cmb = lsqvce.Combine(3);
    src = @(k) deal(randn(rs,20,3), randn(rs,20,2), randn(rs,20,1), ones(20,1));
    verifyError(tc, @() cmb.addSystem(5, src, 1), ...
        'lsqvce:combine:badAL');       % declared 5, returned 2
end

function test_correlated_batch_refused(tc)
% categories: contract, robustness | tags: diagnostics
    % Combine does not support matrix weights inside an eliminated system: the
    % generalised (P Qvv P)_ii with a'Qxx a spanning both blocks was not
    % validated, and shipping unverified algebra is exactly the failure mode
    % this toolbox exists to prevent. Refuse, and say how to work around it.
    rs = RandStream('mt19937ar','Seed', 19);
    n = 25;
    cmb = lsqvce.Combine(3);
    Pk = eye(n) + 0.3*diag(ones(n-1,1),1) + 0.3*diag(ones(n-1,1),-1);
    src = @(k) deal(randn(rs,n,3), randn(rs,n,2), randn(rs,n,1), Pk);
    cmb.addSystem(2, src, 1);
    sol = cmb.solve();
    verifyError(tc, @() cmb.diagnose(sol), ...
        'lsqvce:combine:correlatedBatch');
end

%% ================================================ STRUCTURE
function test_system_with_no_locals(tc)
% categories: precision, robustness
    % uL = 0 is legal: a data set that sees ONLY common parameters.
    S = arrowProblem(21, 4, [0 5], [50 60]);
    [xr, ~, ~, ~, ~] = denseRef(S);

    cmb = buildCombine(S, 'qr');
    sol = cmb.solve();
    % tol: independent dense reference, same-precision algebra; eps*cond headroom (house bounds, HANDOVER 6.1.2)
    verifyEqual(tc, sol.xC, xr(1:S.uc), 'RelTol', 1e-8);
    verifyEmpty(tc, sol.xL{1});
end

function test_multiple_batches_per_system(tc)
% categories: vectorisation
    % The result must not depend on how a system's data is chopped into
    % batches -- that is the whole premise of a streaming accumulator.
    S = arrowProblem(23, 4, [6 4], [90 80]);

    xs = [];
    for nb = [1 3 9]
        cmb = lsqvce.Combine(S.uc, 'method', 'qr');
        for i = 1:S.nsys
            st = S.sys{i};
            e = round(linspace(0, st.n, nb+1));
            src = @(k) deal(st.AC(e(k)+1:e(k+1), :), ...
                            st.AL(e(k)+1:e(k+1), :), ...
                            st.l(e(k)+1:e(k+1)), ...
                            st.p(e(k)+1:e(k+1)));
            cmb.addSystem(st.ul, src, nb);
        end
        xs = [xs, cmb.solve().xC]; %#ok<AGROW>
    end
    verifyLessThan(tc, max(max(abs(xs - xs(:,1)))), 1e-9, ...
        'the solution must not depend on the batch partition');
end

function test_qr_beats_normal_when_ill_conditioned(tc)
% categories: robustness
    % Forming N squares the condition number. Measured, at cond(A) = 2.4e8 the
    % normal-equation path has NO correct digits while Householder still has
    % eight. Verify the ORDERING here rather than a specific tolerance.
    muteWarn(tc, 'lsqvce:normalEq:negativeVtv');   % same regime as the precision sibling that flaked on CI: cond 2.4e8 on the normal path CAN round vPv negative (WP-47 r1e)
    rs = RandStream('mt19937ar','Seed', 31);
    uc = 5;  ul = 4;  n = 80;  nsys = 3;
    eps_ = 1e-7;

    Afull = [];  pfull = [];  lfull = [];
    sysd = cell(nsys,1);
    utot = uc + nsys*ul;
    xt = randn(rs, utot, 1);
    off = uc;
    for i = 1:nsys
        AC = randn(rs, n, uc);
        AC(:,2) = AC(:,1) + eps_*randn(rs, n, 1);   % near-collinear
        AL = randn(rs, n, ul);
        p  = ones(n,1);
        sysd{i} = struct('AC',AC,'AL',AL,'p',p,'ul',ul,'n',n,'off',off);
        Ab = zeros(n, utot);
        Ab(:,1:uc) = AC;  Ab(:, off+(1:ul)) = AL;
        Afull = [Afull; Ab]; %#ok<AGROW>
        pfull = [pfull; p];  %#ok<AGROW>
        off = off + ul;
    end
    lfull = Afull*xt + 0.01*randn(rs, size(Afull,1), 1);
    r0 = 0;
    for i = 1:nsys
        sysd{i}.l = lfull(r0+(1:n));  r0 = r0 + n;
    end

    % reference via SVD (does not form N)
    W = sqrt(pfull);
    [U,Sv,V] = svd(W.*Afull, 'econ');
    xr = V * ((U' * (W.*lfull)) ./ diag(Sv));

    e = zeros(2,1);
    ms = {'qr','normal'};
    for j = 1:2
        cmb = lsqvce.Combine(uc, 'method', ms{j});
        for i = 1:nsys
            st = sysd{i};
            cmb.addSystem(ul, @(k) deal(st.AC, st.AL, st.l, st.p), 1);
        end
        sol = cmb.solve();
        e(j) = norm(sol.xC - xr(1:uc)) / norm(xr(1:uc));
    end

    verifyLessThan(tc, e(1), e(2), ...
        sprintf(['QR (%.2e) should beat normal equations (%.2e) on an ' ...
                 'ill-conditioned problem'], e(1), e(2)));
end

function test_scaling_claim(tc)
% categories: contract
    % Not a numerical test -- a STRUCTURAL one. The running memory must be
    % O(uC^2 + max(uL)^2), never O(u_total^2).
    S = arrowProblem(41, 3, [4 4 4 4 4], [40 40 40 40 40]);
    cmb = buildCombine(S, 'qr');
    sol = cmb.solve();

    verifySize(tc, sol.QCC, [S.uc S.uc], ...
        'the reduced covariance must be uC x uC regardless of the local count');
    verifyEqual(tc, sol.uTotal, S.utot);
    verifyEqual(tc, sol.uL, sum(S.ul));
end

%% ============================== SENSOR x DAY HIERARCHY (BatchIndex)
function test_sensorday_batching_matches_flat(tc)
% categories: precision
    % THE CLAIM: with PER-SYSTEM locals, a (sensor, day) batch is just another
    % batch inside its system's elimination block. Verify the sensor x day
    % structure gives IDENTICAL results to one batch per system.
    rs = RandStream('mt19937ar','Seed', 501);
    uc = 5;  ul = [3 4];  nSens = [2 3];  nDay = 4;  nk = 20;

    utot = uc + sum(ul);
    xt = randn(rs, utot, 1);

    % build every (system, sensor, day) block, and the concatenation
    S = cell(2,1);
    off = uc;
    for i = 1:2
        blk = {};
        for s = 1:nSens(i)
            for d = 1:nDay
                AC = randn(rs, nk, uc);
                AL = randn(rs, nk, ul(i));
                sg = (0.3 + 0.4*s) * (0.5 + rand(rs, nk, 1));
                p  = 1 ./ sg.^2;
                l  = AC*xt(1:uc) + AL*xt(off+(1:ul(i))) + sg .* randn(rs,nk,1);
                blk{end+1} = struct('AC',AC,'AL',AL,'l',l,'p',p); %#ok<AGROW>
            end
        end
        S{i} = blk;
        off = off + ul(i);
    end

    % (a) sensor x day batches
    cmbA = lsqvce.Combine(uc, 'method', 'qr');
    for i = 1:2
        srcA = @(k) deal(S{i}{k}.AC, S{i}{k}.AL, S{i}{k}.l, S{i}{k}.p);
        cmbA.addSystem(ul(i), srcA, numel(S{i}));
    end
    solA = cmbA.solve();

    % (b) one batch per system (everything concatenated)
    cmbB = lsqvce.Combine(uc, 'method', 'qr');
    for i = 1:2
        AC = cell2mat(cellfun(@(b) b.AC, S{i}(:), 'uni', 0));
        AL = cell2mat(cellfun(@(b) b.AL, S{i}(:), 'uni', 0));
        l  = cell2mat(cellfun(@(b) b.l,  S{i}(:), 'uni', 0));
        p  = cell2mat(cellfun(@(b) b.p,  S{i}(:), 'uni', 0));
        cmbB.addSystem(ul(i), @(k) deal(AC, AL, l, p), 1);
    end
    solB = cmbB.solve();

    % tol: invariance claim (same solve, different bookkeeping); eps-level accumulation-order difference only
    verifyEqual(tc, solA.xC, solB.xC, 'RelTol', 1e-9, ...
        'sensor x day batching must give the same answer as one batch');
    for i = 1:2
        % tol: invariance claim (same solve, different bookkeeping); eps-level accumulation-order difference only
        verifyEqual(tc, solA.xL{i}, solB.xL{i}, 'RelTol', 1e-8);
    end
end

function test_batch_order_does_not_matter(tc)
% categories: vectorisation
    % QR accumulation is order-independent. You may feed days and sensors in
    % any order and interleave them freely. Verified to 3e-16 in Python.
    rs = RandStream('mt19937ar','Seed', 511);
    uc = 4;  ul = 3;  nB = 8;  nk = 25;
    blk = cell(nB,1);
    for k = 1:nB
        blk{k} = struct('AC', randn(rs,nk,uc), 'AL', randn(rs,nk,ul), ...
                        'l', randn(rs,nk,1),   'p', ones(nk,1));
    end

    xs = [];
    for trial = 1:3
        ord = randperm(RandStream('mt19937ar','Seed', 600+trial), nB);
        cmb = lsqvce.Combine(uc, 'method', 'qr');
        b = blk(ord);
        cmb.addSystem(ul, @(k) deal(b{k}.AC, b{k}.AL, b{k}.l, b{k}.p), nB);
        xs = [xs, cmb.solve().xC]; %#ok<AGROW>
    end
    verifyLessThan(tc, max(max(abs(xs - xs(:,1)))), 1e-10, ...
        'the solution must not depend on the order of the batches');
end

%% ===================================================== BatchIndex
function test_batchindex_decode(tc)
% categories: contract
    bx = lsqvce.BatchIndex();
    bx.addSystem('GRACE-C', {'KBR','LRI'}, 1:3);
    for s = {'KBR','LRI'}
        for d = 1:3
            bx.setCount('GRACE-C', s{1}, d, 100);
        end
    end
    bx.lock();

    verifyEqual(tc, bx.nBatch(1), 6);
    verifyEqual(tc, bx.nSystem(), 1);

    seen = cell(6,1);
    for k = 1:6
        info = bx.decode(1, k);
        verifyEqual(tc, info.system, 'GRACE-C');
        verifyEqual(tc, info.n, 100);
        seen{k} = sprintf('%s-%d', info.sensor, info.day);
    end
    verifyEqual(tc, numel(unique(seen)), 6, ...
        'every (sensor, day) must map to a distinct batch');
end

function test_batchindex_missing_sensor_day(tc)
% categories: robustness
    % A sensor that was down on a given day is simply not registered. No empty
    % batches, and no requirement that every sensor cover every day.
    bx = lsqvce.BatchIndex();
    bx.addSystem('S', {'A','B'}, 1:4);
    for d = 1:4, bx.setCount('S', 'A', d, 50); end
    for d = [1 2 4], bx.setCount('S', 'B', d, 50); end     % B down on day 3
    bx.lock();
    verifyEqual(tc, bx.nBatch(1), 7);      % 4 + 3, not 8
end

function test_batchindex_rejects_double_registration(tc)
% categories: contract
    % Registering the same (sensor, day) twice would DOUBLE-COUNT those
    % observations -- in the solve and in the diagnostics.
    bx = lsqvce.BatchIndex();
    bx.addSystem('S', {'A'}, 1:2);
    bx.setCount('S', 'A', 1, 50);
    verifyError(tc, @() bx.setCount('S', 'A', 1, 50), ...
        'lsqvce:batchIndex:duplicateBatch');
end

function test_batchindex_requires_lock(tc)
% categories: contract
    bx = lsqvce.BatchIndex();
    bx.addSystem('S', {'A'}, 1);
    bx.setCount('S', 'A', 1, 10);
    verifyError(tc, @() bx.decode(1,1), 'lsqvce:batchIndex:notLocked');
end

function test_batchindex_unknown_sensor(tc)
% categories: contract
    bx = lsqvce.BatchIndex();
    bx.addSystem('S', {'KBR'}, 1);
    verifyError(tc, @() bx.setCount('S', 'LRI', 1, 10), ...
        'lsqvce:batchIndex:noSuchSensor');
end

function test_sensorgroups_length(tc)
% categories: contract
    bx = lsqvce.BatchIndex();
    bx.addSystem('S', {'A','B'}, 1:2);
    bx.setCount('S','A',1,30);  bx.setCount('S','A',2,40);
    bx.setCount('S','B',1,50);  bx.setCount('S','B',2,60);
    bx.lock();

    g = bx.sensorGroups(1);
    verifyEqual(tc, numel(g), 4);
    for k = 1:4
        info = bx.decode(1,k);
        verifyEqual(tc, numel(g{k}), info.n, ...
            'the group vector must have one entry per observation');
        verifyEqual(tc, numel(unique(g{k})), 1, ...
            'all observations of a batch belong to the SAME sensor');
    end
end

%% ============================ PER-SENSOR VARIANCE COMPONENTS
function test_varcomp_recovers_sensor_noise(tc)
% categories: precision
    % Two sensors observe the same arcs. Sensor 2 is truly 3x noisier, but we
    % CLAIM sigma = 1 for both. Helmert must find the factor.
    %
    % A SINGLE realisation is not a test: with r_g of a few hundred the
    % sampling sd of a sigma RATIO is ~6%. Verified unbiased over 30
    % realisations (3.0 -> 3.032). Here we allow 3 sampling sd.
    rs = RandStream('mt19937ar','Seed', 521);
    uc = 5;  ul = 3;  nSens = 2;  nBatchPer = 8;  nk = 40;
    sigTrue = [1.0, 3.0];

    utot = uc + ul;
    xt = randn(rs, utot, 1);
    blk = cell(nSens*nBatchPer, 1);
    grp = cell(nSens*nBatchPer, 1);
    q = 0;
    for s = 1:nSens
        for b = 1:nBatchPer
            q = q + 1;
            AC = randn(rs, nk, uc);
            AL = randn(rs, nk, ul);
            sg = sigTrue(s) * ones(nk,1);
            l  = AC*xt(1:uc) + AL*xt(uc+1:end) + sg .* randn(rs,nk,1);
            blk{q} = struct('AC',AC,'AL',AL,'l',l,'p',ones(nk,1));
            grp{q} = repmat(s, nk, 1);
        end
    end

    cmb = lsqvce.Combine(uc, 'method', 'normal');
    cmb.addSystem(ul, @(k) deal(blk{k}.AC, blk{k}.AL, blk{k}.l, blk{k}.p), ...
                  numel(blk));
    sol = cmb.solve();

    VC = cmb.varComp(sol, @(i,k) grp{k}, 2, 'verbose', false);
    ratio = VC.sigma(2) / VC.sigma(1);

    verifyTrue(tc, VC.converged);
    % tol: recovery of a planted value; estimator sampling scatter on this fixture, measured stable across seeds (not an eps claim)
    verifyEqual(tc, ratio, 3.0, 'RelTol', 0.20, ...
        sprintf('recovered %.3f, truth 3.0', ratio));
end

function test_varcomp_no_false_correction(tc)
% categories: precision
    % With CORRECT weights the estimator must not invent a correction beyond
    % its own resolution limit.
    rs = RandStream('mt19937ar','Seed', 531);
    uc = 4;  ul = 3;  nB = 16;  nk = 40;
    utot = uc + ul;
    xt = randn(rs, utot, 1);
    blk = cell(nB,1);  grp = cell(nB,1);
    for k = 1:nB
        AC = randn(rs, nk, uc);  AL = randn(rs, nk, ul);
        l  = AC*xt(1:uc) + AL*xt(uc+1:end) + randn(rs, nk, 1);
        blk{k} = struct('AC',AC,'AL',AL,'l',l,'p',ones(nk,1));
        grp{k} = repmat(1 + mod(k,2), nk, 1);       % two groups, same noise
    end

    cmb = lsqvce.Combine(uc, 'method', 'normal');
    cmb.addSystem(ul, @(k) deal(blk{k}.AC, blk{k}.AL, blk{k}.l, blk{k}.p), nB);
    sol = cmb.solve();
    VC = cmb.varComp(sol, @(i,k) grp{k}, 2, 'verbose', false);

    ratio = VC.sigma(2)/VC.sigma(1);
    tol = 3 * sqrt(sum(VC.resolution.^2));
    % tol: 3-sigma of the estimator's own reported resolution/SE -- self-scaling statistical bar
    verifyEqual(tc, ratio, 1, 'AbsTol', max(tol, 0.10), ...
        'varComp invented a weight correction on correct weights');
end

function test_varcomp_reports_resolution(tc)
% categories: precision, robustness
    rs = RandStream('mt19937ar','Seed', 541);
    uc = 4;  ul = 2;  nB = 10;  nk = 30;
    xt = randn(rs, uc+ul, 1);
    blk = cell(nB,1);
    for k = 1:nB
        AC = randn(rs,nk,uc);  AL = randn(rs,nk,ul);
        l = AC*xt(1:uc) + AL*xt(uc+1:end) + randn(rs,nk,1);
        blk{k} = struct('AC',AC,'AL',AL,'l',l,'p',ones(nk,1));
    end
    cmb = lsqvce.Combine(uc, 'method', 'normal');
    cmb.addSystem(ul, @(k) deal(blk{k}.AC, blk{k}.AL, blk{k}.l, blk{k}.p), nB);
    sol = cmb.solve();
    VC = cmb.varComp(sol, @(i,k) ones(nk,1), 1, 'verbose', false);

    verifyTrue(tc, all(isfinite(VC.resolution)));
    % tol: closed-form resolution 0.5*sqrt(2/r) recomputed from the SAME r; house exact-identity bound
    verifyEqual(tc, VC.resolution(1), 0.5*sqrt(2/VC.redund(1)), 'RelTol', 1e-9);
end

%% =================================================== ReducedContribution
function test_reduced_contribution_contract(tc)
% categories: contract
% WP-43 (audit check 16): ReducedContribution was the ONE public unit with no
% direct test -- both documented rejection ids were unexercised. Asserts the
% help-text promises: constructor rejections with the documented ids (and only
% those), column-vector normalisation of n, and immutability of every property.
    N = eye(3);  n = (1:3).';  ltl = 2.5;  dof = 7;

    % documented rejections, each with its id
    tc.verifyError(@() lsqvce.ReducedContribution(ones(3,2), n, ltl, dof), ...
                   'lsqvce:reduced:notSquare');
    tc.verifyError(@() lsqvce.ReducedContribution(N, (1:4).', ltl, dof), ...
                   'lsqvce:reduced:sizeMismatch');

    % shapes and the row-input normalisation promise (n stored as a column)
    red = lsqvce.ReducedContribution(N, 1:3, ltl, dof);
    tc.verifyEqual(red.N,   N);
    tc.verifyEqual(red.n,   n);                 % 1x3 in -> 3x1 stored
    tc.verifyEqual(red.ltl, ltl);
    tc.verifyEqual(red.dof, dof);

    % immutable by design: every property refuses assignment
    for p = {'N', 'n', 'ltl', 'dof'}
        ok = false;
        try
            red.(p{1}) = 0;
        catch e
            ok = strcmp(e.identifier, 'MATLAB:class:SetProhibited');
        end
        tc.verifyTrue(ok, sprintf('property %s must be immutable', p{1}));
    end
end
