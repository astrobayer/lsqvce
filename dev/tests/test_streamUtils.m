function tests = test_streamUtils
%TEST_STREAMUTILS  Direct tests for the functions moved into +lsqvce this session:
%
%     lsqvce.kernel.whitenessCorr          whiteness statistic
%     lsqvce.engine.collectResiduals    stack a residual field over arcs
%     lsqvce.engine.makeLoader  / lsqvce.engine.covModel  the error information of a day (R6 WP-90)
%
%   Reference properties/tolerances established in validation/v9verify.py (14).
%
%   Run:
%     results = runtests('tests/test_streamUtils');  disp(table(results))
%
%   See also LSQVCE.ENGINE.MAKELOADER, LSQVCE.ENGINE.COVMODEL, LSQVCE.ENGINE.COLLECTRESIDUALS
%
%   Updated 2026-10-17 (CEST, UTC+02:00), WP-148: a day's excluded epoch is
%   deleted from a correlated law (normals against A_o' C_oo^-1 A_o). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-138 (code review of WP-135):
%   lsqvce.util.priorStream makes and recognises the prior pseudo-streams.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-133: a zero epoch weight is
%   excluded and whitens to ~0 in place (the nugget rotation). With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    pathBefore = path;   % R7-45: a standalone run (run_lsqvce / runtests) leaves the path as found
    tc.addTeardown(@() path(pathBefore));
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(root);  addpath(fullfile(root, 'examples'));   % fixture_grace (moved tests, R5 M4 4c)
end

%% ===================================================== whitenessCorr
function testMaxLagCorrWhiteVsCorrelated(tc)
% categories: contract
    rs = RandStream('mt19937ar','Seed',9);
    w  = randn(rs, 4000, 1);
    tc.verifyLessThan(lsqvce.kernel.whitenessCorr(w, 5), 3/sqrt(4000));
    ar1 = filter(1, [1 -0.8], randn(rs, 4000, 1));
    mc  = lsqvce.kernel.whitenessCorr(ar1, 5);
    tc.verifyGreaterThan(mc, 3/sqrt(4000));
    tc.verifyGreaterThanOrEqual(mc, 0);
    tc.verifyLessThanOrEqual(mc, 1);
end

function testMaxLagCorrShortSeries(tc)
% categories: contract
    tc.verifyEqual(lsqvce.kernel.whitenessCorr(7, 5), 0);     % 1 element -> 0
end

%% ===================================================== collectResiduals
function testCollectResidualsFilters(tc)
% categories: robustness
    res = makeRes();
    wC = lsqvce.engine.collectResiduals(res, 'C', 'pos', 'w');
    tc.verifyEqual(numel(wC), 18);                  % 2 arcs x 3 epochs x 3 axes (C/pos)
    vAll = lsqvce.engine.collectResiduals(res, [], 'pos', 'v');
    tc.verifyEqual(numel(vAll), 36);                % C+D pooled (4 arcs x 9)
    empty = lsqvce.engine.collectResiduals(res, 'C', 'nope');
    tc.verifyEmpty(empty);
    wDef = lsqvce.engine.collectResiduals(res, 'D', 'pos'); % default field 'w'
    tc.verifyEqual(numel(wDef), 18);
end

function testCollectResidualsRejectsMissingField(tc)
% categories: contract
% WP-43: the documented rejection lsqvce:collectResiduals:badRes was never
% exercised directly (audit check 16 found the gap). A residual structure
% whose matching arcs lack the requested field must fire that id -- and only
% that id -- naming the offending field.
    res = makeRes();
    % try/catch capture (WP-39 pattern): verifyError's output-argument form
    % errored on the function-based TestCase in CI (R2026a hosted runner) --
    % round-1 red, both new err-capturing tests. Captured lesson: assert the
    % id AND the message through the caught exception, never via the output.
    caught = MException.empty;
    try
        lsqvce.engine.collectResiduals(res, 'C', 'pos', 'zz');
    catch e
        caught = e;
    end
    tc.verifyNotEmpty(caught, 'the absent field must be rejected');
    tc.verifyEqual(caught.identifier, 'lsqvce:collectResiduals:badRes');
    tc.verifySubstring(caught.message, 'zz');   % the message names the field
end

function testStreamAndArcKeysAreDefinedOnce(tc)
% categories: contract
    % R7-43 (WP-115 part 4): the stream key and the arc key each have ONE
    % definition (lsqvce.util.streamKey / arcKey); they were composed inline
    % in eight places and parsed back by a regex. The formats are the ones
    % out.arcKeys / R.noise.stream always had.
    tc.verifyEqual(lsqvce.util.streamKey(1, 'kbr'), 'x1|kbr');
    tc.verifyEqual(lsqvce.util.streamKey("C", "pos"), 'C|pos');
    [k, s] = lsqvce.util.streamKey(7);
    tc.verifyEqual({k, s}, {'x7', 'x7'}, 'without obsType the key is the satKey');
    tc.verifyTrue(isvarname(s), 'satKey is always a legal field name');
    tc.verifyEqual(lsqvce.util.arcKey(71, 'pos', 12, 3), 's71_pos_d12_a3');
    tc.verifyEqual(lsqvce.util.arcKey(1, 'acc_x', 20260901, 2), 's1_acc_x_d20260901_a2', ...
        'a YYYYMMDD day prints exactly; an underscore in the obsType is carried, not parsed');
end

function testPriorStreamIsNamedOnce(tc)
% categories: contract, metamorphic
    % WP-138 (code review of WP-135): the prior pseudo-streams' prefix was
    % spelled out at nine places (and fit parsed a stream key for '|prior_');
    % lsqvce.util.priorStream makes the obsType and recognises it.
    tc.verifyEqual(lsqvce.util.priorStream.of('gravity field'), 'prior_gravityField');
    tc.verifyEqual(lsqvce.util.priorStream.of("datum"), 'prior_datum');
    tc.verifyTrue(lsqvce.util.priorStream.is('prior_x'));
    tc.verifyEqual(lsqvce.util.priorStream.is({'prior_x', 'kbr', 'priorx'}), [true false false]);
    % metamorphic: every name it makes, it recognises; the array call equals
    % the per-name calls
    names = {'a', 'gravity field', '1st', 'x_y'};
    made = cellfun(@lsqvce.util.priorStream.of, names, 'UniformOutput', false);
    tc.verifyTrue(all(lsqvce.util.priorStream.is(made)));
    tc.verifyEqual(lsqvce.util.priorStream.is(made), cellfun(@lsqvce.util.priorStream.is, made));
end

%% =================================================================== helpers
function res = makeRes()
% Minimal res: C has 2 pos arcs (3 epochs each), D has 2 pos arcs.
    res.sys(1).sysID = 'C';
    res.sys(1).obs   = [arc('pos',1,1,3), arc('pos',1,2,3)];
    res.sys(2).sysID = 'D';
    res.sys(2).obs   = [arc('pos',1,1,3), arc('pos',1,2,3)];
end

function a = arc(obsType, day, arcId, nE)
    a = struct('obsType',obsType, 'day',day, 'arc',arcId, ...
               'v', (1:3*nE).', 'w', (1:3*nE).');
end

%% ================= moved from test_adjustVCE (R5 M4 4c): loader / covModel
function test_covModel_composes_cep_with_the_temporal_law(tc)
% categories: contract, precision
    % ENGINE tier (R6 WP-90): a TEMPORAL arma law plus cm.Cep -- what the
    % loader sets from the day -- whitens a 3-axis arc built as
    % D (R_t kron I_3) D to unit-variance, temporally white output; the
    % axis count and the filter memory (p epochs = p*d rows) are reported.
    d = 3;  nE = 400;  dt = 5;
    Ld = [1 0 0; 0.3 1 0; 0.2 0.4 2.5];  Dep = Ld*Ld.';
    ar = [1 -0.5 0.1];  sig = 0.9;
    rs = RandStream('mt19937ar','Seed',11);
    e = zeros(d*nE,1);
    W = zeros(d,nE);
    for ax = 1:d, W(ax,:) = filter(1,ar,randn(rs,nE,1)).'; end
    W = sig * (chol(Dep,'lower') * W);
    e(:) = W(:);                                  % epoch-major
    t = reshape(repmat((0:nE-1)*dt, d, 1), [], 1);
    Cep = repmat(reshape(Dep,1,d,d), nE, 1, 1);
    cm = struct('type','arma','ar',ar,'sigma',sig,'Cep',Cep);
    [~, Wfcn, info] = lsqvce.engine.covModel(cm, t);
    verifyEqual(tc, info.nAxis, d);
    verifyEqual(tc, info.filterMemory, 2 * d);
    verifyEqual(tc, size(info.Cep), [nE d d]);
    A = randn(rs, d*nE, 2);
    [Aw,~,ew] = Wfcn(1:d*nE, A, zeros(d*nE,0), e);
    % whitened residual variance ~ 1 (unit innovation after full whitening)
    % tol: statistical identity on a finite sample; scatter ~ 1/sqrt(dof) of the fixture, bar sized 4/sqrt(Z)-style (HANDOVER 6.1.2)
    verifyEqual(tc, var(ew), 1, 'RelTol', 0.25);
    verifySize(tc, Aw, [d*nE, 2]);
end

function test_covOverride_bySys(tc)
% categories: precision
    % A bySys override for the merged pair system (sysID 1) must rescale the
    % whitener. bySys is keyed by makeValidName(string(sysID)) -> 'x1'.
    [cfg, L] = fixture_grace();
    Lcd = L{1};                                  % the merged pair system, sysID 1
    baseOpt = struct('sysID', 1);
    arcBase = Lcd.closure(cfg.days(1), 'rr', 1, baseOpt);

    optOv = struct('sysID', 1);
    optOv.covOverride.bySys.x1.rr = struct('sigma', 10);
    arcOv = Lcd.closure(cfg.days(1), 'rr', 1, optOv);

    % same design/observations, different whitener scaling -> whitened norm shrinks
    % tol: 0 -- an EXACT identity claim (same stored value, no arithmetic)
    verifyEqual(tc, arcOv.A, arcBase.A, 'AbsTol', 0);
    wB = arcBase.Wfcn(1:numel(arcBase.l), arcBase.A, arcBase.B, arcBase.l);
    wO = arcOv.Wfcn(1:numel(arcOv.l),     arcOv.A,   arcOv.B,   arcOv.l);
    verifyGreaterThan(tc, norm(wB(:)), norm(wO(:)));   % bigger sigma -> smaller w
end

function test_covOverride_flat_backcompat(tc)
% categories: contract
    % The old flat opt.covOverride.(obsType) form must still work.
    [cfg, L] = fixture_grace();
    Lcd = L{1};
    opt = struct('sysID',1, 'covOverride', struct('rr', struct('sigma', 5)));
    arc = Lcd.closure(cfg.days(1), 'rr', 1, opt);
    verifyTrue(tc, isa(arc.Wfcn, 'function_handle'));
end

%% ===================================================== SCALE (VCE) recovers weights


%% ==================================== R5-D8 round 3: validity-aware expansion
function test_validity_columns_expand_day_rows_only_where_valid(tc)
% categories: contract, precision
    % A 'day' template row with a validity [start end] (as [MJD sec] pairs or
    % fractional MJD, in the columns named by ValidityColumns) is instantiated
    % ONLY on the days it overlaps; nl is derived accordingly; buildDay's info
    % says which rows are present on the day (info.dayRows / info.arcRows).
    days = 100:102;
    T = {'day', 'bias',     [], [];                          % always (empty validity)
         'day', 'manoeuvre', [101 0], [101 86390];          % day 101 only
         'day', 'ramp',      100.5, 101.6;                  % days 100 and 101 (fractional MJD)
         'arc', 'ic',        [], []};                       % always, per arc
    seen = containers.Map('KeyType', 'double', 'ValueType', 'any');
    bd = @(day, ot, common, local, po, info) recordInfo(seen, day, info);
    L = lsqvce.engine.makeLoader(bd, 'ObservationTypes', {'obs'}, ...
            'Covariance', struct('obs', struct('type', 'diag', 'sigma', 1)), ...
            'Arcs', @(d, ot) 1:2, 'Days', days, 'LocalParams', T, 'Scoped', true, ...
            'NumCommon', 3, 'ValidityColumns', [3 4], 'SystemId', 1);
    plan = L.plan;
    tc.verifyEqual(plan.dayCount, [2 3 1], 'bias + ramp | bias + manoeuvre + ramp | bias');
    tc.verifyEqual(plan.dayRowsOf{1}, [1 3]);
    tc.verifyEqual(plan.dayRowsOf{2}, [1 2 3]);
    tc.verifyEqual(plan.dayRowsOf{3}, 1);
    % nl = 2 arcs/day x 1 arc row x 3 days + (2 + 3 + 1) day rows
    tc.verifyEqual(plan.nl, 6 + 6);
    % the builder is told the rows present on the day
    L.day(101, 'obs');
    tc.verifyTrue(isKey(seen, 101));
    tc.verifyEqual(seen(101).dayRows, [1 2 3]);
    tc.verifyEqual(seen(101).arcRows, 4);
    tc.verifyEqual(seen(101).dayCount, 3);
    % without validity columns every row is valid always (the v4 behaviour)
    L0 = lsqvce.engine.makeLoader(bd, 'ObservationTypes', {'obs'}, ...
            'Covariance', struct('obs', struct('type', 'diag', 'sigma', 1)), ...
            'Arcs', @(d, ot) 1:2, 'Days', days, 'LocalParams', T, 'Scoped', true, ...
            'NumCommon', 3, 'SystemId', 1);
    tc.verifyEqual(L0.plan.dayCount, [3 3 3]);
    % refusals by name
    tc.verifyError(@() lsqvce.engine.makeLoader(bd, 'ObservationTypes', {'obs'}, ...
            'Covariance', struct('obs', struct('type', 'diag', 'sigma', 1)), ...
            'Arcs', @(d, ot) 1:2, 'Days', days, 'LocalParams', T, 'Scoped', true, ...
            'NumCommon', 3, 'ValidityColumns', [3 9], 'SystemId', 1), 'lsqvce:makeLoader:badValidity');
end

function D = recordInfo(seen, day, info)
    seen(day) = info;
    D = struct('A', {{}}, 'B', {{}}, 'l', {{}}, 't', {{}});      % v6: a day of no cells
end

function test_prior_group_pattern_builds_L_from_the_name_column(tc)
% categories: contract, precision
    % R5-D8 round 6: Priors(g).group is a regexp on the name column; the
    % loader builds L with one unit row per matched template row per
    % instance -- identical to the explicit L the user could have written
    % from lsqvce.engine.localColumns.
    days = 100:101;
    T = {'day', 'bias',  [], [];  'day', 'xacc', [], [];  'day', 'yacc', [], [];  'arc', 'ic', [], []};
    bd = @(day, ot, common, local, po, info) struct('A', {{}}, 'B', {{}}, 'l', {{}}, 't', {{}});   % v6: no cells
    mk = @(P) lsqvce.engine.makeLoader(bd, 'ObservationTypes', {'obs'}, ...
            'Covariance', struct('obs', struct('type', 'diag', 'sigma', 1)), ...
            'Arcs', @(d, ot) 1:2, 'Days', days, 'LocalParams', T, 'Scoped', true, ...
            'NumCommon', 2, 'SystemId', 1, 'Priors', P);
    L = mk(struct('group', 'acc$', 'sigma', 0.01, 'name', 'acc', 'estimate', "scale"));
    g = L.priors(1);
    loc = lsqvce.engine.localColumns(L.plan);
    want = zeros(0, L.plan.nl);
    for c = 1:L.plan.nl
        if ismember(loc(c).row, [2 3]), r = zeros(1, L.plan.nl);  r(c) = 1;  want(end+1, :) = r; end %#ok<AGROW>
    end
    % tol: exact -- unit rows
    tc.verifyEqual(g.L, want, 'xacc and yacc on both days: 4 rows');
    tc.verifyEqual(size(g.L, 1), 4);
    tc.verifyEqual(g.obsType, 'prior_acc');
    tc.verifyEqual(g.estimate, 'scale');
    tc.verifyTrue(any(strcmp(L.obsTypes, 'prior_acc')), 'the group is a stream of its own');
    % explicit L: used as is; both: refused; no match: refused
    L2 = mk(struct('L', want, 'sigma', 0.01, 'name', 'acc'));
    tc.verifyEqual(L2.priors(1).L, want);
    % R6 E-14: a scalar x0 is every row's (the data.priors rule; the GRAVFIRE
    % template emits column 10 as a scalar per group and failed on the bundle)
    L3 = mk(struct('group', 'acc$', 'x0', 0.5, 'sigma', 0.01, 'name', 'acc'));
    tc.verifyEqual(L3.priors(1).x0, 0.5 * ones(4, 1));      % tol: exact -- a repmat
    tc.verifyError(@() mk(struct('group', 'acc$', 'L', want, 'sigma', 1)), 'lsqvce:makeLoader:badPrior');
    tc.verifyError(@() mk(struct('group', '^zzz$', 'sigma', 1)), 'lsqvce:makeLoader:badPrior');
end

%% ============================ R6 WP-90: the error information travels with the day
function test_day_cep_composes_with_the_temporal_law_through_the_loader(tc)
% categories: contract, precision
    % A 3-axis day returns D.Cep (nEpoch x 3 x 3) and declares only the
    % TEMPORAL law (arma). The loader derives nAxis, slices Cep per arc and
    % covModel composes C = D (R_t kron I_3) D: the whitened normals of the
    % arc must equal the dense decorrelation of that covariance -- the kernel
    % identity of test_vceMath, now proven at the SEAM the user's data enters.
    rs = RandStream('mt19937ar', 'Seed', 90);
    d = 3;  nE = 40;  m = nE * d;
    a = [1 -1.2 0.4];  s = 0.8;
    r = s^2 * lsqvce.kernel.armaACV(a, 1, nE);
    Rt = toeplitz(r / r(1));
    Cep = zeros(nE, d, d);  Db = zeros(m);
    for k = 1:nE
        Q = 0.2 * randn(rs, d, d) + diag([1 1 2.5]);  Ck = Q * Q.';
        Cep(k, :, :) = Ck;
        ib = (k-1)*d + (1:d);  Db(ib, ib) = chol(Ck, 'lower');
    end
    C = Db * (kron(Rt, eye(d)) * r(1)) * Db.';
    D = struct('A', randn(rs, m, 4), 'B', randn(rs, m, 2), 'l', randn(rs, m, 1), ...
               't', reshape(repmat((0:nE-1)*5, d, 1), [], 1), 'Cep', Cep);   % NO nAxis: derived
    cov.pos = struct('type', 'arma', 'ar', a, 'sigma', s);               % the temporal law only
    Dc = arcCells(D);                                                    % v6: one cell = one arc
    L = lsqvce.engine.makeLoader(@(day, ot) Dc, 'ObservationTypes', {'pos'}, 'Covariance', cov, ...
                                 'Arcs', 1, 'SystemId', 1, 'NumLocal', 2, 'NumCommon', 4);
    arc = L.closure(1, 'pos', 1, struct('sysID', 1));
    tc.verifyEqual(arc.nAxis, d, 'the axis count is derived from Cep');
    tc.verifyEqual(size(arc.Cep), [nE d d], 'the per-epoch part travels with the arc');
    tc.verifyEqual(arc.filterMemory, 2 * d, 'the filter memory is p epochs = p*d rows');
    [Aw, Bw, lw] = arc.Wfcn(1:m, arc.A, arc.B, arc.l, struct());
    W = lsqvce.kernel.decorrFactor(C);
    Ad = W * D.A;  Bd = W * D.B;  ld = W * D.l;
    % tol: rotation-invariant quadratic forms of the same whitening; eps*cond
    % headroom -- the bound test_vceMath uses for the kernel identity
    tc.verifyEqual(Aw.'*Aw, Ad.'*Ad, 'RelTol', 1e-8);
    tc.verifyEqual(Aw.'*lw, Ad.'*ld, 'RelTol', 1e-8);     % tol: same source
    tc.verifyEqual(Bw.'*Bw, Bd.'*Bd, 'RelTol', 1e-8);     % tol: same source
    % Cfcn is the TEMPORAL covariance over epochs (what the joint estimator forms per axis)
    % tol: the same armaACV arithmetic on the same doubles -- a few ulps
    tc.verifyEqual(arc.Cfcn(1:nE, struct()), s^2 * toeplitz(lsqvce.kernel.armaACV(a, 1, nE)), 'RelTol', 1e-12);
end

function test_epoch_weights_are_an_additive_nugget(tc)
% categories: contract, precision
    % R7 WP-115 part 3b (R7-10, user decision (b)): the joint estimator's
    % per-epoch robust factors w reach the loader as opt.epochWeights and
    % compose as an ADDITIVE nugget on the temporal law,
    %   C = D ((R_t + diag(q .* Sigma_ee)) kron I_d) D,  q = 1/w - 1,
    % not as the old fold D / sqrt(w), which scaled the neighbours'
    % covariances too. The arc's whitened normals must equal the dense
    % decorrelation of THAT covariance -- a zero weight (floored 1e-12)
    % included -- and the arc must announce itself as whitened whole.
    rs = RandStream('mt19937ar', 'Seed', 1153);
    d = 3;  nE = 40;  m = nE * d;
    a = [1 -0.9];  s = 0.8;
    r = s^2 * lsqvce.kernel.armaACV(a, 1, nE);
    Cep = zeros(nE, d, d);  Db = zeros(m);
    for k = 1:nE
        Q = 0.2 * randn(rs, d, d) + diag([1 1 2.5]);  Ck = Q * Q.';
        Cep(k, :, :) = Ck;
        ib = (k-1)*d + (1:d);  Db(ib, ib) = chol(Ck, 'lower');
    end
    w = ones(nE, 1);  w([3 11 12 30]) = [0.4 0 0.05 0.9];
    q = 1 ./ max(w, 1e-12) - 1;
    C = Db * kron(toeplitz(r) + diag(q * r(1)), eye(d)) * Db.';
    D = struct('A', randn(rs, m, 4), 'B', randn(rs, m, 2), 'l', randn(rs, m, 1), ...
               't', reshape(repmat((0:nE-1)*5, d, 1), [], 1), 'Cep', Cep);
    cov.pos = struct('type', 'arma', 'ar', a, 'sigma', s);
    L = lsqvce.engine.makeLoader(@(day, ot) arcCells(D), 'ObservationTypes', {'pos'}, 'Covariance', cov, ...
                                 'Arcs', 1, 'SystemId', 1, 'NumLocal', 2, 'NumCommon', 4);
    ew = containers.Map({lsqvce.util.rowWeightKey(1, 1, 'pos', 1)}, {w});
    arc = L.closure(1, 'pos', 1, struct('sysID', 1, 'epochWeights', ew));
    tc.verifyTrue(arc.wholeArc, 'a nugget arc is whitened in ONE batch');
    % tol: 1/max(w,1e-12) - 1 recomputed on the same doubles -- exact
    tc.verifyEqual(arc.nugget, q, 'the relative nugget travels with the arc');
    tc.verifyEqual(arc.Cep, Cep, 'Cep stays the DECLARED error information (no fold)');
    [Aw, Bw, lw] = arc.Wfcn(1:m, arc.A, arc.B, arc.l, struct());
    W = lsqvce.kernel.decorrFactor(C);
    Ad = W * D.A;  Bd = W * D.B;  ld = W * D.l;
    % tol: rotation-invariant quadratic forms of two square roots of the same
    % inverse (G F vs the dense factor); the mirror reproduces A'C^-1 A to
    % 2e-15 with a 1e12 nugget, 1e-8 is the composition tests' bound
    tc.verifyEqual(Aw.'*Aw, Ad.'*Ad, 'RelTol', 1e-8);
    tc.verifyEqual(Aw.'*lw, Ad.'*ld, 'RelTol', 1e-8);     % tol: same source
    tc.verifyEqual(Bw.'*Bw, Bd.'*Bd, 'RelTol', 1e-8);     % tol: same source
    % WP-133 (review 3 finding 14): the zero-weight epoch 11 is excluded -- its
    % d rows -- and its whitened rows are ~0 IN PLACE (G F spread the deleted
    % direction over the next epochs; the rotation puts it back on its rows)
    ex = false(m, 1);  ex(31:33) = true;
    tc.verifyEqual(arc.excluded, ex, 'the rejected epoch''s rows are marked excluded');
    % tol: the floor keeps sqrt(1e-12) ~ 1e-6 of the epoch's weight (the
    % mirror wp133_nugget_rotation_mirror.py: <= 4e-7 of the largest row)
    rowNorm = vecnorm([Aw Bw lw], 2, 2);
    tc.verifyLessThan(max(rowNorm(ex)) / max(rowNorm), 1e-5, 'the deleted rows whiten to ~0 in place');
    tc.verifyError(@() arc.Wfcn(1:m-d, arc.A(1:m-d, :), arc.B(1:m-d, :), arc.l(1:m-d), struct()), ...
        'lsqvce:covModel:nuggetNeedsWholeArc');
    % the batch iterator honours wholeArc: BatchSize 7 accumulates the same
    % normals as one batch (a causal-filter overlap would be wrong here)
    sols = cell(1, 2);  bss = [7 1e9];
    for kb = 1:2
        acc = lsqvce.engine.BlockAccumulator(4, 'Backend', 'qr', 'BatchSize', bss(kb), ...
                                             'UserOpt', struct('epochWeights', ew));
        acc.addSystem(1, 2, L.closure, 1);
        sols{kb} = acc.solve();
    end
    % tol: one QR of the same whitened rows, merged in two orders
    tc.verifyEqual(sols{1}.x, sols{2}.x, 'RelTol', 1e-10, 'BatchSize does not touch a nugget arc');
    tc.verifyEqual(sols{1}.dof, m - 3 - 6, 'the 3 rejected rows leave the dof (m - 3 - nc - nl)');
end

function test_excluded_epoch_is_deleted_from_a_correlated_law(tc)
% categories: reference, precision
    % WP-148: a day's .excluded rows (zero in A, B, l) kept their epoch in the
    % law's covariance, a phantom observation that conditioned its neighbours:
    % the normals were A_o' inv(C)_oo A_o, 2e-3 off the deletion at phi 0.9
    % (measured). An epoch whose d rows are all excluded is now deleted by the
    % nugget of an epoch weight 0: the normals are A_o' C_oo^-1 A_o, with C_oo
    % the covariance of the KEPT rows on the untouched grid.
    rs = RandStream('mt19937ar', 'Seed', 148);
    d = 3;  nE = 40;  m = nE * d;  e = [9 20];
    a = [1 -0.9];  s = 0.8;
    r = s^2 * lsqvce.kernel.armaACV(a, 1, nE);
    Cep = zeros(nE, d, d);  Db = zeros(m);
    for k = 1:nE
        Q = 0.2 * randn(rs, d, d) + diag([1 1 2.5]);  Ck = Q * Q.';
        Cep(k, :, :) = Ck;
        ib = (k-1)*d + (1:d);  Db(ib, ib) = chol(Ck, 'lower');
    end
    C = Db * kron(toeplitz(r), eye(d)) * Db.';
    ex = false(m, 1);  for k = e, ex((k-1)*d + (1:d)) = true; end
    A = randn(rs, m, 4);  B = randn(rs, m, 2);  l = randn(rs, m, 1);
    A(ex, :) = 0;  B(ex, :) = 0;  l(ex) = 0;                      % weight 0: zero rows
    D = struct('A', A, 'B', B, 'l', l, 't', reshape(repmat((0:nE-1)*5, d, 1), [], 1), ...
               'Cep', Cep);
    Dc = arcCells(D);  Dc.excluded = {ex};                        % one arc: its excluded rows
    cov.pos = struct('type', 'arma', 'ar', a, 'sigma', s);
    L = lsqvce.engine.makeLoader(@(day, ot) Dc, 'ObservationTypes', {'pos'}, 'Covariance', cov, ...
                                 'Arcs', 1, 'SystemId', 1, 'NumLocal', 2, 'NumCommon', 4);
    arc = L.closure(1, 'pos', 1, struct('sysID', 1));
    tc.verifyTrue(arc.wholeArc, 'the deletion is the nugget: the arc is whitened whole');
    tc.verifyEqual(arc.excluded, ex, 'the excluded rows stay marked');
    [Aw, Bw, lw] = arc.Wfcn(1:m, arc.A, arc.B, arc.l, struct());
    k = ~ex;  Co = C(k, k);
    % tol: 1e-8 relative -- the nugget tests' bound (two square roots of the
    % same inverse; the scalar mirror of this case agrees to 3e-15)
    tc.verifyEqual(Aw.'*Aw, A(k,:).' / Co * A(k,:), 'RelTol', 1e-8, 'A_o'' C_oo^-1 A_o');
    tc.verifyEqual(Aw.'*lw, A(k,:).' / Co * l(k), 'RelTol', 1e-8);     % tol: same source
    tc.verifyEqual(Bw.'*Bw, B(k,:).' / Co * B(k,:), 'RelTol', 1e-8);   % tol: same source
    acc = lsqvce.engine.BlockAccumulator(4, 'Backend', 'qr');
    acc.addSystem(1, 2, L.closure, 1);
    tc.verifyEqual(acc.solve().dof, nnz(k) - 4 - 2, 'the excluded epochs leave the dof');
    % a 'diag' law needs no nugget: its epochs are independent, a zero row is a deletion
    covD.pos = struct('type', 'diag', 'sigma', s);
    LD = lsqvce.engine.makeLoader(@(day, ot) Dc, 'ObservationTypes', {'pos'}, 'Covariance', covD, ...
                                  'Arcs', 1, 'SystemId', 1, 'NumLocal', 2, 'NumCommon', 4);
    arcD = LD.closure(1, 'pos', 1, struct('sysID', 1));
    tc.verifyFalse(lsqvce.util.getfielddef(arcD, 'wholeArc', false), 'a diag law keeps its batches');
end

function test_excluded_rows_the_nugget_misses_are_deleted_exactly(tc)
% categories: reference, precision
    % WP-153 (review 12.84, C4). WP-148 deletes an epoch whose d rows are all
    % excluded, through the epoch nugget; the rest stayed zero rows in the
    % covariance and conditioned their neighbours -- every excluded row of a
    % 'var' law (no nugget there; normals 1.1e-2 off for one epoch of a
    % VAR(1)) and an epoch excluded on SOME of its axes (x 0.35 formal sd off
    % for one axis at phi 0.9; a 'diag' law with a FULL Cep too, through the
    % epoch's own covariance). covModel's deleteRows deletes them exactly.
    % THE ORACLE: A_o' C_oo^-1 A_o, C_oo the covariance of the kept rows; and
    % through fit, x and the dof of that solve, the verdict's sum(r) closing.
    rs = RandStream('mt19937ar', 'Seed', 153);
    dt = 5;  nE = 40;  d = 2;  m = nE * d;
    t = reshape(repmat((0:nE-1) * dt, d, 1), [], 1);
    cases = struct('name', {}, 'spec', {}, 'Cep', {}, 'C', {}, 'ex', {});
    cmv = struct('type', 'var', 'ar', [0.9 0.05; 0 0.8], 'sigma', eye(2));
    Cv = lsqvce.engine.covModel(cmv, t, struct('nAxis', d));
    ex = false(m, 1);  ex([19 20 41]) = true;        % epoch 10 whole, epoch 21 axis 1
    cases(end+1) = struct('name', 'var', 'spec', cmv, 'Cep', [], 'C', Cv(1:m, struct('nAxis', d)), 'ex', ex);
    cma = struct('type', 'arma', 'ar', [1 -0.9], 'sigma', 0.8);
    Ct = lsqvce.engine.covModel(cma, t(1:d:end), struct());
    Cep = zeros(nE, d, d);  Db = zeros(m);
    for k = 1:nE
        Q = 0.2 * randn(rs, d, d) + diag([1 2]);  Ck = Q * Q.';
        Cep(k, :, :) = Ck;  ib = (k-1)*d + (1:d);  Db(ib, ib) = chol(Ck, 'lower');
    end
    ex = false(m, 1);  ex([20 41]) = true;           % one axis of epochs 10 and 21
    cases(end+1) = struct('name', 'arma + full Cep', 'spec', cma, 'Cep', Cep, ...
                          'C', Db * kron(Ct(1:nE, struct()), eye(d)) * Db.', 'ex', ex);
    cases(end+1) = struct('name', 'diag + full Cep', 'spec', struct('type', 'diag', 'sigma', 1), ...
                          'Cep', Cep, 'C', Db * Db.', 'ex', ex);
    for c = cases
        A = randn(rs, m, 3);  l = chol(c.C, 'lower') * randn(rs, m, 1) + A * [1; 2; 3];
        A(c.ex, :) = 0;  l(c.ex) = 0;                % weight 0: zero rows
        D = struct('A', {{A}}, 'B', {{zeros(m, 0)}}, 'l', {{l}}, 't', {{t}}, 'excluded', {{c.ex}});
        if isempty(c.Cep), D.nAxis = d; else, D.Cep = {c.Cep}; end
        L = lsqvce.engine.makeLoader(@(day, ot) D, 'ObservationTypes', {'o'}, ...
                'Covariance', struct('o', c.spec), 'Arcs', 1, 'SystemId', 1, 'NumLocal', 0, 'NumCommon', 3);
        arc = L.closure(1, 'o', 1, struct('sysID', 1));
        tc.verifyTrue(arc.wholeArc, [c.name ': a deletion whitens the arc whole']);
        [Aw, ~, lw] = arc.Wfcn(1:m, arc.A, arc.B, arc.l, struct());
        k = ~c.ex;  Co = c.C(k, k);
        % tol: 1e-8 relative -- the nugget tests' bound (two square roots of
        % the same inverse; measured 1e-14 through fit)
        tc.verifyEqual(Aw.'*Aw, A(k,:).' / Co * A(k,:), 'RelTol', 1e-8, [c.name ': A_o'' C_oo^-1 A_o']);
        tc.verifyEqual(Aw.'*lw, A(k,:).' / Co * l(k), 'RelTol', 1e-8);     % tol: same source
        % tol: 1e-12 of the row scale -- a deleted row whitens to ~0 in place (WP-133)
        tc.verifyLessThan(max(abs(lw(c.ex))), 1e-12 * max(abs(lw)), [c.name ': the deleted rows are ~0 in place']);
        data = struct('commonParams', {{'common','c1';'common','c2';'common','c3'}}, 'days', 1);
        data.systems = struct('buildDay', @(day, ot) D, 'covariance', struct('o', c.spec), ...
                              'arcs', 1, 'id', 1, 'name', 'c4');
        R = lsqvce.fit(data, lsqvce.Options(Verbose=0));
        Ne = A(k,:).' / Co * A(k,:);  xe = Ne \ (A(k,:).' / Co * l(k));
        % tol: 1e-9 formal sd -- measured 6e-15 to 3e-14; the zero rows put
        % x 0.035-0.35 sd off
        tc.verifyLessThan(max(abs(R.x - xe) ./ sqrt(diag(inv(Ne)))), 1e-9, [c.name ': fit''s x is the deletion''s']);
        tc.verifyEqual(R.sol.dof, nnz(k) - 3, [c.name ': the excluded rows leave the dof']);
        % tol: the verdict's own sum(r) integrity line is 1e-6
        tc.verifyLessThan(R.diag.sumRRelErr, 1e-9, [c.name ': sum(r) counts the live rows']);
    end
    % the loader's field, refused by name when it does not fit the arc
    tc.verifyError(@() lsqvce.engine.covModel(struct('type', 'diag', 'sigma', 1, 'deleteRows', true(3, 1)), ...
                   (1:5).' * dt, struct()), 'lsqvce:covModel:badDeleteRows');
end

function test_cep_of_one_axis_is_standard_deviations(tc)
% categories: contract
    % WP-153 (review 12.84, C3): MATLAB stores an nEpoch x 1 x 1 array as
    % nEpoch x 1, so a d = 1 per-epoch VARIANCE array cannot be told from the
    % documented nEpoch x 1 SIGMAS -- a variance of 4 whitens 1 to 0.25. The
    % contract is stated (covModel, the API reference): a d = 1 Cep is
    % standard deviations; this pins it.
    Cv = 4 * ones(5, 1, 1);
    tc.verifyEqual(size(Cv), [5 1], 'MATLAB drops the trailing singletons');
    [~, W] = lsqvce.engine.covModel(struct('type', 'diag', 'sigma', 1, 'Cep', Cv), (1:5).', struct());
    [~, ~, lw] = W(1:5, zeros(5, 0), zeros(5, 0), ones(5, 1));
    % tol: one division
    tc.verifyEqual(lw, 0.25 * ones(5, 1), 'AbsTol', 1e-15, 'read as sigma 4');
end

function test_nugget_and_ma_startup_respect_max_dense_arc(tc)
% categories: robustness, contract
    % WP-153 (review 12.84, C5). Two allocations the maxDenseArc limit did not
    % reach: (1) a nugget whitener holds three m x k arrays, k ~ m under a
    % redescending influence's partial weights (180 GB at 86 400 epochs), and
    % it read the law's variances through a dense k x k Ct(z), whose guard
    % then blamed the noise sweep; (2) an MA part's exact start-up is a dense
    % min(filter memory, m)^2 block -- an MA root near the unit circle asked
    % for up to 65 536^2 doubles. Both are refused above the limit, by name.
    t = (0:29).' * 5;
    cm = struct('type', 'arma', 'ar', [1 -0.5], 'sigma', 1, 'Cep', ones(30, 1), ...
                'nugget', 0.5 * ones(30, 1), 'maxDenseArc', 10);
    [~, W] = lsqvce.engine.covModel(cm, t, struct());
    err = caught(@() W(1:30, ones(30, 1), zeros(30, 0), ones(30, 1), struct()));
    tc.verifyEqual(err.identifier, 'lsqvce:covModel:denseArcTooLong');
    tc.verifySubstring(err.message, 'robust nuggets', 'the message names the nugget, not the sweep');
    cm.maxDenseArc = 30;                              % 30 x 30 epoch products fit
    [~, W] = lsqvce.engine.covModel(cm, t, struct());
    [Aw, ~, ~] = W(1:30, ones(30, 1), zeros(30, 0), ones(30, 1), struct());
    tc.verifyTrue(all(isfinite(Aw)), 'within the limit the nugget whitens');
    t = (0:99).' * 5;
    ma = struct('type', 'arma', 'ar', [1 -0.5], 'ma', [1 0.95], 'sigma', 1, 'maxDenseArc', 50);
    err = caught(@() lsqvce.engine.covModel(ma, t, struct()));
    tc.verifyEqual(err.identifier, 'lsqvce:covModel:denseArcTooLong');
    tc.verifySubstring(err.message, 'MA part', 'the message names the MA start-up');
    ma.ma = [1 0.1];                                  % decays within 17 samples
    tc.verifyWarningFree(@() lsqvce.engine.covModel(ma, t, struct()), 'a benign MA part is not refused');
end

function err = caught(fh)
%CAUGHT  The exception fh() throws (identifier '' and message '' when none).
    err = struct('identifier', '', 'message', '');
    try
        fh();
    catch e
        err = struct('identifier', e.identifier, 'message', e.message);
    end
end

function test_nugget_refusals_by_name(tc)
% categories: contract
    % covModel's cm.nugget composes with the per-epoch Cep only, one finite
    % non-negative value per epoch (R7 WP-115 part 3b)
    t = (0:9).' * 5;
    law = struct('type', 'arma', 'ar', [1 -0.5], 'sigma', 1);
    tc.verifyError(@() lsqvce.engine.covModel(setfield(law, 'nugget', ones(10, 1)), t), ...
        'lsqvce:covModel:badNugget');                               % no Cep
    lc = setfield(law, 'Cep', ones(10, 1));                          %#ok<SFLD>
    tc.verifyError(@() lsqvce.engine.covModel(setfield(lc, 'nugget', ones(9, 1)), t), ...
        'lsqvce:covModel:badNugget');                               % wrong length
    tc.verifyError(@() lsqvce.engine.covModel(setfield(lc, 'nugget', -ones(10, 1)), t), ...
        'lsqvce:covModel:badNugget');                               % negative
    [~, ~, info] = lsqvce.engine.covModel(setfield(lc, 'nugget', zeros(10, 1)), t);
    tc.verifyFalse(isfield(info, 'wholeArc'), 'an all-zero nugget is the plain composition');
end

function test_scalar_stream_cep_is_a_per_epoch_sigma(tc)
% categories: contract, precision
    % A scalar stream returns D.Cep as nEpoch x 1: the per-epoch standard
    % deviation (finding E-03: this information could not travel with the
    % day before). With a white law of level sigma the whitened rows are
    % l ./ (sigma * Cep); with an arma law the composition equals the dense
    % reference D R_t D.
    rs = RandStream('mt19937ar', 'Seed', 91);
    m = 50;  sig = 0.3;
    cep = 0.5 + rand(rs, m, 1);
    D = struct('A', randn(rs, m, 2), 'B', zeros(m, 0), 'l', randn(rs, m, 1), 't', (0:m-1).' * 5, 'Cep', cep);
    Dc = arcCells(D);                                                    % v6: one cell = one arc
    mk = @(cm) lsqvce.engine.makeLoader(@(day, ot) Dc, 'ObservationTypes', {'rr'}, ...
             'Covariance', struct('rr', cm), 'Arcs', 1, 'SystemId', 1, 'NumLocal', 0, 'NumCommon', 2);
    arc = mk(struct('type', 'diag', 'sigma', sig)).closure(1, 'rr', 1, struct('sysID', 1));
    tc.verifyFalse(isfield(arc, 'nAxis'), 'a scalar stream stays scalar');
    [Aw, ~, lw] = arc.Wfcn(1:m, arc.A, arc.B, arc.l, struct());
    % tol: two divisions vs one product-and-division of the same doubles -- a few ulps
    tc.verifyEqual(lw, D.l ./ (sig * cep), 'RelTol', 1e-14);
    tc.verifyEqual(Aw, D.A ./ (sig * cep), 'RelTol', 1e-14);   % tol: same source
    a = [1 -0.7];  r = sig^2 * lsqvce.kernel.armaACV(a, 1, m);
    arc = mk(struct('type', 'arma', 'ar', a, 'sigma', sig)).closure(1, 'rr', 1, struct('sysID', 1));
    tc.verifyEqual(size(arc.Cep), [m 1], 'the per-epoch sigma travels with the arc');
    C = diag(cep) * toeplitz(r) * diag(cep);
    W = lsqvce.kernel.decorrFactor(C);
    [Aw, ~, lw] = arc.Wfcn(1:m, arc.A, arc.B, arc.l, struct());
    % tol: rotation-invariant quadratic forms of the same whitening; eps*cond headroom
    tc.verifyEqual(Aw.'*Aw, (W*D.A).'*(W*D.A), 'RelTol', 1e-8);
    tc.verifyEqual(Aw.'*lw, (W*D.A).'*(W*D.l), 'RelTol', 1e-8);   % tol: same source
end

function test_cep_refusals_by_name(tc)
% categories: contract, robustness
    % Every documented refusal of the WP-90 day contract fires with its own
    % id and only that id; the retired 'separable' type names its
    % replacement; a multi-axis day WITHOUT Cep is composed with D = I and
    % said ONCE per run (the omission is distinguishable from a choice).
    rs = RandStream('mt19937ar', 'Seed', 92);
    nE = 10;  d = 3;  m = nE * d;
    base = arcCells(struct('A', randn(rs, m, 2), 'B', zeros(m, 0), 'l', randn(rs, m, 1), ...
                  't', reshape(repmat((0:nE-1)*5, d, 1), [], 1)));    % v6: one cell = one arc
    mk = @(D, cm) lsqvce.engine.makeLoader(@(day, ot) D, 'ObservationTypes', {'pos'}, ...
             'Covariance', struct('pos', cm), 'Arcs', 1, 'SystemId', 1, 'NumLocal', 0, 'NumCommon', 2);
    read = @(D, cm) mk(D, cm).closure(1, 'pos', 1, struct('sysID', 1));
    arm = struct('type', 'arma', 'ar', [1 -0.5], 'sigma', 1);
    D = base;  D.Cep = {ones(nE, d)};  D.nAxis = 2;                     % nAxis disagrees with Cep
    tc.verifyError(@() read(D, arm), 'lsqvce:makeLoader:axisConflict');
    D = base;  D.Cep = {ones(nE - 1, d)};                               % Cep does not cover the arc
    tc.verifyError(@() read(D, arm), 'lsqvce:makeLoader:badDayCep');
    D = base;  D.Cep = {ones(1, m)};                                    % a row is not a per-epoch sigma
    tc.verifyError(@() read(D, arm), 'lsqvce:makeLoader:badDayCep');
    D = base;  D.Cep = {ones(nE, d, 2)};                                % d x d only
    tc.verifyError(@() read(D, arm), 'lsqvce:makeLoader:badDayCep');
    D = base;  D.Cep = ones(nE, d);                                     % v6: Cep must be a cell like A
    tc.verifyError(@() read(D, arm), 'lsqvce:makeLoader:badDayCep');
    D = base;  D.nAxis = 4;                                             % nAxis must divide the rows
    tc.verifyError(@() read(D, arm), 'lsqvce:makeLoader:axisConflict');
    D = base;  D.Cep = {ones(nE, d)};
    cm = arm;  cm.d = 2;                                                % the spec's d disagrees
    tc.verifyError(@() read(D, cm), 'lsqvce:makeLoader:axisConflict');
    cm = arm;  cm.Cep = ones(nE, d);                                    % Cep on the spec AND the day
    tc.verifyError(@() read(D, cm), 'lsqvce:makeLoader:cepTwice');
    vr = struct('type', 'var', 'ar', 0.3 * eye(d), 'sigma', eye(d));    % a VAR has its own cross-axis structure
    tc.verifyError(@() read(D, vr), 'lsqvce:covModel:cepOnVar');
    caught = MException.empty;                                          % the retired type
    try
        read(D, struct('type', 'separable', 'ar', [1 -0.5], 'sigma', 1, 'd', d));
    catch e
        caught = e;
    end
    tc.verifyNotEmpty(caught, 'the retired type is refused');
    tc.verifyEqual(caught.identifier, 'lsqvce:covModel:unknownType');
    tc.verifySubstring(caught.message, 'Cep');
    % a multi-axis day without Cep: D = I, reported once per run
    w = warning('off', 'lsqvce:makeLoader:noCep');
    tc.addTeardown(@() warning(w));
    lsqvce.util.log.inventory('reset');
    D = base;  D.nAxis = d;
    L = mk(D, arm);
    arc = L.closure(1, 'pos', 1, struct('sysID', 1));
    tc.verifyEqual(arc.nAxis, d);
    tc.verifyEqual(arc.Cep, ones(nE, d), 'composed with D = I');
    L.closure(1, 'pos', 1, struct('sysID', 1));                          % a second read: silent
    inv = lsqvce.util.log.inventory();
    tc.verifyEqual(sum(strcmp(inv, 'lsqvce:makeLoader:noCep')), 1, 'said once per run');
end

function test_cep_is_the_arcs_own_cell(tc)
% categories: contract, robustness
    % Finding E-01 (v5): the day's Cep was sliced with the arc ID instead of
    % the arc's POSITION. v6 (R6 WP-102) removes the slice: an arc's Cep is
    % its own cell. A day generated with arc 3's rows first and arc 1's
    % second, cut into the schedule [1 3], must whiten arc 3 with ITS
    % epochs' sigmas. Scalar stream, white law of level 1: whitened l = l ./ Cep.
    rs = RandStream('mt19937ar', 'Seed', 93);
    n1 = 8;  n3 = 12;  m = n1 + n3;
    cep = [10 * ones(n3, 1); 0.1 * ones(n1, 1)];        % rows of arc 3 first, then arc 1
    D = struct('A', randn(rs, m, 2), 'B', zeros(m, 0), 'l', randn(rs, m, 1), 't', (0:m-1).' * 5, ...
               'Cep', cep, 'arcIdx', {{1:n3, n3 + (1:n1)}}, 'arcId', [3 1]);
    Dc = arcCells(D, [1 3]);                             % cells in schedule order
    tc.verifyEqual(Dc.Cep{2}, cep(1:n3), 'arc 3 is the second cell, with its own sigmas');
    L = lsqvce.engine.makeLoader(@(day, ot) Dc, 'ObservationTypes', {'rr'}, ...
             'Covariance', struct('rr', struct('type', 'diag', 'sigma', 1)), ...
             'Arcs', @(day, ot) [1 3], 'SystemId', 1, 'NumLocal', 0, 'NumCommon', 2);
    a3 = L.closure(1, 'rr', 3, struct('sysID', 1));
    [~, ~, w3] = a3.Wfcn(1:n3, a3.A, a3.B, a3.l, struct());
    % tol: a division by the same double (10, then by sigma 1) -- exact
    tc.verifyEqual(w3, D.l(1:n3) / 10, 'AbsTol', 0);
    a1 = L.closure(1, 'rr', 1, struct('sysID', 1));
    [~, ~, w1] = a1.Wfcn(1:n1, a1.A, a1.B, a1.l, struct());
    % tol: the same division by the same double -- a few ulps at most
    tc.verifyEqual(w1, D.l(n3 + (1:n1)) / 0.1, 'RelTol', 1e-15);
end

function test_psd_law_composes_on_the_epoch_grid(tc)
% categories: precision
    % A PSD law on a 3-axis day: the temporal covariance must be built with
    % the EPOCH spacing (5 s), not the row spacing an epoch-major t shows
    % (0, 0, 5, 0, 0, 5, ...) -- the round-3b class of bug. Reference:
    % psd2acv at dt = 5 over nE epochs.
    rs = RandStream('mt19937ar', 'Seed', 94);
    d = 3;  nE = 64;  m = nE * d;  dt = 5;
    f = linspace(0, 0.5/dt, 257).';
    H = 1 ./ (exp(-1i*2*pi*f*dt*(0:1)) * [1 -0.6].');
    S = 0.7^2 * dt * 2 * abs(H).^2;
    D = struct('A', randn(rs, m, 2), 'B', zeros(m, 0), 'l', randn(rs, m, 1), ...
               't', reshape(repmat((0:nE-1)*dt, d, 1), [], 1), 'Cep', ones(nE, d));
    Dc = arcCells(D);                                                    % v6: one cell = one arc
    L = lsqvce.engine.makeLoader(@(day, ot) Dc, 'ObservationTypes', {'pos'}, ...
             'Covariance', struct('pos', struct('type', 'psd', 'f', f, 'S', S)), ...
             'Arcs', 1, 'SystemId', 1, 'NumLocal', 0, 'NumCommon', 2);
    arc = L.closure(1, 'pos', 1, struct('sysID', 1, 'psdMode', 'toeplitz'));
    tc.verifyEqual(arc.nAxis, d);
    Cref = toeplitz(lsqvce.kernel.psd2acv(f, S, dt, nE));
    % tol: the same psd2acv quadrature on the same grid; covModel's 1e-8
    % ridge on r(1) (psdRidge default) is the only difference -> 1e-7 relative
    tc.verifyEqual(arc.Cfcn(1:nE, struct()), Cref, 'RelTol', 1e-7);
    % and the dense composition whitens: A' inv(C) A of the row covariance
    Cfull = kron(Cref, eye(d));
    [Aw, ~, ~] = arc.Wfcn(1:m, arc.A, arc.B, arc.l, struct('psdMode', 'toeplitz'));
    % tol: quadratic forms of the same whitening (dense per axis); eps*cond headroom
    tc.verifyEqual(Aw.'*Aw, D.A.' * (Cfull \ D.A), 'RelTol', 1e-6);
end
