function tests = test_reliability
%TEST_RELIABILITY  Error handling, edge cases and degenerate inputs.
%
%   A numerical library is only as good as its behaviour when the inputs are
%   wrong. These tests assert that every failure mode is either handled
%   gracefully or reported with a specific, actionable error identifier --
%   never a silent wrong answer, and never a bare MATLAB index error.
%
%   Run:  results = runtests('tests/test_reliability');  disp(table(results))
%
%   See also TEST_PRECISION, TEST_SPEED, TEST_LSQACCUM
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    tc.TestData.cfg = makeCfg();
end

%% ============================================ input validation: primitives
function testPSD2ACVRejectsBadInput(tc)
% categories: contract
    f = linspace(0, 0.1, 100).';   S = ones(100,1);
    tc.verifyError(@() lsqvce.kernel.psd2acv(f, S(1:50), 5, 10), 'MATLAB:psd2acv:incorrectNumel');
    tc.verifyError(@() lsqvce.kernel.psd2acv(f, -S,      5, 10), 'MATLAB:psd2acv:expectedNonnegative');
    tc.verifyError(@() lsqvce.kernel.psd2acv(f, S,       0, 10), 'MATLAB:psd2acv:expectedPositive');
    tc.verifyError(@() lsqvce.kernel.psd2acv(f, S,       5,  0), 'MATLAB:psd2acv:expectedPositive');
    tc.verifyError(@() lsqvce.kernel.psd2acv(flipud(f), S, 5, 10), ...
                   'MATLAB:psd2acv:expectedIncreasing');
end

function testPSD2ACVWarnsAboveNyquist(tc)
% categories: contract
% Power above Nyquist is aliased in the data but not in the integral, so the
% resulting ACV would not match the sampled process.
    dt = 5;
    f  = linspace(0, 2/dt, 200).';        % 2x above Nyquist
    S  = ones(200,1);
    tc.verifyWarning(@() lsqvce.kernel.psd2acv(f, S, dt, 10), 'lsqvce:psd2acv:aboveNyquist');
end

function testPSD2ACVZeroLagIdentity(tc)
% categories: precision
    f = linspace(0, 0.1, 2048).';
    S = exp(-f/0.01);
    r = lsqvce.kernel.psd2acv(f, S, 5, 3);
    % tol: r(0)=trapz(f,S) is the quadrature identity (same sum, same order);
    % house round-trip bound 1e-12 (HANDOVER 6.1.2)
    tc.verifyEqual(r(1), trapz(f, S), 'RelTol', 1e-12);
end

function testARMAACVRejectsNonStationary(tc)
% categories: contract
% A non-stationary AR has no autocovariance. Better to say so than to return
% numbers that look plausible.
    tc.verifyError(@() lsqvce.kernel.armaACV([1 -1.5], 1, 10), 'lsqvce:armaACV:nonStationary');
    tc.verifyError(@() lsqvce.kernel.armaACV([1 -1.0], 1, 10), 'lsqvce:armaACV:nonStationary');
    tc.verifyError(@() lsqvce.kernel.armaACV([0  0.5], 1, 10), 'lsqvce:armaACV:zeroLeadingAR');
end

function testLevinsonRejectsBadInput(tc)
% categories: contract
    r = lsqvce.kernel.armaACV([1 -0.9], 1, 5);
    tc.verifyError(@() lsqvce.kernel.levinsonDurbin(r, 10), 'lsqvce:levinson:orderTooHigh');
    tc.verifyError(@() lsqvce.kernel.levinsonDurbin([-1 0.5 0.2], 2), 'lsqvce:levinson:notPD');
end

function testLevinsonWarnsOnNonPDSequence(tc)
% categories: contract
% An autocovariance that is not positive definite (from a PSD with zeros, say)
% must be caught, not silently turned into an unstable filter.
    r = [1 0.9 0.9 0.9 -0.9];             % not a valid ACV (Toeplitz has a negative eigenvalue)
    [a, ~] = tc.verifyWarning(@() lsqvce.kernel.levinsonDurbin(r, 4), 'lsqvce:levinson:notPD');
    tc.verifyLessThan(numel(a)-1, 4);     % truncated below the requested order
end

function testArmaWhitenRejectsSizeMismatch(tc)
% categories: contract
    a = [1 -0.9];  r = lsqvce.kernel.armaACV(a, 1, 10);
    A = zeros(20, 3);   l = zeros(15, 1);
    tc.verifyError(@() lsqvce.kernel.armaWhiten(a, 1, 0.5, r, A, zeros(20,2), l, true), ...
                   'lsqvce:armaWhiten:sizeMismatch');
    tc.verifyError(@() lsqvce.kernel.armaWhiten(a, 1, 0, r, A, zeros(20,2), zeros(20,1), true), ...
                   'lsqvce:armaWhiten:badSigma');
end

function testDecorrFactorRejectsNonCovariance(tc)
% categories: contract
    C = -eye(5);
    tc.verifyError(@() lsqvce.kernel.decorrFactor(C), 'lsqvce:decorrFactor:notPD');
end

%% ============================================ input validation: covariance
function testCovModelRejectsUnknownType(tc)
% categories: contract
    t = (0:99).'*5;
    tc.verifyError(@() lsqvce.engine.covModel(struct('type','bogus'), t, struct()), ...
                   'lsqvce:covModel:unknownType');
    tc.verifyError(@() lsqvce.engine.covModel(struct('nope',1), t, struct()), ...
                   'lsqvce:covModel:noType');
end

function testCovModelRejectsBadCep(tc)
% categories: contract
    % R7-28 / 12.14 E-10: badCep survived WP-90 and was recorded as pinned,
    % yet no suite fired it. Three of its sites: rows vs epochs x axes, a
    % non-square per-epoch covariance, and an axis count stated twice.
    t = (0:5).'*5;                                   % 6 rows
    tc.verifyError(@() lsqvce.engine.covModel(struct('type','diag','sigma',1, ...
                   'Cep',ones(4,2)), t, struct()), 'lsqvce:covModel:badCep');   % 8 rows declared
    tc.verifyError(@() lsqvce.engine.covModel(struct('type','diag','sigma',1, ...
                   'Cep',ones(3,2,3)), t, struct()), 'lsqvce:covModel:badCep');  % d x d' not square
    tc.verifyError(@() lsqvce.engine.covModel(struct('type','diag','sigma',1, ...
                   'Cep',ones(3,2),'d',3), t, struct()), 'lsqvce:covModel:badCep'); % d stated twice
end

function testCovModelRejectsBadSigma(tc)
% categories: contract
    t = (0:99).'*5;
    tc.verifyError(@() lsqvce.engine.covModel(struct('type','diag','sigma',0), t, struct()), ...
                   'lsqvce:covModel:badSigma');
    tc.verifyError(@() lsqvce.engine.covModel(struct('type','diag','sigma',-1), t, struct()), ...
                   'lsqvce:covModel:badSigma');
    tc.verifyError(@() lsqvce.engine.covModel(struct('type','diag','sigma',ones(5,1)), ...
                                           t, struct()), ...
                   'lsqvce:covModel:sigmaLength');
end

function testCovModelRejectsImpossibleWhitenMode(tc)
% categories: contract
% 'acf' and 'custom' have no filter form; asking for one must fail loudly
% rather than quietly falling back to the dense path.
    t  = (0:99).'*5;
    cm = struct('type','acf', 'acf', @(dt,o) exp(-abs(dt)/50), 'sigma',1);
    tc.verifyError(@() lsqvce.engine.covModel(cm, t, struct('whitenMode','filter')), ...
                   'lsqvce:covModel:noFilterForm');
    tc.verifyError(@() lsqvce.engine.covModel(cm, t, struct('whitenMode','nonsense')), ...
                   'lsqvce:covModel:badWhitenMode');
end

function testCovModelWarnsOnNegativePSD(tc)
% categories: contract
% A raw periodogram difference can go negative. It is not physical, and it
% breaks the Levinson recursion downstream.
    t  = (0:99).'*5;
    f  = linspace(0, 0.1, 200).';
    S  = ones(200,1);   S(50:60) = -0.1;
    cm = struct('type','psd', 'f',f, 'S',S);
    tc.verifyWarning(@() lsqvce.engine.covModel(cm, t, struct('arOrder',5)), ...
                     'lsqvce:covModel:negativePSD');
end

function testCovModelFilterMemoryIsReported(tc)
% categories: robustness
% The accumulator relies on filterMemory to set the batch overlap. If it comes
% back zero for a filter model, batching would silently corrupt the answer.
    t  = (0:999).'*5;
    cm = struct('type','arma', 'ar',[1 -1.6 0.9 -0.25], 'ma',1, 'sigma',0.7);
    [~, Wf, info] = lsqvce.engine.covModel(cm, t, struct());
    tc.verifyNotEmpty(Wf);
    tc.verifyEqual(info.filterMemory, 3);

    cmD = struct('type','diag','sigma',1);
    [~, ~, infoD] = lsqvce.engine.covModel(cmD, t, struct());
    tc.verifyEqual(infoD.filterMemory, 0);      % white noise needs no overlap
end

function testCovModelMaPartExtendsTheFilterMemory(tc)
% categories: contract, precision
% WP-125: with an MA part the whitener 1/B(z) is IIR; the reported memory is
% p plus the samples after which its impulse response is below eps of its
% peak, and a non-invertible MA part is refused by name.
    t  = (0:999).'*5;
    b  = [1 0.8];
    cm = struct('type','arma', 'ar',[1 -0.5], 'ma',b, 'sigma',0.7);
    [~, ~, info] = lsqvce.engine.covModel(cm, t, struct());
    N  = info.filterMemory - 1;                 % p = 1
    tc.verifyGreaterThan(N, 1, 'the MA part adds to the AR order');
    h  = filter(1, b, [1; zeros(N - 1, 1)]);
    % tol: the definition of the memory -- the tail of 1/B(z) below eps of its peak
    tc.verifyLessThanOrEqual(abs(h(N)), eps * max(abs(h)));
    tc.verifyError(@() lsqvce.engine.covModel(setfield(cm, 'ma', [1 1.2]), t, struct()), ...
        'lsqvce:covModel:maNotInvertible');                         %#ok<SFLD>
end

function testCovModelRefusesADenseArcAboveTheLimit(tc)
% categories: contract, robustness
% WP-126: a law without a filter form is factorised over the whole arc
% (WP-125), so its length is bounded by the declared cm.maxDenseArc (default
% 20 000 epochs, the batch the dense factor never exceeded before). The
% boundary is exact: m == limit runs, m == limit + 1 is refused by name; a
% filter law is not bounded; a limit that is not a positive number is refused.
    cm = struct('type','acf', 'acf',@(dt,o) exp(-abs(dt)/50), 'sigma',1);
    t  = (0:99).'*5;                                            % m = 100
    [~, ~, info] = lsqvce.engine.covModel(setfield(cm, 'maxDenseArc', 100), t, struct()); %#ok<SFLD>
    tc.verifyTrue(info.wholeArc, 'a dense law is whitened whole');
    tc.verifyError(@() lsqvce.engine.covModel(setfield(cm, 'maxDenseArc', 99), t, struct()), ...
        'lsqvce:covModel:denseArcTooLong');                         %#ok<SFLD>
    tc.verifyError(@() lsqvce.engine.covModel(cm, (0:20000).'*5, struct()), ...
        'lsqvce:covModel:denseArcTooLong');                         % the default, 20 000
    tc.verifyError(@() lsqvce.engine.covModel(setfield(cm, 'maxDenseArc', 0), t, struct()), ...
        'lsqvce:covModel:denseArcTooLong');                         %#ok<SFLD>
    ar = struct('type','arma', 'ar',[1 -0.5], 'sigma',0.7, 'maxDenseArc',10);
    [~, W, ~] = lsqvce.engine.covModel(ar, t, struct());
    tc.verifyNotEmpty(W, 'a filter law is not bounded by the dense limit');
end

%% ============================================ input validation: accumulator
function testAccumulatorRejectsBadConstruction(tc)
% categories: contract
    % WP-39 c2: the ctor's scalar constraint and enum checks moved into
    % the typed arguments signature -- expected ids follow the migration
    % (the b1-r3 rule); the CondWarn<=CondAbort cross-check stays a body
    % check with its lsqvce: id.
    tc.verifyError(@() lsqvce.engine.BlockAccumulator(0),  'MATLAB:validators:mustBePositive');
    tc.verifyError(@() lsqvce.engine.BlockAccumulator(-3), 'MATLAB:validators:mustBePositive');
    tc.verifyError(@() lsqvce.engine.BlockAccumulator(5, 'Backend','bogus'), ...
                   'MATLAB:validators:mustBeMember');
    tc.verifyError(@() lsqvce.engine.BlockAccumulator(5, 'Reg','bogus'), ...
                   'MATLAB:validators:mustBeMember');
    tc.verifyError(@() lsqvce.engine.BlockAccumulator(5, 'CondWarn',1e15, 'CondAbort',1e10), ...
                   'lsqvce:accum:badThresholds');
end

function testAccumulatorRejectsBadSystem(tc)
% categories: contract
    cfg = tc.TestData.cfg;
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Verbose',false);

    tc.verifyError(@() acc.addSystem(1, cfg.nl, 'notafunction', 1), ...
                   'lsqvce:accum:badLoader');
    tc.verifyError(@() acc.addSystem(1, cfg.nl, binder(cfg,1), []), ...
                   'lsqvce:accum:noDays');

    acc.addSystem(1, cfg.nl, binder(cfg,1), 1);
    tc.verifyError(@() acc.addSystem(1, cfg.nl, binder(cfg,1), 1), ...
                   'lsqvce:accum:duplicateSystem');
end

function testAccumulatorDetectsWrongArcDimensions(tc)
% categories: contract
% The loader is the most likely place for a bug. Every dimension is checked,
% with an error message that names the day, obs type and arc.
    cfg = tc.TestData.cfg;

    function out = badNc(d, ot, a, uo)
        out = lsqvce.engine.synthLoader(cfg, 1, d, ot, a, uo);
        if ~ischar(d), out.A = out.A(:, 1:end-1); end     % wrong nc
    end
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Verbose',false);
    tc.verifyError(@() acc.addSystem(1, cfg.nl, @badNc, 1), 'lsqvce:accum:arcBadNc');

    function out = badRows(d, ot, a, uo)
        out = lsqvce.engine.synthLoader(cfg, 1, d, ot, a, uo);
        if ~ischar(d), out.l = out.l(1:end-5); end        % A and l disagree
    end
    acc2 = lsqvce.engine.BlockAccumulator(cfg.nc, 'Verbose',false);
    tc.verifyError(@() acc2.addSystem(1, cfg.nl, @badRows, 1), ...
                   'lsqvce:accum:arcRowMismatch');
end

function testAccumulatorRejectsNaNObservations(tc)
% categories: contract, robustness
% NaN in l would propagate into the normal matrix and poison everything. The
% loader must remove bad epochs; the accumulator cannot know which are valid.
    cfg = tc.TestData.cfg;
    function out = withNaN(d, ot, a, uo)
        out = lsqvce.engine.synthLoader(cfg, 1, d, ot, a, uo);
        if ~ischar(d), out.l(17) = NaN; end
    end
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Verbose',false);
    tc.verifyError(@() acc.addSystem(1, cfg.nl, @withNaN, 1), ...
                   'lsqvce:accum:arcNotFinite');
end

function testAccumulatorRejectsBadMeta(tc)
% categories: contract, robustness
    cfg = tc.TestData.cfg;
    function out = badMeta(d, ot, a, uo)                  %#ok<INUSD>
        if ischar(d), out = struct('obsTypes', {{'a','b'}}, 'arcs', {{1}});
        else,         out = lsqvce.engine.synthLoader(cfg, 1, d, ot, a, uo);
        end
    end
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Verbose',false);
    tc.verifyError(@() acc.addSystem(1, cfg.nl, @badMeta, 1), ...
                   'lsqvce:accum:metaMismatch');
end

%% ============================================ degenerate data
function testEmptyArcsAreSkipped(tc)
% categories: robustness
% Data gaps are normal. An empty arc must be skipped silently, not crash.
    cfg = tc.TestData.cfg;
    function out = gappy(d, ot, a, uo)
        out = lsqvce.engine.synthLoader(cfg, 1, d, ot, a, uo);
        if ~ischar(d) && a == 2
            out.A = zeros(0, cfg.nc);
            out.B = zeros(0, cfg.nl);
            out.l = [];
        end
    end
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Verbose',false);
    acc.addSystem(1, cfg.nl, @gappy, 1:2);
    sol = acc.solve();
    tc.verifyTrue(all(isfinite(sol.x)));
end

function testAllArcsEmptyIsAnError(tc)
% categories: contract, robustness
% A system that contributes nothing is a bug in the loader, not a valid state.
    cfg = tc.TestData.cfg;
    function out = allEmpty(d, ot, a, uo)
        out = lsqvce.engine.synthLoader(cfg, 1, d, ot, a, uo);
        if ~ischar(d)
            out.A = zeros(0, cfg.nc);  out.B = zeros(0, cfg.nl);  out.l = [];
        end
    end
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Verbose',false);
    tc.verifyError(@() acc.addSystem(1, cfg.nl, @allEmpty, 1:2), 'lsqvce:accum:noArcs');
end

function testNoRedundancyIsAnError(tc)
% categories: contract
% Fewer observations than parameters: sigma0 is undefined. Say so.
    cfg = tc.TestData.cfg;
    cfg.mArc = 2;   cfg.arcsPerDay = 1;   cfg.obsTypes = {'kbr'};
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Verbose',false);
    tc.verifyError(@() acc.addSystem(1, cfg.nl, binder(cfg,1), 1), ...
                   'lsqvce:accum:localNotRedundant');
end

function testSolveWithoutSystemsIsAnError(tc)
% categories: contract, robustness
    acc = lsqvce.engine.BlockAccumulator(5, 'Verbose',false);
    tc.verifyError(@() acc.solve(), 'lsqvce:accum:noSystems');
end

%% ============================================ rank deficiency
function testExactRankDeficiencyIsRescued(tc)
% categories: robustness
% Two identical local columns -> Nll is exactly singular. With 'tsvd' the run
% must complete; the pinv fallback must not be needed.
    muteWarn(tc, 'lsqvce:accum:pinvFallback');   % provoked on purpose here (WP-47, gate rule 7)
    cfg = tc.TestData.cfg;
    cfg.collinear = true;   cfg.collinearEps = 0;

    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','qr', 'Reg','tsvd', 'Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), 1:2); end
    sol = acc.solve();

    tc.verifyTrue(all(isfinite(sol.x)));
    tc.verifyGreaterThanOrEqual(acc.diag(1).rankDef, 1);
    tc.verifyEqual(acc.diag(1).regUsed, 'tsvd');
end

function testPinvFallbackNeverCrashes(tc)
% categories: contract, robustness
% With Reg = 'none' and a singular local block, pre-elimination must fall back
% to pinv and warn -- never abort. A day of processing should not be lost to
% one bad arc.
    cfg = tc.TestData.cfg;
    cfg.collinear = true;   cfg.collinearEps = 0;

    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','qr', 'Reg','none', 'Verbose',false);
    tc.verifyWarning(@() acc.addSystem(1, cfg.nl, binder(cfg,1), 1:2), ...
                     'lsqvce:accum:pinvFallback');
    tc.verifyEqual(acc.diag(1).regUsed, 'pinv-fallback');
end

function testNullDirectionIdentifiesTheCulprit(tc)
% categories: precision
% The diagnostic that actually helps: which combination of local parameters is
% unconstrained. For two collinear columns it must be their difference.
    muteWarn(tc, 'lsqvce:accum:pinvFallback');   % provoked on purpose here (WP-47, gate rule 7)
    cfg = tc.TestData.cfg;
    cfg.collinear = true;   cfg.collinearEps = 1e-12;

    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','qr', 'Reg','tsvd', 'Verbose',false);
    acc.addSystem(1, cfg.nl, binder(cfg,1), 1:2);

    v = acc.diag(1).nullDir;
    v = v / max(abs(v));

    % the last two local columns are collinear, so the null direction must be
    % concentrated on them, with opposite signs
    tc.verifyGreaterThan(abs(v(end)) + abs(v(end-1)), 1.5);
    tc.verifyLessThan(abs(v(end) + v(end-1)), 0.3);
    tc.verifyLessThan(max(abs(v(1:end-2))), 0.3);
end

function testDiagnosticsFlagIllConditioning(tc)
% categories: robustness
    muteWarn(tc, 'lsqvce:accum:pinvFallback');   % provoked on purpose here (WP-47, gate rule 7)
    cfg = tc.TestData.cfg;
    cfg.collinear = true;   cfg.collinearEps = 1e-9;

    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), 1:2); end

    tc.verifyGreaterThan(acc.diag(1).condNll, 1e8);
    tc.verifyNotEqual(acc.diag(1).status, 'ok');
end

%% ============================================ regularisation branches
function testAutoRegDispatchesOnStatus(tc)
% categories: robustness
    muteWarn(tc, 'lsqvce:accum:pinvFallback');   % provoked on purpose here (WP-47, gate rule 7)
    cfg = tc.TestData.cfg;

    a1 = lsqvce.engine.BlockAccumulator(cfg.nc, 'Reg','auto', 'Verbose',false);
    a1.addSystem(1, cfg.nl, binder(cfg,1), 1:2);
    tc.verifyEqual(a1.diag(1).regUsed, 'none');          % well conditioned

    cfg.collinear = true;   cfg.collinearEps = 1e-11;
    a2 = lsqvce.engine.BlockAccumulator(cfg.nc, 'Reg','auto', 'Verbose',false);
    a2.addSystem(1, cfg.nl, binder(cfg,1), 1:2);
    tc.verifyNotEqual(a2.diag(1).regUsed, 'none');       % ill conditioned
end

function testAllRegularisationBranchesRun(tc)
% categories: robustness
    muteWarn(tc, 'lsqvce:accum:pinvFallback');   % provoked on purpose here (WP-47, gate rule 7)
    cfg = tc.TestData.cfg;
    cfg.collinear = true;   cfg.collinearEps = 1e-10;

    for r = {'none','tikhonov','tsvd','lm','auto'}
        acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','qr', 'Reg',r{1}, ...
                              'RegLambda',1e-4, 'Verbose',false);
        acc.addSystem(1, cfg.nl, binder(cfg,1), 1:2);
        sol = acc.solve();
        tc.verifyTrue(all(isfinite(sol.x)), sprintf('Reg = %s produced NaN', r{1}));
    end
end

function testRegMatrixDimensionIsChecked(tc)
% categories: contract
    Nll = eye(4);   nl = ones(4,1);
    dg  = lsqvce.engine.localDiagnostics(Nll, struct('TruncTol',1e-10));
    opt = struct('Reg','tikhonov', 'RegLambda',1, 'RegMatrix',eye(3), ...
                 'TruncTol',1e-10);
    tc.verifyError(@() lsqvce.engine.applyRegularisation(Nll, nl, dg, opt), ...
                   'lsqvce:applyReg:badRegMatrix');
end

function testGCVWarnsAtGridBoundary(tc)
% categories: contract
% GCV under-regularises when the residuals are still correlated. It must say so
% rather than silently returning lambda ~ 0.
    N = diag([1 1 1 1]);          % perfectly conditioned -> no regularisation wanted
    n = ones(4,1);
    tc.verifyWarning(@() lsqvce.engine.gcvLambda(N, n, eye(4)), 'lsqvce:gcv:atLowerBound');
end

%% ============================================ loader contract
function testUserOptIsForwardedEverywhere(tc)
% categories: contract
    cfg  = tc.TestData.cfg;
    seen = 0;

    function out = probe(d, ot, a, uo)
        tc.verifyTrue(isfield(uo,'myFlag') && uo.myFlag == 77);
        tc.verifyTrue(isfield(uo,'sysID'));
        tc.verifyTrue(isfield(uo,'nc') && uo.nc == cfg.nc);
        tc.verifyTrue(isfield(uo,'nl') && uo.nl == cfg.nl);
        if ~ischar(d)
            tc.verifyTrue(isfield(uo,'day') && isfield(uo,'arc'));
            tc.verifyTrue(isfield(uo,'obsType'));
        end
        seen = seen + 1;
        out  = lsqvce.engine.synthLoader(cfg, 1, d, ot, a, uo);
    end

    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Verbose',false, 'UserOpt',struct('myFlag',77));
    acc.addSystem(1, cfg.nl, @probe, 1);
    tc.verifyGreaterThan(seen, 1);
end

function testNonDeterministicLoaderIsDetected(tc)
% categories: contract
% The residual pass re-runs the loader. If it returns different data, sigma0
% will not reproduce -- and that MUST be reported, because silently wrong
% residuals are worse than none.
    cfg = tc.TestData.cfg;
    call = 0;

    function out = drifting(d, ot, a, uo)
        out = lsqvce.engine.synthLoader(cfg, 1, d, ot, a, uo);
        if ~ischar(d)
            call = call + 1;
            if call > 4                     % perturb only on the second pass
                out.l = out.l + 0.5;
            end
        end
    end

    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Verbose',false);
    acc.addSystem(1, cfg.nl, @drifting, 1);
    sol = acc.solve();

    tc.verifyWarning(@() acc.residuals(sol, 'Verbose',false), ...
                     'lsqvce:residuals:sigma0Mismatch');
end

function testKeepLocalFalseBlocksBackSubstitution(tc)
% categories: contract, robustness
    cfg = tc.TestData.cfg;
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'KeepLocal',false, 'Verbose',false);
    acc.addSystem(1, cfg.nl, binder(cfg,1), 1);
    sol = acc.solve();
    tc.verifyTrue(all(isfinite(sol.x)));                    % solving still works
    tc.verifyError(@() acc.backSubstitute(1, sol.x), 'lsqvce:accum:noStore');
end

function testUnknownSystemIsAnError(tc)
% categories: contract
    cfg = tc.TestData.cfg;
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Verbose',false);
    acc.addSystem(1, cfg.nl, binder(cfg,1), 1);
    sol = acc.solve();
    tc.verifyError(@() acc.backSubstitute(99, sol.x), 'lsqvce:accum:unknownSystem');
end

%% ============================================ reproducibility
function testResultIsIndependentOfSystemOrder(tc)
% categories: precision
% Systems are independent. Adding them in a different order must change only
% the floating-point rounding.
    cfg = tc.TestData.cfg;

    a1 = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false);
    for s = [1 2 3], a1.addSystem(s, cfg.nl, binder(cfg,s), 1:2); end
    s1 = a1.solve();

    a2 = lsqvce.engine.BlockAccumulator(cfg.nc, 'Backend','qr', 'Verbose',false);
    for s = [3 1 2], a2.addSystem(s, cfg.nl, binder(cfg,s), 1:2); end
    s2 = a2.solve();

    % tol: two backends on one well-conditioned problem; difference ~ cond*eps, bound with release-drift headroom
    tc.verifyEqual(s2.x,      s1.x,      'RelTol', 1e-10);
    tc.verifyEqual(s2.sigma0, s1.sigma0, 'RelTol', 1e-10);
end

function testReportTableIsWellFormed(tc)
% categories: contract | tags: diagnostics
    cfg = tc.TestData.cfg;
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Verbose',false);
    for s = cfg.sats, acc.addSystem(s, cfg.nl, binder(cfg,s), 1:2); end
    tbl = acc.report();
    tc.verifyEqual(height(tbl), numel(cfg.sats));
    tc.verifyTrue(all(ismember({'sysID','condNll','status','regUsed'}, ...
                               tbl.Properties.VariableNames)));
end

%% ============================================ helpers
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
    f = @(d,ot,a,uo) lsqvce.engine.synthLoader(cfg, s, d, ot, a, uo);
end

%% ===================================================== non-finite propagation
function testSignalKernelsRejectNonFiniteAtTheDoor(tc)
% categories: robustness
% WP-43 (T-04): non-finite handling was asserted only on the accumulator
% observation path. For the signal kernels the contract is REJECTION AT THE
% DOOR (validateattributes 'finite') -- a NaN PSD sample or ACV value must
% never yield a silently finite, wrong autocovariance or whitener. This test
% pins that guard so a refactor cannot drop it unnoticed.
    f = linspace(0, 0.1, 64).';
    S = ones(64, 1);

    Snan = S;  Snan(10) = NaN;
    tc.verifyError(@() lsqvce.kernel.psd2acv(f, Snan, 5, 8), ...
                   'MATLAB:psd2acv:expectedFinite');
    fInf = f;  fInf(end) = Inf;
    tc.verifyError(@() lsqvce.kernel.psd2acv(fInf, S, 5, 8), ...
                   'MATLAB:psd2acv:expectedFinite');

    r = 0.5.^(0:19);
    rNan = r;  rNan(3) = NaN;
    tc.verifyError(@() lsqvce.kernel.levinsonDurbin(rNan, 4), ...
                   'MATLAB:levinsonDurbin:expectedFinite');

    tc.verifyError(@() lsqvce.kernel.armaACV([1 NaN], 1, 8), ...
                   'MATLAB:armaACV:expectedFinite');
    tc.verifyError(@() lsqvce.kernel.armaACV([1 -0.5], [1 Inf], 8), ...
                   'MATLAB:armaACV:expectedFinite');
end
