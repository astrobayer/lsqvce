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
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-138 (code review of WP-135): the
%   accumulator judges id kinds before duplicates (67 after 'C'); the GCV
%   errors under lsqvce:gcv:*. With assistance from Claude (Anthropic), model
%   Claude Opus 5.5.
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-137 (code review of WP-135): the
%   GCV test on deterministic data, no muted id, its tolerance at its source;
%   snoop's re-solve keeps the solve's Reg. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-14 (CEST, UTC+02:00), WP-135 (third review findings 11,
%   12, 16 and lows): gcvLambda is the data's GCV (against direct solves), a
%   null-space L, its refusals; the accumulator refuses repeated days and
%   mixed ids; static arc ids are checked; the template refuses robust
%   factors; a zero local block is wholly rank deficient. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-134 (third review): a filter law's dense form is guarded; per-row sigma, uneven grid.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-130: a rank-deficient local
%   block reduced on the factor (x, dof, sigma0, B y, Cyy, the redundancy
%   identity, the downdate). With assistance from Claude (Anthropic), model
%   Claude Opus 5.5.
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

function testPSD2ACVRefusesAGridWithoutZero(tc)
% categories: contract, precision
% WP-153 (review 12.84, C1). The integral runs over [f(1), f(end)]: a grid
% that starts above 0 left the band [0, f(1)] out of r -- an AR(0.99) law on
% 512 bins read r(0) 33.4 instead of 50.25 (-34 %), silently. It is refused
% when the band, estimated flat as f(1) S(1), exceeds 1e-6 of r(0); the
% habitual start just above 0 passes when its band is negligible. The
% lags are formed in blocks (memory O(n + nf)): a long call returns the same
% leading lags as a short one.
    dt = 5;  phi = 0.99;
    f = linspace(0, 0.5/dt, 513).';
    S = dt * 2 ./ abs(1 - phi * exp(-1i*2*pi*f*dt)).^2;
    tc.verifyError(@() lsqvce.kernel.psd2acv(f(2:end), S(2:end), dt, 4), 'lsqvce:psd2acv:noZeroFrequency');
    tc.verifyError(@() lsqvce.engine.covModel(struct('type', 'psd', 'f', f(2:end), 'S', S(2:end)), ...
                   (0:99).' * dt, struct('psdMode', 'toeplitz')), 'lsqvce:psd2acv:noZeroFrequency');
    f9 = [1e-12; f(2:end)];                         % a negligible band below the grid
    % tol: the omitted band is 1e-12 Hz x S(0) = 2e-9 of r(0) at phi 0.99 -- r(0) to 1e-5
    tc.verifyEqual(lsqvce.kernel.psd2acv(f9, S, dt, 1), lsqvce.kernel.psd2acv(f, S, dt, 1), 'RelTol', 1e-5, ...
        'a negligible band below the grid passes');
    r = lsqvce.kernel.psd2acv(f, S, dt, 2e4);       % 2e4 x 512 products: five blocks
    % tol: 1e-3 relative -- the piecewise-linear interpolant of the AR(0.99)
    % spectrum on 513 bins against the closed form 1 / (1 - phi^2)
    tc.verifyEqual(r(1), 1 / (1 - phi^2), 'RelTol', 1e-3);
    % tol: exact -- the same products in the same order, blocked or not
    tc.verifyEqual(r(1:100), lsqvce.kernel.psd2acv(f, S, dt, 100), 'AbsTol', 0, ...
        'the blocked lags equal the unblocked ones');
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

function testCofactorOfAnIllConditionedLocalBlockIsExact(tc)
% categories: precision, robustness
% WP-128: the QR backend eliminates a local block exactly on its factor up to
% cond(B) ~ 4.5e12; the diagnostics cofactor was rebuilt from the SQUARED
% block Nll = B'B and, from cond(B) ~ 7e7, cut its eigenvalues at 1e-10 max --
% one redundancy unit lost per cut direction while the solve kept them all.
% A quadratic local polynomial in raw seconds over 1e4 s (cond(B) = 1.3e8) is
% enough (dev/validation/wp128_local_block_mirror.py: sum h 6 instead of 7).
% Integrity identity: sum(r_i) = m - u, through acc.redundancy.
    rs = RandStream('mt19937ar', 'Seed', 128);
    m = 200;  nc = 4;  t = linspace(0, 1e4, m).';
    A = randn(rs, m, nc);  B = [ones(m, 1), t, t.^2];
    l = A*randn(rs, nc, 1) + B*[1; 1e-3; 1e-7] + 1e-3*randn(rs, m, 1);
    function out = poly(d, ot, a, uo)                          %#ok<INUSD>
        if ischar(d) || isstring(d), out = struct('obsTypes', {{'x'}}, 'arcs', {{1}}); return, end
        out = struct('A', A, 'B', B, 'l', l);
    end
    acc = lsqvce.engine.BlockAccumulator(nc, 'Backend', 'qr', 'KeepLocal', true, 'Verbose', false);
    acc.addSystem(1, 3, @poly, 1);
    sol = acc.solve();
    red = acc.redundancy(sol, 'Verbose', false);
    % tol: 1e-8 -- the identity holds to rounding on the factor; the cut lost a whole unit
    tc.verifyLessThan(red.check, 1e-8, 'sum(r_i) = m - u with the exact local cofactor');
    % and the local covariance from the factor equals the dense reference
    Cyy = acc.localCovariance(1, sol.Cxx, sol.s0sq);
    [~, Rf] = qr([B, A], 0);  Ri = Rf \ eye(nc + 3);  Qf = Ri * Ri.';   % the reference, on a factor
    % tol: 1e-6 relative -- two factorisations of a design at cond 1.3e8:
    % eps * cond ~ 3e-8 on these entries, with headroom
    tc.verifyEqual(diag(Cyy), sol.s0sq * diag(Qf(1:3, 1:3)), 'RelTol', 1e-6);
end

function testRankDeficientLocalBlockIsReducedOnTheFactor(tc)
% categories: precision, robustness
% WP-130: one duplicated local column (B = [1 t t^2 1], t over 1e4 s) fails
% R11's pivots, and the fallback squared the block: pinv(Nll, 1e-10 norm) cut
% the weak but DETERMINED t^2 direction with the null one (eig(Nll)/max = 2e-31,
% 1e-16, 1e-9, 1), dof counted mObs - nl, and y / Cyy were solved with
% backslash on the singular Nll_r, not with the pinv: x 31 formal sd off,
% sigma0 35x (dev/validation/wp130_pinv_fallback_mirror.py). Now the rank is
% decided on the factor (SVD of R11, QRSystem.RankTol) and every reader of the
% block uses the solve's pseudo-inverse. Reference: the design without the
% duplicate, whose locals are the same estimable functions.
    muteWarn(tc, 'lsqvce:accum:pinvFallback');   % provoked on purpose here (WP-47, gate rule 7)
    rs = RandStream('mt19937ar', 'Seed', 130);
    m = 200;  nc = 4;  t = linspace(0, 1e4, m).';
    A = randn(rs, m, nc);  B = [ones(m, 1), t, t.^2, ones(m, 1)];
    l = A*randn(rs, nc, 1) + B(:, 1:3)*[1; 1e-3; 1e-7] + 1e-2*randn(rs, m, 1);
    function out = dupPoly(d, ot, a, uo)                          %#ok<INUSD>
        if ischar(d) || isstring(d), out = struct('obsTypes', {{'x'}}, 'arcs', {{1}}); return, end
        out = struct('A', A, 'B', B, 'l', l);
    end
    acc = lsqvce.engine.BlockAccumulator(nc, 'Backend', 'qr', 'KeepLocal', true, 'Verbose', false);
    acc.addSystem(1, 4, @dupPoly, 1);
    tc.assertEqual(acc.diag(1).regUsed, 'pinv-fallback', 'the fixture reaches the fallback');
    sol = acc.solve();
    Dl = [A, B(:, 1:3)];  [~, Rf] = qr(Dl, 0);  xr = Rf \ (Rf.' \ (Dl.' * l));
    vr = l - Dl*xr;  dofR = m - nc - 3;  Ri = Rf \ eye(nc + 3);  Cr = (vr.'*vr / dofR) * (Ri * Ri.');
    tc.verifyEqual(sol.dof, dofR, 'dof counts the rank, not nl');
    % tol: 1e-6 formal sd -- two factorisations of a design at cond 1.3e8
    % (eps * cond ~ 3e-8 relative); the cut route was 31 sd off
    tc.verifyLessThan(max(abs(sol.x - xr(1:nc)) ./ sqrt(diag(Cr(1:nc, 1:nc)))), 1e-6, 'x');
    % tol: 1e-8 relative -- the same residual sum over the same dof
    tc.verifyEqual(sol.s0sq, vr.'*vr / dofR, 'RelTol', 1e-8, 'sigma0^2');
    y = acc.backSubstitute(1, sol.x);
    % tol: 1e-8 relative -- the estimable local function B*y, two routes
    tc.verifyEqual(B*y, B(:, 1:3)*xr(nc+1:end), 'RelTol', 1e-8, 'B y is the estimable local signal');
    % tol: 1e-10 absolute -- the minimum-norm y splits the duplicate equally
    tc.verifyEqual(y(1), y(4), 'AbsTol', 1e-10 * max(abs(y)), 'minimum-norm y');
    Cyy = acc.localCovariance(1, sol.Cxx, sol.s0sq);
    % tol: 1e-6 relative -- the covariance of the estimable local signal B y
    % from the kept pseudo-inverse vs the reference's (was 2.8e14-scale garbage)
    tc.verifyEqual(diag(B*Cyy*B.'), diag(B(:, 1:3)*Cr(nc+1:end, nc+1:end)*B(:, 1:3).'), 'RelTol', 1e-6);
    red = acc.redundancy(sol, 'Verbose', false);
    % tol: 1e-8 -- sum(r_i) = m - rank with the solve's pseudo-inverse (identityViolated before)
    tc.verifyLessThan(red.check, 1e-8, 'the redundancy identity counts the rank');
    % snoop's downdate on the kept basis: sum h over the survivors = the rank
    Q = lsqvce.engine.BlockCofactor.fromAccumulator(acc);
    [Q, ok] = Q.downdate([A(57, :), B(57, :)], 1);
    tc.assertTrue(ok);
    keep = (1:m).' ~= 57;
    % tol: 1e-8 -- the identity on the downdated factor
    tc.verifyEqual(sum(Q.leverage(1, A(keep, :), B(keep, :))), nc + 3, 'AbsTol', 1e-8, ...
        'sum h over the survivors = nc + rank after a downdate');
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
    % WP-127: the design too -- a NaN in A of a system WITHOUT locals passed
    % (the only check sat in the local block) and came out as x = NaN
    rs = RandStream('mt19937ar', 'Seed', 41);
    A  = randn(rs, 30, 3);  A(7, 2) = NaN;
    tc.verifyError(@() lsqvce.fit(struct('A', A, 'l', randn(rs, 30, 1)), Verbose=false), ...
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
    tc.verifyError(@() lsqvce.engine.applyRegularisation(Nll, nl, dg, opt, 10, 4), ...
                   'lsqvce:applyReg:badRegMatrix');
end

function testAccumulatorRefusesMixedIdsAndRepeatedDays(tc)
% categories: contract
    % WP-135 (review 3 findings 11, 12): a repeated day was accumulated twice,
    % silently; a character id next to a number made [store.sysID] a char row
    % (3 became char(3), robust factors keyed under it were never read).
    cfg = tc.TestData.cfg;
    acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'Verbose', false);
    tc.verifyError(@() acc.addSystem(1, cfg.nl, binder(cfg, 1), [1 2 2]), 'lsqvce:accum:duplicateDay');
    acc.addSystem('C', cfg.nl, binder(cfg, 1), 1:2);
    tc.verifyError(@() acc.addSystem(3, cfg.nl, binder(cfg, 1), 1:2), 'lsqvce:accum:mixedIdTypes');
    specs = struct('sysID', 3, 'nl', cfg.nl, 'dataFcn', binder(cfg, 1), 'dayList', 1:2);
    tc.verifyError(@() acc.addSystems(specs), 'lsqvce:accum:mixedIdTypes');
    % WP-138: the kinds are judged before the duplicate test -- 67 is char(67)
    % = 'C', and [obj.diag.sysID] == 67 / [specs.sysID] called it a duplicate
    tc.verifyError(@() acc.addSystem(67, cfg.nl, binder(cfg, 1), 1:2), 'lsqvce:accum:mixedIdTypes');
    specs.sysID = 67;
    tc.verifyError(@() acc.addSystems(specs), 'lsqvce:accum:mixedIdTypes');
end

function testStaticArcIdsAreChecked(tc)
% categories: contract
    % WP-135 (review 3, lane 4): arcs = [1 1 2] read cell 1 twice and cell 2
    % never (find(ids == arcID, 1)), silently -- a resolver's list was checked.
    tc.verifyError(@() staticArcsMeta([1 1 2]), 'lsqvce:makeLoader:badArcs');
    m = staticArcsMeta([3 1 2]);                      % unique ids, any order
    tc.verifyEqual(m.arcs{1}, [3 1 2]);
end

function m = staticArcsMeta(ids)
    L = lsqvce.engine.makeLoader(@(d, ot) error('test:reached', 'no read'), ...
            'ObservationTypes', {'kbr'}, 'Covariance', struct('kbr', struct('type', 'diag', 'sigma', 1)), ...
            'Arcs', ids, 'NumLocal', 0, 'NumCommon', 2, 'Days', 1);
    m = L.meta(1);
end

function testRobustFactorsAreRefusedByTheTemplate(tc)
% categories: contract
    % WP-135 (review 3, lane 4): synthLoader -- the loader template -- applied
    % no robust factors while the driver reported them applied.
    cfg = tc.TestData.cfg;
    key = lsqvce.util.rowWeightKey(1, 1, cfg.obsTypes{1}, 1);
    opt = struct('sysID', 1, 'rowWeights', containers.Map({key}, {ones(cfg.mArc, 1)}));
    tc.verifyError(@() lsqvce.engine.synthLoader(cfg, 1, 1, cfg.obsTypes{1}, 1, opt), ...
                   'lsqvce:synthLoader:robustFactors');
    opt.rowWeights = containers.Map({'1|99|x|1'}, {1});           % another arc's: not this one's
    out = lsqvce.engine.synthLoader(cfg, 1, 1, cfg.obsTypes{1}, 1, opt);
    tc.verifyNumElements(out.l, cfg.mArc);
end

function testAZeroLocalBlockIsWhollyRankDeficient(tc)
% categories: robustness
    % WP-135 (review 3, lane 1): an all-zero block read rankDef 0 (ev < tol * 0
    % is never true) while its status said rank-deficient.
    dg = lsqvce.engine.localDiagnostics(zeros(3), struct('TruncTol', 1e-10));
    tc.verifyEqual(dg.rankDef, 3);
    tc.verifyEqual(dg.status, 'rank-deficient');
end

function testGCVWarnsAtGridBoundary(tc)
% categories: contract
% GCV under-regularises when the residuals are still correlated. It must say so
% rather than silently returning lambda ~ 0. WP-135: the functional is the
% data's -- an exactly fitting block (l'l = n'N^-1 n, no residual) wants no
% regularisation; a block that explains nothing (n = 0) wants the most.
    N = diag([1 1 1 1]);  n = ones(4,1);
    tc.verifyWarning(@() lsqvce.engine.gcvLambda(N, n, eye(4), 10, 4), 'lsqvce:gcv:atLowerBound');
    tc.verifyWarning(@() lsqvce.engine.gcvLambda(N, zeros(4,1), eye(4), 10, 4), ...
                     'lsqvce:gcv:atUpperBound');
end

function testGcvLambdaRefusesWhatItCannotChoose(tc)
% categories: contract, robustness
% WP-135: m and l'l are the data's (m > nl, l'l >= 0); an L that penalises
% nothing warns and returns 0; a direction L leaves free that the data do not
% see is refused -- no lambda regularises it.
    N = eye(4);  n = ones(4,1);
    tc.verifyError(@() lsqvce.engine.gcvLambda(N, n, [], 4, 4), 'lsqvce:gcv:badData');
    tc.verifyError(@() lsqvce.engine.gcvLambda(N, n, [], 10, -1), 'lsqvce:gcv:badData');
    lam = tc.verifyWarning(@() lsqvce.engine.gcvLambda(N, n, zeros(4), 10, 4), 'lsqvce:gcv:degenerate');
    tc.verifyEqual(lam, 0);
    tc.verifyError(@() lsqvce.engine.gcvLambda(diag([1 1 1 0]), [1;1;1;0], [eye(3), zeros(3,1)], 10, 4), ...
                   'lsqvce:gcv:nullUnobserved');
end

function testGcvLambdaIsTheGcvOfTheData(tc)
% categories: precision, reference
% WP-135 (review 3 finding 16): the curve is Golub-Heath-Wahba's
% ||B y - l||^2 / (m - tr H)^2 of the local problem -- the code minimised
% ||f.*b||^2/(sum f)^2 of N y = n, and on a cond 1.5e10 polynomial block
% |y - y_true| was 2770 at its lambda (1.0 at GCV's). A RegMatrix with a
% null space (first difference) gets a grid below |N| (it started 1.5e4 |N|
% above it). Reference: direct solves of (N + lambda L'L) on the same grid.
% WP-137 (code review of WP-135): deterministic data whose troughs are
% interior for both L (grid index 51 of 80, the mirror on the same numbers),
% so no edge warning can fire -- the two muted ids could not, the input was
% fixed (HANDOVER 6.1.3: never mute an id that cannot fire); now asserted.
    m = 500;  t = linspace(0, 1, m).';
    B = t .^ (0:5);
    l = B * [1; -2; 0.5; 3; -1; 0.7] + 0.4 * sin(2*pi*23*t) + 0.4 * cos(2*pi*61*t);
    N = B.'*B;  n = B.'*l;
    for Lc = {eye(6), diff(eye(6))}
        L = Lc{1};
        [lam, g, grid] = tc.verifyWarningFree(@() lsqvce.engine.gcvLambda(N, n, L, m, l.'*l));
        gd = zeros(size(grid));
        for k = 1:numel(grid)
            M = N + grid(k) * (L.'*L);
            r = B * (M \ n) - l;
            gd(k) = (r.'*r) / (m - trace(M \ N))^2;
        end
        [~, kd] = min(gd);
        % the eigenbasis curve against direct solves of cond(N) = 1.5e7: worst
        % case cond*eps = 3e-9 relative, measured 9.4e-13 on these data
        % (dev/validation/wp135_gcv_mirror.py, HANDOVER 6.1.4)
        % tol: 1e-9 -- 1000x the measured for another BLAS, under the worst case
        tc.verifyEqual(g, gd, 'RelTol', 1e-9);
        tc.verifyEqual(lam, grid(kd));
        tc.verifyLessThan(grid(1), 1e-6 * norm(N));
        tc.verifyTrue(kd > 1 && kd < numel(grid), 'the trough is interior (no edge warning)');
    end
end

function testSnoopReSolveKeepsTheSolvesRegularisation(tc)
% categories: metamorphic, contract
% WP-137 (code review of WP-135, its lane-1 low "snoop's re-diagnosis loses
% options and Reg"): snoop's clean re-solve built its accumulator with the
% default options, so a fit whose local block was regularised (Reg,
% RegLambda) was re-solved WITHOUT it -- x after snoop jumped to the
% unregularised estimate. diagSol.accOpt carries the solve's options since
% WP-135; nothing tested it. Reference: the same regularised accumulation
% with the rejected rows at weight 0 (a zero row is a no-op in the QR).
    rs = RandStream('mt19937ar', 'Seed', 137);
    nc = 2;  nl = 3;  nB = 4;  mB = 60;
    blk = cell(nB, 1);
    for k = 1:nB
        Ak = [randn(rs, mB, nc), randn(rs, mB, nl)];
        blk{k} = {Ak, Ak * [1; -1; 0.5; 2; -0.3] + 0.1 * randn(rs, mB, 1)};
    end
    blk{2}{2}(17) = blk{2}{2}(17) + 5;                   % a 50-sigma blunder, row mB + 17
    lam = 0.5 * nB * mB;                                 % a ridge of half |N_ll|: y moves by O(1)
    function acc = regAcc(p)
        acc = lsqvce.engine.BlockAccumulator(nc, 'Backend', 'qr', 'Reg', 'tikhonov', ...
                                             'RegLambda', lam, 'Verbose', false);
        L = lsqvce.engine.batchSystem(nc, @(k) deal(blk{k}{1}, blk{k}{2}, p{k}), nB, ...
                                      'LocalColumns', nc + (1:nl));
        acc.addSystem(L.sysID, nl, L.closure, 1:nB);
    end
    acc = regAcc(repmat({ones(mB, 1)}, nB, 1));
    [D, dsol, bf] = lsqvce.engine.verdict.fromAccumulator(acc, acc.solve());
    [solS, ~, REJ] = lsqvce.engine.verdict.snoop(dsol, D, bf, 'verbose', false);
    tc.assertTrue(ismember(mB + 17, REJ.idx), 'the fixture: the blunder is rejected');
    p = repmat({ones(mB, 1)}, nB, 1);
    for i = REJ.idx(:).'
        k = ceil(i / mB);  p{k}(i - (k - 1) * mB) = 0;
    end
    accR = regAcc(p);
    solR = accR.solve();
    % tol: two QR accumulations of the same rows under the same ridge (the
    % zeroed rows are no-op reflectors) -- eps-level; the unregularised
    % re-solve was off by O(lam / |N_ll|) = O(1) relative in y
    tc.verifyEqual(solS.x, [solR.x; accR.backSubstitute(1, solR.x)], 'RelTol', 1e-10);
end

function testGcvLambdaScalesWithTheData(tc)
% categories: metamorphic
% WP-135: B -> c B, l -> c l scales N, n and l'l by c^2: the GCV curve scales
% by c^2 (m - tr H is scale-free) and lambda by c^2 -- the same grid index.
    t = linspace(0, 1, 200).';  B = t .^ (0:3);  l = B * [1; -2; 0.5; 3] + 0.1 * cos(40 * t);
    [lam1, ~, g1] = lsqvce.engine.gcvLambda(B.'*B, B.'*l, [], 200, l.'*l);
    c2 = 9;
    [lam2, ~, g2] = lsqvce.engine.gcvLambda(c2 * (B.'*B), c2 * (B.'*l), [], 200, c2 * (l.'*l));
    [~, k1] = min(g1);  [~, k2] = min(g2);
    tc.verifyEqual(k2, k1);
    % tol: log10 / logspace of a grid end scaled by 9 round to a few ulp
    tc.verifyEqual(lam2, c2 * lam1, 'RelTol', 1e-12);
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

function testFilterLawDenseFormIsGuardedToo(tc)
% categories: contract, robustness
% WP-134 (review 3 finding 8, user decision): the joint noise estimator's
% sweep builds a FILTER law's dense Cfcn per arc, whatever the whitener --
% past maxDenseArc a long arc ran out of memory instead of being refused. The
% dense form is guarded at the size it is asked for, the filter is not.
    ar = struct('type','arma', 'ar',[1 -0.5], 'sigma',0.7, 'maxDenseArc',10);
    t  = (0:99).'*5;
    [Cf, Wf] = lsqvce.engine.covModel(ar, t, struct());
    tc.verifyNotEmpty(Wf, 'the filter is not bounded');
    tc.verifySize(Cf(1:10, struct()), [10 10]);
    tc.verifyError(@() Cf(1:11, struct()), 'lsqvce:covModel:denseArcTooLong');
end

function testCovModelSigmaPerRowAndUnevenGrid(tc)
% categories: contract, precision
% WP-134 (review 3 lows): 'acf' averaged a per-row sigma silently -- it is
% C = D R D now; 'arma' refused a vector sigma with an unnamed || error --
% by name now; a psd / arma law on an uneven grid (a removed epoch) takes the
% gap for one step -- reported.
    t = (0:19).'*5;  s = linspace(0.5, 2, 20).';
    acf = @(dt, o) exp(-abs(dt)/10);
    Cf = lsqvce.engine.covModel(struct('type','acf', 'acf',acf, 'sigma',s), t, struct());
    R  = acf(t - t.', []);
    % tol: 4 eps relative -- the same two products per entry
    tc.verifyEqual(Cf(1:20, struct()), (s .* s.') .* R, 'RelTol', 4*eps);
    tc.verifyError(@() lsqvce.engine.covModel(struct('type','acf', 'acf',acf, 'sigma',s(1:5)), t, struct()), ...
                   'lsqvce:covModel:sigmaLength');
    tc.verifyError(@() lsqvce.engine.covModel(struct('type','arma', 'ar',[1 -0.5], 'sigma',s), t, struct()), ...
                   'lsqvce:covModel:badSigma');
    muteWarn(tc, 'lsqvce:covModel:unevenGrid');
    lsqvce.util.log.inventory('reset');
    lsqvce.engine.covModel(struct('type','arma', 'ar',[1 -0.5], 'sigma',1), t([1:9 11:20]), struct());
    tc.verifyTrue(any(strcmp(lsqvce.util.log.inventory(), 'lsqvce:covModel:unevenGrid')));
    lsqvce.util.log.inventory('reset');
end
