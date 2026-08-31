function out = iterateNoiseModel(accFcn, cov0, opt)
%LSQVCE.ENGINE.ITERATENOISEMODEL  Iterate the covariance model until it is RIGHT, not assumed.
%
%   out = LSQVCE.ENGINE.ITERATENOISEMODEL(accFcn, cov0, 'Name', value, ...)
%
%   THE PROBLEM THIS CLOSES
%     The covariance model is normally assumed, then only SELF-CONSISTENCY is
%     checked (sigma0 -> 1). But the PSD you would estimate from residuals
%     comes from a solution that already assumed a PSD -- circular. This
%     helper breaks the circle by iterating:
%
%        solve -> residuals in the ORIGINAL metric -> empirical ACV ->
%        re-fit AR (order on ACV fidelity) -> re-solve -> check
%
%     and declaring victory only when BOTH hold:
%        sigma0 ~ 1          the noise LEVEL is right
%        acf(w) ~ delta      the noise COLOUR is right
%
%     Converges in 2 iterations from a completely wrong start (white assumed,
%     AR(2) true) -- measured in validation/v4verify.py, which also verifies
%     that the converged model matches the TRUE ACV, not merely a
%     self-consistent one.
%
%   INPUT
%     accFcn  @(cov) -> acc   builds AND accumulates a fresh accumulator with
%                             the given covariance specs (not yet solved).
%                             cov is a struct with one field per obsType,
%                             each a LSQVCE.ENGINE.COVMODEL spec.
%     cov0    the initial spec struct (your best a-priori guess; white is a
%             legitimate start)
%
%   OPTIONS
%     'UpdateTypes' cellstr of obsTypes whose model is re-estimated
%                   (default: all fields of cov0)
%     'MaxIter'     default 4 (2 is typical)
%     'SigmaTol'    |sigma0 - 1| tolerance (default 0.05)
%     'WhiteFactor' whiteness band = WhiteFactor/sqrt(n) on max |acf(w)| over
%                   lags 1..WhiteLags (default 3)
%     'WhiteLags'   default 5
%     'MaxLag'      lags of the empirical ACV used for the AR fit
%                   (default 128)
%     'OrderGrid'   AR orders to try (default [2 5 10 20 30 40 60])
%     'AcvTol'      target relative ACV fidelity of the fit (default 5e-2 --
%                   the EMPIRICAL acv carries sampling noise; do not chase it
%                   below its own scatter)
%     'Detrend'     subtract the per-arc mean before the ACV (default true)
%     'Verbose'     default true
%
%   OUTPUT
%     out.cov        the converged covariance specs ('arma' type, ar + sigma)
%     out.sol        the final solution under the converged model
%     out.res        the final LSQVCE.ENGINE.RESIDUALS output
%     out.acc        the final accumulator
%     out.converged  logical
%     out.nIter
%     out.history    struct array per iteration: .sigma0, .maxAcf (per
%                    obsType), .order, .sigmaInn, .acvRelErr
%
%   WHAT IT DOES NOT DO
%     It re-estimates the noise COLOUR and LEVEL per obsType, jointly for all
%     systems. It does not estimate per-arc variance components -- that is
%     LSQVCE.ENGINE.VCE, and the two compose: run iterateNoiseModel first to get the colour
%     right, then VCE for per-arc scale factors. It also assumes the residual
%     noise is stationary within an obsType; a manoeuvre-contaminated arc
%     should be cut in the loader, not modelled.
%
%   THE RESIDUAL-ACV BIAS, MEASURED
%     v = (I - H) e, so the empirical ACV of the residuals is biased by the
%     projector. The bias is O(u/m): measured < 2% relative at u/m = 1.6e-3
%     (v4verify.py). At GRACE-FO redundancy (m >> u) it is negligible. If
%     your problem has u/m > 0.05, distrust the second decimal of the fitted
%     model.
%
%   EXAMPLE
%     cov0.pos = struct('type','diag','sigma',1);      % start white
%     cov0.kbr = struct('type','diag','sigma',1);
%     accFcn   = @(cov) buildAccumulator(cfg, cov);    % you write this
%     it = lsqvce.engine.iterateNoiseModel(accFcn, cov0);
%     lsqvce.util.log.msg('%safter %d iterations, sigma0 = %.4f\n', ...
%             string(it.converged), it.nIter, it.sol.sigma0);
%
%   ERRORS
%     lsqvce:iterateNoiseModel:badAccFcn   the accumulator handle is not callable
%     lsqvce:iterateNoiseModel:badCov0     the starting covariance spec is malformed
%     lsqvce:iterateNoiseModel:noObsType   the stream has no observation type
%     lsqvce:iterateNoiseModel:tooShort    the series is too short to fit the
%                                          requested model order
%
%   See also LSQVCE.ENGINE.SELECTARORDER, LSQVCE.ENGINE.CHECKPSD, LSQVCE.ENGINE.VCE, LSQVCE.ENGINE.RESIDUALS
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (a) (F-07): typed
%   `arguments` signature; option names are UpperCamelCase (UpdateTypes,
%   MaxIter, SigmaTol, WhiteFactor, WhiteLags, MaxLag, OrderGrid, AcvTol,
%   Detrend, Verbose). NOTE (premise corrected round a1-2, measured on
%   CI): NV NAMES still match case-insensitively (MATLAB semantics) --
%   what tightened are the typed constraints and unknown-name rejection
%   (release note, HANDOVER §9 WP-39). With assistance from Claude
%   (Anthropic), model Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-40 (fixes F-13): the
%   whiteness statistic (maxAcf) pools per arc -- no lagged product
%   crosses an arc boundary, mirroring pooledACV. With assistance from
%   Claude (Anthropic), model Claude Fable 5.
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.
arguments
    accFcn (1,1) function_handle
    cov0   (1,1) struct
    opt.UpdateTypes cell = {}
    opt.MaxIter  (1,1) double {mustBeInteger, mustBePositive} = 4
    opt.SigmaTol (1,1) double {mustBePositive} = 0.05
    opt.WhiteFactor (1,1) double {mustBePositive} = 3
    opt.WhiteLags (1,1) double {mustBeInteger, mustBePositive} = 5
    opt.MaxLag   (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(opt.MaxLag, 8)} = 128
    opt.OrderGrid (1,:) double = [2 5 10 20 30 40 60]
    opt.AcvTol   (1,1) double {mustBePositive} = 5e-2
    opt.Detrend  (1,1) logical = true
    opt.Verbose  (1,1) logical = true
end

    narginchk(2, inf);
    if ~isa(accFcn, 'function_handle')
        error('lsqvce:iterateNoiseModel:badAccFcn', ...
              'accFcn must be a function handle @(cov) -> accumulator.');
    end
    if ~isstruct(cov0) || isempty(fieldnames(cov0))
        error('lsqvce:iterateNoiseModel:badCov0', ...
              'cov0 must be a struct with one covariance spec per obsType.');
    end

    % WP-39 stage (a): typed `arguments` signature; names UpperCamelCase
    % (release note in HANDOVER §9 WP-39). UpdateTypes defaults to all
    % obsTypes in cov0 -- the sentinel [] keeps that data-dependent default.
    o = struct('updateTypes', {opt.UpdateTypes}, 'maxIter', opt.MaxIter, ...
               'sigmaTol', opt.SigmaTol, 'whiteFactor', opt.WhiteFactor, ...
               'whiteLags', opt.WhiteLags, 'maxLag', opt.MaxLag, ...
               'orderGrid', opt.OrderGrid, 'acvTol', opt.AcvTol, ...
               'detrend', opt.Detrend, 'verbose', opt.Verbose);
    if isempty(o.updateTypes), o.updateTypes = fieldnames(cov0); end

    cov  = cov0;
    hist = struct([]);
    conv = false;

    for it = 1:o.maxIter
        acc = accFcn(cov);
        sol = acc.solve();
        res = acc.residuals(sol, 'Collect', true, 'Verbose', false);

        % ---- diagnose the CURRENT model -------------------------------
        H = struct('sigma0', sol.sigma0, 'obsType', struct([]));
        allWhite = true;
        for q = 1:numel(o.updateTypes)
            ot = o.updateTypes{q};
            % WP-40 (fixes F-13): the whiteness statistic pools PER ARC --
            % no lagged product crosses an arc boundary -- mirroring
            % pooledACV. Concatenation put spurious correlation into the
            % boundary products (a bias for many short arcs).
            wArcs = gatherResidualArcs(res, ot, 'w');
            n  = sum(cellfun(@numel, wArcs));
            if n < 10 * o.whiteLags
                error('lsqvce:iterateNoiseModel:tooShort', ...
                    'obsType %s has only %d whitened residuals.', ot, n);
            end
            mx = maxAcf(wArcs, o.whiteLags);
            band = o.whiteFactor / sqrt(n);
            e = struct('obsType', ot, 'maxAcf', mx, 'band', band, ...
                       'white', mx < band, 'order', NaN, 'sigmaInn', NaN, ...
                       'acvRelErr', NaN);
            allWhite = allWhite && e.white;
            if isempty(H.obsType), H.obsType = e;
            else,                  H.obsType(end+1) = e; end       %#ok<AGROW>
        end

        levelOK = abs(sol.sigma0 - 1) < o.sigmaTol;

        if lsqvce.util.verbosity.progress(o.verbose)
            lsqvce.util.log.msg('[psdIt %d ] sigma0 = %.4f   max|acf(w)| = %s\n', it, ...
                sol.sigma0, strjoin(compose('%s %.4f', ...
                string({H.obsType.obsType}).', [H.obsType.maxAcf].'), '  '));
        end

        if levelOK && allWhite
            conv = true;
            if isempty(hist), hist = H; else, hist(end+1) = H; end
            break
        end

        % ---- re-estimate the model from RAW-metric residuals ----------
        for q = 1:numel(o.updateTypes)
            ot = o.updateTypes{q};
            r  = pooledACV(res, ot, o.maxLag, o.detrend);
            [a, s2, relErr, ord] = fitAROnACV(r, o.orderGrid, o.acvTol);
            cov.(ot) = struct('type','arma', 'ar', a(:).', 'sigma', sqrt(s2));
            H.obsType(q).order     = ord;
            H.obsType(q).sigmaInn  = sqrt(s2);
            H.obsType(q).acvRelErr = relErr;
            if lsqvce.util.verbosity.progress(o.verbose)
                lsqvce.util.log.msg('           %s -> AR(%d), sigma = %.4g, acvRelErr = %.2e\n', ...
                        ot, ord, sqrt(s2), relErr);
            end
        end
        if isempty(hist), hist = H; else, hist(end+1) = H; end     %#ok<AGROW>
    end

    if ~conv
        lsqvce.util.log.warn('lsqvce:iterateNoiseModel:notConverged', ...
            ['Not converged in %d iterations (sigma0 = %.4f). Either the ' ...
             'noise is non-stationary within an obsType (cut the offending ' ...
             'arcs in the loader), the AR order grid tops out too low, or ' ...
             'the spectrum has sharp lines an AR cannot represent -- see ' ...
             'lsqvce.engine.selectArOrder guard 2.'], o.maxIter, sol.sigma0);
    end

    out = struct('cov', cov, 'sol', sol, 'res', res, 'acc', acc, ...
                 'converged', conv, 'nIter', it, 'history', hist);

    if lsqvce.util.verbosity.progress(o.verbose)
        if conv, s = 'converged'; else, s = 'STOPPED (not converged)'; end
        lsqvce.util.log.msg('%s%s after %d iteration(s)\n', lsqvce.util.verbosity.banner('PSD iteration'), s, it);
    end
end

%==========================================================================

%==========================================================================
function mx = maxAcf(wArcs, nLag)
%MAXACF  max |acf| over lags 1..nLag, POOLED PER ARC (WP-40, fixes F-13).
%   Arcs are independent realisations: the per-arc mean is removed and no
%   lagged product crosses an arc boundary -- the same discipline as
%   pooledACV. Numerators and the denominator pool over arcs.
    nA  = numel(wArcs);
    num = zeros(nLag, 1);
    den = 0;
    for a = 1:nA
        w = wArcs{a}(:);
        w = w - mean(w);
        den = den + w.'*w;
        for k = 1:min(nLag, numel(w)-1)
            num(k) = num(k) + w(1:end-k).'*w(1+k:end);
        end
    end
    mx = max(abs(num)) / den;
end

function C = gatherResidualArcs(res, obsType, fld)
%GATHERRESIDUALARCS  One obsType's residuals as a per-arc cell (WP-40).
%   Residuals per arc, for statistics that must not
%   cross arc boundaries (maxAcf).
    C = {};
    for k = 1:numel(res.sys)
        obs = res.sys(k).obs;
        for j = 1:numel(obs)
            if strcmp(obs(j).obsType, obsType)
                C{end+1} = obs(j).(fld)(:);                      %#ok<AGROW>
            end
        end
    end
    if isempty(C)
        error('lsqvce:iterateNoiseModel:noObsType', ...
              'No residuals found for obsType %s.', obsType);
    end
end

%==========================================================================
function r = pooledACV(res, obsType, maxLag, detrend)
%POOLEDACV  Biased sample ACV of the RAW residuals, pooled over arcs.
%   Arcs are independent realisations: lagged products never cross an arc
%   boundary. The biased (1/n) normalisation is deliberate -- it keeps the
%   ACV sequence positive semi-definite, which Levinson requires.
    num = zeros(maxLag + 1, 1);
    den = 0;
    for k = 1:numel(res.sys)
        obs = res.sys(k).obs;
        for j = 1:numel(obs)
            if ~strcmp(obs(j).obsType, obsType), continue, end
            v = obs(j).v(:);
            if detrend, v = v - mean(v); end
            n = numel(v);
            L = min(maxLag, n - 1);
            for tau = 0:L
                num(tau+1) = num(tau+1) + v(1:n-tau).'*v(1+tau:n);
            end
            den = den + n;
        end
    end
    if den == 0
        error('lsqvce:iterateNoiseModel:noObsType', ...
              'No residuals for obsType ''%s''.', obsType);
    end
    r = num / den;
end

%==========================================================================
function [a, s2, relErr, ord] = fitAROnACV(r, grid, tol)
%FITARONACV  AR fit on the EMPIRICAL ACV, order chosen on ACV fidelity --
%   the same criterion as LSQVCE.ENGINE.SELECTARORDER, on measured data.
    r    = r(:);
    L    = numel(r) - 1;
    grid = grid(grid < L);
    if isempty(grid), grid = min(2, L - 1); end
    best = [];
    for pp = grid(:).'
        [a, s2] = lsqvce.kernel.levinsonDurbin(r, pp);
        rf = s2 * lsqvce.kernel.armaACV(a, 1, L + 1);
        e  = norm(rf(:) - r) / norm(r);
        if isempty(best) || e < best.e
            best = struct('a', a, 's2', s2, 'e', e, 'p', numel(a) - 1);
        end
        if e < tol, break, end
    end
    a = best.a;  s2 = best.s2;  relErr = best.e;  ord = best.p;
end
