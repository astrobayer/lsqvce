function info = selectArOrder(cm, dt, opt)
%LSQVCE.ENGINE.SELECTARORDER  Choose the AR order automatically, on the RIGHT criterion.
%
%   info = LSQVCE.ENGINE.SELECTARORDER(cm, dt, 'Name', value, ...)
%
%   Selects the AR order by ACV FIDELITY -- how well the fitted model reproduces
%   the autocovariance implied by your PSD -- not by AIC or BIC.
%
%   WHY NOT AIC/BIC
%     They target PREDICTION error. We want DECORRELATION fidelity. These are
%     not the same objective, and the difference bites: AIC will happily select
%     a low order that predicts well on average while leaving structure in the
%     tail of the autocorrelation -- which is exactly what corrupts the normal
%     equations. AIC and BIC are still REPORTED here, for information, but they
%     do not drive the choice.
%
%   INPUT
%     cm   struct with .f (Hz) and .S (one-sided PSD)
%     dt   sampling interval [s]
%
%   OPTIONS
%     'orderGrid'  orders to try (default [2 5 10 20 30 40 60 80 120 160])
%     'acvTol'     target relative ACV error (default 1e-3)
%     'checkLags'  lags used in the comparison (default 512)
%     'nData'      sample size assumed for AIC/BIC (default 10000)
%     'MaxMemory'  warn if filterMemory exceeds this fraction of BatchSize
%     'batchSize'  for that warning (default 20000)
%     'verbose'    default true
%     'plot'       default false
%
%   OUTPUT
%     info.p           selected order
%     info.acvRelErr   achieved ACV fidelity
%     info.sigma       innovation std
%     info.a           the AR polynomial
%     info.verdict     'good' | 'marginal' | 'bad'
%     info.curve       struct array: p, acvRelErr, aic, bic, maxRoot
%     info.warnings    cellstr of actionable warnings
%
%   THE GUARDS -- each catches a real failure mode
%     1. NEAR-UNIT ROOT (|root| > 0.999)
%        The fit is chasing a pole at f = 0. The PSD probably has an unmodelled
%        trend. Raising the order makes this WORSE, not better -- detrend the
%        residuals before estimating the PSD.
%
%     2. PLATEAU
%        The error stops improving before reaching tolerance. AR simply cannot
%        represent this spectrum -- usually sharp spectral lines (1/rev, 2/rev).
%        Use psdMode = 'arma' (LSQVCE.KERNEL.AR2ARMA) or fall back to 'toeplitz'.
%
%     3. filterMemory vs BatchSize
%        A large p means a large batch overlap. If p approaches BatchSize/10 you
%        are recomputing most rows. Raise BatchSize (free on the filter path) or
%        accept a looser AcvTol.
%
%   RULES OF THUMB (5 s sampling)
%     white / near-white ........ p = 0..2
%     smooth coloured floor ..... p = 20..40
%     floor + orbital harmonics . p = 60..100, or an ARMA
%   p scales with the correlation length in SAMPLES, so halving dt roughly
%   doubles the p you need.
%
%   EXAMPLE
%     info = lsqvce.engine.selectArOrder(cfg.cov.pos, 5, 'acvTol', 1e-3, 'plot', true);
%     assert(strcmp(info.verdict,'good'), strjoin(info.warnings, ' | '));
%     opt.arOrder = info.p;                       % feed it to the accumulator
%
%   ERRORS
%     lsqvce:selectAR:badModel  the candidate model spec is unrecognised
%     lsqvce:selectAR:noFit     no order in the search range produced a fit
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (b) b2 (F-07):
%   typed `arguments` signature, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   See also LSQVCE.ENGINE.CHECKPSD, LSQVCE.KERNEL.LEVINSONDURBIN, LSQVCE.KERNEL.AR2ARMA, LSQVCE.KERNEL.CROSSCORR
%
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

arguments
    cm
    dt
    opt.OrderGrid (1,:) double = [2 5 10 20 30 40 60 80 120 160]
    opt.AcvTol (1,1) double {mustBePositive} = 1e-3
    opt.CheckLags (1,1) double {mustBeGreaterThanOrEqual(opt.CheckLags,8)} = 512
    opt.NData (1,1) double {mustBePositive} = 10000
    opt.BatchSize (1,1) double {mustBePositive} = 20000
    opt.Verbose (1,1) logical = true
    opt.Plot (1,1) logical = false
end

    narginchk(2,inf);

    % WP-39 stage (b) b2: typed `arguments` signature (canonical
    % UpperCamelCase; legacy casings keep matching -- MATLAB NV
    % semantics; enum VALUES case-sensitive lowercase).
    o = struct('orderGrid', opt.OrderGrid, 'acvTol', opt.AcvTol, 'checkLags', opt.CheckLags, 'nData', opt.NData, 'batchSize', opt.BatchSize, 'verbose', opt.Verbose, 'plot', opt.Plot);

    if ~isstruct(cm) || ~isfield(cm,'f') || ~isfield(cm,'S')
        error('lsqvce:selectAR:badModel', 'cm must have fields f and S.');
    end

    f  = cm.f(:);
    S  = cm.S(:);
    S  = max(S, 1e-10*max(S));
    nL = o.checkLags;

    r_ref = lsqvce.kernel.psd2acv(f, S, dt, nL);

    orderGrid = o.orderGrid(:).';
    curve  = struct('p',{},'acvRelErr',{},'aic',{},'bic',{},'maxRoot',{},'sigma',{});
    best   = [];
    ws     = warning('off','lsqvce:levinson:notPD');
    cleanup = onCleanup(@() warning(ws));

    for k = 1:numel(orderGrid)
        pk = orderGrid(k);
        if pk + 1 > numel(f)*4, continue, end

        [a, s2] = lsqvce.kernel.levinsonDurbin(lsqvce.kernel.psd2acv(f, S, dt, pk+1), pk);
        r_fit   = s2 * lsqvce.kernel.armaACV(a, 1, nL);

        e   = norm(r_fit - r_ref) / norm(r_ref);
        rts = abs(roots(a));
        if isempty(rts), rts = 0; end

        rec = struct('p', numel(a)-1, 'acvRelErr', e, ...
                     'aic', o.nData*log(s2) + 2*pk, ...
                     'bic', o.nData*log(s2) + pk*log(o.nData), ...
                     'maxRoot', max(rts), 'sigma', sqrt(s2));
        rec.a = a;
        curve(end+1) = rmfield(rec,'a');                          %#ok<AGROW>

        if isempty(best) || e < best.acvRelErr
            best = rec;
        end
        if e < o.acvTol
            break
        end
    end

    if isempty(best)
        error('lsqvce:selectAR:noFit', 'No AR order could be fitted.');
    end

    % ---- verdict ----------------------------------------------------------
    if     best.acvRelErr < o.acvTol,      verdict = 'good';
    elseif best.acvRelErr < 10*o.acvTol,   verdict = 'marginal';
    else,                                  verdict = 'bad';
    end

    % ---- guards -----------------------------------------------------------
    warns = {};

    if best.maxRoot > 0.999
        warns{end+1} = sprintf(['near-unit AR root (|root| = %.6f): the PSD ' ...
            'probably has an unmodelled trend or a pole at f = 0. Detrend the ' ...
            'residuals BEFORE estimating the PSD -- raising the order will make ' ...
            'this worse, not better.'], best.maxRoot);                %#ok<AGROW>
    end

    if numel(curve) >= 2 && best.acvRelErr > o.acvTol
        impr = curve(end).acvRelErr / max(curve(end-1).acvRelErr, realmin);
        if impr > 0.9
            warns{end+1} = ['the fit has plateaued: AR cannot represent this ' ...
                'spectrum (sharp spectral lines?). Use psdMode = ''arma'' ' ...
                '(lsqvce.kernel.ar2arma) or fall back to psdMode = ''toeplitz''.'];  %#ok<AGROW>
        end
    end

    if best.p > o.batchSize/10
        warns{end+1} = sprintf(['filterMemory = %d is more than 10%% of ' ...
            'BatchSize = %d: the batch overlap will dominate the cost. Raise ' ...
            'BatchSize (free on the filter path) or loosen AcvTol.'], ...
            best.p, o.batchSize);                                     %#ok<AGROW>
    end

    if strcmp(verdict,'bad')
        warns{end+1} = sprintf(['best AR(%d) reproduces the ACV to only %.2e ' ...
            '(target %.0e). Do NOT use the filter path at this order.'], ...
            best.p, best.acvRelErr, o.acvTol);                        %#ok<AGROW>
    end

    info = struct('p', best.p, 'acvRelErr', best.acvRelErr, 'sigma', best.sigma, ...
                  'a', best.a, 'verdict', verdict, 'curve', curve, ...
                  'warnings', {warns});

    if lsqvce.util.verbosity.progress(o.verbose)
        lsqvce.util.log.msg('%sp = %d   acvRelErr = %.3e   verdict = %s\n', lsqvce.util.verbosity.banner('AR order selection'), ...
                info.p, info.acvRelErr, verdict);
        for k = 1:numel(warns)
            lsqvce.util.log.msg('           ! %s\n', warns{k});
        end
    end

    if o.plot
        e = [curve.acvRelErr];  pp = [curve.p];
        figure('Name','selectArOrder');
        subplot(1,2,1);
        semilogy(pp, e, 'o-', 'LineWidth',1.2); grid on; hold on
        yline(o.acvTol, 'r--', 'target');
        xline(info.p, 'g--', 'selected');
        xlabel('AR order p'); ylabel('ACV relative error');
        title(sprintf('selected p = %d (%s)', info.p, verdict));
        subplot(1,2,2);
        plot(pp, [curve.aic], 'o-', pp, [curve.bic], 's-', 'LineWidth',1.2);
        grid on; legend('AIC','BIC','Location','best');
        xlabel('AR order p'); ylabel('criterion');
        title('AIC / BIC (reported, NOT used to select)');
    end
end
