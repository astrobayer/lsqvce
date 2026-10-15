function chk = checkPSD(cm, dt, opt, doPlot)
%LSQVCE.ENGINE.CHECKPSD  Verify that the AR/ARMA reduction reproduces the PSD.
%
%   chk = LSQVCE.ENGINE.CHECKPSD(cm, dt, opt)        cm.f, cm.S one-sided PSD
%   chk = LSQVCE.ENGINE.CHECKPSD(cm, dt, opt, true)  also plots the comparison
%
%   RUN THIS ONCE PER OBSERVATION TYPE before trusting the filter path. The
%   AR whitening filter is exact for the AR model it was fitted with -- but
%   only useful if that model actually reproduces YOUR noise. This function
%   closes that loop.
%
%   INPUT
%     cm      struct  .f (nf x 1) [Hz], .S (nf x 1) one-sided PSD
%     dt      scalar  sampling interval [s]
%     opt     struct  .arOrder    (default 30)
%                     .psdFloor   (default 1e-10)
%                     .checkLags  lags used for the ACV comparison (default 512)
%     doPlot  logical (default false)
%
%   OUTPUT
%     chk.acvRelErr  relative error of the refitted autocovariance  <- the one
%                    that matters
%     chk.psdRelErr  median relative error in the spectral domain
%     chk.arOrder    order actually fitted (may be < opt.arOrder if truncated)
%     chk.sigma      innovation std
%     chk.stable     all AR roots strictly inside the unit circle
%     chk.a          the fitted AR polynomial
%     chk.S_ar       the AR model's PSD on the same grid
%     chk.verdict    'good' | 'marginal' | 'bad'
%
%   HOW TO READ IT
%     acvRelErr < 1e-3   good      -- use the filter path
%     1e-3 .. 1e-2       marginal  -- raise arOrder, or accept a small bias
%     > 1e-2             bad       -- raise arOrder a lot, or use 'toeplitz'
%
%     A high acvRelErr with a LOW psdRelErr means the AR fits the spectrum in
%     the log-domain but misses the low-frequency power that dominates the
%     covariance. That is the dangerous case: the whitening looks right on a
%     log-log plot and is wrong where it counts. Always trust acvRelErr.
%
%   IF THE FIT IS BAD
%     1. Raise arOrder. Correlation length in samples is a decent first guess.
%     2. If the PSD has sharp lines (1/rev, 2/rev), a pure AR needs a very
%        high order to resolve them. Use psdMode = 'arma' (LSQVCE.KERNEL.AR2ARMA).
%     3. If chk.stable is false, the PSD has an unmodelled trend or a pole at
%        f = 0. Detrend the residuals before estimating the PSD.
%     4. If nothing works, use psdMode = 'toeplitz' and pay the O(b^3).
%
%   EXAMPLE
%     dt = 5;
%     f  = linspace(0, 0.5/dt, 4096).';
%     H  = 1 ./ (exp(-1i*2*pi*f*dt*(0:3)) * [1 -1.6 0.9 -0.25].');
%     cm = struct('type','psd','f',f,'S', 0.7^2*dt*2*abs(H).^2);
%
%     chk = lsqvce.engine.checkPSD(cm, dt, struct('arOrder',20), true);
%     assert(strcmp(chk.verdict,'good'), 'raise arOrder');
%
%   ACCURACY
%     The output is itself the measurement: acvRelErr compares the AR model's
%     ACV (LSQVCE.KERNEL.ARMAACV) with the PSD's own (LSQVCE.KERNEL.PSD2ACV)
%     over checkLags lags, kernels of the HANDOVER 6.1.4 row "Closed-form
%     ARMA(1,1)/AR ACV algebra" whose ACCURACY blocks state theirs. On a known
%     AR(3) spectrum an AR(20) refit gives acvRelErr < 1e-3, verdict 'good',
%     and an AR(1) refit is rated 'bad' (test_precision
%     (testCheckPSDVerdicts)). WP-139.
%
%   ERRORS
%     lsqvce:checkPSD:badModel      the model struct is not a recognisable
%                                   spectral model
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:checkPSD:badFit                     [warning] AR(X) reproduces the autocovariance to only X. Do NOT use the filter...
%     lsqvce:checkPSD:marginalFit                [warning] AR(X) reproduces the autocovariance to only X (target < 1e-3)
%     lsqvce:checkPSD:unstable                   [warning] The fitted AR(X) is not stationary (max |root| = X)
%
%   See also LSQVCE.ENGINE.COVMODEL, LSQVCE.KERNEL.LEVINSONDURBIN, LSQVCE.KERNEL.AR2ARMA
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-139: the help states its
%   accuracy, every number traced to the test that pins it. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    narginchk(2,4);
    if nargin < 3 || isempty(opt),    opt    = struct(); end
    if nargin < 4 || isempty(doPlot), doPlot = false;    end

    if ~isstruct(cm) || ~isfield(cm,'f') || ~isfield(cm,'S')
        error('lsqvce:checkPSD:badModel', 'cm must be a struct with fields f and S.');
    end
    validateattributes(dt, {'double'}, {'scalar','finite','positive'}, mfilename, 'dt', 2);

    p  = lsqvce.util.getfielddef(opt, 'arOrder',   30);
    fl = lsqvce.util.getfielddef(opt, 'psdFloor',  1e-10);
    nL = lsqvce.util.getfielddef(opt, 'checkLags', 512);

    f = cm.f(:);
    S = cm.S(:);
    S = max(S, fl*max(S));

    % --- reference ACV straight from the PSD -------------------------------
    r_psd = lsqvce.kernel.psd2acv(f, S, dt, nL);

    % --- AR fit, then its ACV ----------------------------------------------
    [a, s2] = lsqvce.kernel.levinsonDurbin(lsqvce.kernel.psd2acv(f, S, dt, p+1), p);
    r_ar    = s2 * lsqvce.kernel.armaACV(a, 1, nL);

    chk          = struct();
    chk.a        = a;
    chk.arOrder  = numel(a) - 1;
    chk.sigma    = sqrt(s2);
    chk.stable   = all(abs(roots(a)) < 1);
    chk.acvRelErr = norm(r_ar - r_psd) / norm(r_psd);

    % --- spectral-domain comparison ----------------------------------------
    H        = arFreqResp(a, f, dt);
    S_ar     = s2 * dt * 2 * abs(H).^2;
    chk.S_ar = S_ar;
    chk.psdRelErr = median(abs(S_ar - S) ./ S);

    % --- verdict ------------------------------------------------------------
    if     chk.acvRelErr < 1e-3, chk.verdict = 'good';
    elseif chk.acvRelErr < 1e-2, chk.verdict = 'marginal';
    else,                        chk.verdict = 'bad';
    end

    if ~chk.stable
        lsqvce.util.log.warn('lsqvce:checkPSD:unstable', ...
            ['The fitted AR(%d) is not stationary (max |root| = %.6f). The PSD ' ...
             'probably has an unmodelled trend or a pole at f = 0. Detrend the ' ...
             'residuals before estimating it.'], chk.arOrder, max(abs(roots(a))));
    end

    switch chk.verdict
    case 'marginal'
        lsqvce.util.log.warn('lsqvce:checkPSD:marginalFit', ...
            ['AR(%d) reproduces the autocovariance to only %.2e (target < 1e-3). ' ...
             'Raise arOrder, or accept a small bias in the weighting.'], ...
             chk.arOrder, chk.acvRelErr);
    case 'bad'
        lsqvce.util.log.warn('lsqvce:checkPSD:badFit', ...
            ['AR(%d) reproduces the autocovariance to only %.2e. Do NOT use the ' ...
             'filter path with this order. Raise arOrder substantially, try ' ...
             'psdMode = ''arma'', or fall back to psdMode = ''toeplitz''.'], ...
             chk.arOrder, chk.acvRelErr);
    end

    % --- plot ---------------------------------------------------------------
    if doPlot
        figure('Name', sprintf('checkPSD: AR(%d), %s', chk.arOrder, chk.verdict));

        subplot(1,2,1);
        loglog(f, S, 'k-', f, S_ar, 'r--', 'LineWidth', 1.2); grid on
        xlabel('frequency [Hz]'); ylabel('PSD [unit^2/Hz]');
        legend('input','AR fit','Location','best');
        title(sprintf('PSD, median rel err %.2e', chk.psdRelErr));

        subplot(1,2,2);
        lags = (0:nL-1)*dt;
        plot(lags, r_psd, 'k-', lags, r_ar, 'r--', 'LineWidth', 1.2); grid on
        xlabel('lag [s]'); ylabel('autocovariance');
        legend('input','AR fit','Location','best');
        title(sprintf('ACV, rel err %.2e  (%s)', chk.acvRelErr, chk.verdict));
    end
end

%==========================================================================
function H = arFreqResp(a, f, dt)
%ARFREQRESP  Frequency response of 1/A(z) at frequencies f. No toolbox needed.
    p = numel(a) - 1;
    z = exp(-1i*2*pi*f(:)*dt*(0:p));      % nf x (p+1)
    H = 1 ./ (z * a(:));
end
