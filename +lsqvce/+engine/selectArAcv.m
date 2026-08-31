function [a, s2, p, sel] = selectArAcv(acv, nData, opts)
%LSQVCE.ENGINE.SELECTARACV  Automatic AR-order selection by ACV (decorrelation) fidelity.
%   [a, s2, p, sel] = lsqvce.engine.selectArAcv(acv, nData)
%   [a, s2, p, sel] = lsqvce.engine.selectArAcv(acv, nData, opts)
%
%   Chooses the AR order WITHOUT a candidate grid. Because Levinson-Durbin is
%   recursive, one recursion to pMax yields the innovation variance of every
%   lower order (via the reflection coefficients). This sweeps orders and takes
%   the SMALLEST whose refitted autocovariance is within AcvTol -- the same
%   decorrelation objective a whitened-residual whiteness test uses. AICc, BIC
%   and the PACF band are computed and REPORTED (in sel) but do NOT drive the
%   choice: they target prediction and under-fit the correlation tail that
%   corrupts the normal equations. This is the ACV-domain counterpart of
%   LSQVCE.ENGINE.SELECTARORDER (which is PSD-driven).
%
%   INPUT
%     acv     (L+1)x1   double   Sample autocovariance r(0..L); r(1)=r(0)=
%                                variance. Required. L = numel(acv)-1.
%     nData   1x1       numeric  Effective sample size, for AICc/BIC/PACF
%                                reporting only. Required; [] -> 1e4.
%     opts    1x1       struct   Options ([] or omitted -> all defaults):
%       .pMax     1x1 integer   Largest order to consider.   (default min(L,120))
%       .AcvTol   1x1 double    Target relative ACV error.   (default 1e-3)
%       .RootTol  1x1 double    Near-unit-root guard on |root|.(default 0.999)
%       .Plateau  1x1 double    Relative gain per order below which the curve
%                               counts as a knee.            (default 0.02)
%
%   OUTPUT
%     a    1x(p+1)   double   AR polynomial [1 a1 ... ap] (a(1)=1).
%     s2   1x1       double   Innovation variance.
%     p    1x1       double   Selected order (0 if the series is white).
%     sel  1x1       struct   Report:
%            .p         1x1        selected order
%            .acvRelErr 1x1        achieved relative ACV error at p
%            .sigma     1x1        sqrt(s2)
%            .verdict   1xk char   'good' | 'marginal' | 'bad' | 'white'
%            .warnings  1xW cell   actionable warning strings (may be empty)
%            .aic       (pMax+1)x1 Akaike IC per order (0..pMax)
%            .aicc      (pMax+1)x1 corrected AIC per order
%            .bic       (pMax+1)x1 Schwarz IC per order
%            .pacfBand  1x1        ~95% PACF band 2/sqrt(nData)
%            .refl      pMaxx1     reflection coefficients (partial autocorr.)
%            .maxRoot   1x1        max |root| of a
%            .curve     (pMax+1)x1 relative ACV error per order (inf where unvisited)
%
%   GUARDS (each a real failure mode)
%     near-unit root (|root|>RootTol)  an unmodelled trend; detrend first.
%     plateau (AcvTol never reached)   AR cannot fit sharp spectral lines
%                                      (1/rev, 2/rev); use LSQVCE.KERNEL.AR2ARMA.
%
%   EXAMPLE
%     r = 0.7^2 * lsqvce.kernel.armaACV([1 -1.6 0.9 -0.25], 1, 200);
%     [a, s2, p, sel] = lsqvce.engine.selectArAcv(r, 5000);
%     assert(p == 3 && strcmp(sel.verdict, 'good'));
%
%   See also LSQVCE.ENGINE.SELECTARORDER, LSQVCE.KERNEL.LEVINSONDURBIN, LSQVCE.KERNEL.ARMAACV,
%            LSQVCE.KERNEL.AR2ARMA, LSQVCE.KERNEL.RESIDUALACV
%
%   Created 2026-07-15 13:40 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.
    narginchk(2, 3);
    if nargin < 3, opts = struct(); end
    if isempty(nData), nData = 1e4; end
    validateattributes(acv, {'numeric'}, {'vector','finite','nonempty'}, mfilename, 'acv', 1);
    validateattributes(nData, {'numeric'}, {'scalar','positive'}, mfilename, 'nData', 2);
    gf   = @lsqvce.util.getfielddef;
    acv  = double(acv(:));  nLag = numel(acv);  r0 = acv(1);
    pMax = min(nLag - 1, gf(opts, 'pMax',    120));
    tolA = gf(opts, 'AcvTol',  1e-3);
    tolR = gf(opts, 'RootTol', 0.999);
    plat = gf(opts, 'Plateau', 0.02);

    if r0 <= 0 || pMax < 1                            % white / degenerate
        [a, s2, p, sel] = whiteReturn(r0);  return
    end

    % --- one recursion: reflection ladder -> sigma^2_p for every order -----
    [~, ~, refl] = lsqvce.kernel.levinsonDurbin(acv, pMax);
    refl = refl(:);  pMax = numel(refl);             % truncates if not PD
    if pMax < 1
        [a, s2, p, sel] = whiteReturn(r0);  return
    end
    s2Lad = r0 * [1; cumprod(1 - refl.^2)];          % sigma^2_p, index p+1
    kk    = (0:pMax).';
    aic   = nData*log(s2Lad) + 2*kk;
    aicc  = aic + 2*kk.*(kk+1) ./ max(nData - kk - 1, 1);
    bic   = nData*log(s2Lad) + kk*log(max(nData, 2));
    pacf  = 2 / sqrt(nData);                          % ~95% PACF band

    % --- ACV-fidelity sweep: smallest order reaching AcvTol, else the knee -
    % p = 0 (the white model, variance r0) is the first candidate, so genuinely
    % white input returns order 0 rather than an AR(1) with a zero coefficient.
    normAcv = max(norm(acv), realmin);
    err = inf(pMax + 1, 1);
    err(1) = norm([r0; zeros(nLag-1, 1)] - acv) / normAcv;   % p = 0
    prev = err(1);  knee = 0;  pStar = -1;
    if err(1) <= tolA, pStar = 0; end
    for pk = 1:pMax
        if pStar >= 0, break, end
        [ak, s2k] = lsqvce.kernel.levinsonDurbin(acv, pk);
        rFit = s2k * lsqvce.kernel.armaACV(ak, 1, nLag);    % SAME length as acv
        err(pk+1) = norm(rFit - acv) / normAcv;
        if knee == 0 && (prev - err(pk+1)) <= plat*prev, knee = pk - 1; end
        prev = err(pk+1);
        if err(pk+1) <= tolA, pStar = pk; end
    end
    if pStar < 0
        if knee >= 1, pStar = knee; else, [~, ix] = min(err); pStar = ix - 1; end
    end
    if pStar == 0
        a = 1;  s2 = r0;  p = 0;
    else
        [a, s2] = lsqvce.kernel.levinsonDurbin(acv, pStar);  p = pStar;
    end

    % --- guards + verdict --------------------------------------------------
    w = {};  mroot = max([abs(roots(a)); 0]);
    if mroot > tolR
        w{end+1} = sprintf(...
            'near-unit root |root|=%.4f: an unmodelled trend -- detrend first.', mroot);
    end
    if err(pStar+1) > tolA
        w{end+1} = ['AcvTol not reached (plateau): sharp spectral lines -- ' ...
                    'use lsqvce.kernel.ar2arma (ARMA) or a Toeplitz whitener.'];
    end
    verdict = 'good';
    if err(pStar+1) > 10*tolA, verdict = 'marginal'; end
    if mroot > tolR,           verdict = 'bad';      end
    if pStar == 0,             verdict = 'white';    end

    sel = struct('p',pStar, 'acvRelErr',err(pStar+1), 'sigma',sqrt(s2), ...
                 'verdict',verdict, 'warnings',{w}, 'aic',aic, 'aicc',aicc, ...
                 'bic',bic, 'pacfBand',pacf, 'refl',refl, 'maxRoot',mroot, ...
                 'curve',err);
end

%==========================================================================
function [a, s2, p, sel] = whiteReturn(r0)
%WHITERETURN  The order-0 (white) result and its report, used for degenerate or
%   genuinely white input (r0 <= 0, pMax < 1, or a non-PD ACV that truncates to
%   order 0).
%   INPUT
%     r0   1x1   double   lag-0 autocovariance (process variance); may be <= 0
%   OUTPUT
%     a    1x1   double   the trivial AR polynomial, 1
%     s2   1x1   double   innovation variance = max(r0, realmin)
%     p    1x1   double   selected order, 0
%     sel  1x1   struct   report with .verdict='white' and empty IC/curve fields
    a = 1;  s2 = max(r0, realmin);  p = 0;
    sel = struct('p',0, 'acvRelErr',0, 'sigma',sqrt(s2), 'verdict','white', ...
                 'warnings',{{}}, 'aic',[], 'aicc',[], 'bic',[], ...
                 'pacfBand',[], 'refl',[], 'maxRoot',0, 'curve',[]);
end
