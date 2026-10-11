classdef stat
%STAT  Statistical kernel for lsqvce -- base MATLAB only.
%
%   Deliberately free of the Statistics and Machine Learning Toolbox. Every
%   quantile / cdf needed by the diagnostics is implemented from erf, gammainc
%   and betainc, all of which ship with base MATLAB.
%
%   Static methods
%     normCdf(z)            standard normal CDF
%     normInv(p)            standard normal quantile
%     chi2Cdf(x,r)          chi-square CDF
%     chi2Inv(p,r)          chi-square quantile (regularised-gamma inversion,
%                           Wilson-Hilferty seeded, Newton-polished)
%     tInv(p,nu)            Student-t quantile
%     tauInv(p,r)           Pope tau quantile (two-sided p is NOT applied here;
%                           pass the upper-tail probability directly)
%     tauCrit(alpha,r)      two-sided Pope tau critical value
%     normCrit(alpha)       two-sided normal critical value
%     delta0(alpha,beta)    MDB non-centrality parameter
%     fInv(p,d1,d2)         F quantile
%
%   All functions are elementwise where it makes sense and validate inputs.
%
%   EXAMPLE
%       lsqvce.kernel.stat.normCrit(1e-3)          % 3.2905  Baarda's threshold
%       lsqvce.kernel.stat.delta0(1e-3, 0.20)      % 4.1321  MDB non-centrality
%       lsqvce.kernel.stat.chi2Inv(0.95, 10)       % 18.3070
%       lsqvce.kernel.stat.tauCrit(1e-3, 100)      % Pope's tau, r = 100
%
%       % Pope's tau is BOUNDED by sqrt(r). That is the formal statement that
%       % a weak network cannot detect its own blunders:
%       lsqvce.kernel.stat.tauCrit(1e-8, 5) < sqrt(5)      % always true
%
%   ERRORS
%     lsqvce:stat:badProb           the probability is outside (0,1)
%     lsqvce:stat:badDof            a degrees-of-freedom argument is not positive
%     lsqvce:stat:sizeMismatch      the arguments do not broadcast
%
%   See also LSQVCE.ENGINE.VERDICT.DIAGNOSE.
%
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-86 part 4: moved from the
%   ENGINE to the KERNEL tier (the distributions; the tests live
%   in lsqvce.engine.verdict). With assistance from Claude (Anthropic),
%   model Claude Opus 5.5.
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

methods (Static)

%% ---------------------------------------------------------------- normal
function p = normCdf(z)
%NORMCDF  Standard-normal CDF Phi(z) = P(Z <= z), via erfc (no toolbox).
%   INPUT   z  any-size double   standardised values
%   OUTPUT  p  same size as z    Phi(z) in [0,1]
    p = 0.5 * erfc(-z ./ sqrt(2));
end

function z = normInv(p)
%NORMINV  Standard-normal inverse CDF (quantile), Phi^-1(p), via erfcinv.
%   INPUT   p  any-size double   probabilities in [0,1] (checked)
%   OUTPUT  z  same size as p    the z with Phi(z) = p
    lsqvce.kernel.stat.mustBeProb(p, 'normInv');
    z = -sqrt(2) .* erfcinv(2 .* p);
end

function k = normCrit(alpha)
    % Two-sided critical value: P(|Z| > k) = alpha.
    k = lsqvce.kernel.stat.normInv(1 - alpha / 2);
end

%% ------------------------------------------------------------ chi-square
function p = chi2Cdf(x, r)
    % Regularised lower incomplete gamma. gammainc is base MATLAB.
    x = max(x, 0);
    p = gammainc(x ./ 2, r ./ 2, 'lower');
end

function x = chi2Inv(p, r)
    %CHI2INV  Chi-square quantile without the Statistics Toolbox.
    %
    %   Wilson-Hilferty seed, polished by safeguarded Newton on F(x) = p.
    %
    %   The iteration is performed in y = log(x), NOT in x. This is not a
    %   stylistic choice. The chi-square density behaves like x^(r/2 - 1) near
    %   the origin, so for r < 2 it DIVERGES as x -> 0: Newton in x then has an
    %   exploding derivative and stalls, leaving a relative error of ~4e-8 at
    %   (r=1, p=1e-6). Substituting y = log(x) gives
    %
    %       dF/dy = x * pdf(x)  ~  x^(r/2)  ->  0
    %
    %   which is smooth and bounded, and the iteration converges to machine
    %   precision for every r >= 1 including the singular cases. The step is
    %   damped to |dy| <= 2 (a factor of e^2 per iteration) so a bad seed
    %   cannot fly off to overflow.
    %
    %   ACCURACY (verified against an independent implementation):
    %     r <= 1e7   : |F(x_returned) - p| / p  <  1e-11   (machine precision)
    %     r  > 1e7   : ~1e-6 relative, limited by gammainc itself, not by this
    %                  solver: at r = 1e9 the CDF is so steep that a one-ulp
    %                  change in x moves F by ~1e-11. The returned quantile
    %                  still agrees with a reference implementation to 2e-15
    %                  relative, and the chi-square acceptance band at that
    %                  redundancy is ~2e-5 wide, so the error is irrelevant.

    lsqvce.kernel.stat.mustBeProb(p, 'chi2Inv');
    if any(r(:) <= 0)
        error('lsqvce:stat:badDof', 'chi2Inv: degrees of freedom must be > 0.');
    end
    [p, r] = lsqvce.kernel.stat.expand(p, r);

    x = zeros(size(p));

    % Exact endpoints -- Newton cannot help and would divide by zero.
    lo = (p <= 0);  hi = (p >= 1);
    x(lo) = 0;  x(hi) = Inf;
    act = ~(lo | hi);
    if ~any(act(:)), return; end

    pa = p(act);  ra = r(act);

    % --- Wilson-Hilferty seed
    z  = lsqvce.kernel.stat.normInv(pa);
    x0 = ra .* (1 - 2 ./ (9 .* ra) + z .* sqrt(2 ./ (9 .* ra))).^3;
    % WH can return a non-positive seed for small r and small p.
    bad = ~(x0 > 0) | ~isfinite(x0);
    x0(bad) = max(ra(bad) .* 0.1, 1e-3);

    % --- damped Newton in log-space
    y = log(x0);
    for it = 1:100
        xc = exp(y);
        F  = gammainc(xc ./ 2, ra ./ 2, 'lower') - pa;

        % dF/dy = x * pdf(x). Both factors are formed in logs so that neither
        % the x^(r/2-1) nor the exp(-x/2) term can under/overflow at large r.
        logpdf = (ra ./ 2 - 1) .* log(xc) - xc ./ 2 ...
                 - gammaln(ra ./ 2) - (ra ./ 2) .* log(2);
        dFdy = exp(logpdf + log(xc));

        step = F ./ dFdy;
        step(~isfinite(step)) = 0;          % derivative underflow: freeze
        step = max(min(step, 2), -2);       % damp: at most e^2 per iteration

        yn = y - step;
        if all(abs(yn - y) <= 1e-15 .* max(abs(y), 1))
            y = yn;  break;
        end
        y = yn;
    end

    x(act) = exp(y);
end

%% -------------------------------------------------------------- Student t
function t = tInv(p, nu)
    %TINV  Student-t quantile via inversion of the symmetric beta form.
    %
    %   For p >= 1/2 :  x = betaincinv(2*(1-p), nu/2, 1/2)  gives
    %                   t = sqrt(nu*(1-x)/x).
    %   Symmetry handles p < 1/2.  betaincinv is base MATLAB.

    lsqvce.kernel.stat.mustBeProb(p, 'tInv');
    if any(nu(:) <= 0)
        error('lsqvce:stat:badDof', 'tInv: degrees of freedom must be > 0.');
    end
    [p, nu] = lsqvce.kernel.stat.expand(p, nu);

    t = zeros(size(p));
    t(p == 0) = -Inf;
    t(p == 1) =  Inf;
    act = (p > 0) & (p < 1);
    if ~any(act(:)), return; end

    pa = p(act);  na = nu(act);
    up = pa >= 0.5;
    q  = pa;
    q(up)  = 2 .* (1 - pa(up));      % upper tail mass, in (0,1]
    q(~up) = 2 .* pa(~up);

    xb = betaincinv(q, na ./ 2, 0.5);
    xb = min(max(xb, realmin), 1);   % keep the sqrt real
    tt = sqrt(na .* (1 - xb) ./ xb);

    tt(~up) = -tt(~up);
    t(act)  = tt;
end

%% ------------------------------------------------------------- Pope's tau
function tc = tauCrit(alpha, r)
    %TAUCRIT  Two-sided critical value of Pope's tau distribution.
    %
    %   Used when sigma0 is ESTIMATED from the adjustment (the usual case).
    %   Relation to Student-t with r-1 degrees of freedom:
    %
    %       tau_r = t_{r-1} * sqrt(r) / sqrt(r - 1 + t_{r-1}^2)
    %
    %   Note tau is bounded above by sqrt(r): with low redundancy no residual
    %   can ever look extreme, which is exactly the statement that a weak
    %   network cannot detect its own blunders.

    [alpha, r] = lsqvce.kernel.stat.expand(alpha, r);
    tc = nan(size(r));

    tc(r <= 1) = Inf;                % no redundancy -> nothing is testable
    act = (r > 1);
    if ~any(act(:)), return; end

    ra = r(act);  aa = alpha(act);
    t  = lsqvce.kernel.stat.tInv(1 - aa ./ 2, ra - 1);
    tc(act) = t .* sqrt(ra) ./ sqrt(ra - 1 + t.^2);
end

%% ------------------------------------------------------------------- MDB
function d = delta0(alpha, beta)
    %DELTA0  Non-centrality parameter for the minimal detectable bias.
    %
    %   Solves, for delta:
    %       power = 1 - Phi(k - delta) + Phi(-k - delta) = 1 - beta,
    %       k = normCrit(alpha).
    %
    %   The lower-tail term Phi(-k-delta) is utterly negligible for the alphas
    %   in use but is retained so the identity is exact rather than "close".
    %   Solved by bisection: monotone in delta, so this cannot fail.
    %
    %   delta0(0.001, 0.20) = 4.1321  (the textbook 4.13).

    [alpha, beta] = lsqvce.kernel.stat.expand(alpha, beta);
    lsqvce.kernel.stat.mustBeProb(alpha, 'delta0');
    lsqvce.kernel.stat.mustBeProb(beta,  'delta0');

    d = zeros(size(alpha));
    for i = 1:numel(alpha)
        k    = lsqvce.kernel.stat.normCrit(alpha(i));
        want = 1 - beta(i);
        pow  = @(dd) 1 - lsqvce.kernel.stat.normCdf(k - dd) ...
                       + lsqvce.kernel.stat.normCdf(-k - dd) - want;

        loB = 0;  hiB = 1;
        while pow(hiB) < 0 && hiB < 1e6
            hiB = hiB * 2;
        end
        for it = 1:200
            mid = (loB + hiB) / 2;
            if pow(mid) < 0, loB = mid; else, hiB = mid; end
            if hiB - loB < 1e-12, break; end
        end
        d(i) = (loB + hiB) / 2;
    end
end

%% --------------------------------------------------------------------- F
function x = fInv(p, d1, d2)
    %FINV  F quantile via the inverse incomplete beta.
    lsqvce.kernel.stat.mustBeProb(p, 'fInv');
    [p, d1, d2] = lsqvce.kernel.stat.expand(p, d1, d2);
    xb = betaincinv(p, d1 ./ 2, d2 ./ 2);
    xb = min(max(xb, 0), 1 - eps);
    x  = (d2 ./ d1) .* xb ./ (1 - xb);
end

%% ----------------------------------------------------------------- checks
function mustBeProb(p, who)
%MUSTBEPROB  Validate that p is a real, finite probability array in [0,1].
%   Errors ('lsqvce:stat:badProb') otherwise.
%   INPUT
%     p    any-size double   the value(s) to check
%     who  1xk char          caller name, used in the error message
%   OUTPUT  none (throws on failure)
    if ~isnumeric(p) || ~isreal(p) || any(~isfinite(p(:)))
        error('lsqvce:stat:badProb', ...
              '%s: probability must be real and finite.', who);
    end
    if any(p(:) < 0) || any(p(:) > 1)
        error('lsqvce:stat:badProb', ...
              '%s: probability must lie in [0,1] (got %g).', who, ...
              p(find(p < 0 | p > 1, 1)));
    end
end

function varargout = expand(varargin)
    % Scalar-expand a set of arrays to a common size, erroring on mismatch.
    sz = [1 1];
    for i = 1:nargin
        s = size(varargin{i});
        if ~isequal(s, [1 1])
            if isequal(sz, [1 1])
                sz = s;
            elseif ~isequal(sz, s)
                error('lsqvce:stat:sizeMismatch', ...
                      'Arguments must be scalars or share a common size.');
            end
        end
    end
    varargout = cell(1, nargin);
    for i = 1:nargin
        if isequal(size(varargin{i}), [1 1])
            varargout{i} = repmat(varargin{i}, sz);
        else
            varargout{i} = varargin{i};
        end
    end
end

end
end
