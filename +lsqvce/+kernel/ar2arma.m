function [a, b, s2] = ar2arma(aHi, s2Hi, q, opt)
%LSQVCE.KERNEL.AR2ARMA  Reduce a high-order AR to a compact ARMA(p,q) (Durbin's method).
%
%   [a, b, s2] = LSQVCE.KERNEL.AR2ARMA(aHi, s2Hi, q, opt)
%
%   WHEN TO USE THIS
%     Only when the PSD has sharp spectral LINES that a pure AR would need a
%     very high order to represent -- for GRACE-FO that means the 1/rev and
%     2/rev harmonics in KBR or accelerometer residuals. A narrow spectral peak
%     is a near-unit-root pole pair; AR needs p ~ 100+ to place it accurately,
%     while ARMA can place it with a handful of coefficients.
%
%     For a smooth coloured floor with no lines, do NOT use this: a plain AR
%     of moderate order is simpler, always stationary, and just as accurate.
%
%   INPUT
%     aHi   1 x (P+1)  high-order AR polynomial (P >> p), from LSQVCE.KERNEL.LEVINSONDURBIN
%     s2Hi  scalar     its innovation variance
%     q     scalar     desired MA order, >= 1
%     opt   struct     .armaP  desired AR order (default 10)
%                      .armaN  impulse-response length used in the fit (default 4096)
%
%   OUTPUT
%     a, b  the ARMA(p,q) polynomials, both normalised to a(1) = b(1) = 1
%     s2    innovation variance (carried over from the AR fit)
%
%   METHOD
%     Durbin's two-stage procedure. The long AR is a nonparametric handle on
%     the spectrum; its impulse response h = 1/A_hi(z) is then fitted by a
%     compact rational form B(z)/A(z) using a Hankel least-squares solve.
%
%   STABILITY IS NOT GUARANTEED
%     The Hankel fit can produce a non-stationary AR part, or a non-invertible
%     MA part (a root of B(z) on or outside the unit circle, which makes the
%     whitener A(z)/B(z) unstable). Both are detected here and the function
%     falls back to the input AR with a warning rather than returning a model
%     whose autocovariance does not exist or that cannot whiten. If you see
%     either warning, lower armaP or q, or simply use the pure AR. An
%     over-parameterised fit (armaP above the true order) is where the MA part
%     goes non-invertible.
%
%   ALWAYS VERIFY
%     ARMA reduction is the least robust step in this toolbox. Check it:
%         [a,b,s2] = lsqvce.kernel.ar2arma(aHi, s2Hi, q, opt);
%         r_arma   = s2 * lsqvce.kernel.armaACV(a, b, 512);
%         r_psd    = lsqvce.kernel.psd2acv(f, S, dt, 512);
%         relErr   = norm(r_arma - r_psd) / norm(r_psd);   % want < 1e-3
%     If relErr is worse than the plain AR it replaced, keep the AR.
%
%   EXAMPLE
%     opt = struct('armaP',10, 'armaN',4096);
%     [aHi, s2] = lsqvce.kernel.levinsonDurbin(r, 100);
%     [a, b, s2] = lsqvce.kernel.ar2arma(aHi, s2, 2, opt);
%     numel(a)-1, numel(b)-1        % 10, 2  instead of 100
%
%   ERRORS
%     lsqvce:ar2arma:tooShort    opt.armaN is too short to fit the requested
%                                ARMA(p,q); at least 4*(p+q+1) is needed
%     lsqvce:ar2arma:degenerate  the reduction produced a degenerate MA part
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:ar2arma:degenerateMA                [warning] The MA fit is degenerate; falling back to the AR model
%     lsqvce:ar2arma:unstable                    [warning] The ARMA reduction produced a non-stationary AR part (max |root| = X)
%     lsqvce:ar2arma:nonInvertibleMA             [warning] The ARMA reduction produced a non-invertible MA part
%
%   See also LSQVCE.KERNEL.LEVINSONDURBIN, LSQVCE.KERNEL.ARMAACV, LSQVCE.ENGINE.CHECKPSD
%
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-126: the Hankel fit runs the
%   recursion forwards (it returned reciprocal roots); a non-invertible MA part
%   falls back to the AR (lsqvce:ar2arma:nonInvertibleMA). With assistance from
%   Claude (Anthropic), model Claude Opus 5.5.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+kernel (R4 decision R4-D1,
%   WP-64a): the KERNEL tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    narginchk(3,4);
    if nargin < 4 || isempty(opt), opt = struct(); end

    aHi = aHi(:).';
    validateattributes(aHi,  {'double'}, {'row','finite','nonempty'}, mfilename, 'aHi', 1);
    validateattributes(s2Hi, {'double'}, {'scalar','finite','positive'}, mfilename, 's2Hi', 2);
    validateattributes(q,    {'numeric'}, {'scalar','integer','positive'}, mfilename, 'q', 3);

    p = lsqvce.util.getfielddef(opt, 'armaP', 10);
    n = lsqvce.util.getfielddef(opt, 'armaN', 4096);

    validateattributes(p, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'opt.armaP');
    validateattributes(n, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'opt.armaN');

    p = double(p);  q = double(q);  n = double(n);

    if n < 4*(p+q+1)
        error('lsqvce:ar2arma:tooShort', ...
            ['armaN = %d is too short to fit ARMA(%d,%d). Use at least %d.'], ...
             n, p, q, 4*(p+q+1));
    end

    % --- stage 1: impulse response of the long AR --------------------------
    h = filter(1, aHi, [1; zeros(n-1,1)]);

    % --- stage 2: Hankel least squares for the AR part ---------------------
    H  = hankel(h(q+1:end-p), h(end-p:end));
    if size(H,1) < p
        error('lsqvce:ar2arma:degenerate', ...
              'The Hankel system is underdetermined (%d rows, %d unknowns).', ...
              size(H,1), p);
    end
    % Row i of H is h(q+i .. q+i+p). For n > q the ARMA impulse response obeys
    % h(n) + a_1 h(n-1) + ... + a_p h(n-p) = 0: the LAST column is the target
    % and a_k multiplies column p+1-k. WP-126: the fit used to regress the
    % FIRST column on the later ones -- the recursion run backwards, whose
    % roots are the reciprocals of the true ones (0.94 -> 1.064), so the
    % stationarity guard below fired on almost every input and the 'arma'
    % route silently returned the long AR.
    aa = -fliplr(H(:,1:end-1)) \ H(:,end);
    a  = [1; aa(:)].';

    % --- MA part from the filtered impulse response ------------------------
    hb = filter(a, 1, h);
    if abs(hb(1)) < eps
        lsqvce.util.log.warn('lsqvce:ar2arma:degenerateMA', ...
                'The MA fit is degenerate; falling back to the AR model.');
        a = aHi;  b = 1;  s2 = s2Hi;
        return
    end
    b  = [1; hb(2:q+1)].';
    b  = b / b(1);
    s2 = s2Hi;

    % --- stationarity guard -------------------------------------------------
    rts = roots(a);
    if ~isempty(rts) && max(abs(rts)) >= 1
        lsqvce.util.log.warn('lsqvce:ar2arma:unstable', ...
            ['The ARMA reduction produced a non-stationary AR part (max |root| ' ...
             '= %.6f). Falling back to the input AR(%d). Lower armaP or q, or ' ...
             'just use the pure AR.'], max(abs(rts)), numel(aHi)-1);
        a = aHi;  b = 1;  s2 = s2Hi;
        return
    end

    % --- invertibility guard (WP-126) -----------------------------------------
    % The whitener is A(z)/B(z): a root of B(z) on or outside the unit circle
    % makes 1/B(z) unstable (LSQVCE.ENGINE.COVMODEL refuses such a law). The
    % MA part is a fit, not a declaration, so the reduction falls back to the
    % input AR -- as for the AR part above -- instead of failing the run.
    rtb = roots(b);
    if ~isempty(rtb) && max(abs(rtb)) >= 1
        lsqvce.util.log.warn('lsqvce:ar2arma:nonInvertibleMA', ...
            ['The ARMA reduction produced a non-invertible MA part (max |root| ' ...
             '= %.6f). Falling back to the input AR(%d). Lower q or armaP, or ' ...
             'just use the pure AR.'], max(abs(rtb)), numel(aHi)-1);
        a = aHi;  b = 1;  s2 = s2Hi;
    end
end
