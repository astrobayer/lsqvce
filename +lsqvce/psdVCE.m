function out = psdVCE(omegaT, sT, Cprior, varargin)
%LSQVCE.PSDVCE  Per-frequency PSD amplitude by variance component estimation.
%
%   out = LSQVCE.PSDVCE(omegaT, sT, Cprior)
%   out = LSQVCE.PSDVCE(omegaT, sT, Cprior, 'Name', value, ...)
%
%   The core of Ellmer's (2018) non-parametric stochastic model: instead of
%   fitting a parametric AR/ARMA process, the AMPLITUDE OF THE PSD AT EVERY
%   FREQUENCY is estimated by VCE, treating each lag of the covariance function
%   as its own variance component (eqs. 6.5.19-6.5.27).
%
%       Shat^j = alpha2_j * S^j,     alpha2_j = OmegaTilde_j / sTilde_j   (6.5.19-20)
%       OmegaTilde_j = sum_m OmegaTilde^m_j,  sTilde_j = sum_m sTilde^m_j (6.5.21)
%       OmegaTilde^m = X * Omega^m,   sTilde^m = X * s^m                  (6.5.25-26)
%
%   with X the type-I DCT (LSQVCE.DCTI). The per-arc TEMPORAL vectors are the
%   weighted per-lag quantities of eqs. (6.5.22)-(6.5.23),
%
%       Omega^m_n = sigma2_m * C^n_xx * (ehat' Sigma^-T V_n Sigma^-1 ehat)
%       s^m_n     = sigma2_m * C^n_xx * trace(R V_n)
%
%   i.e. the RAW quantities from LSQVCE.REDUNDANCYLAGS scaled by the CURRENT arc
%   factor and lag covariance. This function takes the ALREADY-WEIGHTED temporal
%   vectors (see 'omegaT', 'sT' below) and does the spectral half.
%
%   WHY THE PSD AND NOT THE COVARIANCE FUNCTION (eq. 6.5.19 discussion)
%     A VCE factor is estimated as a square, so a covariance value estimated
%     directly in the time domain could NEVER change sign. The PSD is positive by
%     definition, so estimating it is unrestricted -- and the inverse DCT of a
%     positive PSD is free to produce NEGATIVE covariance lags, which real noise
%     processes require. This is the whole reason the estimation happens in the
%     spectral domain.
%
%   *** DEGENERATE START WARNING (verified numerically; read this) ***
%     Ellmer's literal white start (eq. 6.5.62), Cprior = [1 0 0 ... 0], is a
%     DEGENERATE FIXED POINT of his own weighting: with C^n = 0 for every n > 0,
%     eqs. (6.5.22)-(6.5.23) annihilate every lag but n = 0, so those lags can
%     never be updated. The iteration then reports convergence (|alpha2-1| ~ 1e-13)
%     while the recovered PSD is WRONG. Numerically: starting AT the true PSD gives
%     corr(Shat,Strue) = 0.9999 (the estimator is correct and holds the truth), a
%     non-degenerate flat-ACV start converges to corr 0.945, but the delta start
%     lands at corr -0.19 and calls it converged. This function DETECTS such a
%     start and warns ('lsqvce:psdVCE:degenerateStart'); pass 'AllowDegenerate'
%     true to proceed anyway (e.g. to reproduce ITSG bit-exactly, knowingly).
%
%   INPUT
%     omegaT   nLag x M double  per-arc TEMPORAL Omega^m_n vectors, ALREADY
%                               weighted by sigma2_m * C^n (eq. 6.5.22), one
%                               COLUMN per arc. Rows are lags n = 0..nLag-1.
%     sT       nLag x M double  the matching weighted redundancies s^m_n
%                               (eq. 6.5.23), same shape.
%     Cprior   nLag x 1 double  the CURRENT discretised covariance function C^n_xx
%                               the weighting used. Its DCT gives the a-priori
%                               S^j that alpha2_j scales (eq. 6.5.19).
%
%   OPTIONS (Name/Value)
%     'Norm'            char    DCT normalisation, 'ortho' (default) | 'ellmer'
%                               (see LSQVCE.DCTI).
%     'ZeroFreqRule'    char    how to treat the zero-frequency amplitude, which
%                               is ill-defined for a zero-mean process (Ellmer
%                               §6.5.6, "Estimation of the Zero-Frequency"):
%                               'lowest' (default) replace Shat^0 by the lowest
%                               non-zero frequency amplitude; 'keep' leave as is.
%     'MinAmplitude'    1x1     FLOOR applied to Shat (max(Shat, MinAmplitude);
%                               negative estimates land ON the floor, they are
%                               never abs()-folded -- WP-35/F-04, see the
%                               policy comment at the update line) to keep the
%                               PSD strictly
%                               positive (default 1e-12). A PSD must be >= 0.
%     'Smooth'          1x1     half-width of a moving average applied to the
%                               estimated amplitudes (default 0 = off). The
%                               per-frequency estimate alpha2_j = Omega_j/s_j is
%                               a RATIO from finite residuals, so it is noisy;
%                               smoothing across neighbouring frequencies trades
%                               spectral resolution for variance, the standard
%                               lag-window / spectral-averaging idea (and the
%                               approach favoured by the Delft school). Use a
%                               half-width of 1-3 for a mildly noisy spectrum.
%                               NOTE this does NOT fix a rank-deficient
%                               covariance -- that is a RESOLUTION problem, cured
%                               by using enough frequencies (nf >= m/2 + 1), not
%                               by smoothing.
%     'AllowDegenerate' logical proceed without warning on a degenerate start
%                               (default false).
%
%   OUTPUT
%     out   1x1 struct
%       .Shat      nLag x 1  the updated PSD amplitudes (strictly positive)
%       .alpha2    nLag x 1  the per-frequency variance factors alpha2_j. The
%                            iteration has converged when these are all ~ 1.
%       .Cnew      nLag x 1  the updated covariance function, inverse-DCT of Shat.
%                            MAY CONTAIN NEGATIVE LAGS -- that is correct and is
%                            the point of estimating in the spectral domain.
%       .Sprior    nLag x 1  the a-priori amplitudes X*Cprior
%       .maxDev    1x1       max|alpha2 - 1|, the convergence measure
%       .degenerate 1x1 logical  true if the start had no support beyond lag 0
%
%   EXAMPLE
%     % per arc: r = lsqvce.redundancyLags(...); weight by sigma2_m and Cprior
%     omegaT(:,m) = sig2(m) * Cprior .* r.omega;
%     sT(:,m)     = sig2(m) * Cprior .* r.s;
%     p = lsqvce.psdVCE(omegaT, sT, Cprior);
%     Cprior = p.Cnew;                       % next iteration's covariance
%     p.maxDev                               % -> 0 as the model settles
%
%   ERRORS
%     lsqvce:psdVCE:arcMismatch     the arc list disagrees with the residuals
%     lsqvce:psdVCE:shapeMismatch   the spectral shape has the wrong length
%
%   See also LSQVCE.DCTI, LSQVCE.REDUNDANCYLAGS, LSQVCE.TOEPLITZCOFACTOR,
%            LSQVCE.ARCVARIANCEFACTORS
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-35 (fixes F-04): negative
%   amplitude estimates are FLOORED at MinAmplitude, never abs()-folded --
%   abs injected ~2x the power at zero-power bins (v10review.py T6); bins
%   with genuine positive power are unchanged. With assistance from Claude
%   (Anthropic), model Claude Fable 5.
%   Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.

    p = inputParser;  p.FunctionName = 'lsqvce.psdVCE';
    p.addParameter('Norm'           , 'ortho');
    p.addParameter('ZeroFreqRule'   , 'lowest', ...
                   @(s) any(strcmpi(s, {'lowest','keep'})));
    p.addParameter('MinAmplitude'   , 1e-12, @(v) isscalar(v) && v >= 0);
    p.addParameter('Smooth'         , 0, @(v) isscalar(v) && v >= 0 && v == floor(v));
    p.addParameter('AllowDegenerate', false, ...
                   @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
    p.parse(varargin{:});
    o = p.Results;

    Cprior = Cprior(:);
    nLag   = numel(Cprior);
    if size(omegaT,1) ~= nLag || size(sT,1) ~= nLag
        error('lsqvce:psdVCE:shapeMismatch', ...
             ['omegaT (%d x %d) and sT (%d x %d) must have one ROW per lag to ' ...
              'match Cprior (%d lags).'], size(omegaT,1), size(omegaT,2), ...
              size(sT,1), size(sT,2), nLag);
    end
    if size(omegaT,2) ~= size(sT,2)
        error('lsqvce:psdVCE:arcMismatch', ...
              'omegaT has %d arcs but sT has %d.', size(omegaT,2), size(sT,2));
    end

    % ---- degenerate-start detection (see the header warning) -------------
    degenerate = all(abs(Cprior(2:end)) <= 0);
    if degenerate && ~logical(o.AllowDegenerate)
        lsqvce.util.log.warn('lsqvce:psdVCE:degenerateStart', ...
           ['The prior covariance has NO support beyond lag 0 (a delta / pure ' ...
            'white start). The weighting sigma^2*C^n then annihilates every lag ' ...
            'n > 0, so those lags can NEVER be updated: the iteration will ' ...
            'report convergence while the PSD stays wrong. Start from a ' ...
            'covariance with some correlation (e.g. an empirical residual ACV, ' ...
            'or C = [1 c c ...] with small c), or pass ''AllowDegenerate'',true ' ...
            'if you are reproducing this behaviour deliberately.']);
    end

    % ---- spectral transform of the per-arc temporal vectors  (6.5.25-26) --
    %   X acts column-wise, so all arcs transform in one call.
    OmTilde = lsqvce.dctI(omegaT, 'Norm', o.Norm);      % nLag x M
    sTilde  = lsqvce.dctI(sT    , 'Norm', o.Norm);      % nLag x M

    % ---- accumulate over arcs                                    (6.5.21) --
    Omj = sum(OmTilde, 2);
    sj  = sum(sTilde , 2);

    % ---- the per-frequency variance factor                       (6.5.20) --
    alpha2 = ones(nLag, 1);
    ok     = isfinite(Omj) & isfinite(sj) & (abs(sj) > 0);
    alpha2(ok) = Omj(ok) ./ sj(ok);
    alpha2(~isfinite(alpha2)) = 1;

    % ---- updated amplitudes                                      (6.5.19) --
    Sprior = lsqvce.dctI(Cprior, 'Norm', o.Norm);
    Shat   = alpha2 .* Sprior;
    % POLICY (WP-35, fixes F-04): FLOOR, never abs. A PSD is non-negative by
    % definition; a negative estimate at a (near-)zero-power bin is NOISE
    % around zero, and abs() FOLDS that noise -- E|x| = 2 E[x+] for zero-mean
    % x, i.e. abs injects ~2x the power a floor admits (measured:
    % validation/v10review.py T6). Bins with genuine positive power are
    % untouched by either policy. The comment above this line always
    % promised the floor; the code now keeps the promise.
    Shat   = max(Shat, o.MinAmplitude);
    if o.Smooth > 0
        Shat = smoothSpectrum(Shat, o.Smooth);
        Shat = max(Shat, o.MinAmplitude);
    end

    % ---- zero-frequency handling (Ellmer §6.5.6) --------------------------
    %   For a zero-mean process the theoretical S^0 is well defined but the
    %   discretised estimate can be ~0, making eq. (6.5.19) ill-conditioned.
    if strcmpi(o.ZeroFreqRule, 'lowest') && nLag >= 2
        Shat(1) = Shat(2);
    end

    % ---- back to the covariance function ---------------------------------
    %   Inverse DCT. For 'ortho' the transform is self-inverse; for 'ellmer'
    %   the explicit inverse is applied. Negative lags here are EXPECTED.
    Cnew = lsqvce.dctI(Shat, 'Norm', o.Norm, 'Inverse', true);
    % NOTE: no ridge is applied to Cnew here. The structural singularity of a
    % PSD-derived Toeplitz covariance is handled where the FULL-LENGTH ACV is
    % built, in LSQVCE.COVMODEL's 'toeplitz' branch (option 'psdRidge'); doing it
    % on this truncated nLag-vector would not help, because truncating an ACV is
    % itself not a valid covariance.

    out = struct('Shat', Shat, 'alpha2', alpha2, 'Cnew', Cnew, ...
                 'Sprior', Sprior, 'maxDev', max(abs(alpha2 - 1)), ...
                 'degenerate', degenerate);
end

function S = smoothSpectrum(S, h)
%SMOOTHSPECTRUM  Moving average of half-width h, with reflected edges.
%   Reflection keeps the endpoints (DC and Nyquist) from being pulled toward
%   zero, which a zero-padded average would do.
%   INPUT   S  nx1 double   amplitudes;  h 1x1 int  half-width
%   OUTPUT  S  nx1 double   smoothed amplitudes, same length
    n = numel(S);
    if h < 1 || n < 3, return, end
    h = min(h, floor((n-1)/2));
    idx = [h+1:-1:2, 1:n, n-1:-1:n-h];      % reflected padding
    Sp  = S(idx);
    w   = ones(2*h+1,1) / (2*h+1);
    Sc  = conv(Sp, w, 'same');
    S   = Sc(h+1 : h+n);
end
