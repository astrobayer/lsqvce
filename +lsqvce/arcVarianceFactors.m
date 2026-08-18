function out = arcVarianceFactors(omegaT, sT, sigma2, varargin)
%LSQVCE.ARCVARIANCEFACTORS  Relative per-arc variance factors within a group.
%
%   out = LSQVCE.ARCVARIANCEFACTORS(omegaT, sT, sigma2)
%   out = LSQVCE.ARCVARIANCEFACTORS(omegaT, sT, sigma2, 'Name', value, ...)
%
%   Estimates Ellmer's (2018) arc-wise variance factors: the RELATIVE weight of
%   each short arc with respect to the other arcs of the SAME observation type
%   (eqs. 6.5.27-6.5.30). Where the PSD amplitudes are found by estimating one
%   frequency across all arcs, the arc factors are the mirror image -- the
%   CUMULATIVE amplitude over all frequencies for one arc:
%
%       alpha2_m = sum_j OmegaTilde^m_j / sum_j sTilde^m_j            (6.5.27)
%       gamma2_m = alpha2_m * sigma2_m                                (6.5.28)
%       gammaBar = mean over arcs of gamma_m                          (6.5.29)
%       sigma2hat_m = gamma2_m / gammaBar^2                           (6.5.30)
%
%   WHY THE gamma DETOUR (eq. 6.5.28 discussion) -- do not "simplify" this
%     The estimate is NOT simply sigma2hat_m = alpha2_m * sigma2_m. The overall
%     magnitude of the covariance matrix has ALREADY been adjusted by the
%     per-frequency PSD amplitudes (LSQVCE.PSDVCE), so taking alpha2_m * sigma2_m
%     directly would count the same variance component TWICE. The intermediate
%     gamma2_m is therefore normalised by the mean over arcs, leaving the arc
%     factors carrying only the RELATIVE arc-to-arc weighting while the
%     covariance function carries the absolute magnitude and structure.
%
%   THE MEAN: ROBUST BY DEFAULT
%     'robust' (default, what GROOPS does, eq. 6.5.57): the mean of the CENTRAL
%       HALF of the sorted gamma_m,  gammaBar = (2/M) * sum_{m=M/4}^{3M/4} gamma_m,
%       ignoring outliers at either tail. Recommended: a handful of bad arcs
%       (attitude-hold periods, manoeuvres) otherwise drag the normalisation.
%     'plain' (eq. 6.5.29): the arithmetic mean of all gamma_m. This is the only
%       variant for which mean(gamma_m/gammaBar) == 1 EXACTLY, i.e. the mean of
%       the arc factors sqrt(sigma2hat_m) is exactly 1. The robust variant gives
%       up that exact property in exchange for outlier resistance -- with planted
%       outliers the plain mean was 1.69 against a robust 1.24 and a median 1.24.
%
%   INPUT
%     omegaT   nLag x M double  per-arc TEMPORAL Omega^m_n vectors, already
%                               weighted by sigma2_m * C^n (eq. 6.5.22); one
%                               COLUMN per arc, rows are lags.
%     sT       nLag x M double  the matching weighted redundancies s^m_n.
%     sigma2   M x 1 double     the CURRENT arc variance factors sigma2_m that
%                               were used in the weighting (eq. 6.5.28 needs them).
%
%   OPTIONS (Name/Value)
%     'Norm'      char    DCT normalisation, 'ortho' (default) | 'ellmer'
%     'MeanType'  char    'robust' (default) | 'plain'   (see above)
%     'TrimFrac'  1x1     fraction trimmed from EACH tail for 'robust'
%                         (default 0.25, i.e. Ellmer's central half)
%
%   OUTPUT
%     out   1x1 struct
%       .sigma2   M x 1   the updated arc variance factors sigma2hat_m
%       .alpha2   M x 1   the raw cumulative factors alpha2_m (eq. 6.5.27)
%       .gamma2   M x 1   the intermediate gamma2_m (eq. 6.5.28)
%       .gammaBar 1x1     the normalising mean actually used
%       .meanType char    which mean was applied
%
%   EXAMPLE
%     a = lsqvce.arcVarianceFactors(omegaT, sT, sig2);
%     sig2 = a.sigma2;                 % feed back into the next iteration
%     mean(sqrt(a.sigma2))             % == 1 exactly only for 'plain'
%
%   ERRORS
%     lsqvce:arcVarianceFactors:arcMismatch  the arc list disagrees with the
%                                   accumulated arcs
%
%   See also LSQVCE.PSDVCE, LSQVCE.REDUNDANCYLAGS, LSQVCE.DCTI
%
%   Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.

    p = inputParser;  p.FunctionName = 'lsqvce.arcVarianceFactors';
    p.addParameter('Norm'    , 'ortho');
    p.addParameter('MeanType', 'robust', ...
                   @(s) any(strcmpi(s, {'robust','plain'})));
    p.addParameter('TrimFrac', 0.25, @(v) isscalar(v) && v >= 0 && v < 0.5);
    p.parse(varargin{:});
    o = p.Results;

    sigma2 = sigma2(:);
    M      = numel(sigma2);
    if size(omegaT,2) ~= M || size(sT,2) ~= M
        error('lsqvce:arcVarianceFactors:arcMismatch', ...
             ['omegaT has %d and sT has %d columns, but sigma2 lists %d arcs; ' ...
              'one COLUMN per arc is required.'], size(omegaT,2), size(sT,2), M);
    end

    % ---- spectral transform, then cumulative over frequencies    (6.5.27) --
    OmTilde = lsqvce.dctI(omegaT, 'Norm', o.Norm);
    sTilde  = lsqvce.dctI(sT    , 'Norm', o.Norm);
    num     = sum(OmTilde, 1).';                 % M x 1
    den     = sum(sTilde , 1).';                 % M x 1

    alpha2      = ones(M,1);
    ok          = isfinite(num) & isfinite(den) & (abs(den) > 0);
    alpha2(ok)  = num(ok) ./ den(ok);
    alpha2(~isfinite(alpha2)) = 1;

    % ---- gamma detour, then normalise                       (6.5.28-6.5.30) --
    gamma2 = alpha2 .* sigma2;
    gamma  = sqrt(abs(gamma2));                  % magnitudes for the mean
    gammaBar = meanOf(gamma, o.MeanType, o.TrimFrac);
    if ~(isfinite(gammaBar) && gammaBar > 0)
        lsqvce.util.log.warn('lsqvce:arcVarianceFactors:badMean', ...
           ['The normalising mean of the arc factors is %g, which cannot scale ' ...
            'anything. Leaving the arc factors unchanged this pass.'], gammaBar);
        out = struct('sigma2', sigma2, 'alpha2', alpha2, 'gamma2', gamma2, ...
                     'gammaBar', gammaBar, 'meanType', lower(o.MeanType));
        return
    end
    sigma2Hat = gamma2 / gammaBar^2;

    out = struct('sigma2', sigma2Hat, 'alpha2', alpha2, 'gamma2', gamma2, ...
                 'gammaBar', gammaBar, 'meanType', lower(o.MeanType));
end

%==========================================================================
function gb = meanOf(gamma, meanType, trimFrac)
%MEANOF  Plain arithmetic mean, or Ellmer's robust central-half mean (6.5.57).
    if strcmpi(meanType, 'plain')
        gb = mean(gamma);
        return
    end
    M = numel(gamma);
    g = sort(gamma(:));
    lo = floor(M*trimFrac) + 1;
    hi = M - floor(M*trimFrac);
    if hi < lo                      % too few arcs to trim: fall back
        gb = mean(gamma);
    else
        gb = mean(g(lo:hi));
    end
end
