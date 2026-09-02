function out = arcVarianceFactors(omegaT, sT, sigma2, opt)
%LSQVCE.ENGINE.ARCVARIANCEFACTORS  Relative per-arc variance factors within a group.
%
%   out = LSQVCE.ENGINE.ARCVARIANCEFACTORS(omegaT, sT, sigma2)
%   out = LSQVCE.ENGINE.ARCVARIANCEFACTORS(omegaT, sT, sigma2, 'Name', value, ...)
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
%     per-frequency PSD amplitudes (LSQVCE.ENGINE.PSDVCE), so taking alpha2_m * sigma2_m
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
%       .clamped  M x 1 logical  true where sigma2hat_m fell below ClampEps
%                         (negative or degenerate alpha2) and was clamped to
%                         it (WP-35, fixes F-05); all-false on healthy input
%       .meanType char    which mean was applied
%
%   EXAMPLE
%     a = lsqvce.engine.arcVarianceFactors(omegaT, sT, sig2);
%     sig2 = a.sigma2;                 % feed back into the next iteration
%
%   OPTIONS (Name/Value, in addition to 'MeanType'/'TrimFrac')
%     'ClampEps' 1x1 double > 0, default 1e-8.  Floor for the NORMALISED arc
%                factors (sigma2hat is already divided by gammaBar^2, so
%                healthy factors are O(1)): any sigma2hat_m below it --
%                negative gamma2 from a negative alpha2, or a degenerate
%                ratio -- is clamped to ClampEps, warned about with the arc
%                indices (id lsqvce:arcVarianceFactors:clamped) and flagged
%                in .clamped. A negative variance factor must never leave
%                this function: it weights a covariance downstream.
%     mean(sqrt(a.sigma2))             % == 1 exactly only for 'plain'
%
%   ERRORS
%     lsqvce:arcVarianceFactors:arcMismatch  the arc list disagrees with the
%                                   accumulated arcs
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:arcVarianceFactors:badMean          [warning] The normalising mean of the arc factors is X, which cannot scale anyt...
%
%   See also LSQVCE.ENGINE.PSDVCE, LSQVCE.KERNEL.REDUNDANCYLAGS, LSQVCE.KERNEL.DCTI
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (b) b2 (F-07):
%   typed `arguments` signature, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-35 (fixes F-05): negative /
%   degenerate sigma2hat_m are clamped at ClampEps (documented above), the
%   arcs are warned about, and .clamped reports the mask; healthy inputs
%   are byte-identical. With assistance from Claude (Anthropic), model
%   Claude Fable 5.
%   Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    % WP-39 stage (b) b2: typed `arguments` signature (canonical
    % UpperCamelCase; legacy casings keep matching -- MATLAB NV
    % semantics; enum VALUES case-sensitive lowercase).
arguments
    omegaT
    sT
    sigma2
    opt.Norm (1,1) string = "ortho"
    opt.MeanType (1,1) string {mustBeMember(opt.MeanType, ["robust","plain"])} = "robust"
    opt.TrimFrac (1,1) double {mustBeGreaterThanOrEqual(opt.TrimFrac,0), mustBeLessThan(opt.TrimFrac,0.5)} = 0.25
    opt.ClampEps (1,1) double {mustBePositive} = 1e-8
end

    o = struct('Norm', char(opt.Norm), 'MeanType', char(opt.MeanType), 'TrimFrac', opt.TrimFrac, 'ClampEps', opt.ClampEps);

    sigma2 = sigma2(:);
    M      = numel(sigma2);
    if size(omegaT,2) ~= M || size(sT,2) ~= M
        error('lsqvce:arcVarianceFactors:arcMismatch', ...
             ['omegaT has %d and sT has %d columns, but sigma2 lists %d arcs; ' ...
              'one COLUMN per arc is required.'], size(omegaT,2), size(sT,2), M);
    end

    % ---- spectral transform, then cumulative over frequencies    (6.5.27) --
    OmTilde = lsqvce.kernel.dctI(omegaT, 'Norm', o.Norm);
    sTilde  = lsqvce.kernel.dctI(sT    , 'Norm', o.Norm);
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
                     'gammaBar', gammaBar, 'meanType', lower(o.MeanType), ...
                     'clamped', false(size(sigma2)));
        return
    end
    sigma2Hat = gamma2 / gammaBar^2;

    % ---- clamp (WP-35, fixes F-05) ---------------------------------------
    %   gamma2 = alpha2 .* sigma2 can be NEGATIVE (a negative alpha2 on a
    %   pathological arc), and a negative variance factor must never leave
    %   this function: downstream it weights a covariance. Floor at
    %   ClampEps -- an epsilon on the NORMALISED scale (sigma2Hat is already
    %   divided by gammaBar^2, so healthy factors are O(1)) -- warn with the
    %   arc indices and the count, and report the mask. On healthy inputs
    %   the mask is all-false and sigma2Hat is untouched (the no-op
    %   regression test pins that).
    clamped = sigma2Hat < o.ClampEps;
    if any(clamped)
        idx = find(clamped);
        show = idx(1:min(10, numel(idx)));
        lsqvce.util.log.warnOnce('lsqvce:arcVarianceFactors:clamped', ...
           ['%d of %d arc variance factors were < ClampEps = %g (negative ' ...
            'or degenerate alpha2) and were clamped to it. Arc indices: ' ...
            '%s%s. Inspect those arcs -- a persistent clamp usually means ' ...
            'a data problem, not a weighting one.'], ...
            numel(idx), numel(sigma2Hat), o.ClampEps, ...
            mat2str(show(:).'), ternary(numel(idx) > 10, ' ...', ''));
        sigma2Hat(clamped) = o.ClampEps;
    end

    out = struct('sigma2', sigma2Hat, 'alpha2', alpha2, 'gamma2', gamma2, ...
                 'gammaBar', gammaBar, 'meanType', lower(o.MeanType), ...
                 'clamped', clamped);
end

function s = ternary(cond, a, b)
%TERNARY  Inline conditional for the warning text.
    if cond, s = a; else, s = b; end
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
