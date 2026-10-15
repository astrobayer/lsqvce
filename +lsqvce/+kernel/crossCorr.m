function c = crossCorr(w, d, opt)
%LSQVCE.KERNEL.CROSSCORR  Cross-axis / temporal correlation of the whitened residuals.
%
%   c = LSQVCE.KERNEL.CROSSCORR(w, d, 'Name', value, ...)
%
%   THE DIAGNOSTIC THAT TELLS YOU WHETHER YOUR SEPARABLE MODEL IS SAFE.
%
%   Under a CORRECT separable covariance the whitened residuals are white in
%   time AND across axes, so R(0) = I and R(tau>0) = 0. Departures from that
%   diagnose exactly WHICH part of the model is wrong -- and the four cases call
%   for four different fixes, which is why a single "is it white" test is not
%   enough.
%
%   INPUT
%     w    m x 1     whitened residuals, EPOCH-MAJOR: [x1 y1 z1 x2 y2 z2 ...];
%                    whole epochs (m a multiple of d), at least two
%     d    scalar    number of axes (3 for a 3-D position)
%
%   OPTIONS
%     'maxLag'  default 20 epochs
%     'alpha'   significance level (default 0.05 -> band 1.96/sqrt(nEpoch))
%     'plot'    default false
%
%   OUTPUT
%     c.R        (maxLag+1) x d x d   empirical cross-covariance R(tau)
%     c.Rcorr    (maxLag+1) x d x d   the same as a correlation,
%                R(tau)_ij / sqrt(R(0)_ii R(0)_jj) -- what the verdict tests
%     c.band     significance band per element of c.Rcorr
%     c.verdict  'ok' | 'cross-axis-zero-lag' | 'ar-order-low' | 'needs-var'
%     c.advice   what to do about it
%
%   READING THE VERDICT
%     'ok'
%         R(0) = I, R(tau>0) = 0. The separable model is adequate. Done.
%
%     'cross-axis-zero-lag'   R(0) has significant OFF-diagonals
%         The axes are correlated AT THE SAME EPOCH. Your per-epoch covariance
%         is wrong -- almost always because you supplied a diagonal sigma when
%         your kinematic PPP actually gives you the full 3x3 (the vertical is
%         2-3x worse than horizontal AND correlated with it). Fix D, not the
%         temporal model. This is the most common failure and the easiest fix.
%
%     'ar-order-low'          R(tau>0) has a nonzero DIAGONAL
%         Leftover TEMPORAL correlation within each axis. Raise the AR order
%         (or estimate the law: Estimate = "noise").
%
%     'needs-var'             R(tau>0) has significant OFF-diagonals
%         Genuine cross-axis TEMPORAL correlation: axis i at epoch k correlates
%         with axis j at epoch k-tau. Separability FAILS. You need a vector-AR
%         (VAR) model with matrix coefficients. This toolbox does not implement
%         that yet -- but this test tells you empirically whether you need it,
%         rather than assuming either way.
%
%   RUN THIS BEFORE TRUSTING VCE OR OUTLIER DETECTION
%     Both assume the weighting is right. If R(0) is far from I, both are biased.
%
%   ACCURACY
%     c.R is a plain mean of lag products and c.Rcorr its normalisation by
%     the lag-0 standard deviations: no approximation of their own. The
%     verdict is statistical, its band z(1 - alpha/2)/sqrt(nEpoch) that of
%     white residuals: on simulated series of 4000-6000 epochs white noise
%     reads 'ok' at sigma 1 and 2 and each of the three departures is named
%     (test_vceMath (testCrossCorrWhiteIsOK, testCrossCorrVerdictIsScaleFree,
%     testCrossCorrDetectsCrossAxis, testCrossCorrDetectsARTooLow,
%     testCrossCorrDetectsNeedsVAR)). WP-139.
%
%   EXAMPLE
%     res = acc.residuals(sol);
%     w   = vertcat(res.sys(1).obs.w);
%     c   = lsqvce.kernel.crossCorr(w, 3, 'plot', true);
%
%     lsqvce.util.log.msg('%s: %s\n', c.verdict, c.advice);
%     disp(squeeze(c.Rcorr(1,:,:)))      % R(0) as a correlation -- ~ the identity
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (b) b2 (F-07):
%   typed `arguments` signature, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   ERRORS
%     lsqvce:crossCorr:badLength   numel(w) is not a whole number of epochs
%                                  of d axes, or fewer than two epochs (WP-134:
%                                  a partial epoch was dropped silently, and
%                                  one epoch raised an unnamed index error)
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:crossCorr:shortSeries               [warning] Only X epochs for a max lag of X. The cross-covariance estimates will...
%
%   See also LSQVCE.ENGINE.RESIDUALS, LSQVCE.ENGINE.COVMODEL
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-139: the help states its
%   accuracy, every number traced to the test that pins it. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-134: a partial epoch or fewer
%   than two epochs is refused by name (lsqvce:crossCorr:badLength). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-131: the verdict tests the
%   correlation (c.Rcorr) against the band, not the covariance -- it no longer
%   depends on the residuals' variance. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-118 part 1 (R7-34): the
%   deleted SELECTARORDER references removed. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+kernel (R4 decision R4-D1,
%   WP-64a): the KERNEL tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

arguments
    w
    d
    opt.MaxLag (1,1) double {mustBeInteger, mustBePositive} = 20
    opt.Alpha (1,1) double {mustBePositive, mustBeLessThan(opt.Alpha,1)} = 0.05
    opt.Plot (1,1) logical = false
end

    narginchk(2,inf);

    % WP-39 stage (b) b2: typed `arguments` signature (canonical
    % UpperCamelCase; legacy casings keep matching -- MATLAB NV
    % semantics; enum VALUES case-sensitive lowercase).
    o = struct('maxLag', opt.MaxLag, 'alpha', opt.Alpha, 'plot', opt.Plot);

    validateattributes(w, {'double'}, {'vector','finite','nonempty'}, mfilename, 'w', 1);
    validateattributes(d, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'd', 2);

    w  = w(:);
    d  = double(d);
    nE = floor(numel(w)/d);
    if mod(numel(w), d) ~= 0 || nE < 2
        error('lsqvce:crossCorr:badLength', ...
              ['w has %d values: %d whole epochs of %d axes are needed (rows epoch-major), ' ...
               'and at least two of them.'], numel(w), nE, d);
    end
    if nE < 10*o.maxLag
        lsqvce.util.log.warn('lsqvce:crossCorr:shortSeries', ...
            ['Only %d epochs for a max lag of %d. The cross-covariance estimates ' ...
             'will be noisy.'], nE, o.maxLag);
    end

    W = reshape(w(1:nE*d), d, nE).';           % nEpoch x d, epoch-major

    mx = min(o.maxLag, nE-2);
    R  = zeros(mx+1, d, d);
    for tau = 0:mx
        R(tau+1,:,:) = (W(1:nE-tau,:).' * W(1+tau:end,:)) / (nE - tau);
    end

    zcrit = norminvApprox(1 - o.alpha/2);
    band  = zcrit / sqrt(nE);

    % WP-131: the band z/sqrt(nE) is that of a CORRELATION; the verdict tested
    % the covariance R, so its false-alarm rate and power scaled with the
    % residuals' variance (white noise at sigma 2 failed every time, a real
    % AR(1) leftover at sigma 0.5 was never seen; dev/validation/wp131_kernel_mirror.py).
    % The verdict reads Rn = R(tau)_ij / sqrt(R0_ii R0_jj); c.R stays the covariance.
    sd = sqrt(max(diag(reshape(R(1,:,:), d, d)), realmin));
    Rn = R ./ reshape(sd * sd.', [1 d d]);

    R0     = reshape(Rn(1,:,:), d, d);
    cross0 = max(max(abs(R0 - diag(diag(R0)))));

    diagT  = 0;  crossT = 0;
    for t = 2:mx+1
        Rt     = reshape(Rn(t,:,:), d, d);
        diagT  = max(diagT,  max(abs(diag(Rt))));
        crossT = max(crossT, max(max(abs(Rt - diag(diag(Rt))))));
    end

    % ---- verdict, in priority order --------------------------------------
    verdict = 'ok';
    advice  = 'the separable model is adequate';

    if cross0 > 3*band
        verdict = 'cross-axis-zero-lag';
        advice  = ['the axes are correlated at the same epoch -- supply the ' ...
                   'full d x d per-epoch covariance instead of a diagonal sigma'];
    elseif crossT > 3*band
        verdict = 'needs-var';
        advice  = ['genuine cross-axis temporal correlation -- separability ' ...
                   'fails, a vector-AR (VAR) model is required'];
    elseif diagT > 3*band
        verdict = 'ar-order-low';
        advice  = 'leftover temporal correlation -- raise the AR order';
    end

    c = struct('R', R, 'Rcorr', Rn, 'band', band, 'verdict', verdict, 'advice', advice, ...
               'cross0', cross0, 'diagT', diagT, 'crossT', crossT, ...
               'nEpoch', nE, 'd', d);

    if ~strcmp(verdict, 'ok')
        wid = sprintf('lsqvce:crossCorr:%s', strrep(verdict,'-','_'));
        lsqvce.util.log.warn(wid, '%s', advice);
    end

    if o.plot
        figure('Name', sprintf('crossCorr: %s', verdict));
        for i = 1:d
            for j = 1:d
                subplot(d, d, (i-1)*d + j);
                stem(0:mx, squeeze(Rn(:,i,j)), '.'); hold on
                yline([-1 1]*band, 'r--');
                grid on;  ylim([-1 1]*max(0.3, 1.2*max(abs(squeeze(Rn(:,i,j))))));
                title(sprintf('R_{%d%d}', i, j));
                if i == d, xlabel('lag [epochs]'); end
            end
        end
        sgtitle(sprintf('%s -- %s', verdict, advice));
    end
end

%==========================================================================
function z = norminvApprox(p)
%NORMINVAPPROX  Inverse normal CDF, no Statistics Toolbox needed.
    z = sqrt(2) * erfinv(2*p - 1);
end
