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
%     w    m x 1     whitened residuals, EPOCH-MAJOR: [x1 y1 z1 x2 y2 z2 ...]
%     d    scalar    number of axes (3 for a 3-D position)
%
%   OPTIONS
%     'maxLag'  default 20 epochs
%     'alpha'   significance level (default 0.05 -> band 1.96/sqrt(nEpoch))
%     'plot'    default false
%
%   OUTPUT
%     c.R        (maxLag+1) x d x d   empirical cross-covariance R(tau)
%     c.band     significance band per element
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
%         Leftover TEMPORAL correlation within each axis. Raise the AR order.
%         See LSQVCE.ENGINE.SELECTARORDER.
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
%   EXAMPLE
%     res = acc.residuals(sol);
%     w   = vertcat(res.sys(1).obs.w);
%     c   = lsqvce.kernel.crossCorr(w, 3, 'plot', true);
%
%     lsqvce.util.log.msg('%s: %s\n', c.verdict, c.advice);
%     disp(squeeze(c.R(1,:,:)))          % R(0) -- should be the identity
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (b) b2 (F-07):
%   typed `arguments` signature, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:crossCorr:shortSeries               [warning] Only X epochs for a max lag of X. The cross-covariance estimates will...
%
%   See also LSQVCE.ENGINE.SELECTARORDER, LSQVCE.ENGINE.RESIDUALS, LSQVCE.ENGINE.COVMODEL
%
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

    R0     = squeeze(R(1,:,:));
    cross0 = max(max(abs(R0 - diag(diag(R0)))));

    diagT  = 0;  crossT = 0;
    for t = 2:mx+1
        Rt     = squeeze(R(t,:,:));
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

    c = struct('R', R, 'band', band, 'verdict', verdict, 'advice', advice, ...
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
                stem(0:mx, squeeze(R(:,i,j)), '.'); hold on
                yline([-1 1]*band, 'r--');
                grid on;  ylim([-1 1]*max(0.3, 1.2*max(abs(squeeze(R(:,i,j))))));
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
