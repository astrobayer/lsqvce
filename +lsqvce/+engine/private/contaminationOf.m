function D = contaminationOf(D, fhat, nScored, madAvailable)
%CONTAMINATIONOF  diagnose's contamination gauge: fhat's null law and the ladder.
%   D = contaminationOf(D, fhat, nScored, madAvailable)
%
%   fhat is the fraction of the scored residuals beyond 2.5 robust sigmas.
%   This unit turns it into the verdict's contamination level: it states the
%   null law of fhat on clean data at the sample's size, the band a clean
%   sample stays in, whether a collapse of fhat (the robust scale captured by
%   the blunders) could be told at this size, and the level -- 'clean',
%   'moderate', 'severe', 'refuse', 'breakdown', 'indeterminate' or
%   'unknown'. diagnose reports it (D.contamLevel), snoop acts on it
%   (verdict>guardContamination), LSQVCE.REPORT prints it.
%
%   WP-138 (code review of WP-135): lifted out of verdict.m. The verdict
%   class stays one (WP-86 part 4, user decision), but its file changed for
%   many reasons, and this one -- the gauge's statistics, recalibrated in
%   WP-135 and completed in WP-137 -- has its own file now: a change to the
%   law is a change here, against its mirror, and nowhere else.
%
%   INPUT
%     D             1x1 struct   diagnose's verdict so far; reads .backstop
%                                (s0 / sigma0) and .backstopTol (its clean band)
%     fhat          1x1 double   the fraction beyond 2.5 robust sigmas (NaN
%                                when no robust scale exists)
%     nScored       1x1 double   the residuals that were scored
%     madAvailable  1x1 logical  a robust scale exists (else 'unknown')
%   OUTPUT
%     D   1x1 struct   D with .fhatFloor (P(|Z| > 2.5)), .fhatTol (the 99.75 %
%                      point of the null), .collapseDetectable (P(fhat = 0)
%                      <= 1 % at this n), .fhatExcess (fhat - floor, >= 0)
%                      and .contamLevel (char, the level)
%
%   ACCURACY
%     The null law against Monte Carlo on clean data, n = 20..2000: the
%     false-alarm rate of fhatTol 0.06-0.25 % (claimed 0.25 %), P(fhat = 0)
%     within 0.01 at every n. Authority: HANDOVER 6.1.4, the row of the third
%     review's mirrors (dev/validation/wp135_contamination_mirror.py, faithful
%     to this file since WP-137). The normal limit's quantile zHi is
%     normInv(1 - aHi) to 16 digits (test_lsqvce_ext, RelTol 1e-12).
%
%   ERRORS
%     none of its own
%
%   See also LSQVCE.ENGINE.VERDICT, LSQVCE.REPORT, LSQVCE.KERNEL.STAT
%
%   Created 2026-10-15 (CEST, UTC+02:00), WP-138 (code review of WP-135):
%   verdict>contaminationOf, >contamLaw, >contamNull, >fhatCollapsed and
%   >contamCdf moved here verbatim. With assistance from Claude (Anthropic),
%   model Claude Opus 5.5.
    % ---- fhat's NULL DISTRIBUTION -------------------------------------------
    % fhat counts residuals beyond 2.5 robust sigmas, so on PERFECTLY CLEAN data it
    % does not go to zero -- it goes to the tail mass of the normal:
    %       P(|Z| > 2.5) = 0.0124
    % and it fluctuates around that. A FIXED threshold of 0.02 sits INSIDE that
    % null distribution and cries wolf: it escalated 87 of 300 clean datasets (29%)
    % before this was corrected.
    %
    % WP-135 (review 3 finding 19; user decision: an n-dependent recalibration,
    % measured -- dev/validation/wp135_contamination_mirror.py). The band was
    % calibrated at n ~ 300 (P + 4 * 1.30 * binomial sd) and failed at small n:
    % n = 30 clean data raised severe / refuse in 3.1 % (claimed 0.25 %), and
    % "fhat <= 0.005" -- the breakdown fingerprint -- is what HALF of all clean
    % samples of n = 30 show (P(fhat = 0) = 59 %), so a sigma0 declared 3.5x too
    % small read "breakdown" in 50 % and snoop refused. Measured on clean data,
    % n = 20..2000:
    %   * the mean of fhat is P (1 + 27/n): a small sample's MAD is biased low;
    %   * around it the count K = n fhat is over-dispersed by a CONSTANT 1.30^2
    %     (the 1.30 measured at n ~ 300 was the same thing);
    % so K is negative binomial, mean mu = n P (1 + 27/n), variance 1.69 mu
    % (P(fhat = 0) to within 0.01 at every n). Both tails come from that law,
    % n = the scored residuals (contamNull):
    %   fhatTol    = the 99.75 % point (false alarms <= 0.25 %, the rate the
    %                k = 4 band always claimed; measured 0.06-0.25 %)
    %   collapsed  = P(K <= n fhat) <= 1 %: fhat implausibly LOW for clean data.
    %                Below n ~ 450 even fhat = 0 is not; there a robust scale
    %                captured by the blunders cannot be told from a sigma0
    %                declared too small, and the level says so: 'indeterminate'.
    law = contamLaw();
    D.fhatFloor = law.p;
    [D.fhatTol, pZero] = contamNull(nScored, law);
    D.collapseDetectable = pZero <= law.aLo;

    % The EXCESS is what actually matters: how much contamination is there beyond
    % what a clean normal sample would produce anyway?
    if isfinite(fhat)
        D.fhatExcess = max(fhat - law.p, 0);
    else
        D.fhatExcess = NaN;
    end

    % Escalation. fhat collapsed below its own null together with a large
    % backstop is the fingerprint of TOTAL breakdown. Within the null band the
    % data look clean -- a uniformly inflated residual is the global test's
    % finding, not contamination (WP-135: the band now wins over the backstop) --
    % unless the backstop is large and the collapse is undetectable at this n.
    bigBackstop = isfinite(D.backstop) && D.backstop > max(3, 2 * D.backstopTol);
    if ~madAvailable
        D.contamLevel = 'unknown';
    elseif bigBackstop && fhatCollapsed(fhat, nScored, law)
        % fhat says zero but the a-priori backstop says the residuals are enormous:
        % the robust scale has been captured by the blunders.
        D.contamLevel = 'breakdown';
    elseif fhat <= D.fhatTol
        if bigBackstop && ~D.collapseDetectable
            D.contamLevel = 'indeterminate';   % captured scale or sigma0 too small?
        else
            D.contamLevel = 'clean';           % within the null distribution
        end
    elseif fhat >= 0.25
        D.contamLevel = 'refuse';
    elseif fhat >= 0.10
        D.contamLevel = 'severe';
    else
        D.contamLevel = 'moderate';
    end
end

function law = contamLaw()
%CONTAMLAW  The ONE statement of the contamination ladder's null law (WP-135,
%   measured in dev/validation/wp135_contamination_mirror.py; WP-137: named --
%   its numbers were spelled out in four functions). On clean data the count
%   K = n fhat is negative binomial, mean mu = n p (1 + bias / n), variance
%   phi mu.
%   OUTPUT law 1x1 struct  .p .bias .phi .normalMu .aHi .zHi .aLo
    law = struct( ...
        'p',        2 * lsqvce.kernel.stat.normCdf(-2.5), ... % P(|Z| > 2.5) = 0.012419
        'bias',     27, ...                    % measured: the small-sample MAD's excess
        'phi',      1.69, ...                  % measured: the over-dispersion, 1.30^2
        'normalMu', 200, ...                   % above this mean: the normal limit
        'aHi',      0.0025, ...                % fhatTol: the false-alarm rate
        'zHi',      2.807033768343811, ...     % normInv(1 - aHi), to 16 digits
        'aLo',      0.01);                     % collapsed: P(K <= n fhat) at most this
end

function [tol, pZero] = contamNull(nS, law)
%CONTAMNULL  fhat's null law on clean data (contamLaw): tol = the (1 - aHi)
%   point of K / nS; pZero = P(K = 0).
%   INPUT  nS 1x1 the scored residuals; law 1x1 struct (contamLaw)
    if nS < 1, tol = Inf;  pZero = 1;  return, end
    [mu, cdf, pZero] = contamCdf(nS, law, []);
    if mu > law.normalMu                          % the normal limit
        tol = (mu + law.zHi * sqrt(law.phi * mu)) / nS;
    else
        tol = (find(cdf >= 1 - law.aHi, 1) - 1) / nS;
    end
end

function tf = fhatCollapsed(fhat, nS, law)
%FHATCOLLAPSED  P(K <= nS fhat) <= aLo under contamLaw: fhat is implausibly
%   low for clean data (WP-135; it was fhat <= 0.005 at any n).
    tf = false;
    if ~isfinite(fhat) || nS < 1, return, end
    k = round(fhat * nS);
    [mu, cdf] = contamCdf(nS, law, k);
    if mu > law.normalMu
        tf = lsqvce.kernel.stat.normCdf((k + 0.5 - mu) / sqrt(law.phi * mu)) <= law.aLo;
    else
        tf = cdf(k + 1) <= law.aLo;
    end
end

function [mu, cdf, pZero] = contamCdf(nS, law, kMax)
%CONTAMCDF  P(K <= k), k = 0..kMax, of contamLaw's negative binomial (log-space
%   recursion); kMax = [] runs to 12 sd above the mean. Empty above
%   law.normalMu (the normal limit); pZero = P(K = 0) at every mu, from the
%   closed form (WP-135: contamNull read cdf(1) of the empty table and every
%   run of more than ~16 000 scored residuals stopped in diagnose).
    mu = nS * law.p * (1 + law.bias / nS);
    r = mu / (law.phi - 1);  q = mu / (r + mu);  % 1.69 - 1 is 0.69 exactly
    pZero = exp(r * log1p(-q));
    cdf = [];
    if mu > law.normalMu, return, end
    if isempty(kMax), kMax = ceil(mu + 12 * sqrt(law.phi * mu) + 20); end
    k = (1:kMax).';
    lp = r * log1p(-q) + [0; cumsum(log((r + k - 1) ./ k * q))];
    cdf = cumsum(exp(lp));
end
