function [du, dt] = truncationDeficit(P, X, excl, del, tauBar, nL)
%LSQVCE.KERNEL.TRUNCATIONDEFICIT  The redundancy a robust sample loses to its truncation, lag by lag.
%
%   [du, dt] = LSQVCE.KERNEL.TRUNCATIONDEFICIT(P, X, excl, del, tauBar, nL)
%
%   The robust noise estimator deletes the epochs it rejects, and the epochs
%   it KEEPS are a Gaussian sample truncated at the rejection point in their
%   own w direction (LSQVCE.KERNEL.WSNOOP). For a kept epoch j that direction
%   is R e_j / sqrt(R_jj), so the expected outer product of etilde = P e (what
%   the per-lag Omega reads) loses
%       (1 - tau) R D^-1 R,     D = diag(R) over the kept epochs,
%   and the estimator's expectation of Omega_n is s_n - du_n instead of s_n.
%   This returns du in LSQVCE.KERNEL.REDUNDANCYLAGS' convention for s (lag 0,
%   then twice each sub-diagonal sum) and dt, the deficit of the whole law's
%   Helmert pair, trace(R D^-1 R T) = the kept count (R T R = R for this R).
%   Subtracting them from the redundancy makes every lag -- the law's shape,
%   its level and sigma0 -- consistent to first order in the rejection rate
%   for ANY law (WP-156, design F): a uniform division by tau is exact for the
%   total only, and over-corrected a red stream's marginal variance by up to
%   1 / tau - 1 = 2.7 % (HANDOVER 12.91). The arc is d single-axis series
%   sharing the temporal law (the separable model); du and dt are summed
%   over the axes.
%
%   INPUT
%     P       m x m   inverse temporal law (LSQVCE.KERNEL.LAWINVERSE)
%     X       1 x d   cell, X{k} m x n_k: the parameters of axis k estimated
%                     with this law (the arc's local design, descaled)
%     excl    m x d   logical, rows excluded from the data
%     del     m x 1   logical, the epochs the robust pass rejected (deleted)
%     tauBar  1x1     1 - tau, tau the truncation factor of one axis at the
%                     rejection point (0: no deficit)
%     nL      1x1     lags wanted (<= m)
%
%   OUTPUT
%     du   nL x 1  the deficit of s, summed over the axes
%     dt   1x1     the deficit of the whole law's trace, tauBar x kept rows
%
%   COST. O(m^2) to form R per axis and O(m^2 log m) for the lag sums (one
%   FFT per column, in blocks of 512 columns).
%
%   EXAMPLE
%     Sig = toeplitz(0.9 .^ (0:127));  P = lsqvce.kernel.lawInverse(Sig);
%     du = lsqvce.kernel.truncationDeficit(P, {ones(128, 1)}, false(128, 1), ...
%                                          false(128, 1), 1 - 0.9733, 128);
%
%   ACCURACY
%     Exact linear algebra on P: du equals the lag sums of the dense (1 -
%     tau) R D^-1 R to 1e-12 and dt the kept count to 1e-10 (test_precision
%     (testTruncationDeficitMatchesDenseForm)); with it the per-lag Omega of
%     4000 snooped AR(1) arcs, phi 0.9, matches s - du at lags 0 and 1 within
%     four Monte Carlo standard errors, where tau s misses lag 1 by 3 %
%     (testDeficitMakesEveryLagConsistent). In the full loop (20 seeds,
%     paired against Robust "none"): a red stream's variance -0.10 % at phi
%     0.9 where the uniform division read +2.4 % (HANDOVER 12.91).
%
%   ERRORS
%     lsqvce:truncationDeficit:badSize   X, excl or del do not match P
%
%   See also LSQVCE.KERNEL.WSNOOP, LSQVCE.KERNEL.LAWINVERSE,
%            LSQVCE.KERNEL.REDUNDANCYLAGS
%
%   Created 2026-10-18 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- WP-156 round 2, design F.

    m = size(P, 1);
    d = size(excl, 2);
    if ~isequal(size(P), [m m]) || size(excl, 1) ~= m || numel(del) ~= m || numel(X) ~= d
        error('lsqvce:truncationDeficit:badSize', ...
              'P must be square and excl / del / X match its %d rows (got excl %d x %d, del %d, %d designs).', ...
              m, size(excl, 1), d, numel(del), numel(X));
    end
    du = zeros(nL, 1);  dt = 0;
    if tauBar == 0, return, end
    excl = logical(excl);  del = logical(del(:));
    gone = all(excl, 2) | del;
    PD = P;
    for j = find(gone).'
        pj = PD(:, j);
        if pj(j) > 0, PD = PD - (pj * pj.') / pj(j); end
        PD(:, j) = 0;  PD(j, :) = 0;
    end
    for k = 1:d
        part = find(excl(:, k) & ~gone);
        Xk = [full(sparse(part, 1:numel(part), 1, m, numel(part))), X{k}];
        Xk(gone, :) = 0;
        Xk = Xk(:, any(Xk ~= 0, 1));
        R = PD;
        if ~isempty(Xk)
            G = PD * Xk;
            [U, S] = eig((Xk.' * G + G.' * Xk) / 2);
            s = diag(S);
            keep = s > 1e3 * eps * max([s; 0]);
            T = G * (U(:, keep) ./ sqrt(s(keep)).');
            R = R - T * T.';
        end
        kept = find(~gone & ~excl(:, k) & diag(R) > 0);
        u = zeros(nL, 1);
        sq = sqrt(diag(R));
        for c0 = 1:512:numel(kept)
            cc = kept(c0:min(c0 + 511, numel(kept)));
            u = u + columnLagSums(R(:, cc) ./ sq(cc).', nL);
        end
        du = du + tauBar * [u(1); 2 * u(2:end)];
        dt = dt + tauBar * numel(kept);
    end
end
