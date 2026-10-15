function [lam, gcv, grid] = gcvLambda(N, n, L, m, ltl)
%LSQVCE.ENGINE.GCVLAMBDA  Generalised cross-validation choice of the Tikhonov weight.
%
%   lam = LSQVCE.ENGINE.GCVLAMBDA(N, n, L, m, ltl) minimises the GCV functional
%   (Golub, Heath & Wahba 1979) of the local least-squares problem
%       min ||B y - l||^2 + lambda ||L y||^2,   N = B'B, n = B'l, ltl = l'l,
%   over a log grid:
%       G(lambda) = ||B y_lambda - l||^2 / (m - tr H_lambda)^2.
%   The commons enter at their a-priori values (x = 0 of the linearised
%   problem): lambda is chosen per system, before x is known.
%
%   [lam, gcv, grid] = LSQVCE.ENGINE.GCVLAMBDA(...) also returns the GCV curve, so you
%   can check that the minimum is a real trough and not a grid endpoint.
%
%   INPUT
%     N    nl x nl  local normal matrix B'B (symmetric positive semidefinite)
%     n    nl x 1   local right-hand side B'l
%     L    q x nl   regularisation matrix ([] = identity); it may have a null
%                   space (a first difference, a subset): those directions
%                   are fitted freely at every lambda
%     m    1x1      rows of the local problem (whitened observations)
%     ltl  1x1      l'l of the whitened observations
%
%   OUTPUT
%     lam   scalar   selected weight
%     gcv   ng x 1   GCV functional over the grid
%     grid  ng x 1   the lambda grid (80 points, log-spaced over the
%                    eigenvalues of the penalised part, 1e-12 * max .. max)
%
%   METHOD
%     The standard form (Elden 1977): the eigenvectors of L'L split nl into
%     null(L) (eigenvalues below 1e3 nl eps |L'L|) and its complement W. The
%     null(L) part is unpenalised and eliminated like a common (a Schur step
%     on N, n and l'l); on W, L'L = diag(e), and one symmetric eig of
%     diag(e)^-1/2 N_W diag(e)^-1/2 = Q diag(d) Q' diagonalises the rest. With
%     b = Q' diag(e)^-1/2 n_W and f = lambda ./ (d + lambda):
%         ||B y - l||^2 = rss_perp + sum (b.^2 ./ d) .* f.^2,
%         tr H          = nl - sum f,
%     rss_perp the unregularised residual. Every lambda costs O(nl).
%
%   ACCURACY
%     The curve against the dense GCV (direct solves of N + lambda L'L):
%     worst case ~cond(N) eps relative, measured 9e-13 at cond(N) = 1.5e7 and
%     6e-10 at 1.5e10; the same grid point for L = I, a first difference and a
%     subset. Authority: HANDOVER 6.1.4, the row of the third review's mirrors
%     (dev/validation/wp135_gcv_mirror.py); the test pins 1e-9 at 1.5e7.
%
%   WHEN GCV FAILS
%     GCV assumes the residuals are white. After only approximate whitening
%     they are not, and GCV then systematically UNDER-regularises -- it reads
%     the leftover correlation as signal. Symptoms: lam lands on the smallest
%     grid point, or the GCV curve is monotone with no interior minimum. Both
%     are detected and warned about here.
%
%     For production, do not call GCV per day. Fix lambda once from an L-curve
%     on a representative day and hold it constant, otherwise the amount of
%     regularisation drifts between days and the estimated common parameters
%     inherit that drift as a spurious time signal.
%
%   EXAMPLE
%     B = [ones(50,1), (1:50).'/50];  l = B*[1; 2] + 0.1*randn(50,1);
%     [lam, gcv, grid] = lsqvce.engine.gcvLambda(B.'*B, B.'*l, [], 50, l.'*l);
%     loglog(grid, gcv); xline(lam,'r--'); xlabel('\lambda'); ylabel('GCV');
%
%   ERRORS
%     lsqvce:gcv:badData          m is not an integer > nl, or ltl is not a
%                                 finite scalar >= 0
%     lsqvce:gcv:nullUnobserved   the data do not determine a direction
%                                 of null(L): no lambda regularises it
%     (WP-138: one segment, gcv, beside the warnings; both were raised as
%     lsqvce:gcvLambda:* in 6.22.0 - 6.23.0)
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:gcv:atLowerBound                    [warning] GCV selected the smallest lambda on the grid (X), i.e
%     lsqvce:gcv:atUpperBound                    [warning] GCV selected the largest lambda on the grid (X)
%     lsqvce:gcv:degenerate                      [warning] L penalises no direction, or the data none; lambda = 0
%
%   See also LSQVCE.ENGINE.APPLYREGULARISATION
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-138 (code review of WP-135): its
%   errors under the one segment of its warnings, lsqvce:gcv:badData and
%   lsqvce:gcv:nullUnobserved (they were lsqvce:gcvLambda:*). With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-137: ACCURACY states the
%   measured and the worst-case error and cites its HANDOVER 6.1.4 row. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-14 (CEST, UTC+02:00), WP-135 (review 3 finding 16): the
%   GCV of the local least-squares problem (it minimised ||f.*b||^2/(sum
%   f)^2 of N y = n, as if N were the design: |y - y_true| 9 to 2770 times the
%   GCV choice's), with m and l'l from the pre-elimination; the standard form
%   for a RegMatrix with a null space (the eps ridge put the whole grid above
%   1e3 |N| and corrupted the finite eigenvalues); the QZ fallback is gone
%   (one symmetric eig). With assistance from Claude (Anthropic), model
%   Claude Opus 5.5.
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    narginchk(5,5);
    nlp = size(N,1);
    if isempty(L), L = eye(nlp); end

    validateattributes(N, {'double'}, {'2d','square','finite'}, mfilename, 'N', 1);
    validateattributes(n, {'double'}, {'vector','finite','numel',nlp}, mfilename, 'n', 2);
    if ~(isnumeric(m) && isscalar(m) && isfinite(m) && m == round(m) && m > nlp) || ...
            ~(isnumeric(ltl) && isscalar(ltl) && isfinite(ltl) && ltl >= 0)
        error('lsqvce:gcv:badData', ...
              'm must be an integer > nl (%d) and ltl a finite scalar >= 0.', nlp);
    end
    m = double(m);  ltl = double(ltl);

    N   = (N + N.')/2;
    LtL = L.'*L;
    [U, e] = eig((LtL + LtL.')/2, 'vector');
    e = real(e);  U = real(U);
    emax = max([e; 0]);
    pen = e > 1e3 * nlp * eps * emax;          % tol: eig(L'L) is exact to ~nl eps |L'L|
    if emax <= 0 || ~any(pen)
        lsqvce.util.log.warn('lsqvce:gcv:degenerate', ...
            'L penalises no direction, or the data none; lambda = 0.');
        lam = 0;  gcv = [];  grid = [];
        return
    end
    W = U(:, pen);  Z = U(:, ~pen);  sc = 1 ./ sqrt(e(pen));
    Nw = W.'*N*W;  nw = W.'*n(:);
    if ~isempty(Z)                              % null(L): fitted freely -- a Schur step
        [cz, p] = chol(Z.'*N*Z);
        if p ~= 0
            error('lsqvce:gcv:nullUnobserved', ...
                  ['The data do not determine a direction L leaves unpenalised: no ' ...
                   'lambda regularises it. Penalise it (extend L) or remove it.']);
        end
        G  = cz.' \ (Z.'*N*W);                % cz' \ Nzw
        g  = cz.' \ (Z.'*n(:));
        Nw = Nw - G.'*G;  nw = nw - G.'*g;  ltl = ltl - g.'*g;
    end
    S = sc .* Nw .* sc.';
    [Q, d] = eig((S + S.')/2, 'vector');
    d = real(d);  b = real(Q).' * (sc .* nw);
    dmax = max(d);
    if ~(dmax > 0)
        lsqvce.util.log.warn('lsqvce:gcv:degenerate', ...
            'L penalises no direction, or the data none; lambda = 0.');
        lam = 0;  gcv = [];  grid = [];
        return
    end
    ok = d > nlp * eps * dmax;                  % identified by the data
    q  = zeros(size(d));  q(ok) = b(ok).^2 ./ d(ok);
    rssPerp = max(ltl - sum(q), 0);             % the unregularised residual

    grid = logspace(log10(dmax*1e-12), log10(dmax), 80).';
    gcv  = zeros(size(grid));
    for k = 1:numel(grid)
        f = ones(size(d));                      % an unidentified direction: fully damped
        f(ok) = grid(k) ./ (d(ok) + grid(k));   % residual filter factors
        gcv(k) = (rssPerp + sum(q .* f.^2)) / max((m - nlp + sum(f))^2, realmin);
    end

    [~, k] = min(gcv);
    lam    = grid(k);

    if k == 1
        lsqvce.util.log.warn('lsqvce:gcv:atLowerBound', ...
            ['GCV selected the smallest lambda on the grid (%.3e), i.e. no ' ...
             'regularisation. This usually means the whitened residuals are ' ...
             'still correlated, so GCV reads leftover colour as signal. Fix ' ...
             'lambda from an L-curve instead.'], lam);
    elseif k == numel(grid)
        lsqvce.util.log.warn('lsqvce:gcv:atUpperBound', ...
            ['GCV selected the largest lambda on the grid (%.3e). The local ' ...
             'block carries almost no information -- consider removing the ' ...
             'local parameters instead of damping them to zero.'], lam);
    end
end
