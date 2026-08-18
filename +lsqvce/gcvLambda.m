function [lam, gcv, grid] = gcvLambda(N, n, L)
%LSQVCE.GCVLAMBDA  Generalised cross-validation choice of the Tikhonov weight.
%
%   lam = LSQVCE.GCVLAMBDA(N, n, L) minimises the GCV functional over a log grid,
%   for the regularised problem (N + lambda*L'*L) y = n.
%
%   [lam, gcv, grid] = LSQVCE.GCVLAMBDA(...) also returns the GCV curve, so you
%   can check that the minimum is a real trough and not a grid endpoint.
%
%   INPUT
%     N   nl x nl  local normal matrix (symmetric positive semidefinite)
%     n   nl x 1   local right-hand side
%     L   nl x nl  regularisation matrix (default identity)
%
%   OUTPUT
%     lam   scalar   selected weight
%     gcv   ng x 1   GCV functional over the grid
%     grid  ng x 1   the lambda grid (80 points, log-spaced)
%
%   METHOD
%     Uses the generalised eigendecomposition of (N, L'*L), which diagonalises
%     both matrices simultaneously, so every lambda on the grid costs O(nl)
%     instead of a fresh factorisation.
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
%     [lam, gcv, grid] = lsqvce.gcvLambda(Nll, nl, eye(size(Nll,1)));
%     loglog(grid, gcv); xline(lam,'r--'); xlabel('\lambda'); ylabel('GCV');
%
%   See also LSQVCE.APPLYREGULARISATION
%
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.

    narginchk(2,3);
    nlp = size(N,1);
    if nargin < 3 || isempty(L), L = eye(nlp); end

    validateattributes(N, {'double'}, {'2d','square','finite'}, mfilename, 'N', 1);
    validateattributes(n, {'double'}, {'vector','finite','numel',nlp}, mfilename, 'n', 2);

    N   = (N + N.')/2;
    LtL = L.'*L;
    LtL = (LtL + LtL.')/2 + eps*eye(nlp);       % guarantee PD for the chol form

    try
        [V, D] = eig(N, LtL, 'chol');
    catch
        [V, D] = eig(N, LtL);                    % fall back to the QZ form
    end
    d = real(diag(D));
    b = real(V.') * n(:);

    dmax = max(d);
    if ~isfinite(dmax) || dmax <= 0
        lsqvce.util.log.warn('lsqvce:gcv:degenerate', ...
            'The local block has no positive generalised eigenvalue; lambda = 0.');
        lam = 0;  gcv = [];  grid = [];
        return
    end

    grid = logspace(log10(dmax*1e-12), log10(dmax), 80).';
    gcv  = zeros(size(grid));
    for k = 1:numel(grid)
        f      = grid(k) ./ (d + grid(k));       % residual filter factors
        gcv(k) = norm(f.*b)^2 / max(sum(f)^2, realmin);
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
