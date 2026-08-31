function [W, info] = decorrFactor(C, tol)
%LSQVCE.KERNEL.DECORRFACTOR  Decorrelation factor W with W'*W = inv(C).
%
%   W = LSQVCE.KERNEL.DECORRFACTOR(C) returns the upper-triangular W = chol(C)^-1, so
%   that whitening a batch is Aw = W*A and the weighted normals A'*inv(C)*A
%   become Aw'*Aw.
%
%   [W, info] = LSQVCE.KERNEL.DECORRFACTOR(C) also returns
%     info.method   'chol' | 'eig-floored'
%     info.rcond    reciprocal condition estimate of C
%     info.floored  number of eigenvalues floored (0 if method is 'chol')
%
%   W = LSQVCE.KERNEL.DECORRFACTOR(C, tol) sets the relative eigenvalue floor used only
%   by the rescue path (default 1e-10).
%
%   INPUT
%     C    n x n   covariance of one batch, symmetric positive definite
%     tol  scalar  relative eigenvalue floor for the rescue path, default 1e-10
%
%   COST
%     O(n^3) and O(n^2) memory. This is the expensive path. If the noise is
%     AR or ARMA, use the filter form (LSQVCE.KERNEL.ARMAWHITEN) instead: O(n*p) and no
%     matrix at all, with identical normal equations. Only use this factor for
%     covariances with no rational spectral representation.
%
%   BATCH SIZE
%     At n = 5000 the matrix is 200 MB and the Cholesky is ~40 GFLOP. That,
%     not memory, is the practical ceiling on BatchSize for the dense path.
%
%   NOT POSITIVE DEFINITE
%     If the Cholesky fails, a factor is built from the SYMMETRIC eigen-
%     decomposition with the eigenvalues floored at tol*lambda_max, and a
%     warning is issued. This rescues the run but the weighting is no longer
%     the one you specified -- an indefinite covariance usually means the PSD
%     had negative values, or the batch straddles a data gap with an
%     inconsistent model.
%
%     The floor must NOT be eps: flooring at eps*lambda_max leaves C with
%     condition number ~1/eps, so the resulting W satisfies W*C*W' = I only to
%     ~1e-1. The default 1e-10 gives ~1e-7, which is usable. Tighten tol only
%     if you understand why the covariance was indefinite.
%
%   OUTPUT
%     W     m x m  upper triangular, W'*W = inv(C)
%     info  1x1 struct  .method ('chol' | 'eig-floored'), .rcond, .floored
%
%   EXAMPLE
%     r = lsqvce.kernel.armaACV([1 -0.9], 1, 100);
%     W = lsqvce.kernel.decorrFactor(toeplitz(r));
%     norm(W.'*W - inv(toeplitz(r))) / norm(inv(toeplitz(r)))    % ~ 1e-14
%
%   ERRORS
%     lsqvce:decorrFactor:notPD     C is not positive definite even after the
%                                   eigenvalue floor; it cannot be a covariance
%
%   See also LSQVCE.KERNEL.ARMAWHITEN, LSQVCE.ENGINE.COVMODEL
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-37: diagnostics (incl. the
%   condest estimate) are computed only when info is requested
%   (nargout > 1) -- the per-batch hot path pays nothing. With assistance
%   from Claude (Anthropic), model Claude Fable 5.
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+kernel (R4 decision R4-D1,
%   WP-64a): the KERNEL tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    narginchk(1,2);
    if nargin < 2 || isempty(tol), tol = 1e-10; end
    validateattributes(C, {'double','single'}, {'2d','square','finite','nonempty'}, ...
                       mfilename, 'C', 1);
    validateattributes(tol, {'double'}, {'scalar','positive','<',1}, mfilename, 'tol', 2);

    n = size(C,1);
    C = (C + C.')/2;                     % enforce symmetry

    % WP-37 (hot-loop toll, F-09 rider): the diagnostics -- above all the
    % condest() estimate -- are computed ONLY when the caller asks for info
    % (nargout > 1). The accumulator calls this per batch with one output;
    % condest on every healthy batch was pure overhead.
    wantInfo = nargout > 1;
    if wantInfo
        info = struct('method','chol', 'rcond',NaN, 'floored',0);
    end

    [L, fail] = chol(C, 'lower');
    if fail == 0
        W = L \ eye(n);                  % triangular solve, never inv()
        if wantInfo
            info.rcond = 1/(condest(C) + realmin);
        end
        return
    end

    % --- rescue path -------------------------------------------------------
    % Symmetric eigendecomposition: V is orthonormal, so
    %     C = V*diag(d)*V'   ->   inv(C) = V*diag(1./d)*V'
    %     W = diag(1./sqrt(d))*V'   satisfies   W'*W = inv(C)   exactly.
    % Using the general eig() here would give a non-orthonormal V and a wrong
    % factor.
    [V, d] = eig(C, 'vector');
    d      = real(d);
    dmax   = max(d);
    if dmax <= 0
        error('lsqvce:decorrFactor:notPD', ...
            ['The batch covariance has no positive eigenvalue (max = %g). ' ...
             'It cannot be a covariance.'], dmax);
    end

    fl  = tol*dmax;
    nfl = sum(d < fl);
    lsqvce.util.log.warn('lsqvce:decorrFactor:notPD', ...
        ['The batch covariance is not positive definite (min eigenvalue ' ...
         '%.3e; %d of %d eigenvalues below the floor %.3e). Flooring and ' ...
         'continuing, but the applied weighting is no longer the model you ' ...
         'specified. Check the PSD for negative values.'], min(d), nfl, n, fl);

    d(d < fl)    = fl;
    % W = diag(1./sqrt(d)) * V.'  -- use column broadcasting to avoid forming the
    % dense n x n diagonal and the O(n^3) product; this scales each ROW of V.' by
    % 1/sqrt(d), which is exactly diag(1./sqrt(d))*V.', at O(n^2).
    W            = (1./sqrt(d)) .* V.';
    if wantInfo
        info.method  = 'eig-floored';
        info.floored = nfl;
        info.rcond   = min(d)/dmax;
    end
end
