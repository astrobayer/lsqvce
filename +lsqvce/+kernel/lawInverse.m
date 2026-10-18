function P = lawInverse(Sig)
%LSQVCE.KERNEL.LAWINVERSE  The explicit inverse of an arc's covariance law.
%
%   P = LSQVCE.KERNEL.LAWINVERSE(Sig)
%
%   Returns P = inv(Sig) for a symmetric positive definite N x N covariance.
%   A stationary law -- every covModel ACF / ARMA / PSD law, a Toeplitz Sig --
%   is inverted through its structure: the first column x of inv(Sig) from the
%   Levinson-Durbin recursion (LSQVCE.KERNEL.LEVINSONDURBIN), then the rest by
%   the Trench recursion of the Gohberg-Semencul form
%       P(i+1, j+1) = P(i, j) + (x(i+1) x(j+1) - z(i+1) z(j+1)) / x(1),
%   z = [0; x(N:-1:2)] -- O(N^2) time instead of the O(N^3) of a dense
%   inverse. Any other Sig goes through its Cholesky factor. The robust noise
%   estimator needs P itself (WP-156): Baarda's w per epoch reads its
%   diagonal, an epoch is deleted from it exactly by a rank-one update, and
%   the truncation deficit reads its columns.
%
%   INPUT
%     Sig  N x N double  symmetric positive definite covariance
%
%   OUTPUT
%     P    N x N double  inv(Sig), symmetric
%
%   STRUCTURE TEST AND GUARD
%     Sig is taken as Toeplitz when every diagonal is constant to 1e-12 of
%     Sig(1,1). The structured result is then checked on three columns
%     (first, middle, last): if max |Sig * P(:, j) - e_j| exceeds 1e-8 --
%     a law so ill-conditioned that the recursion lost digits -- the
%     Cholesky path is taken instead, so the fast path never returns a worse
%     inverse than the dense one.
%
%   EXAMPLE
%     Sig = toeplitz(0.95 .^ (0:999));
%     P = lsqvce.kernel.lawInverse(Sig);
%     norm(Sig * P - eye(1000), 'fro')        % ~1e-12
%
%   ACCURACY
%     The structured path equals the dense inverse to 1e-11 relative in norm
%     (AR(1) phi 0.95 and 0.995, N 500: measured 6e-14 and 1.5e-12; a PSD law
%     from psd2acv to 1e-10), Sig -> c Sig gives P / c to 1e-12, and a
%     non-Toeplitz Sig (a per-row variance) takes the Cholesky path, equal to
%     inv(Sig) to 1e-12 (test_precision (testLawInverseMatchesDenseInverse,
%     testLawInverseTakesTheDensePathOffToeplitz,
%     testRobustKernelIsScaleEquivariant)). Measured: 4e-14 at
%     N = 1000-4000, phi 0.95; 0.06 / 0.13 / 0.25 s at N = 2000 / 3000 /
%     4000 against 0.08 / 0.25 / 0.48 s for the dense inverse, structure test
%     and guard included (test_speedBudget testLawInverseBeatsTheDenseInverse:
%     < 0.85 of it at N = 3000; HANDOVER 12.91).
%
%   ERRORS
%     lsqvce:lawInverse:notSquare   Sig is not a square matrix
%     lsqvce:lawInverse:notPD       Sig is not positive definite
%
%   See also LSQVCE.KERNEL.LEVINSONDURBIN, LSQVCE.KERNEL.WSNOOP,
%            LSQVCE.KERNEL.TRUNCATIONDEFICIT, LSQVCE.KERNEL.REDUNDANCYLAGS
%
%   Created 2026-10-18 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- WP-156 round 2, design F.

    [N, N2] = size(Sig);
    if N ~= N2
        error('lsqvce:lawInverse:notSquare', 'Sig must be square (got %d x %d).', N, N2);
    end
    Sig = (Sig + Sig.') / 2;
    if N == 0, P = zeros(0); return, end
    if isdiag(Sig)                                     % a 'diag' / 'var' law: no coupling
        dg = diag(Sig);
        if any(~(dg > 0))
            error('lsqvce:lawInverse:notPD', 'Sig is not positive definite (a diagonal entry <= 0).');
        end
        P = diag(1 ./ dg);
        return
    end
    if N >= 3 && isToeplitz(Sig)
        P = trenchInverse(Sig(:, 1));
        if ~isempty(P) && guardOk(Sig, P), return, end
    end
    [R, pd] = chol(Sig);
    if pd > 0
        error('lsqvce:lawInverse:notPD', 'Sig is not positive definite (chol failed at pivot %d).', pd);
    end
    X = R \ eye(N);
    P = X * X.';
    P = (P + P.') / 2;
end

function tf = isToeplitz(Sig)
%ISTOEPLITZ  Every diagonal constant to 1e-12 of the variance.
    tf = max(max(abs(Sig(2:end, 2:end) - Sig(1:end-1, 1:end-1)))) <= 1e-12 * abs(Sig(1, 1));
end

function P = trenchInverse(c)
%TRENCHINVERSE  inv(toeplitz(c)) by Levinson-Durbin + the Trench recursion;
%   [] when the recursion finds the law not positive definite.
    N = numel(c);
    lw = warning('off', 'lsqvce:levinson:notPD');  cw = onCleanup(@() warning(lw));
    try
        [a, s2] = lsqvce.kernel.levinsonDurbin(c(:), N - 1);
    catch ME
        % the recursion finds the law not positive definite: the Cholesky path
        % decides (and refuses it, lsqvce:lawInverse:notPD); anything else is
        % a fault and propagates
        if ~strcmp(ME.identifier, 'lsqvce:levinson:notPD'), rethrow(ME); end
        P = [];  return
    end
    if numel(a) < N || ~(s2 > 0), P = [];  return, end   % truncated: not PD
    x = a(:) / s2;                                     % inv(Sig) e_1
    z = [0; x(N:-1:2)];
    P = zeros(N);
    P(:, 1) = x;  P(1, :) = x.';
    for j = 1:N-1
        P(2:N, j+1) = P(1:N-1, j) + (x(2:N) * x(j+1) - z(2:N) * z(j+1)) / x(1);
    end
    P = (P + P.') / 2;
end

function ok = guardOk(Sig, P)
%GUARDOK  The structured inverse checked on three columns.
    N = size(Sig, 1);
    j = unique([1, ceil(N / 2), N]);
    E = zeros(N, numel(j));  E(sub2ind(size(E), j, 1:numel(j))) = 1;
    ok = max(max(abs(Sig * P(:, j) - E))) <= 1e-8;
end
