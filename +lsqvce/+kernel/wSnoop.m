function [z, q, del] = wSnoop(P, E, X, excl, thr, scale)
%LSQVCE.KERNEL.WSNOOP  Baarda's w per epoch of one arc, with iterative data snooping.
%
%   [z, q, del] = LSQVCE.KERNEL.WSNOOP(P, E, X, excl, thr, scale)
%
%   The robust noise estimator's test statistic (WP-156, design F). For an
%   observation i of an adjustment with residual v, R = P - P X (X' P X)^-1
%   X' P (P the inverse law, X the parameters estimated with it),
%       w_i = (R v)_i / sqrt(R_ii)
%   is Baarda's w for CORRELATED observations, and deleting observation i
%   removes exactly w_i^2 from v'Pv. Rejecting on w therefore truncates the
%   very quantity a deletion removes, which is what makes the white
%   truncation factor exact on a coloured stream; the raw residual it
%   replaces correlates with w_i only by sqrt(c_i / Sigma_ii) (HANDOVER
%   12.91). The arc is d single-axis series that share the temporal law
%   (the separable model), so axis k has its own residual E(:, k) and its
%   own parameter design X{k}, and the epoch statistic is the sum of its
%   axes' w^2.
%
%   SNOOPING. The largest |z| beyond THR is deleted (all its axes), every w
%   of the arc recomputed, and so on until none is left (Baarda's iterative
%   data snooping). Deleting one epoch changes its neighbours' w on a
%   correlated stream, so a simultaneous rule can delete a pair that each
%   pass restores; snooping removes one epoch at a time, the largest first.
%   It starts from the excluded rows only -- the caller rebuilds the set
%   from scratch every pass and never inherits it.
%
%   INPUT
%     P      m x m   inverse temporal law of the arc (LSQVCE.KERNEL.LAWINVERSE)
%     E      m x d   residual per axis, descaled by the epoch's Cep
%     X      1 x d   cell, X{k} m x n_k: the parameters of axis k estimated
%                    with this law (the arc's local design, descaled), m x 0 for none
%     excl   m x d   logical, rows excluded from the data (not observations)
%     thr    1x1     |z| beyond which an epoch is snooped; Inf: none
%     scale  1x1     the stream's scale s (w units, the previous pass's)
%
%   OUTPUT
%     z    m x 1  the epoch statistic in units of s: w / s for one live axis,
%                 sign(sum w) sqrt(2) erfcinv(Q(d_e / 2, q / (2 s^2))) for d_e
%                 > 1 live axes (the chi^2 tail mapped to a normal deviate, as
%                 R7-08); +-Inf on a snooped epoch, NaN on an epoch excluded on
%                 every axis
%     q    m x 1  the epoch's sum of w^2 over its live axes (unscaled; NaN as z)
%     del  m x 1  logical, the snooped epochs
%
%   COST. One rank-one update of P per snooped epoch (O(m^2)) and per axis a
%   matrix-vector product; the parameter part is O(m n_k^2).
%
%   EXAMPLE
%     Sig = toeplitz(0.9 .^ (0:255));  P = lsqvce.kernel.lawInverse(Sig);
%     e = chol(Sig, 'lower') * randn(256, 1);  e(100) = e(100) + 3;   % a blunder
%     [z, ~, del] = lsqvce.kernel.wSnoop(P, e, {ones(256, 1)}, false(256, 1), 3, 1);
%     find(del)                                % 100
%
%   ACCURACY
%     Deleting an epoch removes w_i^2 from v'Pv to 1e-10 relative, local
%     parameters re-estimated, and a two-axis epoch the sum of its axes' w^2
%     (test_precision (testDeletionRemovesExactlyWSquared)); on 40 AR(1)
%     arcs, phi 0.98, every planted 2-sigma blunder is snooped and the clean
%     epochs are rejected at the Gaussian 0.27 % within four Poisson sd
%     (testSnoopingFindsSmallBlundersOnARedStream); the statistic and the set
%     are invariant to the level (testRobustKernelIsScaleEquivariant).
%     The controlled study (dev/validation/wp156_selection_measure.m):
%     the level with tau within 0.15 % for phi <= 0.98.
%
%   ERRORS
%     lsqvce:wSnoop:badSize   E, X or excl do not match P
%
%   See also LSQVCE.KERNEL.LAWINVERSE, LSQVCE.KERNEL.TRUNCATIONDEFICIT,
%            LSQVCE.KERNEL.ROBUSTWEIGHTS, LSQVCE.ENGINE.OUTLIERS
%
%   Created 2026-10-18 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- WP-156 round 2, design F.

    [m, d] = size(E);
    if ~isequal(size(P), [m m]) || ~isequal(size(excl), [m d]) || numel(X) ~= d
        error('lsqvce:wSnoop:badSize', ...
              'P must be %d x %d, excl %d x %d and X hold %d designs (got P %d x %d, excl %d x %d, %d designs).', ...
              m, m, m, d, d, size(P, 1), size(P, 2), size(excl, 1), size(excl, 2), numel(X));
    end
    excl = logical(excl);
    del = false(m, 1);
    gone = all(excl, 2);                              % excluded on every axis
    PD = P;
    for j = find(gone).'
        PD = deleteEpoch(PD, j);
    end
    W = axisW(PD, E, X, excl, del);
    [z, q] = epochStat(W, excl, del, scale);
    while true
        a = abs(z);  a(del | ~isfinite(a)) = 0;
        [mx, j] = max(a);
        if isempty(j) || ~(mx > thr), break, end
        del(j) = true;
        PD = deleteEpoch(PD, j);
        W = axisW(PD, E, X, excl, del);
        [z, q] = epochStat(W, excl, del, scale);
    end
    s = sign(sum(E(del, :) .* ~excl(del, :), 2));  s(s == 0) = 1;
    z(del) = s * Inf;
end

function PD = deleteEpoch(PD, j)
%DELETEEPOCH  Epoch j out of the inverse law, exactly: the inverse of the law
%   restricted to the other epochs, row and column j zero.
    pj = PD(:, j);
    if pj(j) > 0
        PD = PD - (pj * pj.') / pj(j);
    end
    PD(:, j) = 0;  PD(j, :) = 0;
end

function W = axisW(PD, E, X, excl, del)
%AXISW  w of every live row, per axis: R = PD - PD Xk (Xk' PD Xk)^-1 Xk' PD,
%   Xk = [unit columns of the rows excluded on this axis only, X{k}].
    [m, d] = size(E);
    W = nan(m, d);
    gone = all(excl, 2) | del;
    dP = diag(PD);
    for k = 1:d
        part = find(excl(:, k) & ~gone);
        Xk = [full(sparse(part, 1:numel(part), 1, m, numel(part))), X{k}];
        Xk(gone, :) = 0;
        Xk = Xk(:, any(Xk ~= 0, 1));
        Pe = PD * E(:, k);
        dR = dP;
        if ~isempty(Xk)
            G = PD * Xk;
            % N = Xk' PD Xk by its eigenvectors: a column the live rows do not
            % determine (an arc parameter seen only on deleted rows) drops out
            [U, S] = eig((Xk.' * G + G.' * Xk) / 2);
            s = diag(S);
            keep = s > 1e3 * eps * max([s; 0]);
            T = G * (U(:, keep) ./ sqrt(s(keep)).');  % T T' = PD Xk N^+ Xk' PD
            Pe = Pe - T * (T.' * E(:, k));
            dR = dR - sum(T.^2, 2);
        end
        live = ~gone & ~excl(:, k) & dR > 0;
        W(live, k) = Pe(live) ./ sqrt(dR(live));
    end
end

function [z, q] = epochStat(W, excl, del, scale)
%EPOCHSTAT  One statistic per epoch in units of the stream's scale.
    live = ~excl & ~isnan(W);
    live(del, :) = false;
    de = sum(live, 2);
    W0 = W;  W0(~live) = 0;
    q = sum(W0.^2, 2);
    z = nan(size(q));
    one = de == 1;
    z(one) = sum(W0(one, :), 2) / scale;
    multi = de > 1;
    if any(multi)
        tail = gammainc(q(multi) / (2 * scale^2), de(multi) / 2, 'upper');
        sg = sign(sum(W0(multi, :), 2));  sg(sg == 0) = 1;
        z(multi) = sg .* min(sqrt(2) * erfcinv(tail), 40);
    end
    q(de == 0) = NaN;
end
