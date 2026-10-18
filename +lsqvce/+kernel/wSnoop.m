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
%   COST. One block update for the excluded epochs and one product P E per
%   call (O(m^2)); each snooping step projects the snooped epochs out with the
%   parameters, O(m r^2) for r projected columns (WP-155: a rank-one update
%   of P per deletion was O(m^2)).
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
%   Updated 2026-10-18 (CEST, UTC+02:00), WP-155: the snooped epochs are
%   projected out with the parameters (a deletion is a free parameter for
%   the observation) instead of updating the inverse law per deletion; the
%   same w to 1e-10 (test_precision testSnoopingByProjectionIsTheReducedSystem).
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
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
    PD = deleteRows(P, gone);                         % once, BLAS-3
    PE = PD * E;                                      % once per call: every axis' P e
    PX = cellfun(@(Xk) PD * Xk, X, 'UniformOutput', false);
    W = axisW(PD, PE, PX, E, X, excl, del);
    [z, q] = epochStat(W, excl, del, scale);
    while true
        a = abs(z);  a(del | ~isfinite(a)) = 0;
        [mx, j] = max(a);
        if isempty(j) || ~(mx > thr), break, end
        del(j) = true;
        W = axisW(PD, PE, PX, E, X, excl, del);
        [z, q] = epochStat(W, excl, del, scale);
    end
    s = sign(sum(E(del, :) .* ~excl(del, :), 2));  s(s == 0) = 1;
    z(del) = s * Inf;
end

function PD = deleteRows(P, gone)
%DELETEROWS  The epochs GONE out of the inverse law, exactly (the inverse of
%   the law restricted to the others: P - P_:g P_gg^-1 P_g:), their rows and
%   columns zero -- one block update.
    PD = P;
    g = find(gone);
    if isempty(g), return, end
    Pg = P(:, g);
    PD = P - Pg * (P(g, g) \ Pg.');
    PD(:, g) = 0;  PD(g, :) = 0;
end

function W = axisW(PD, PE, PX, E, X, excl, del)
%AXISW  w of every live row, per axis: R = PD - PD Z (Z' PD Z)^+ Z' PD, Z =
%   [unit columns of the snooped epochs and of the rows excluded on this axis
%   only, X{k}]. WP-155: deleting an observation is projecting out its unit
%   vector (it gains a free parameter), so the snooped epochs join the
%   projection instead of updating PD -- O(m r^2) per step for r projected
%   columns, where the rank-one update of PD cost O(m^2) per deletion. The
%   same R to rounding.
    [m, d] = size(E);
    W = nan(m, d);
    gone = all(excl, 2) | del;
    dP = diag(PD);
    sn = find(del);
    for k = 1:d
        part = find(excl(:, k) & ~gone);
        u = [sn; part];                               % unit columns
        G = [PD(:, u), PX{k}];
        N = [G(u, :); X{k}.' * G];                    % Z' PD Z
        Pe = PE(:, k);
        dR = dP;
        if ~isempty(G)
            % by its eigenvectors: a column the live rows do not determine (an
            % arc parameter seen only on deleted rows) drops out
            [U, S] = eig((N + N.') / 2);
            sv = diag(S);
            keep = sv > 1e3 * eps * max([sv; 0]);
            T = G * (U(:, keep) ./ sqrt(sv(keep)).'); % T T' = PD Z N^+ Z' PD
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
