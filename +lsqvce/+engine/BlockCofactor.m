classdef BlockCofactor
%LSQVCE.ENGINE.BLOCKCOFACTOR  The cofactor Qxx = inv(N) of a block-arrow system, kept as blocks.
%
%   The diagnostics need Qxx of the FULL parameter set [commons | locals of
%   every system] -- a u x u matrix with u = nc + sum(nl), 13 729 at a
%   GRAVFIRE month (1.5 GB per dense copy). They never need it dense. With
%   the reduced commons factor Rc (Rc'Rc = N_red, the prior folded in) and,
%   per system, the local normal block Nll and G = inv(Nll)*Nlc,
%
%       Qcc = inv(N_red)                  Qlc,s = -G_s Qcc
%       Qll,s = inv(Nll_s) + G_s Qcc G_s'    Qll,st = G_s Qcc G_t'  (s ~= t)
%
%   so for a full-width row a = [a_c ; b] (b on the locals),
%
%       a' Qxx a = || Rc' \ (a_c - sum_s G_s' b_s) ||^2 + sum_s || H_s b_s ||^2
%
%   with H_s'H_s = inv(Nll_s). Every question diagnose, snoop, ellipse and
%   report ask of Qxx is answered from these blocks: memory O(nc^2 + sum
%   nl*nc), cost O(nc^2 + nl*nc) per row instead of O(u^2) (WP-86 part 2).
%
%   CONSTRUCTOR
%     Q = LSQVCE.ENGINE.BLOCKCOFACTOR(Rc, blocks)
%       Rc      nc x nc  upper triangular, Rc'*Rc = the reduced common normal
%                        matrix INCLUDING the prior (BlockAccumulator.commonFactor)
%       blocks  1 x S struct  .cols (the system's local columns in the full
%                        vector), .Nll (the REGULARISED local block the solve
%                        used), .Nlc, optionally .truncTol (the relative
%                        eigenvalue cut of the solve's pinv fallback, the
%                        accumulator's TruncTol); [] for a local-free problem
%
%     Q = LSQVCE.ENGINE.BLOCKCOFACTOR.fromAccumulator(acc)
%       the cofactor of a solved lsqvce.engine.BlockAccumulator ('KeepLocal'),
%       its systems in the accumulator's order, local columns after the commons
%
%   METHODS
%     W = halfRows(Q, A)           u x nk, W'*W = A*Qxx*A' for full-width rows
%                                  A (nk x u): the dense code's R' \ A'
%     h = leverage(Q, s, Ac, B)    the same diagonal for rows of system s given
%                                  as their commons part Ac (nk x nc) and local
%                                  part B (nk x nl_s); s = 0 for local-free rows
%     s = systemOf(Q, A)           the system whose locals full-width rows A
%                                  touch (0: none)
%     Y = qxxTimes(Q, A)           u x nk, Qxx*A'
%     d = diagonal(Q)              u x 1, diag(Qxx)
%     C = block(Q, idx)            Qxx(idx, idx)
%     [c, ij] = maxCorrelation(Q)  the largest |corr| off the diagonal, its pair
%     k = condLowerBound(Q)        a LOWER bound on cond(N) (Rayleigh quotients)
%     [Q, ok] = downdate(Q, a, p)  remove one row a (weight p): snoop's step
%
%   A NOTE ON WHAT DOES NOT WORK (moved here from the retired leverage unit,
%   WP-86 part 3). A stochastic (Hutchinson) trace estimator looks attractive
%   because it costs O(nProbe) solves independent of nc. It was tried and
%   REJECTED for the per-observation redundancies: on a realistically
%   ill-conditioned gravity normal matrix the relative error is 20-50% at 32
%   probes, because ||M||_F / trace(M) ~ 1 -- the operator's energy sits in
%   essentially one mode, the worst case for Hutchinson. Neither a sum
%   constraint nor Hutch++ rescued it. The exact triangular-solve route here
%   is both faster (measured 2x against forming inv(N) at nc = 2500) and exact.
%
%   A LOCAL BLOCK THAT IS SINGULAR (the accumulator's pinv fallback) gets the
%   pseudo-inverse half H = S^-1/2 V' over its range -- the same pseudo-inverse
%   the pre-elimination used, so Qxx is the cofactor of the solve actually made.
%
%   EXAMPLE
%       R = lsqvce.fit(data);                          % any fit
%       Q = R.source.Q;                                % built by diagnoseAccum
%       A = R.batchFcn(1);                             % unit-weight rows
%       h = sum(Q.halfRows(A).^2, 1).';                % leverage of batch 1
%       C = Q.block([1 2]);                            % 2 x 2 cofactor of x(1:2)
%
%   ERRORS
%     lsqvce:BlockCofactor:badBlocks       a block's sizes disagree with Rc
%     lsqvce:BlockCofactor:multiSystemRow  downdate (or systemOf) of rows that
%                                          touch the locals of more than one
%                                          system (the block-arrow structure
%                                          would not survive)
%
%   See also LSQVCE.ENGINE.VERDICT.SOURCE, LSQVCE.ENGINE.VERDICT.DIAGNOSE, LSQVCE.ENGINE.VERDICT.SNOOP,
%   LSQVCE.ENGINE.BLOCKACCUMULATOR
%
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-86 part 4: the absorbed
%   diagnostics units are named as lsqvce.engine.verdict methods. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-86 part 3 (R7-24): leverage
%   (per system, without widening the rows), systemOf, fromAccumulator;
%   local-free systems keep their index; the Hutchinson note moved here from
%   the retired leverage. With assistance from Claude (Anthropic), model
%   Claude Opus 5.5.
%   Created 2026-10-10 (CEST, UTC+02:00), R7 WP-86 part 2 (R7-24): the
%   block-native cofactor replacing the dense full-width factor of the
%   diagnostics. With assistance from Claude (Anthropic), model Claude Opus 5.5.

    properties (SetAccess = private)
        Rc          % nc x nc upper triangular, Rc'*Rc = N_red (prior included)
        nc          % number of common parameters
        u           % full parameter count nc + sum(nl)
        sys         % 1 x S struct: cols, Nll, H (H'*H = pinv(Nll)), G = pinv(Nll)*Nlc;
                    % a local-free system keeps its slot (empty cols), so s indexes
                    % the systems in the order given
    end

    methods
        %------------------------------------------------------------------
        function Q = BlockCofactor(Rc, blocks)
            nc = size(Rc, 1);
            Q.Rc = Rc;
            Q.nc = nc;
            Q.sys = struct('cols', {}, 'Nll', {}, 'H', {}, 'G', {});
            u = nc;
            for s = 1:numel(blocks)
                b = blocks(s);
                nl = numel(b.cols);
                if nl == 0                               % local-free: keeps its index
                    Q.sys(end+1) = struct('cols', zeros(1, 0), 'Nll', zeros(0), ...
                                          'H', zeros(0), 'G', zeros(0, nc));
                    continue
                end
                if ~isequal(size(b.Nll), [nl nl]) || ~isequal(size(b.Nlc), [nl nc])
                    error('lsqvce:BlockCofactor:badBlocks', ...
                        ['block %d: Nll is %s and Nlc %s for %d local columns and ' ...
                         '%d commons.'], s, mat2str(size(b.Nll)), mat2str(size(b.Nlc)), nl, nc);
                end
                tol = [];
                if isfield(b, 'truncTol'), tol = b.truncTol; end
                [H, G] = lsqvce.engine.BlockCofactor.localBlock(b.Nll, b.Nlc, tol);
                Q.sys(end+1) = struct('cols', b.cols(:).', 'Nll', b.Nll, 'H', H, 'G', G);
                u = u + nl;
            end
            Q.u = u;
        end

        %------------------------------------------------------------------
        function W = halfRows(Q, A)
        %HALFROWS  W (u x nk) with W'*W = A*Qxx*A': one column per row of A.
            [Z, Bs] = splitRows(Q, A);
            W = zeros(Q.u, size(A, 1));
            W(1:Q.nc, :) = Q.Rc.' \ Z;
            for s = find(~cellfun(@isempty, Bs))
                W(Q.sys(s).cols, :) = Q.sys(s).H * Bs{s};
            end
        end

        %------------------------------------------------------------------
        function h = leverage(Q, s, Ac, B)
        %LEVERAGE  h_i = a_i'*Qxx*a_i for rows [Ac B] of system s (s = 0 or an
        %   empty B: local-free rows), without widening them to u columns:
        %   h = ||Rc' \ (Ac' - G_s'*B')||^2 + ||H_s*B'||^2 per column.
            Z = Ac.';
            if s > 0 && ~isempty(B)
                Bt = B.';
                Z = Z - Q.sys(s).G.' * Bt;
                h = sum((Q.Rc.' \ Z).^2, 1).' + sum((Q.sys(s).H * Bt).^2, 1).';
            else
                h = sum((Q.Rc.' \ Z).^2, 1).';
            end
        end

        %------------------------------------------------------------------
        function s = systemOf(Q, A)
        %SYSTEMOF  The system whose local columns the rows A (full width)
        %   touch: 0 when none (commons or prior rows), an error when several.
            hit = find(arrayfun(@(t) any(any(A(:, t.cols))), Q.sys));
            if numel(hit) > 1
                error('lsqvce:BlockCofactor:multiSystemRow', ...
                      'the rows touch the locals of %d systems; one batch is one system.', numel(hit));
            end
            s = 0;
            if ~isempty(hit), s = hit; end
        end

        %------------------------------------------------------------------
        function Y = qxxTimes(Q, A)
        %QXXTIMES  Qxx*A' (u x nk): the shift of x per unit of each row's bias.
            [Z, Bs] = splitRows(Q, A);
            Yc = Q.Rc \ (Q.Rc.' \ Z);                    % Qcc * z
            Y = zeros(Q.u, size(A, 1));
            Y(1:Q.nc, :) = Yc;
            for s = 1:numel(Q.sys)
                Ys = -Q.sys(s).G * Yc;
                if ~isempty(Bs{s})
                    Ys = Ys + Q.sys(s).H.' * (Q.sys(s).H * Bs{s});
                end
                Y(Q.sys(s).cols, :) = Ys;
            end
        end

        %------------------------------------------------------------------
        function d = diagonal(Q)
        %DIAGONAL  diag(Qxx) (u x 1). One O(nc^3) triangular inverse, then
        %   O(nl*nc^2) per system; nothing u x u is formed.
            Ri = Q.Rc \ eye(Q.nc);                       % Qcc = Ri*Ri'
            d = zeros(Q.u, 1);
            d(1:Q.nc) = sum(Ri.^2, 2);
            for s = 1:numel(Q.sys)
                d(Q.sys(s).cols) = sum(Q.sys(s).H.^2, 1).' + sum((Q.sys(s).G * Ri).^2, 2);
            end
        end

        %------------------------------------------------------------------
        function C = block(Q, idx)
        %BLOCK  Qxx(idx, idx) for a short index list (an ellipse, a pair).
            E = zeros(numel(idx), Q.u);
            E(sub2ind(size(E), 1:numel(idx), idx(:).')) = 1;
            W = halfRows(Q, E);
            C = W.' * W;
        end

        %------------------------------------------------------------------
        function [cmax, pair] = maxCorrelation(Q)
        %MAXCORRELATION  The largest |correlation| between two parameters and
        %   the pair, from column panels Qxx(:, J) of ~64 MB: the u x u
        %   correlation matrix is never held.
            sig = sqrt(max(diagonal(Q), 0));
            pw = max(1, floor(8e6 / Q.u));
            cmax = 0;  pair = [1 1];
            for j0 = 1:pw:Q.u
                J = j0:min(j0 + pw - 1, Q.u);
                E = zeros(numel(J), Q.u);
                E(sub2ind(size(E), 1:numel(J), J)) = 1;
                C = qxxTimes(Q, E) ./ (sig * sig(J).');     % corr(:, J)
                C(sub2ind(size(C), J, 1:numel(J))) = 0;
                C(~isfinite(C)) = 0;                     % a null-space parameter has sigma 0
                [c, k] = max(abs(C(:)));
                if c > cmax
                    [i, jj] = ind2sub(size(C), k);
                    cmax = c;  pair = [i J(jj)];
                end
            end
        end

        %------------------------------------------------------------------
        function k = condLowerBound(Q)
        %CONDLOWERBOUND  A LOWER bound on cond(N) = lambdaMax(N)*lambdaMax(Qxx).
        %   Every Rayleigh quotient and every diagonal entry is <= the largest
        %   eigenvalue, so the product of the best ones found (the diagonals,
        %   then 30 power steps on N and on Qxx) can under-state, never
        %   over-state, the conditioning -- which is how report prints it.
            dN = zeros(Q.u, 1);
            dN(1:Q.nc) = sum(Q.Rc.^2, 1).';
            for s = 1:numel(Q.sys)
                dN(1:Q.nc) = dN(1:Q.nc) + sum(Q.sys(s).G .* (Q.sys(s).Nll * Q.sys(s).G), 1).';
                dN(Q.sys(s).cols) = diag(Q.sys(s).Nll);
            end
            lam = max(dN);
            mu = max(diagonal(Q));
            v = ones(Q.u, 1) / sqrt(Q.u);
            w = v;
            for it = 1:30
                y = normalTimes(Q, v);   lam = max(lam, v.' * y);   v = y / norm(y);
                z = qxxTimes(Q, w.'); mu  = max(mu,  w.' * z);   w = z / norm(z);
            end
            k = lam * mu;
        end

        %------------------------------------------------------------------
        function [Q, ok] = downdate(Q, a, p)
        %DOWNDATE  Remove one row a (1 x u) of weight p from the system.
        %   The Schur complement of a rank-one downdate is a rank-one downdate:
        %   with z = a_c - G_s'*b, c = inv(Nll_s)*b and d = 1 - p*b'*c,
        %       N_red <- N_red - (p/d) z z'      (cholupdate of Rc)
        %       Nll_s <- Nll_s - p b b',   G_s <- G_s - (p/d) c z'
        %   ok = false (Q unchanged) when the removal would make the system
        %   singular: d ~ 0 (the row is the only information on a local) or
        %   the downdated Rc collapses.
            [Z, Bs] = splitRows(Q, a);
            hit = find(~cellfun(@isempty, Bs));
            if numel(hit) > 1
                error('lsqvce:BlockCofactor:multiSystemRow', ...
                      'the row touches the locals of %d systems; one row is one system.', numel(hit));
            end
            tol = 100 * Q.u * eps;
            pt = p;
            if ~isempty(hit)
                S = Q.sys(hit);  b = Bs{hit};
                c = S.H.' * (S.H * b);
                d = 1 - p * (b.' * c);
                if d <= tol
                    ok = false;  return
                end
                pt = p / d;
            end
            [Rn, flag] = cholupdate(Q.Rc, sqrt(pt) * Z, '-');
            dg = abs(diag(Rn));
            ok = flag == 0 && ~any(dg < eps(max(dg)) * Q.u * 100);
            if ~ok, return, end
            Q.Rc = Rn;
            if ~isempty(hit)
                S.Nll = S.Nll - p * (b * b.');
                S.H = lsqvce.engine.BlockCofactor.localBlock(S.Nll, zeros(numel(b), 0), []);
                S.G = S.G - pt * (c * Z.');
                Q.sys(hit) = S;
            end
        end
    end

    methods (Access = private)
        %------------------------------------------------------------------
        function [Z, Bs] = splitRows(Q, A)
        %SPLITROWS  The reduced rows z = a_c - sum G_s'*b_s (nc x nk) and each
        %   system's local part b_s (nl x nk; {} for a system the rows miss).
            Z = A(:, 1:Q.nc).';
            Bs = cell(1, numel(Q.sys));
            for s = 1:numel(Q.sys)
                B = A(:, Q.sys(s).cols).';
                if any(B(:))
                    Z = Z - Q.sys(s).G.' * B;
                    Bs{s} = B;
                end
            end
        end

        %------------------------------------------------------------------
        function Y = normalTimes(Q, V)
        %NORMALTIMES  N*V from the blocks (the power steps of condLowerBound).
            Vc = V(1:Q.nc, :);
            Y = zeros(size(V));
            Y(1:Q.nc, :) = Q.Rc.' * (Q.Rc * Vc);
            for s = 1:numel(Q.sys)
                T = Q.sys(s).Nll * (Q.sys(s).G * Vc + V(Q.sys(s).cols, :));
                Y(1:Q.nc, :) = Y(1:Q.nc, :) + Q.sys(s).G.' * T;
                Y(Q.sys(s).cols, :) = T;
            end
        end
    end

    methods (Static)
        %------------------------------------------------------------------
        function Q = fromAccumulator(acc)
        %FROMACCUMULATOR  The cofactor of a solved accumulator ('KeepLocal'):
        %   the commons factor (acc.commonFactor) and, per system in the
        %   accumulator's order, the REGULARISED local block the solve used,
        %   with the solve's pinv cut. Local columns follow the commons in
        %   system order -- the layout of lsqvce.engine.verdict.source.
            ids = acc.systems();
            nc = acc.nc;
            blocks = struct('cols', cell(1, numel(ids)), 'Nll', [], 'Nlc', [], 'truncTol', []);
            col = nc;
            for k = 1:numel(ids)
                st = acc.getStore(ids(k));
                nl = size(st.Nll_r, 1);
                blocks(k) = struct('cols', col + (1:nl), 'Nll', st.Nll_r, 'Nlc', st.Nlc, ...
                                   'truncTol', acc.opt.TruncTol);   % the solve's pinv cut
                col = col + nl;
            end
            Q = lsqvce.engine.BlockCofactor(acc.commonFactor(), blocks);
        end
    end

    methods (Static, Access = private)
        %------------------------------------------------------------------
        function [H, G] = localBlock(Nll, Nlc, tol)
        %LOCALBLOCK  H with H'*H = pinv(Nll), and G = pinv(Nll)*Nlc. Cholesky
        %   when the block is regular; the eigen pseudo-inverse on its range
        %   when rcond < eps -- the accumulator's own switch to its pinv
        %   fallback, keeping eigenvalues above tol*max ([] = pinv's default
        %   cut nl*eps(max)).
            S = (Nll + Nll.') / 2;
            nl = size(S, 1);
            if rcond(S) >= eps
                [L, f] = chol(S, 'lower');
                if f == 0
                    H = L \ eye(nl);
                    G = L.' \ (L \ Nlc);
                    return
                end
            end
            [V, dv] = eig(S, 'vector');
            dv = real(dv);  V = real(V);
            if isempty(tol)
                cut = nl * eps(max(dv));
            else
                cut = tol * max(dv);
            end
            keep = dv > cut & dv > 0;
            H = zeros(nl);
            H(1:nnz(keep), :) = (V(:, keep) ./ sqrt(dv(keep)).').';
            G = H.' * (H * Nlc);
        end
    end
end
