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
%     Q = LSQVCE.ENGINE.BLOCKCOFACTOR(Rc, blocks, datum)
%       datum   [] | struct .Rm .Dm (BlockAccumulator.datumPart): the common
%                        prior holds a DATUM group (the data and the other
%                        groups do not see its directions), and the PARAMETER
%                        cofactor -- diagonal, block, qxxTimes, maxCorrelation
%                        -- is the minimum-constraint one of the data plus the
%                        other groups, Qcc = inv(Rm'Rm) - Dm*Dm' (R7 WP-114,
%                        R7-05; per group WP-119), independent of the datum's
%                        sigma. leverage / halfRows
%                        keep Rc: on data rows the two agree exactly (a data
%                        row is orthogonal to the datum), and a datum row must
%                        stay absorbed (r = 0) as the solve made it.
%       Rc      nc x nc  upper triangular, Rc'*Rc = the reduced common normal
%                        matrix INCLUDING the prior (BlockAccumulator.commonFactor)
%       blocks  1 x S struct  .cols (the system's local columns in the full
%                        vector), .Nll (the REGULARISED local block the solve
%                        used), .Nlc, optionally .truncTol (the relative
%                        eigenvalue cut of the solve's pinv fallback, the
%                        accumulator's TruncTol), optionally .Rll / .Rlc (the
%                        local factor and coupling of a block eliminated on the
%                        QR factor; used instead of .Nll, WP-128); [] for a
%                        local-free problem
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
%                                  (a block kept on its QR factor: on that factor)
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
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-130: a block reduced by the
%   solve's pseudo-inverse (blocks .Vll) gives H, G on its kept basis and
%   downdates there; an eigen block downdates with its own cut (sys .tol). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-129: downdate of a block kept
%   on its QR factor downdates that factor (cholupdate of Rll), no longer
%   the squared Nll. With assistance from Claude (Anthropic), model Claude
%   Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-128: a local block the solve
%   eliminated on the QR factor (blocks .Rll, .Rlc) gives H and G from that
%   factor, not from the eigen pseudo-inverse of the squared Nll. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-125: leverage() returns the common
%   term hc as its second output. With assistance from Claude (Anthropic), model
%   Claude Opus 5.5.
%   Updated 2026-10-12 (CEST, UTC+02:00), R7 WP-119: the datum part may sit
%   beside regularising prior groups (BlockAccumulator.datumPart per group).
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-114 (R7-05): the optional
%   datum part -- the parameter cofactor of a datum-constrained problem is the
%   minimum-constraint one; downdate updates its factor too. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
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
        sys         % 1 x S struct: cols, Nll, H (H'*H = pinv(Nll)), G = pinv(Nll)*Nlc,
                    % Rll ([] or the local QR factor, positive pivots: downdate's);
                    % a local-free system keeps its slot (empty cols), so s indexes
                    % the systems in the order given
        Rm = []     % [] | nc x nc, Rm'*Rm = N_red + Ki'Ki + c^2 Kd'Kd (the datum re-weighted, R7 WP-114/119)
        Dm = []     % nc x k, Dm*Dm' = K'(KK')^-2 K / c^2: Qcc = inv(Rm'Rm) - Dm*Dm'
    end

    methods
        %------------------------------------------------------------------
        function Q = BlockCofactor(Rc, blocks, datum)
            nc = size(Rc, 1);
            Q.Rc = Rc;
            Q.nc = nc;
            Q.Dm = zeros(nc, 0);
            if nargin > 2 && ~isempty(datum)
                Q.Rm = datum.Rm;  Q.Dm = datum.Dm;
            end
            Q.sys = struct('cols', {}, 'Nll', {}, 'H', {}, 'G', {}, 'Rll', {}, 'Vll', {}, 'tol', {});
            u = nc;
            for s = 1:numel(blocks)
                b = blocks(s);
                nl = numel(b.cols);
                if nl == 0                               % local-free: keeps its index
                    Q.sys(end+1) = struct('cols', zeros(1, 0), 'Nll', zeros(0), ...
                                          'H', zeros(0), 'G', zeros(0, nc), 'Rll', [], ...
                                          'Vll', [], 'tol', []);
                    continue
                end
                if ~isequal(size(b.Nll), [nl nl]) || ~isequal(size(b.Nlc), [nl nc])
                    error('lsqvce:BlockCofactor:badBlocks', ...
                        ['block %d: Nll is %s and Nlc %s for %d local columns and ' ...
                         '%d commons.'], s, mat2str(size(b.Nll)), mat2str(size(b.Nlc)), nl, nc);
                end
                tol = [];  Rl = [];  Vl = [];
                if isfield(b, 'truncTol'), tol = b.truncTol; end
                if isfield(b, 'Rll') && ~isempty(b.Rll)
                    % WP-128: the block the solve eliminated EXACTLY on its QR
                    % factor: H = R11^-T, G = R11^-1 R12 (H'H = inv(Nll)). The
                    % squared Nll took the eigen branch from cond(B) ~ 7e7 and cut
                    % at 1e-10 max -- one redundancy unit lost per cut direction
                    % (a t^2 term over 1e4 s; dev/validation/wp128_local_block_mirror.py)
                    % Kept with POSITIVE pivots (R'R unchanged): downdate's
                    % cholupdate needs a Cholesky factor (WP-129).
                    sg = sign(diag(b.Rll));  sg(sg == 0) = 1;
                    Rl = sg .* b.Rll;
                    H = Rl.' \ eye(size(Rl, 1));  G = Rl \ (sg .* b.Rlc);
                    if isfield(b, 'Vll') && ~isempty(b.Vll)
                        % WP-130: a rank-deficient block reduced by the solve's
                        % pseudo-inverse on the kept basis Vll: H'H = pinv(Nll)
                        Vl = b.Vll;  H = [H * Vl.'; zeros(nl - size(Vl, 2), nl)];  G = Vl * G;
                    end
                else
                    [H, G] = lsqvce.engine.BlockCofactor.localBlock(b.Nll, b.Nlc, tol);
                end
                Q.sys(end+1) = struct('cols', b.cols(:).', 'Nll', b.Nll, 'H', H, 'G', G, 'Rll', Rl, ...
                                      'Vll', Vl, 'tol', tol);
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
        function [h, hc] = leverage(Q, s, Ac, B)
        %LEVERAGE  h_i = a_i'*Qxx*a_i for rows [Ac B] of system s (s = 0 or an
        %   empty B: local-free rows), without widening them to u columns:
        %   h = ||Rc' \ (Ac' - G_s'*B')||^2 + ||H_s*B'||^2 per column.
        %   hc (optional) is the first, COMMON term alone: z_i'*inv(N_red)*z_i
        %   with z_i the row reduced by the elimination, so its sum over a group
        %   is tr(inv(N_red)*N_g^red), the group's share of the common solution
        %   (WP-125: LSQVCE.ENGINE.CONTRIBUTION's Level-1 share).
            Z = Ac.';
            if s > 0 && ~isempty(B)
                Bt = B.';
                Z = Z - Q.sys(s).G.' * Bt;
                hc = sum((Q.Rc.' \ Z).^2, 1).';
                h = hc + sum((Q.sys(s).H * Bt).^2, 1).';
            else
                hc = sum((Q.Rc.' \ Z).^2, 1).';
                h = hc;
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
            F = parFactor(Q);
            Yc = F \ (F.' \ Z) - Q.Dm * (Q.Dm.' * Z);    % Qcc * z (minimum-constraint under a datum)
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
            Ri = parFactor(Q) \ eye(Q.nc);              % Qcc = Ri*Ri' - Dm*Dm'
            d = zeros(Q.u, 1);
            d(1:Q.nc) = sum(Ri.^2, 2) - sum(Q.Dm.^2, 2);
            for s = 1:numel(Q.sys)
                d(Q.sys(s).cols) = sum(Q.sys(s).H.^2, 1).' + sum((Q.sys(s).G * Ri).^2, 2) ...
                                   - sum((Q.sys(s).G * Q.Dm).^2, 2);
            end
        end

        %------------------------------------------------------------------
        function C = block(Q, idx)
        %BLOCK  Qxx(idx, idx) for a short index list (an ellipse, a pair).
            E = zeros(numel(idx), Q.u);
            E(sub2ind(size(E), 1:numel(idx), idx(:).')) = 1;
            Y = qxxTimes(Q, E);                          % the parameter cofactor (datum-aware)
            C = Y(idx(:), :);
            C = (C + C.') / 2;
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
        %   A block kept on its QR factor Rll downdates THAT factor (cholupdate
        %   of Rll with sqrt(p) b); refactoring the squared Nll lost its weak
        %   direction from cond(B) ~ 7e7 -- after one rejection the leverage
        %   sum missed u by a whole unit (WP-129;
        %   dev/validation/wp129_local_downdate_mirror.py).
        %   ok = false (Q unchanged) when the removal would make the system
        %   singular: d ~ 0 (the row is the only information on a local), the
        %   downdated local factor fails, or the downdated Rc collapses.
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
                if ~isempty(S.Rll)
                    bb = b;  if ~isempty(S.Vll), bb = S.Vll.' * b; end   % WP-130: on the kept basis
                    [Rl2, fl] = cholupdate(S.Rll, sqrt(p) * bb, '-');
                    if fl ~= 0, ok = false;  return, end
                end
            end
            [Rn, flag] = cholupdate(Q.Rc, sqrt(pt) * Z, '-');
            dg = abs(diag(Rn));
            ok = flag == 0 && ~any(dg < eps(max(dg)) * Q.u * 100);
            if ~ok, return, end
            if ~isempty(Q.Rm)                            % the datum factor sees the same data
                [Rm2, f2] = cholupdate(Q.Rm, sqrt(pt) * Z, '-');
                if f2 ~= 0, ok = false;  return, end
                Q.Rm = Rm2;
            end
            Q.Rc = Rn;
            if ~isempty(hit)
                S.Nll = S.Nll - p * (b * b.');
                nlb = numel(b);
                if isempty(S.Rll)                        % the block's own cut (WP-130)
                    S.H = lsqvce.engine.BlockCofactor.localBlock(S.Nll, zeros(nlb, 0), S.tol);
                elseif isempty(S.Vll)
                    S.Rll = Rl2;  S.H = Rl2.' \ eye(nlb);         % WP-129: on the factor
                else
                    S.Rll = Rl2;  k = size(Rl2, 1);
                    S.H = [(Rl2.' \ eye(k)) * S.Vll.'; zeros(nlb - k, nlb)];
                end
                S.G = S.G - pt * (c * Z.');
                Q.sys(hit) = S;
            end
        end
    end

    methods (Access = private)
        %------------------------------------------------------------------
        function F = parFactor(Q)
        %PARFACTOR  The commons factor of the PARAMETER cofactor: Rm under a
        %   datum (with Dm subtracted by the caller), else Rc.
            F = Q.Rc;
            if ~isempty(Q.Rm), F = Q.Rm; end
        end

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
            blocks = struct('cols', cell(1, numel(ids)), 'Nll', [], 'Nlc', [], 'truncTol', [], ...
                            'Rll', [], 'Rlc', [], 'Vll', []);
            col = nc;
            for k = 1:numel(ids)
                st = acc.getStore(ids(k));
                nl = size(st.Nll_r, 1);
                blocks(k) = struct('cols', col + (1:nl), 'Nll', st.Nll_r, 'Nlc', st.Nlc, ...
                                   'truncTol', acc.opt.TruncTol, ...   % the solve's pinv cut
                                   'Rll', lsqvce.util.getfielddef(st, 'Rll', []), ...
                                   'Rlc', lsqvce.util.getfielddef(st, 'Rlc', []), ...
                                   'Vll', lsqvce.util.getfielddef(st, 'Vll', []));
                col = col + nl;
            end
            Q = lsqvce.engine.BlockCofactor(acc.commonFactor(), blocks, acc.datumPart());
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
