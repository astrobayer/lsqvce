classdef QRSystem < handle
%LSQVCE.ENGINE.QRSYSTEM  Square-root information accumulator (Householder QR).
%
%   Maintains the upper-triangular factor R of the stacked [A l], so that
%   A'*A = R1'*R1 without A'*A ever being formed. The condition number of R is
%   cond(A), not cond(A)^2 -- roughly twice the flops of LSQVCE.ENGINE.NORMALEQSYSTEM in
%   exchange for double the effective precision.
%
%   This is the default backend, and the right one whenever the local block
%   might be poorly conditioned (near-collinear empirical accelerations, short
%   arcs, weakly separated bias/scale pairs -- i.e. most real problems).
%
%   THE ACCUMULATOR PATH STAYS SQUARE-ROOT (R7 WP-112, review R7-03). A
%   system's local parameters are eliminated ON THE FACTOR: reducedRows
%   re-triangularises R with the local columns first, and the trailing block
%   [R22 r2; 0 rho] is the square root of the reduced contribution
%   (R22'R22 = Ncc - Nlc'Nll^-1 Nlc, rho^2 + r2'r2 = l'l reduced). merge
%   stacks those rows -- no chol, no eig. Until 6.7 the contribution went
%   through N = R1'R1, a normal-form Schur complement and chol in merge:
%   cond(A)^2 after all, 2-3 formal sd error at cond(A) = 1e7-1e9 and
%   full-rank problems refused from 1e8.5 (mirror, HANDOVER 12.45).
%
%   THE RANK DECISION (R7-16) is this class's own: a pivot of the final
%   factor with |R(k,k)| <= RankTol * max|diag R| is refused. RankTol =
%   1e3 eps: a free levelling network leaves pivot ratios <= 5e-16 (2000
%   random networks), a full-rank problem of cond(A) = 1e12 keeps >= 1.9e-12.
%   Before, the decision fell to LAPACK -- chol in merge succeeded or not on
%   rounding, and a free network was refused 61 % / solved with an arbitrary
%   datum 39 %.
%
%   PROPERTIES (read-only)
%     R     (dim+1) x (dim+1)  triangular factor of [A l], compressed
%     m     scalar             number of observations
%     dim   scalar             number of parameters
%
%   MEASURED ADVANTAGE
%     On the test problem with cond(A) = 2.5e7, error in the common
%     parameters versus a monolithic reference solve:
%         normal : 2.1e-06
%         qr     : 2.3e-07        (~9x better)
%     The gap widens as cond(A)^2 approaches 1/eps.
%
%   HOW THE TRIANGULAR FACTOR IS KEPT BOUNDED
%     After each add(), [R ; new rows] is re-triangularised and truncated to
%     dim+1 rows. Row dim+1 carries the accumulated residual norm: it is what
%     makes v'*v available without storing any residuals. Dropping it (or
%     failing to inject it on merge) makes sigma0 collapse silently -- see the
%     comment in merge().
%
%   EXAMPLE
%     acc = lsqvce.engine.QRSystem(10);
%     for k = 1:nBatch
%         [Ab, lb] = getWhitenedBatch(k);
%         acc.add(Ab, lb);
%     end
%     [x, Cxx, s0sq, dof] = acc.solve();
%
%   ERRORS
%     lsqvce:qr:badWidth       a block's column count disagrees with the first
%     lsqvce:qr:sizeMismatch   A and l disagree on their row count
%     lsqvce:qr:badMerge       merging two systems of different widths
%     lsqvce:qr:noData         solve() before anything was added
%     lsqvce:qr:rankDeficient  R is rank deficient and no datum was supplied
%     lsqvce:qr:singularR      a pivot of the factor is below RankTol (relative):
%                              a parameter or combination is not determined
%
%   WARNINGS
%     lsqvce:qr:cancellation   the residual mass was destroyed by
%                              cancellation -- either in a merge() fold of a
%                              pre-eliminated (normal-form) contribution, or
%                              by rotation round-off at extreme dynamic
%                              range -- and vtv/s0sq are returned as NaN,
%                              never as a plausible number (WP-50).
%                              lsqvce.engine.residuals() supersedes them.
%
%   ACCURACY
%     x from Householder/Givens updates: relative error ~ eps*cond(A) (the
%     reason this backend exists; see the measured comparison in
%     test_precision). vtv is carried as the squared residual row of the
%     augmented factor: exact to ~eps*mass through ROTATIONS. Since 6.7
%     (R7 WP-112) a pre-eliminated contribution arrives as rows with its
%     residual row, so this holds through merge() too; only a contribution
%     whose local block was REGULARISED is still folded in normal form and
%     recovers rho2 by subtraction (the WP-50 guard: NaN, never a wrong
%     plausible number).
%     lsqvce:qr:noRedundancy   n = rank, so sigma0 is undefined (the model
%                              reproduces the data exactly by construction)
%
%   See also LSQVCE.ENGINE.NORMALEQSYSTEM, LSQVCE.ENGINE.BLOCKACCUMULATOR
%
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-133: discount(k) takes k rows of
%   weight 0 out of the count m (mExcl keeps them; merge adds them). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-128: reducedRows also returns the
%   local block's coupling R12 and right-hand side r1. With assistance from
%   Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-112 (R7-03, R7-16):
%   reducedRows (local elimination on the factor), merge takes square-root
%   rows, the rank decision is RankTol on the pivots. With assistance from
%   Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-86 part 2 (R7-24):
%   squareRoot, the parameter columns' triangular factor
%   (BlockAccumulator.commonFactor). With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-08-20 (CEST, UTC+02:00), WP-50: honest vtv on both QR
%   loss paths (merge-fold cancellation guard + rotation-round-off
%   guard) -> NaN + lsqvce:qr:cancellation; vtvPoisoned property. With
%   assistance from Claude (Anthropic), model Claude Fable 5.
%
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    properties (SetAccess = private)
        R                       % triangular factor of [A l]
        m   = 0                 % observation count (weight-0 rows out, WP-133)
        mExcl = 0               % rows of weight 0 discounted from m (WP-133)
        dim                     % parameter count
        vtvPoisoned = false     % WP-50: a merge() fold destroyed the residual
                                % mass by cancellation; solve() returns
                                % vtv/s0sq = NaN, never a plausible number
    end

    properties (Constant)
        RankTol = 1e3 * eps     % R7-16: refuse a pivot below RankTol * max|diag R|
    end

    methods
        function o = QRSystem(dim)
        %QRSYSTEM  Construct a dim x dim QR solver (running [R z] triangular factor).
            validateattributes(dim, {'numeric'}, {'scalar','integer','positive'}, ...
                               'lsqvce.engine.QRSystem', 'dim', 1);
            o.dim = double(dim);
            o.R   = zeros(0, o.dim+1);
        end

        function add(o, A, l)
        %ADD  Accumulate one already-whitened batch by re-triangularisation.
            if size(A,2) ~= o.dim
                error('lsqvce:qr:badWidth', ...
                      'The batch has %d columns but the system has %d parameters.', ...
                      size(A,2), o.dim);
            end
            if size(A,1) ~= numel(l)
                error('lsqvce:qr:sizeMismatch', ...
                      'A has %d rows but l has %d elements.', size(A,1), numel(l));
            end
            if isempty(l), return, end

            o.R = triu(qr([o.R ; [A, l(:)]], 0));
            k   = min(size(o.R,1), o.dim+1);
            o.R = o.R(1:k, :);
            o.m = o.m + numel(l);
        end

        function mergeRaw(o, other)
        %MERGERAW  Fold in another RAW square-root system of the same
        %   dimension (v4.1) -- the day-distributed merge, WITHOUT a
        %   normal-form detour.
        %
        %   The other system's compressed factor is a set of pseudo-
        %   observation rows: [R_a ; R_b] re-triangularised has R1'*R1 =
        %   Ra1'*Ra1 + Rb1'*Rb1 and preserves column norms, so N, n AND the
        %   l'*l bookkeeping merge exactly -- including each factor's residual
        %   row, which carries the mass that keeps sigma0 honest. Validated:
        %   validation/v5verify.py tests 2-4 (equality with the fused loop to
        %   ~5e-13, order independence, and the sigma0 guard).
            if ~isa(other, 'lsqvce.engine.QRSystem') || other.dim ~= o.dim
                error('lsqvce:qr:badMerge', ...
                      'mergeRaw needs an lsqvce.engine.QRSystem of dimension %d.', o.dim);
            end
            if isempty(other.R), return, end
            o.R = triu(qr([o.R ; other.R], 0));
            k   = min(size(o.R,1), o.dim+1);
            o.R = o.R(1:k, :);
            o.m = o.m + other.m;
            o.mExcl = o.mExcl + other.mExcl;
        end

        function discount(o, k)
        %DISCOUNT  Take k rows of weight 0 out of the observation count
        %   (WP-133): they were added (their whitened rows are ~0) but are not
        %   observations -- a robust rejection, a zero data weight.
        %   INPUT  k  1x1 nonnegative integer
            o.m = o.m - k;  o.mExcl = o.mExcl + k;
        end

        function [Nll,Nlc,Ncc,nl,nc,ltl,m] = blocks(o, ic, il)        %BLOCKS  Recover the normal-equation blocks from R.
        %   N = R1'*R1 and n = R1'*r -- in normal form, for the local block's
        %   diagnostics, its regularisation and the back-substitution store.
        %   NOT the route of the common solution (R7 WP-112, R7-03: reducedRows
        %   and merge keep that square-root; this method squares cond).
            d  = o.dim;
            R1 = zeros(d, d);
            r  = zeros(d, 1);
            k  = min(size(o.R,1), d);
            R1(1:k, :) = o.R(1:k, 1:d);
            r(1:k)     = o.R(1:k, end);

            N   = R1.'*R1;
            nn  = R1.'*r;
            Nll = N(il,il);   Nlc = N(il,ic);   Ncc = N(ic,ic);
            nl  = nn(il);     nc  = nn(ic);

            ltl = r.'*r;
            if size(o.R,1) > d
                ltl = ltl + o.R(d+1,end)^2;     % the residual row
            end
            m = o.m;
        end

        function [rows, R11, R12, r1] = reducedRows(o, nl)
        %REDUCEDROWS  The square root of the contribution with the LAST nl
        %   columns (the locals) eliminated -- on the factor, never in normal
        %   form (R7 WP-112, R7-03).
        %
        %   R is re-triangularised with the local columns first; that is an
        %   orthogonal transform of R, so Rp'Rp = P'[A l]'[A l]P exactly to
        %   rounding without forming it. With Rp = [R11 R12 r1; 0 R22 r2;
        %   0 0 rho]:
        %     rows = [R22 r2; 0 rho]   (dim-nl)+1 columns: the reduced
        %                              contribution's pseudo-observations
        %     R11                      the local block's factor (its pivots
        %                              decide whether the elimination is sound)
        %     R12, r1                  its coupling to the commons and its
        %                              right-hand side (WP-128: the local block's
        %                              cofactor and solution are formed from
        %                              R11, R12, r1, never from the squared Nll)
        %   nl = 0 returns R itself and an empty R11.
            d  = o.dim;
            nc = d - nl;
            if nl == 0
                rows = o.R;  R11 = zeros(0, 0);  R12 = zeros(0, d);  r1 = zeros(0, 1);
                return
            end
            Rp  = triu(qr(o.R(:, [nc+1:d, 1:nc, d+1]), 0));
            k   = min(size(Rp, 1), d + 1);
            Rp  = Rp(1:k, :);
            kl  = min(nl, k);
            R11 = Rp(1:kl, 1:nl);
            R12 = Rp(1:kl, nl+1:d);  r1 = Rp(1:kl, d+1);
            rows = Rp(kl+1:k, nl+1:end);
        end

        function merge(o, red)
        %MERGE  Fold a pre-eliminated contribution back into the square root.
        %
        %   R7 WP-112: a contribution that carries square-root rows (red.Rows,
        %   from reducedRows) is stacked as it is -- the residual row rides
        %   along, nothing is subtracted, nothing squared. The normal-form
        %   route below remains for a contribution whose local block was
        %   regularised (the regularisation is defined on Nll).
            if ~isempty(red.Rows)
                if size(red.Rows, 2) ~= o.dim + 1
                    error('lsqvce:qr:badMerge', ...
                          'The contribution has %d columns but the system has %d parameters.', ...
                          size(red.Rows, 2) - 1, o.dim);
                end
                o.R = triu(qr([o.R ; red.Rows], 0));
                k   = min(size(o.R,1), o.dim+1);
                o.R = o.R(1:k, :);
                o.m = o.m + red.dof;
                return
            end
        %
        %   The reduced contribution is in normal form (N, n, ltl). Convert it
        %   to pseudo-observation rows: chol(N) = Rc and Rc'*rhs = n.
        %
        %   CRITICAL: the reduced observation vector satisfies
        %       ltl_red = ||rhs||^2 + rho^2
        %   where rho^2 is the residual mass NOT representable in the
        %   triangular factor. It must be appended as an extra pseudo-
        %   observation row. Omitting it makes vtv -- and therefore sigma0 and
        %   every formal error -- collapse towards zero, while still producing
        %   a perfectly plausible-looking parameter vector.
            [Rc, fail] = chol(red.N);
            if fail == 0
                rhs = Rc.' \ red.n;                       % PD: triangular, exact
            else
                % R5 M4 4e (finding H-16, remainder): a SINGLE contribution may
                % legitimately be rank deficient -- one system need not see
                % every common parameter; the SUM over systems does. So the
                % fold must not refuse here, and it must not put a singular
                % factor through `\` (MATLAB:singularMatrix ahead of an honest
                % refusal, the noise of H-16). The square root is taken in the
                % RANGE SPACE of N: r rows of Rc = S_r^1/2 V_r' for the
                % eigenvalues above rank()'s tolerance, and the reduced
                % right-hand side rhs = S_r^-1/2 V_r' n exactly (n lies in the
                % range of N for a consistent contribution). The overall
                % deficiency, if any remains after the sum, is refused by
                % solve()'s pivot test.
                [V, dv] = eig((red.N + red.N.')/2, 'vector');
                dv = real(dv);  V = real(V);
                keep = dv > max(size(red.N)) * eps(max(dv)) & dv > 0;
                Rc  = diag(sqrt(dv(keep))) * V(:, keep).';
                rhs = (V(:, keep).' * red.n) ./ sqrt(dv(keep));
            end
            rho2 = red.ltl - (rhs.'*rhs);

            % WP-50: this subtraction is NORMAL arithmetic inside the QR
            % backend -- the same l'Pl-minus cancellation the other routes
            % guard (and the Schur stage upstream subtracts once more). A
            % destroyed rho2 used to be silently treated as "no residual
            % row", collapsing vtv towards an EXACT ZERO sigma0 -- the wrong
            % plausible number observed in WP-49 (10.11 finding c). Same
            % rule as everywhere: refuse, warn once HERE (the site with the
            % numbers), and let solve() return NaN; lsqvce.engine.residuals()
            % supersedes it. A dof-free contribution legitimately carries
            % rho2 ~ 0 and is exempt.
            if red.dof > 0
                % Threshold 15, MEASURED at three anchors (2026-08-20,
                % pf50 r1/r1b): the 30/1e-4 block fixture folds at
                % digitsLost = 13.4 (ltl_red 7.19e16 / rho2 3.06e3) and its
                % sigma0 was measured GOOD to 8e-5 -- the eps*10^dl bound
                % is ~50x pessimistic here; the WP-49 deep probe sits at
                % ~19 digits with rho2 driven NEGATIVE (destroyed); the
                % first pin at 13 knife-edged the proven-good fixture.
                % 15 leaves 1.6 digits above measured-good and ~4 below
                % measured-destroyed. The normal route's 8 does NOT
                % transfer: its error is kappa^2-amplified, this one's not.
                % R5 M5 5b': a contribution with NO residual mass at all
                % (ltl_red = 0: an exact zero right-hand side) carries
                % rho2 = 0 legitimately -- nothing was there to lose; only a
                % nonzero mass that came out <= 0 or 15 digits smaller was
                % destroyed.
                destroyed = red.ltl > 0 && (rho2 <= 0 || ...
                    log10(max(red.ltl / max(rho2, realmin), 1)) > 15);
                if destroyed
                    lsqvce.util.log.warn('lsqvce:qr:cancellation', ...
                        ['The residual mass of a pre-eliminated contribution ' ...
                         'was destroyed by cancellation in the normal-form ' ...
                         'fold (ltl_red = %.6e, rho^2 = %.6e, dof = %d). ' ...
                         'vtv/s0sq of this system will be NaN, never a ' ...
                         'plausible number; lsqvce.engine.residuals() recomputes ' ...
                         'them honestly and SUPERSEDES them.'], ...
                         red.ltl, rho2, red.dof);
                    o.vtvPoisoned = true;
                end
            end

            blk = [Rc, rhs];
            if rho2 > 0 && ~o.vtvPoisoned
                resRow      = zeros(1, o.dim+1);
                resRow(end) = sqrt(rho2);
                blk         = [blk ; resRow];
            end

            o.R = triu(qr([o.R ; blk], 0));
            k   = min(size(o.R,1), o.dim+1);
            o.R = o.R(1:k, :);
            o.m = o.m + red.dof;
        end

        function [x, Cxx, s0sq, dof, vtv] = solve(o)
        %SOLVE  Triangular back-substitution.
            if o.m == 0
                error('lsqvce:qr:noData', 'No observations were accumulated.');
            end

            d = o.dim;
            if size(o.R,1) < d
                error('lsqvce:qr:rankDeficient', ...
                    ['The triangular factor has only %d of %d rows: the common ' ...
                     'block is rank deficient by %d. Add datum constraints.'], ...
                     size(o.R,1), d, d - size(o.R,1));
            end

            R1 = o.R(1:d, 1:d);
            r  = o.R(1:d, end);

            dR = abs(diag(R1));
            % R5 M4 4b (finding H-16, confirmed regime-dependent on the CI
            % runner): the rank test must NOT sit at the roundoff boundary. A
            % Householder QR of a rank-deficient block leaves a trailing pivot
            % of order eps*|A|, above or below eps*max(dR) depending on the
            % LAPACK build -- refused on Windows, solved on Linux. R7-16 (WP-112):
            % the decision is RankTol (1e3 eps) on the pivots, the class's own
            % threshold -- measured gap: free networks <= 5e-16, full rank at
            % cond(A) = 1e12 >= 1.9e-12 (mirrors/wp112).
            rtol = lsqvce.engine.QRSystem.RankTol * max(dR);
            if min(dR) <= rtol
                [~, k] = min(dR);
                error('lsqvce:qr:singularR', ...
                    ['The triangular factor is singular at parameter %d ' ...
                     '(|R(%d,%d)| = %.3e vs max %.3e, tolerance %.3e). That parameter, ' ...
                     'or a combination involving it, is not determined by the data. ' ...
                     'Add a datum constraint or remove it.'], k, k, k, dR(k), max(dR), rtol);
            end

            x   = R1 \ r;
            vtv = 0;
            if size(o.R,1) > d
                vtv = o.R(d+1,end)^2;
            end
            if o.vtvPoisoned
                vtv = NaN;                 % refused at merge(); warned there
            elseif o.m > d
                % WP-50 rotation guard, MEASURED (pf50 ramp, 2026-08-20):
                % the rotation route is robust far beyond the fold -- at
                % log10(mass/vtv) = 21.3 the error is only 4e-7. The guard
                % at 24 (relerr ~ 2*eps*10^(dl/2) ~ 2e-4) is a backstop for
                % the exact-zero underflow and beyond; healthy use never
                % reaches it.
                mass = o.R(1:min(d+1, size(o.R,1)), end);
                mass = mass.'*mass;
                if mass > 0 && (vtv == 0 || ...
                        log10(max(mass / max(vtv, realmin), 1)) > 24)
                    lsqvce.util.log.warn('lsqvce:qr:cancellation', ...
                        ['The residual row of the triangular factor lost its ' ...
                         'signal to rotation round-off (quadratic mass ' ...
                         '%.6e, vtv %.6e). RETURNING NaN rather than a ' ...
                         'plausible-looking number; lsqvce.engine.residuals() ' ...
                         'recomputes it honestly and SUPERSEDES it.'], ...
                         mass, vtv);
                    vtv = NaN;
                end
            end

            dof = o.m - d;
            if dof <= 0
                error('lsqvce:qr:noRedundancy', ...
                    ['Redundancy is %d (m = %d observations, %d parameters). ' ...
                     'sigma0 is undefined.'], dof, o.m, d);
            end

            s0sq = vtv / dof;
            Ri   = R1 \ eye(d);
            Cxx  = s0sq * (Ri*Ri.');
            Cxx  = (Cxx + Cxx.')/2;
        end

        function R1 = squareRoot(o)
        %SQUAREROOT  The d x d upper-triangular factor R1 of the parameter
        %   columns (A'*A = R1'*R1), zero-padded while fewer than d rows exist.
            d  = o.dim;
            R1 = zeros(d, d);
            k  = min(size(o.R,1), d);
            R1(1:k, :) = o.R(1:k, 1:d);
        end

        function [N, n, ltl, m] = commonBlocks(o)
        %COMMONBLOCKS  Recover the normals from the triangular factor.
            d  = o.dim;
            R1 = zeros(d, d);  r = zeros(d, 1);
            k  = min(size(o.R,1), d);
            R1(1:k,:) = o.R(1:k, 1:d);
            r(1:k)    = o.R(1:k, end);
            N   = R1.'*R1;
            N   = (N + N.')/2;
            n   = R1.'*r;
            ltl = r.'*r;
            if size(o.R,1) > d
                ltl = ltl + o.R(d+1,end)^2;
            end
            m = o.m;
        end
    end
end
