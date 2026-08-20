classdef QRSystem < handle
%LSQVCE.QRSYSTEM  Square-root information accumulator (Householder QR).
%
%   Maintains the upper-triangular factor R of the stacked [A l], so that
%   A'*A = R1'*R1 without A'*A ever being formed. The condition number of R is
%   cond(A), not cond(A)^2 -- roughly twice the flops of LSQVCE.NORMALEQSYSTEM in
%   exchange for double the effective precision.
%
%   This is the default backend, and the right one whenever the local block
%   might be poorly conditioned (near-collinear empirical accelerations, short
%   arcs, weakly separated bias/scale pairs -- i.e. most real problems).
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
%     acc = lsqvce.QRSystem(10);
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
%     lsqvce:qr:singularR      the triangular factor cannot be back-substituted
%
%   WARNINGS
%     lsqvce:qr:cancellation   the residual mass was destroyed by
%                              cancellation -- either in a merge() fold of a
%                              pre-eliminated (normal-form) contribution, or
%                              by rotation round-off at extreme dynamic
%                              range -- and vtv/s0sq are returned as NaN,
%                              never as a plausible number (WP-50).
%                              lsqvce.residuals() supersedes them.
%
%   ACCURACY
%     x from Householder/Givens updates: relative error ~ eps*cond(A) (the
%     reason this backend exists; see the measured comparison in
%     test_precision). vtv is carried as the squared residual row of the
%     augmented factor: exact to ~eps*mass through ROTATIONS, but a merge()
%     of a pre-eliminated contribution recovers it by NORMAL-form
%     subtraction and inherits that route's 8-digit cancellation limit --
%     both guarded since WP-50 (NaN, never a wrong plausible number; the
%     WP-49 deep-regime probe measured an EXACT-ZERO sigma0 before).
%     lsqvce:qr:noRedundancy   n = rank, so sigma0 is undefined (the model
%                              reproduces the data exactly by construction)
%
%   See also LSQVCE.NORMALEQSYSTEM, LSQVCE.BLOCKACCUMULATOR
%
%   Updated 2026-08-20 (CEST, UTC+02:00), WP-50: honest vtv on both QR
%   loss paths (merge-fold cancellation guard + rotation-round-off
%   guard) -> NaN + lsqvce:qr:cancellation; vtvPoisoned property. With
%   assistance from Claude (Anthropic), model Claude Fable 5.
%
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.

    properties (SetAccess = private)
        R                       % triangular factor of [A l]
        m   = 0                 % observation count
        dim                     % parameter count
        vtvPoisoned = false     % WP-50: a merge() fold destroyed the residual
                                % mass by cancellation; solve() returns
                                % vtv/s0sq = NaN, never a plausible number
    end

    methods
        function o = QRSystem(dim)
        %QRSYSTEM  Construct a dim x dim QR solver (running [R z] triangular factor).
            validateattributes(dim, {'numeric'}, {'scalar','integer','positive'}, ...
                               'lsqvce.QRSystem', 'dim', 1);
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
            if ~isa(other, 'lsqvce.QRSystem') || other.dim ~= o.dim
                error('lsqvce:qr:badMerge', ...
                      'mergeRaw needs an lsqvce.QRSystem of dimension %d.', o.dim);
            end
            if isempty(other.R), return, end
            o.R = triu(qr([o.R ; other.R], 0));
            k   = min(size(o.R,1), o.dim+1);
            o.R = o.R(1:k, :);
            o.m = o.m + other.m;
        end

        function [Nll,Nlc,Ncc,nl,nc,ltl,m] = blocks(o, ic, il)        %BLOCKS  Recover the normal-equation blocks from R.
        %   N = R1'*R1 and n = R1'*r. The squaring happens here, but on a
        %   matrix that is already triangular and well scaled, which is a very
        %   different proposition from squaring the raw design matrix.
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

        function merge(o, red)
        %MERGE  Fold a pre-eliminated contribution back into the square root.
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
            if fail > 0
                [V, dv] = eig((red.N + red.N.')/2, 'vector');
                dv = max(real(dv), 0);
                Rc = diag(sqrt(dv)) * real(V).';
            end

            rhs  = Rc.' \ red.n;
            rho2 = red.ltl - (rhs.'*rhs);

            % WP-50: this subtraction is NORMAL arithmetic inside the QR
            % backend -- the same l'Pl-minus cancellation the other routes
            % guard (and the Schur stage upstream subtracts once more). A
            % destroyed rho2 used to be silently treated as "no residual
            % row", collapsing vtv towards an EXACT ZERO sigma0 -- the wrong
            % plausible number observed in WP-49 (10.11 finding c). Same
            % rule as everywhere: refuse, warn once HERE (the site with the
            % numbers), and let solve() return NaN; lsqvce.residuals()
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
                destroyed = rho2 <= 0 || (red.ltl > 0 && ...
                    log10(max(red.ltl / max(rho2, realmin), 1)) > 15);
                if destroyed
                    lsqvce.util.log.warn('lsqvce:qr:cancellation', ...
                        ['The residual mass of a pre-eliminated contribution ' ...
                         'was destroyed by cancellation in the normal-form ' ...
                         'fold (ltl_red = %.6e, rho^2 = %.6e, dof = %d). ' ...
                         'vtv/s0sq of this system will be NaN, never a ' ...
                         'plausible number; lsqvce.residuals() recomputes ' ...
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
            if min(dR) <= eps*max(dR)
                [~, k] = min(dR);
                error('lsqvce:qr:singularR', ...
                    ['The triangular factor is singular at parameter %d ' ...
                     '(|R(%d,%d)| = %.3e vs max %.3e). That parameter, or a ' ...
                     'combination involving it, is not determined by the data. ' ...
                     'Add a datum constraint or remove it.'], k, k, k, dR(k), max(dR));
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
                         'plausible-looking number; lsqvce.residuals() ' ...
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
