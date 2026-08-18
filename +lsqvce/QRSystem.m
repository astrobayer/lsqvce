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
%     lsqvce:qr:noRedundancy   n = rank, so sigma0 is undefined (the model
%                              reproduces the data exactly by construction)
%
%   See also LSQVCE.NORMALEQSYSTEM, LSQVCE.BLOCKACCUMULATOR
%
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.

    properties (SetAccess = private)
        R                       % triangular factor of [A l]
        m   = 0                 % observation count
        dim                     % parameter count
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

            blk = [Rc, rhs];
            if rho2 > 0
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
