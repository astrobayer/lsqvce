classdef NormalEqSystem < handle
%LSQVCE.NORMALEQSYSTEM  Normal-equation accumulator: N = A'*A, n = A'*l.
%
%   The classic approach. Cheap and simple, but it SQUARES the condition
%   number of the design matrix: information below sqrt(eps)*||A|| relative to
%   the leading singular value is destroyed by the product A'*A and cannot be
%   recovered by any solver afterwards.
%
%   Use LSQVCE.QRSYSTEM instead unless you have measured that the conditioning is
%   benign and you need the ~2x speed.
%
%   PROPERTIES (read-only)
%     N     dim x dim  accumulated normal matrix
%     n     dim x 1    accumulated right-hand side
%     ltl   scalar     accumulated l'*l, needed for the residual sum of squares
%     m     scalar     number of observations accumulated
%     dim   scalar     number of parameters
%
%   METHODS
%     add(A, l)             accumulate one (already whitened) batch
%     blocks(ic, il)        split into common/local blocks for pre-elimination
%     merge(red)            fold in a reduced contribution (LSQVCE.REDUCEDCONTRIBUTION)
%     solve()               [x, Cxx, s0sq, dof, vtv]
%
%   THE ltl BOOKKEEPING
%     vtv = l'*l - n'*x is computed from accumulated scalars, never from the
%     residuals themselves (which are not stored). Every path that adds
%     information must also add its share of l'*l, or sigma0 collapses. This
%     is the single easiest thing to get wrong in a batch accumulator.
%
%   EXAMPLE
%     acc = lsqvce.NormalEqSystem(10);
%     for k = 1:nBatch
%         [Ab, lb] = getWhitenedBatch(k);
%         acc.add(Ab, lb);
%     end
%     [x, Cxx, s0sq, dof] = acc.solve();
%
%   ERRORS
%     lsqvce:normalEq:badWidth  a block's column count disagrees with the first
%     lsqvce:normalEq:badMerge  merging two systems of different widths
%     lsqvce:normalEq:noData    solve() before anything was added
%
%   See also LSQVCE.QRSYSTEM, LSQVCE.BLOCKACCUMULATOR
%
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.

    properties (SetAccess = private)
        N                       % dim x dim normal matrix
        n                       % dim x 1   right-hand side
        ltl = 0                 % scalar    l'*l
        m   = 0                 % scalar    observation count
        dim                     % scalar    parameter count
    end

    methods
        function o = NormalEqSystem(dim)
        %NORMALEQSYSTEM  Construct a dim x dim normal-equation solver (N,b accumulation).
            validateattributes(dim, {'numeric'}, {'scalar','integer','positive'}, ...
                               'lsqvce.NormalEqSystem', 'dim', 1);
            o.dim = double(dim);
            o.N   = zeros(o.dim);
            o.n   = zeros(o.dim,1);
        end

        function add(o, A, l)
        %ADD  Accumulate one already-whitened batch.
            if size(A,2) ~= o.dim
                error('lsqvce:normalEq:badWidth', ...
                      'The batch has %d columns but the system has %d parameters.', ...
                      size(A,2), o.dim);
            end
            if size(A,1) ~= numel(l)
                error('lsqvce:normalEq:sizeMismatch', ...
                      'A has %d rows but l has %d elements.', size(A,1), numel(l));
            end
            if isempty(l), return, end

            l = l(:);
            o.N   = o.N   + A.'*A;
            o.n   = o.n   + A.'*l;
            o.ltl = o.ltl + l.'*l;
            o.m   = o.m   + numel(l);
        end

        function [Nll,Nlc,Ncc,nl,nc,ltl,m] = blocks(o, ic, il)
        %BLOCKS  Partition into common (ic) and local (il) index sets.
            Nll = o.N(il,il);   Nlc = o.N(il,ic);   Ncc = o.N(ic,ic);
            nl  = o.n(il);      nc  = o.n(ic);
            ltl = o.ltl;        m   = o.m;
        end

        function merge(o, red)
        %MERGE  Fold in a pre-eliminated contribution.
            o.N   = o.N   + red.N;
            o.n   = o.n   + red.n;
            o.ltl = o.ltl + red.ltl;
            o.m   = o.m   + red.dof;
        end

        function mergeRaw(o, other)
        %MERGERAW  Fold in another RAW system of the same dimension (v4.1).
        %
        %   The reduction behind day-distributed accumulation: partial sums of
        %   N, n, l'l, m add exactly (validation/v5verify.py test 6: the
        %   summation-order gap vs the fused loop is ~1e-14 relative).
            if ~isa(other, 'lsqvce.NormalEqSystem') || other.dim ~= o.dim
                error('lsqvce:normalEq:badMerge', ...
                      'mergeRaw needs an lsqvce.NormalEqSystem of dimension %d.', o.dim);
            end
            o.N   = o.N   + other.N;
            o.n   = o.n   + other.n;
            o.ltl = o.ltl + other.ltl;
            o.m   = o.m   + other.m;
        end

        function [x, Cxx, s0sq, dof, vtv] = solve(o)
        %SOLVE  Cholesky solve with a rank check.
            if o.m == 0
                error('lsqvce:normalEq:noData', 'No observations were accumulated.');
            end

            [R, fail] = chol(o.N);
            if fail > 0
                dg = lsqvce.localDiagnostics(o.N, struct('TruncTol',1e-12));
                error('lsqvce:normalEq:notSPD', ...
                    ['The normal matrix is not positive definite (rank ' ...
                     'deficiency %d, cond %.3e). The common block has a datum ' ...
                     'defect: %d parameter combination(s) are not determined ' ...
                     'by the data. Add datum constraints, or switch to ' ...
                     '''Backend'',''qr'' which will report the singular ' ...
                     'direction.'], dg.rankDef, dg.condNll, dg.rankDef);
            end

            x    = R \ (R.' \ o.n);
            vtv  = o.ltl - o.n.'*x;
            dof  = o.m - o.dim;

            if dof <= 0
                error('lsqvce:normalEq:noRedundancy', ...
                    ['Redundancy is %d (m = %d observations, %d parameters). ' ...
                     'sigma0 is undefined.'], dof, o.m, o.dim);
            end
            if vtv < 0
                lsqvce.util.log.warn('lsqvce:normalEq:negativeVtv', ...
                    ['The weighted sum of squares v''Pv = %.3e is negative, ' ...
                     'which is impossible. This is a cancellation artefact of ' ...
                     'the normal equations -- the condition number is too ' ...
                     'large for this backend. Switch to ''Backend'',''qr''.'], vtv);
                vtv = 0;
            end

            s0sq = vtv / dof;
            Ri   = R \ eye(o.dim);
            Cxx  = s0sq * (Ri*Ri.');
            Cxx  = (Cxx + Cxx.')/2;
        end

        function [N, n, ltl, m] = commonBlocks(o)
        %COMMONBLOCKS  The full accumulated normals (the global system has no
        %   local block, so these ARE the common blocks).
            N = (o.N + o.N.')/2;
            n = o.n;
            ltl = o.ltl;
            m = o.m;
        end
    end
end
