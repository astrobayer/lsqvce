classdef NormalEqSystem < handle
%LSQVCE.ENGINE.NORMALEQSYSTEM  Normal-equation accumulator: N = A'*A, n = A'*l.
%
%   The classic approach. Cheap and simple, but it SQUARES the condition
%   number of the design matrix: information below sqrt(eps)*||A|| relative to
%   the leading singular value is destroyed by the product A'*A and cannot be
%   recovered by any solver afterwards.
%
%   Use LSQVCE.ENGINE.QRSYSTEM instead unless you have measured that the conditioning is
%   benign and you need the ~2x speed.
%
%   PROPERTIES (read-only)
%     N     dim x dim  accumulated normal matrix
%     n     dim x 1    accumulated right-hand side
%     ltl   scalar     accumulated l'*l, needed for the residual sum of squares
%     ltlMass scalar   the l'*l every merged contribution was reduced FROM
%                      (ReducedContribution.ltlMass; = ltl without locals):
%                      the scale of the cancellation guard (WP-151, S2)
%     m     scalar     number of observations accumulated
%     dim   scalar     number of parameters
%
%   METHODS
%     add(A, l)             accumulate one (already whitened) batch
%     discount(k)           take k weight-0 rows out of the count m (WP-133)
%     blocks(ic, il)        split into common/local blocks for pre-elimination
%     merge(red)            fold in a reduced contribution (LSQVCE.ENGINE.REDUCEDCONTRIBUTION)
%     solve()               [x, Cxx, s0sq, dof, vtv]
%
%   THE ltl BOOKKEEPING
%     vtv = l'*l - n'*x is computed from accumulated scalars, never from the
%     residuals themselves (which are not stored). Every path that adds
%     information must also add its share of l'*l, or sigma0 collapses. This
%     is the single easiest thing to get wrong in a batch accumulator.
%
%   EXAMPLE
%     acc = lsqvce.engine.NormalEqSystem(10);
%     for k = 1:nBatch
%         [Ab, lb] = getWhitenedBatch(k);
%         acc.add(Ab, lb);
%     end
%     [x, Cxx, s0sq, dof] = acc.solve();
%
%   ACCURACY
%     x from a Cholesky solve of N = A'A: relative error ~ eps*cond(A)^2
%     (measured against the LAPACK QR reference in
%     test_precision/testQRBeatsNormalOnIllConditioned: agrees at benign
%     conditioning, loses ALL digits at cond(A) = 1e8 where QR keeps 8).
%     vtv is formed as l'Pl - n'x and is cancellation-prone when the fit is
%     tight; since WP-49 a value that is negative or has lost > 8 digits is
%     returned as NaN, never as a plausible-looking number -- the honest
%     value comes from LSQVCE.ENGINE.RESIDUALS' second pass, which supersedes it.
%     The digits are counted against ltlMass, the l'Pl before the local
%     blocks were eliminated (WP-151): the Schur step subtracts first, and a
%     signal the locals absorbed (1e7 sigma) cost 7.8 % of vtv, 1e8 sigma
%     1000 %, while the reduced l'Pl showed nothing (test_precision
%     (testSchurCancellationIsGuarded)).
%
%   ERRORS
%     Accumulation:
%       lsqvce:normalEq:badWidth      a block's column count disagrees with
%                                     the first
%       lsqvce:normalEq:badMerge      merging two systems of different widths
%     Solve:
%       lsqvce:normalEq:noData        solve() before anything was added
%       lsqvce:normalEq:notSPD        the normal matrix is not positive
%                                     definite (rank defect)
%       lsqvce:normalEq:noRedundancy  m <= dim, sigma0 undefined
%
%   WARNINGS
%     lsqvce:normalEq:cancellation    vtv negative OR > 8 digits lost to the
%                                     l'Pl - n'x cancellation -> vtv = NaN.
%                                     (The former 'negativeVtv' identifier in
%                                     this namespace was RETIRED in WP-49 --
%                                     breaking: one failure mode, one name.)
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:normalEq:sizeMismatch               A has X rows but l has X elements.*A; o.n = o.n + A
%
%   See also LSQVCE.ENGINE.QRSYSTEM, LSQVCE.ENGINE.BLOCKACCUMULATOR
%
%   Updated 2026-10-18 (CEST, UTC+02:00), WP-151 (review 12.84, S2): the
%   cancellation guard counts the digits lost against ltlMass, the l'Pl before
%   the local elimination. With assistance from Claude (Anthropic), model
%   Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-133: discount(k) takes k rows of
%   weight 0 out of the count m (mExcl keeps them; merge adds them). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-08-20 (CEST, UTC+02:00), WP-49: a poisoned vtv is returned
%   as NaN (negative OR > 8 digits lost), never clamped to 0 or left as a
%   wrong plausible number; both warnings name the lsqvce.engine.residuals repair
%   path. With assistance from Claude (Anthropic), model Claude Fable 5.
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    properties (SetAccess = private)
        N                       % dim x dim normal matrix
        n                       % dim x 1   right-hand side
        ltl = 0                 % scalar    l'*l
        ltlMass = 0             % scalar    l'*l before any local elimination (WP-151)
        m   = 0                 % scalar    observation count (weight-0 rows out, WP-133)
        mExcl = 0               % scalar    rows of weight 0 discounted from m (WP-133)
        dim                     % scalar    parameter count
    end

    methods
        function o = NormalEqSystem(dim)
        %NORMALEQSYSTEM  Construct a dim x dim normal-equation solver (N,b accumulation).
            validateattributes(dim, {'numeric'}, {'scalar','integer','positive'}, ...
                               'lsqvce.engine.NormalEqSystem', 'dim', 1);
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
            o.ltlMass = o.ltlMass + l.'*l;
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
            o.ltlMass = o.ltlMass + red.ltlMass;    % WP-151: the mass it was reduced from
            o.m   = o.m   + red.dof;
        end

        function mergeRaw(o, other)
        %MERGERAW  Fold in another RAW system of the same dimension (v4.1).
        %
        %   The reduction behind day-distributed accumulation: partial sums of
        %   N, n, l'l, m add exactly (validation/v5verify.py test 6: the
        %   summation-order gap vs the fused loop is ~1e-14 relative).
            if ~isa(other, 'lsqvce.engine.NormalEqSystem') || other.dim ~= o.dim
                error('lsqvce:normalEq:badMerge', ...
                      'mergeRaw needs an lsqvce.engine.NormalEqSystem of dimension %d.', o.dim);
            end
            o.N   = o.N   + other.N;
            o.n   = o.n   + other.n;
            o.ltl = o.ltl + other.ltl;
            o.ltlMass = o.ltlMass + other.ltlMass;
            o.m   = o.m   + other.m;
            o.mExcl = o.mExcl + other.mExcl;
        end

        function discount(o, k)
        %DISCOUNT  Take k rows of weight 0 out of the observation count
        %   (WP-133): they were added (their whitened rows are ~0) but are not
        %   observations -- a robust rejection, a zero data weight.
        %   INPUT  k  1x1 nonnegative integer
            o.m = o.m - k;  o.mExcl = o.mExcl + k;
        end

        function [x, Cxx, s0sq, dof, vtv] = solve(o)
        %SOLVE  Cholesky solve with a rank check.
            if o.m == 0
                error('lsqvce:normalEq:noData', 'No observations were accumulated.');
            end

            [R, fail] = chol(o.N);
            if fail > 0
                dg = lsqvce.engine.localDiagnostics(o.N, struct('TruncTol',1e-12));
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
            % WP-49: the flat path's honesty applied here too. vtv was formed
            % as l'Pl - n'x -- a difference of two nearly equal numbers when
            % the fit is tight. A poisoned value must never leave as a
            % plausible-looking number: NaN gets investigated, a wrong number
            % gets quoted (measured on the flat path: 89% wrong at
            % sigma = 1e-6, factor 1e6 at 1e-9).
            % ONE failure mode, ONE identifier: a negative vtv is just the
            % extreme of the same cancellation (sign at the eps-garbage
            % level is a BLAS coin flip -- WP-47 r1e), so both branches
            % speak as lsqvce:normalEq:cancellation. The former
            % 'negativeVtv' name in this namespace is RETIRED (breaking,
            % WP-49; error identifiers are names and get no aliases).
            % WP-151 (S2): the digits against the mass BEFORE the local
            % elimination -- the Schur step subtracted first (ltlMass >= ltl)
            poisoned = vtv < 0;
            mass = max(o.ltl, o.ltlMass);
            if ~poisoned && mass > 0
                poisoned = log10(max(mass / max(vtv, realmin), 1)) > 8;
            end
            if poisoned
                lsqvce.util.log.warn('lsqvce:normalEq:cancellation', ...
                    ['vtv was formed as l''Pl - n''x and the cancellation ' ...
                     'destroyed it (l''Pl = %.6e before the local elimination, vtv = %.6e%s). RETURNING ' ...
                     'NaN rather than a plausible-looking number. ' ...
                     'lsqvce.engine.residuals() recomputes it honestly from a ' ...
                     'second pass and SUPERSEDES this value; to avoid the ' ...
                     'problem entirely, switch to ''Backend'',''qr''.'], ...
                     mass, vtv, ternaryChar(vtv < 0, ', negative', ''));
                vtv = NaN;
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

function c = ternaryChar(cond, a, b)
%TERNARYCHAR  Inline char pick for the warning text.
    if cond, c = a; else, c = b; end
end
