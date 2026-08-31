classdef ReducedContribution
%LSQVCE.ENGINE.REDUCEDCONTRIBUTION  One system's contribution after pre-elimination.
%
%   The Schur complement of a system, ready to be merged into the global
%   accumulator. Immutable by design: once a system has been reduced, its
%   contribution is fixed.
%
%   PROPERTIES
%     N     nc x nc  Ncc - Nlc'*inv(Nll)*Nlc      reduced normal matrix
%     n     nc x 1   nc  - Nlc'*inv(Nll)*nl       reduced right-hand side
%     ltl   scalar   l'l - nl'*inv(Nll)*nl        reduced residual mass
%     dof   scalar   m - nl                       observations minus locals
%
%   ltl IS NOT OPTIONAL
%     It carries the residual sum of squares of the pre-eliminated system. Both
%     backends must inject it: LSQVCE.ENGINE.NORMALEQSYSTEM adds it to its scalar
%     accumulator, LSQVCE.ENGINE.QRSYSTEM appends sqrt(ltl - ||rhs||^2) as an extra
%     pseudo-observation row. Drop it and sigma0 collapses toward zero while
%     the parameters still look fine.
%
%   ERRORS
%     lsqvce:reduced:notSquare     the contribution block is not square
%     lsqvce:reduced:sizeMismatch  it disagrees with the common block
%
%   See also LSQVCE.ENGINE.BLOCKACCUMULATOR, LSQVCE.ENGINE.QRSYSTEM, LSQVCE.ENGINE.NORMALEQSYSTEM
%
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    properties (SetAccess = immutable)
        N       % reduced normal matrix
        n       % reduced right-hand side
        ltl     % reduced residual mass
        dof     % degrees of freedom contributed
    end

    methods
        function o = ReducedContribution(N, n, ltl, dof)
            narginchk(4,4);
            if size(N,1) ~= size(N,2)
                error('lsqvce:reduced:notSquare', 'N must be square.');
            end
            if numel(n) ~= size(N,1)
                error('lsqvce:reduced:sizeMismatch', ...
                      'N is %dx%d but n has %d elements.', size(N,1), size(N,2), numel(n));
            end
            o.N   = N;
            o.n   = n(:);
            o.ltl = ltl;
            o.dof = dof;
        end
    end
end
