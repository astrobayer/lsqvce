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
%     Rows  k x (nc+1) | []  the SAME contribution as square-root rows
%                    [R22 r2; 0 rho] (LSQVCE.ENGINE.QRSYSTEM.reducedRows):
%                    Rows'Rows = [N n; n' ltl]. Given, QRSystem.merge stacks
%                    them and never touches N (R7 WP-112, R7-03); N, n, ltl
%                    are then derived from them, for the normal backend and
%                    the diagnostics.
%
%   c = LSQVCE.ENGINE.REDUCEDCONTRIBUTION(N, n, ltl, dof)
%   c = LSQVCE.ENGINE.REDUCEDCONTRIBUTION([], [], [], dof, Rows)
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
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-112 (R7-03): the optional
%   square-root form Rows. With assistance from Claude (Anthropic), model
%   Claude Opus 5.5.
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
        Rows = []  % the square-root form, or [] (R7 WP-112)
    end

    methods
        function o = ReducedContribution(N, n, ltl, dof, Rows)
            narginchk(4,5);
            if nargin == 5 && ~isempty(Rows)
                o.Rows = Rows;
                Rc  = Rows(:, 1:end-1);  rv = Rows(:, end);
                N   = Rc.' * Rc;  N = (N + N.')/2;
                n   = Rc.' * rv;
                ltl = rv.' * rv;
            end
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
