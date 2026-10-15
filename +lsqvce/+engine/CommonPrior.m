classdef CommonPrior
%LSQVCE.ENGINE.COMMONPRIOR  Basis-agnostic regularisation of the COMMON parameters.
%
%   Adds a prior  K*x ~ x0  with weight 1/sigma^2, so the reduced system becomes
%
%       (Ncc_red + (1/sigma^2)*K'K) x = ncc_red + (1/sigma^2)*K'*x0
%
%   The toolbox never needs to know what basis you are using. K is just a matrix.
%
%   CONSTRUCTOR
%     pri = LSQVCE.ENGINE.COMMONPRIOR(K)                    x0 = 0, sigma = 1
%     pri = LSQVCE.ENGINE.COMMONPRIOR(K, x0)
%     pri = LSQVCE.ENGINE.COMMONPRIOR(K, x0, sigma)
%     pri = LSQVCE.ENGINE.COMMONPRIOR(K, x0, sigma, group)
%
%   group (nk x 1 positive integers, default ones): the declared group of each
%   row -- lsqvce.fit stacks data.priors(k) as group k. The group is the unit of
%   meaning for the datum test (BlockAccumulator.datumPart, R7 WP-119): a group
%   the data and every OTHER group are blind to is a DATUM and gets the
%   minimum-constraint covariance; a group only partly in the null space is
%   information. One group (the default) is the whole prior, as before 6.11.
%
%   K IS THE SCIENCE -- CHOOSE IT DELIBERATELY
%     K = eye(nc)                     plain damping (rarely what you want)
%     K = diag(1./sigma_l)            Kaula: you build sigma_l = 1e-5./l.^2 from
%                                     your own degree mapping. The toolbox never
%                                     sees "spherical harmonic".
%     K = L                           mascon smoothness: L is the graph Laplacian
%                                     of the mascon adjacency (sparse). Penalises
%                                     differences between neighbouring mascons.
%     K = chol(inv(C_signal))         full signal-covariance prior from a
%                                     background field.
%     K = [zeros(k,n1) eye(k) zeros]  regularise a SUBSET only.
%
%     For a "regularise toward a background model" workflow, set x0 to that
%     model's coefficients rather than zero.
%
%   SIGMA FROM VCE, NOT FROM A TUNING KNOB
%     Treat the prior as a pseudo-observation group with its own variance
%     component. Then lambda = sigma0^2 / sigma_prior^2 falls straight out of the
%     Helmert iteration, and the prior's REDUNDANCY tells you how many degrees of
%     freedom it is actually absorbing. That is far more informative than a bare
%     lambda, and it works identically for spherical harmonics, mascons, or
%     anything else.
%
%   THE DIAGNOSTIC THAT MATTERS
%     [rK, absorbed] = pri.redundancy(Lc)
%
%       absorbed  =  nk - rK  =  trace(Sxx*K'K)/sigma^2
%
%     is the effective number of parameters constrained by the prior.
%       absorbed -> 0   the prior is doing nothing; you do not need it
%       absorbed -> nk  the prior is doing everything; you are recovering x0 and
%                       the data is contributing almost nothing
%     Report this. A lambda without it is uninterpretable.
%
%   THIS BIASES THE LOCAL PARAMETERS TOO
%     Regularising the common block propagates into y_s through back-
%     substitution -- the mirror image of the warning in LSQVCE.ENGINE.APPLYREGULARISATION.
%     There is no free lunch in either direction.
%
%   L2 ONLY
%     This is a normal-equation framework, so the prior is Gaussian / L2. If you
%     want mascon L1 or total-variation regularisation -- which much of the
%     mascon literature actually uses -- that is NOT a normal-equation problem
%     and does not belong here. Use a proper convex solver.
%
%   EXAMPLE
%     % Kaula on a spherical-harmonic field
%     l   = degreeOfEachCoefficient(nc);          % you supply this
%     sig = 1e-5 ./ l.^2;
%     pri = lsqvce.engine.CommonPrior(diag(1./sig), [], 1.0);
%
%     % mascon smoothness
%     L   = masconLaplacian(adjacency);           % sparse, nk x nc
%     pri = lsqvce.engine.CommonPrior(L, [], 1.0);
%
%     sol = acc.solve('Prior', pri);
%     lsqvce.util.log.msg('prior absorbed %.1f of %d DoF\n', sol.priorAbsorbed, pri.nk);
%
%   ACCURACY
%     The prior's normals are the stacked rows (K/sigma) x ~ x0/sigma: the
%     Schur-reduced solve with them equals the dense least-squares solve of
%     the augmented system to 1e-8 relative (test_vceMath
%     (testPriorEqualsAugmentedLeastSquares)). The absorbed DoF equals the
%     dense trace trace((A'A + Nk)^-1 K'K) / sigma^2 to 1e-10 for K = I at
%     sigma 0.01 and 100, where it reads > 0.9 nc and < 0.1 nc (test_vceMath
%     (testPriorAbsorbedDoFBounds); before WP-140 that test pinned only the
%     identity rK = nk - absorbed), and lies in [0, nk] to 1e-8 for an
%     identity, Kaula, first-difference and subset K (test_vceMath
%     (testPriorIsBasisAgnostic)). WP-139, WP-140.
%
%   ERRORS
%     lsqvce:prior:badTarget        the prior target has the wrong length
%     lsqvce:prior:badGroup         group is not nk positive integers
%     lsqvce:prior:zeroRow          a row of K is all zero (it constrains nothing)
%
%   See also LSQVCE.ENGINE.BLOCKACCUMULATOR, LSQVCE.ENGINE.VCE, LSQVCE.ENGINE.APPLYREGULARISATION
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-140: the absorbed DoF of the
%   accuracy statement is pinned against the dense trace trace((A'A + Nk)^-1
%   K'K) / sigma^2 (its test had checked only the identity rK = nk -
%   absorbed). With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-139: the help states its
%   accuracy, every number traced to the test that pins it. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-133: an all-zero row of K is
%   refused (lsqvce:prior:zeroRow). With assistance from Claude (Anthropic),
%   model Claude Opus 5.5.
%
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.
%   2026-10-12 (CEST, UTC+02:00), R7 WP-119: the per-row group (the datum test's
%   unit). With assistance from Claude (Anthropic), model Claude Opus 5.5.

    properties (SetAccess = immutable)
        K            % nk x nc  constraint matrix
        x0           % nk x 1   target (default zero)
        sigma        % scalar   prior standard deviation
        nk           % number of pseudo-observations
        nc           % number of common parameters
        group        % nk x 1   declared group of each row (R7 WP-119)
    end

    methods
        function o = CommonPrior(K, x0, sigma, group)
            narginchk(1,4);
            if nargin < 2 || isempty(x0),    x0    = [];  end
            if nargin < 3 || isempty(sigma), sigma = 1.0; end

            validateattributes(K, {'double'}, {'2d','nonempty','finite'}, ...
                               'lsqvce.engine.CommonPrior', 'K', 1);
            validateattributes(sigma, {'double'}, {'scalar','positive','finite'}, ...
                               'lsqvce.engine.CommonPrior', 'sigma', 3);

            if any(all(K == 0, 2))
                % WP-133 (review 3): a zero row constrains nothing; it passed the
                % blind test as a datum and on the QR backend its re-weighting
                % c = |S| / |Kd| = Inf made Cxx all NaN
                error('lsqvce:prior:zeroRow', ...
                    'Row(s) %s of K are all zero: a prior row must constrain a parameter.', ...
                    mat2str(find(all(K == 0, 2)).'));
            end
            o.K  = K;
            o.nk = size(K,1);
            o.nc = size(K,2);

            if isempty(x0)
                o.x0 = zeros(o.nk, 1);
            else
                x0 = x0(:);
                if numel(x0) ~= o.nk
                    error('lsqvce:prior:badTarget', ...
                        ['x0 has %d elements but K has %d rows. x0 lives in the ' ...
                         'CONSTRAINT space (K*x ~ x0), not the parameter space.'], ...
                         numel(x0), o.nk);
                end
                o.x0 = x0;
            end
            o.sigma = sigma;
            if nargin < 4 || isempty(group)
                o.group = ones(o.nk, 1);
            else
                group = double(group(:));
                if numel(group) ~= o.nk || any(group < 1 | group ~= round(group))
                    error('lsqvce:prior:badGroup', ...
                        'group must be %d positive integers (one per row of K), got %d.', ...
                        o.nk, numel(group));
                end
                o.group = group;
            end
        end

        %------------------------------------------------------------------
        function [N, n] = normals(o)
        %NORMALS  Normal-equation contribution of the prior.
            w  = 1/o.sigma;
            Kw = o.K * w;
            N  = Kw.' * Kw;
            n  = Kw.' * (o.x0 * w);
            N  = (N + N.')/2;
        end

        %------------------------------------------------------------------
        function [rK, absorbed] = redundancy(o, Lc)
        %REDUNDANCY  Redundancy of the prior, and the DoF it absorbs.
        %
        %   rK       = nk - trace(Sxx*K'K)/sigma^2
        %   absorbed = nk - rK   -- the effective number of constrained params
        %
        %   Lc must be chol(N_total, 'lower') INCLUDING the prior's own
        %   contribution, since Sxx is the inverse of the regularised matrix.
            Y        = Lc \ o.K.';                 % nc x nk
            absorbed = sum(Y(:).^2) / o.sigma^2;
            rK       = o.nk - absorbed;
        end

        %------------------------------------------------------------------
        function s = vtv(o, x)
        %VTV  Weighted square sum of the prior's own residual.
            v = (o.K * x(:) - o.x0) / o.sigma;
            s = v.' * v;
        end

        %------------------------------------------------------------------
        function o2 = withSigma(o, sigma)
        %WITHSIGMA  A copy with a new sigma (for the VCE iteration).
            o2 = lsqvce.engine.CommonPrior(o.K, o.x0, sigma, o.group);
        end
    end
end
