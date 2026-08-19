function [t, info] = traceEstimate(M, varargin)
%LSQVCE.TRACEESTIMATE  Stochastic (Hutchinson) trace estimator, Ellmer 2018 6.5.6.
%
%   t = LSQVCE.TRACEESTIMATE(M)
%   t = LSQVCE.TRACEESTIMATE(M, 'Name', value, ...)
%   [t, info] = LSQVCE.TRACEESTIMATE(...)
%
%   Estimates trace(M) WITHOUT forming or touching all of M, by probing it with
%   random vectors. This is the estimator specified in Ellmer (2018, TU Graz,
%   section 6.5.6, eqs. 6.5.39-6.5.43) for the redundancy trace tr(N^-1 N_g) in
%   GRACE variance-component estimation, where forming N^-1 explicitly is
%   prohibitive. It is the ACV-domain workhorse behind the empirical stochastic
%   model: the same estimator computes tr(R Vn) for every lag and the
%   geophysical redundancy, reusing one probe set.
%
%   THE ESTIMATOR (Hutchinson 1990; Ellmer 2018 eqs. 6.5.39-6.5.43)
%     For a P x P matrix X and Z random probe vectors z (entries +1/-1, each
%     with probability 1/2 -- RADEMACHER, eq. 6.5.40),
%
%        trace(X) ~ (1/Z) sum_{z=1..Z} z' X z                        (6.5.39)
%
%     which, stacking the probes as columns of a P x Z matrix and normalising
%     Zbar = Z / sqrt(Z), is the "evil maths trick" of eqs. 6.5.41-6.5.43:
%
%        trace(X) ~ trace( Zbar' X Zbar ).                           (6.5.43)
%
%     The diagonal of Zbar'*X*Zbar carries the per-probe quadratic forms; the
%     off-diagonal cross terms z_i' X z_j (i /= j) are the estimator's noise and
%     vanish in expectation because E[z_i z_j'] = delta_ij I for independent
%     Rademacher probes. The estimator is UNBIASED: E[t] = trace(M) exactly, for
%     any Z >= 1; Z only controls the variance (~ 1/Z).
%
%   MATRIX-FREE FORM (the reason this exists)
%     Ellmer never forms X = U^-T N_g U^-1. He forms Ubar = U^-1 Zbar ONCE (U is
%     the global normal-matrix Cholesky, shared by every arc), then each trace is
%     trace(Zbar' U^-T N_g U^-1 Zbar) = trace(Ubar' N_g Ubar) = sum(sum(Ubar .*
%     (N_g Ubar))). Pass M as a FUNCTION HANDLE @(Zbar)->M*Zbar to exploit this:
%     the estimator then only needs the action of M on the P x Z probe block, not
%     M itself. Supplying the shared Ubar via 'probes' is what makes the whole
%     redundancy sweep O(Z) solves instead of O(P^3).
%
%   INPUT
%     M   P x P double  (symmetric or not), OR
%         function_handle  @(Z) -> M*Z, the ACTION of M on a P x Z probe block
%                          (matrix-free). If a handle, give P via 'size'.
%
%   OPTIONS (Name/Value)
%     'probes'   P x Z double   a FIXED probe block to use (already the RAW +1/-1
%                               Rademacher matrix, NOT normalised -- normalisation
%                               is applied here). Overrides 'numProbes'/'size'/
%                               'seed'. This is how a caller reuses one probe set
%                               (or a pre-solved U^-1*Zbar) across many traces.
%     'numProbes' 1x1 int       Z, number of probes (default 100, Ellmer's value).
%     'size'      1x1 int       P, required when M is a function handle (ignored
%                               when M is an explicit matrix or 'probes' is given).
%     'seed'      1x1 int|[]    RNG seed for reproducible probes (default []: use
%                               the global stream). Uses a private RandStream so
%                               it never disturbs the global RNG state.
%     'distribution' char       'rademacher' (default) | 'gaussian'. Ellmer uses
%                               rademacher; gaussian is offered for comparison.
%
%   OUTPUT
%     t    1x1 double   the trace estimate (unbiased).
%     info 1x1 struct
%          .t          the estimate (same as t)
%          .se         Monte-Carlo standard error of t, sqrt(var(q)/Z), where
%                      q_z = z'Mz are the per-probe quadratic forms. Use this to
%                      judge convergence: the true trace lies within ~2*se of t
%                      with ~95% confidence. NaN if Z < 2.
%          .relse      se / |t|  (relative standard error), NaN if t == 0.
%          .perProbe   1 x Z     the per-probe quadratic forms z'Mz (UN-normalised
%                                by Z), i.e. the summands of eq. 6.5.39.
%          .Z          the number of probes actually used.
%          .P          the matrix dimension.
%          .distribution  the probe distribution used.
%
%   COST
%     Explicit M:      one P x Z matrix product, O(P^2 Z).
%     Matrix-free M:   Z applications of the operator; the whole point is that
%                      the operator (e.g. a triangular solve reusing a shared
%                      factor) is far cheaper than forming M.
%
%   ACCURACY, AND WHEN NOT TO TRUST IT
%     The relative error scales like sqrt(Var/Z), and Var depends on the matrix:
%     for a near-rank-one operator (||M||_F / trace(M) ~ 1) the single-probe
%     variance is large and Z=100 can leave tens of percent of error. The
%     toolbox's EXACT leverage path (LSQVCE.LEVERAGE / LSQVCE.REDUNDANCY) was
%     kept as the default for exactly this reason -- see the note in
%     LSQVCE.LEVERAGE. Ellmer's construction probes the BETTER-conditioned
%     operator U^-T N_g U^-1 (redundancies in [0,1]) and averages over many arcs
%     and iterations, which is why Z=100 suffices there. Check .relse before
%     trusting a single estimate, and prefer the exact path when P is small
%     enough that it is affordable.
%
%   EXAMPLE 1 -- explicit matrix, check against the exact trace
%     A  = randn(500); M = A'*A;
%     te = lsqvce.traceEstimate(M, 'numProbes', 200, 'seed', 1);
%     fprintf('exact %.3f  est %.3f\n', trace(M), te);
%
%   EXAMPLE 2 -- matrix-free, reuse a shared factor (Ellmer eq. 6.5.48)
%     U    = chol(N);                       % N = U'*U, computed once
%     raw  = sign(randn(size(N,1), 100));   % raw Rademacher probes
%     Ubar = U \ raw;                        % ONE triangular solve = U^-1*raw,
%                                            % shared across ALL groups
%     for g = 1:nGroup                      % tr(N^-1 N_g) = tr(Ubar' N_g Ubar),
%         t(g) = lsqvce.traceEstimate(@(Z) Ng{g}*Z, 'probes', Ubar, 'size', size(N,1));
%     end
%     % (U^-1 U^-T = N^-1, so tr(U^-T N_g U^-1) = tr(N^-1 N_g). The estimator
%     %  normalises the supplied probes by 1/sqrt(Z) internally.)
%
%   REFERENCES
%     Hutchinson, M.F. (1990). A stochastic estimator of the trace of the
%       influence matrix for Laplacian smoothing splines. Communications in
%       Statistics - Simulation and Computation, 19(2), 433-450.
%     Ellmer, M. (2018). Contributions to GRACE Gravity Field Recovery.
%       Doctoral thesis, TU Graz. doi:10.3217/978-3-85125-646-8 (section 6.5.6).
%
%   ERRORS
%     lsqvce:traceEstimate:badOperator  the operator handle is not callable
%     lsqvce:traceEstimate:needSize     the operator is a handle, so its size
%                                       must be given explicitly
%     lsqvce:traceEstimate:notSquare    the operator is not square
%     lsqvce:traceEstimate:sizeMismatch the probe dimension disagrees
%
%   See also LSQVCE.REDUNDANCY, LSQVCE.VCE, LSQVCE.LEVERAGE, LSQVCE.DECORRFACTOR
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-36 (fixes F-06): probes from
%   the ONE private generator probeMatrix -- empty seed draws from the
%   package-private session stream (global RNG untouched, as the help
%   always promised); the Rademacher recipe is unified to sign(randn), so
%   seeded realisations differ from pre-WP-36 (unbiasedness unaffected).
%   With assistance from Claude (Anthropic), model Claude Fable 5.
%   Created 2026-07-16 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.

    narginchk(1, inf);

    p = inputParser;  p.FunctionName = 'lsqvce.traceEstimate';
    p.addParameter('probes',       [],   @(v) isnumeric(v) && ismatrix(v));
    p.addParameter('numProbes',    100,  @(v) isscalar(v) && v >= 1 && v == floor(v));
    p.addParameter('size',         [],   @(v) isempty(v) || (isscalar(v) && v >= 1));
    p.addParameter('seed',         [],   @(v) isempty(v) || (isscalar(v) && isnumeric(v)));
    p.addParameter('distribution', 'rademacher', ...
                   @(s) any(strcmpi(s, {'rademacher','gaussian'})));
    p.parse(varargin{:});
    o = p.Results;
    dist = lower(o.distribution);

    isHandle = isa(M, 'function_handle');
    if ~isHandle
        validateattributes(M, {'numeric'}, {'2d'}, mfilename, 'M', 1);
        if size(M,1) ~= size(M,2)
            error('lsqvce:traceEstimate:notSquare', ...
                  'M must be square (got %d x %d).', size(M,1), size(M,2));
        end
    end

    % ---- determine P (matrix dimension) ----------------------------------
    if ~isempty(o.probes)
        P = size(o.probes, 1);
    elseif ~isHandle
        P = size(M, 1);
    elseif ~isempty(o.size)
        P = o.size;
    else
        error('lsqvce:traceEstimate:needSize', ...
              ['When M is a function handle and no ''Probes'' are supplied, ' ...
               'give the matrix dimension via ''Size''.']);
    end
    if isHandle && ~isempty(o.size) && o.size ~= P
        error('lsqvce:traceEstimate:sizeMismatch', ...
              '''Size'' (%d) disagrees with the ''Probes'' row count (%d).', ...
              o.size, P);
    end

    % ---- build (or accept) the RAW probe block, then normalise -----------
    %   Ellmer eq. 6.5.40: entries +1/-1, prob 1/2 (Rademacher).
    %   eq. 6.5.43: Zbar = Z / sqrt(Z) so that trace(Zbar' X Zbar) estimates
    %   trace(X) directly (the 1/Z of eq. 6.5.39 is folded into Zbar).
    if ~isempty(o.probes)
        rawZ = o.probes;
        Z    = size(rawZ, 2);
    else
        Z = o.numProbes;
        % WP-36 (fixes F-06): the ONE generator. Before the fix, an empty
        % seed drew from -- and ADVANCED -- the global stream (the comment
        % claimed otherwise), and this function used the 2*randi-1
        % Rademacher recipe while the other callers used sign(randn): a
        % SEEDED trace therefore reproduces here but its realisation
        % differs from pre-WP-36 values (documented; the estimator is
        % unbiased under either recipe).
        rawZ = probeMatrix(P, Z, dist, o.seed);
    end
    Zbar = rawZ / sqrt(Z);

    % ---- apply the operator: Y = M * Zbar --------------------------------
    if isHandle
        Y = M(Zbar);
        if ~isequal(size(Y), size(Zbar))
            error('lsqvce:traceEstimate:badOperator', ...
                  ['The operator returned a %s block; expected %d x %d ' ...
                   '(= size of the probe block).'], ...
                   mat2str(size(Y)), size(Zbar,1), size(Zbar,2));
        end
    else
        Y = M * Zbar;
    end

    % ---- the estimate: trace(Zbar' M Zbar) = sum of the per-probe forms ---
    %   perProbeBar(z) = Zbar(:,z)' * Y(:,z);  their SUM is trace(Zbar'*M*Zbar).
    perProbeBar = sum(Zbar .* Y, 1);      % 1 x Z, each already divided by Z
    t = sum(perProbeBar);

    % Un-normalised per-probe quadratic forms q_z = z' M z (eq. 6.5.39 summands),
    % recovered as perProbeBar * Z. Their sample scatter gives the MC error.
    perProbe = perProbeBar * Z;

    if nargout > 1
        info = struct();
        info.t        = t;
        info.perProbe = perProbe;
        info.Z        = Z;
        info.P        = P;
        info.distribution = dist;
        if Z >= 2
            % t = mean(perProbe); se of the mean = std(perProbe)/sqrt(Z).
            info.se = std(perProbe) / sqrt(Z);
        else
            info.se = NaN;
        end
        if t ~= 0
            info.relse = info.se / abs(t);
        else
            info.relse = NaN;
        end
    end
end
