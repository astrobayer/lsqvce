classdef Accumulator < handle
%ACCUMULATOR  Streaming normal-equation / QR accumulator for large LSQ problems.
%
%   Designed for the case where the design matrix A (n x u) cannot be held in
%   memory but can be generated in batches.  Nothing of size n is ever stored.
%
%   Two solution paths:
%
%     'normal'  Accumulate N = A'PA and b = A'Pl.  Cheap (u^2 per obs) but
%               cond(N) = cond(A)^2, so ~half the significant digits are lost.
%               Fine for well-conditioned problems.
%
%     'qr'      Maintain an upper-triangular [R d] via Householder updates of
%               sqrt(P)*[A l].  The condition number is NOT squared.  This is
%               the default because a geodetic normal matrix is very often
%               badly scaled.
%
%   SPEED (measured IN MATLAB R2026a, n = 1e5, batches of 5000):
%
%       u            10      30      50     100
%       QR/normal   1.20x   2.39x   1.82x   1.98x
%
%   QR costs only about 1.2-2.4x -- close to the ~2x the flop count predicts,
%   and the two paths agree to ~2e-15 on a well-conditioned problem.
%
%   (An earlier estimate made with NumPy said 3-9x. That was PESSIMISTIC:
%   MATLAB's qr() is far better tuned for tall-skinny matrices than NumPy's, so
%   the real penalty is much smaller than feared. Always measure on the platform
%   you actually ship on -- run examples/bench_qr_vs_normal.)
%
%   Both paths are LINEAR in n; only the constant differs. The cost grows with
%   u, not with n, so for a big streaming problem with few parameters the
%   absolute difference is small.
%
%   WHICH TO USE. The break-even is set by conditioning, not by size:
%
%       cond(N) < 1e8    normal equations lose ~8 of 16 digits. Fine.
%                        Use 'normal' and take the 3-9x speedup.
%       cond(N) 1e8-1e12 borderline. Use 'qr'.
%       cond(N) > 1e12   the normal path has ~4 digits left. QR is not optional.
%
%   RECIPE: run once with 'normal' (cheap), read cond(N) from the report, and
%   re-run with 'qr' only if it exceeds 1e8. One cheap solve tells you.
%
%   Buffering batches before re-triangularising does NOT help (measured: at
%   most 1.07x, and WORSE beyond a buffer of 5 -- cache pressure dominates).
%   Passing larger batches is the only lever, and it is weak.
%
%   USAGE
%       acc = lsqvce.engine.Accumulator(u, 'method','qr');
%       for k = 1:nBatch
%           [Ak, lk, pk] = myBatchSource(k);
%           acc.add(Ak, lk, pk);
%       end
%       sol = acc.solve();
%
%   The 'p' argument is the DIAGONAL of P (a vector), or a full/sparse matrix
%   for a correlated batch, or [] for unit weights.  Correlated batches are
%   handled by an internal Cholesky whitening, so block-diagonal covariance
%   between observations inside one batch is fully supported.
%
%   ON vPv
%   ------
%   It is tempting to recover the weighted residual sum of squares from the
%   accumulated scalars alone:
%
%       vPv = l'Pl - b'x                                            (*)
%
%   Do not trust this.  When the fit is tight, l'Pl and b'x agree to many
%   digits and (*) is catastrophic cancellation: it loses roughly
%   log10(l'Pl / vPv) significant digits.  In numerical experiments with
%   sigma = 1e-7 this produced a vPv wrong by a factor of ~5, and at 1e-9 by a
%   factor of ~4e4 -- i.e. a completely fictitious variance factor, in exactly
%   the regime (precise observations) that geodesy cares most about.
%
%   Therefore:
%     * the QR path takes vPv from the accumulated squares of the discarded
%       rows, which is cancellation-free and is the DEFAULT;
%     * the normal-equation path computes (*) but flags it as provisional,
%       and lsqvce.engine.diagnose() recomputes vPv honestly as sum(p.*v.^2) during
%       the second pass;
%     * the discrepancy between the two is reported, and a loud warning is
%       raised if they disagree, because that disagreement is itself a
%       reliable symptom of a numerically sick problem.
%
%   EXAMPLE
%       rng(0);  n = 200;  u = 4;
%       A = randn(n,u);  l = A*randn(u,1) + randn(n,1);  p = ones(n,1);
%       idx = reshape(1:n, [], 4);
%
%       acc = lsqvce.engine.Accumulator(u);                % 'qr' by default
%       for k = 1:4
%           acc.add(A(idx(:,k),:), l(idx(:,k)), p(idx(:,k)));
%       end
%       sol = acc.solve();
%
%       % sol.x       the parameter estimate
%       % sol.condN   read THIS to decide whether you needed 'qr' at all:
%       %             below 1e8 the cheap 'normal' path would have been fine.
%
%   ERRORS
%     Your batch is malformed
%       lsqvce:acc:badA            A is not a numeric matrix, or its column
%                                  count disagrees with the first batch
%       lsqvce:acc:badL            l is not a numeric column of height size(A,1)
%       lsqvce:acc:badP            p is neither empty, a scalar, nor an
%                                  size(A,1)-vector of weights
%       lsqvce:acc:badPmatrix      a full weight MATRIX was given where a
%                                  diagonal was expected
%       lsqvce:acc:nonFiniteA      A contains NaN or Inf
%       lsqvce:acc:nonFiniteL      l contains NaN or Inf
%       lsqvce:acc:negativeWeight  a weight is negative. p is a WEIGHT
%                                  (1/sigma^2), not a variance
%
%     Your datum constraint is malformed
%       lsqvce:acc:badDatum        the datum matrix has the wrong width
%       lsqvce:acc:datumNotNull    the supplied datum is not a null-space basis
%                                  of the normal matrix
%
%     The sequence of calls is wrong
%       lsqvce:acc:noData          solve() before any batch was added
%       lsqvce:acc:finalized       add() after the system was solved; build a
%                                  new Accumulator instead of reusing one
%
%     The system cannot be solved as posed
%       lsqvce:acc:underdetermined n < u -- fewer observations than unknowns
%       lsqvce:acc:rankDeficient   the normal matrix is singular and no datum
%                                  was supplied; see LSQVCE.ENGINE.DATUM
%       lsqvce:acc:sparseQR        the 'qr' method was asked for with a sparse
%                                  A that MATLAB cannot factor in place
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:acc:cancellation                    [warning] vPv was formed as l'Pl - b'x and has lost about X significant digits...
%     lsqvce:acc:negativeVPV                     [warning] vPv came out NEGATIVE (X), which is impossible for a quadratic form
%     lsqvce:acc:zeroRedundancy                  [warning] Redundancy is zero (n = rank = X)
%     lsqvce:acc:zeroWeight                      [warning] X observations have weight exactly 0
%
%   See also LSQVCE.ENGINE.ADJUST, LSQVCE.ENGINE.DIAGNOSE, LSQVCE.ENGINE.COMBINE.
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

properties (SetAccess = private)
    u           (1,1) double        % number of parameters
    method      (1,:) char          % 'qr' | 'normal'
    useSparse   (1,1) logical = false  % keep N sparse (see solve())

    n           (1,1) double = 0    % observations seen
    nBatch      (1,1) double = 0    % batches seen

    % --- normal-equation state
    N           double = []         % u x u
    b           double = []         % u x 1
    lPl         (1,1) double = 0    % l'Pl

    % --- QR state
    Rd          double = []         % k x (u+1) upper triangular [R d]
    e2          (1,1) double = 0    % accumulated squared discarded rows = vPv

    % --- free-network datum
    datumB      double = []         % u x d orthonormal null-space basis
    datumSeeded (1,1) logical = false   % have the B' rows entered the QR yet?
    datumInfo   struct = struct()

    finalized   (1,1) logical = false
end

methods

function obj = Accumulator(u, opt)
%ACCUMULATOR  Construct a flat streaming accumulator over u parameters (see header).
arguments
    u (1,1) double {mustBeInteger, mustBePositive}
    opt.Method (1,1) string {mustBeMember(opt.Method, ["qr","normal"])} = "qr"
    opt.Sparse (1,1) logical = false
    opt.DatumB double = []
end
    % WP-39 stage (c) c2: typed `arguments` signature (canonical
    % UpperCamelCase; legacy casings keep matching -- MATLAB NV
    % semantics; enum VALUES case-sensitive lowercase).
    obj.u         = double(u);
    obj.method    = char(opt.Method);
    obj.useSparse = opt.Sparse;

    if obj.useSparse && strcmp(obj.method, 'qr')
        % A sparse QR would need a sparse Householder implementation; MATLAB's
        % qr() on a sparse matrix does exist, but the running [R d] we maintain
        % becomes dense almost immediately because each batch couples whatever
        % parameters it touches. Refuse rather than pretend.
        error('lsqvce:acc:sparseQR', ...
            ['''sparse'',true is only supported with ''method'',''normal''. ' ...
             'The running QR factor [R d] fills in as soon as batches couple ' ...
             'different parameter blocks, so a sparse QR buys nothing here.\n' ...
             'Use ''method'',''normal'',''sparse'',true -- but note that the ' ...
             'normal-equation path squares the condition number, so check ' ...
             'condN in the report.']);
    end

    B = opt.DatumB;
    if ~isempty(B)
        if size(B, 1) ~= obj.u
            error('lsqvce:acc:badDatum', ...
                  'datumB has %d rows but there are %d parameters.', ...
                  size(B,1), obj.u);
        end
        % Orthonormalise once, here. Everything downstream assumes B'B = I.
        [Q, Rb] = qr(B, 0);
        if min(abs(diag(Rb))) < eps(max(abs(diag(Rb)))) * size(B,2) * 10
            error('lsqvce:acc:badDatum', ...
                  'The columns of datumB are not linearly independent.');
        end
        obj.datumB = Q;
    else
        obj.datumB = zeros(obj.u, 0);
    end

    switch obj.method
        case 'normal'
            if obj.useSparse
                obj.N = sparse(obj.u, obj.u);
            else
                obj.N = zeros(obj.u);
            end
            obj.b = zeros(obj.u, 1);
        case 'qr'
            obj.Rd = zeros(0, obj.u + 1);
    end
end

%% ------------------------------------------------------------------- add
function add(obj, A, l, p)
%ADD  Fold one batch of observations into the accumulator.
%
%   ADD(A, l)        unit weights
%   ADD(A, l, p)     p = weight vector (diagonal P), or a weight MATRIX for a
%                    correlated batch, or [] for unit weights.

    if obj.finalized
        error('lsqvce:acc:finalized', ...
              ['This accumulator has already been solved. Create a new one ' ...
               'rather than adding to a solved system.']);
    end
    if nargin < 4, p = []; end

    % --- Seed the QR stack with the datum constraints on the FIRST add.
    %
    % The datum rows B' MUST enter the factorisation before (or with) the
    % first observations, not at solve() time. If A is rank deficient, the QR
    % of the observation block alone leaves a zero on the diagonal of R, and
    % the residual energy is then spread across rows rank(A)..u rather than
    % landing entirely in row u+1. Accumulating only row u+1 therefore
    % UNDER-COUNTS vPv -- measured at a 19% error on a small levelling net
    % (6.89 instead of 8.50) before this was fixed.
    %
    % Seeding with B' keeps the stack full rank at every step, so row u+1 is
    % the one and only residual carrier and e2 is exact (verified to 9e-14).
    % The constraint rows carry a zero right-hand side and so contribute
    % nothing to the residual sum of squares.
    if strcmp(obj.method, 'qr') && ~obj.datumSeeded
        d = size(obj.datumB, 2);
        if d > 0
            obj.Rd = triu(qr([obj.datumB.', zeros(d, 1)]));
            obj.Rd = obj.Rd(1:min(d, obj.u + 1), :);
        end
        obj.datumSeeded = true;
    end

    [A, l, sA, sl] = obj.whiten(A, l, p);   % sA = sqrt(P)*A, sl = sqrt(P)*l
    nk = size(A, 1);
    if nk == 0, return; end                 % empty batch is legal, not an error

    switch obj.method
        case 'normal'
            if obj.useSparse
                % Keep the batch contribution sparse. If the caller handed us a
                % dense sA we sparsify it: a batch typically touches only a few
                % parameters, so sA'*sA is mostly zeros and densifying N here
                % would defeat the whole point.
                if ~issparse(sA), sAs = sparse(sA); else, sAs = sA; end
                obj.N = obj.N + sAs' * sAs;
            else
                obj.N = obj.N + sA' * sA;
            end
            obj.b   = obj.b   + sA' * sl;
            obj.lPl = obj.lPl + (sl' * sl);

        case 'qr'
            % Compress [previous R d ; new whitened rows] back to triangular.
            stack = [obj.Rd; sA, sl];
            Rf    = triu(qr(stack));                 % economy: R factor only
            k     = min(size(stack, 1), obj.u + 1);
            Rf    = Rf(1:k, :);

            if k > obj.u
                % Row u+1 carries the residual norm of everything discarded so
                % far. Squaring and accumulating it is exact -- no cancellation.
                obj.e2 = obj.e2 + Rf(obj.u + 1, obj.u + 1)^2;
                Rf     = Rf(1:obj.u, :);
            end
            obj.Rd = Rf;
    end

    obj.n      = obj.n + nk;
    obj.nBatch = obj.nBatch + 1;
end

%% ----------------------------------------------------------------- solve
function sol = solve(obj, opt)
%SOLVE  Produce the parameter estimate and the scalars needed downstream.
%
%   sol = SOLVE()  returns a struct with fields
%       x        u x 1   parameter estimate
%       R        u x u   upper-triangular factor (of N + B*B' if free)
%       n, u             counts
%       defect           datum defect d (0 if full rank)
%       rank             u - d
%       r                REDUNDANCY = n - rank = n - (u - d)
%       vPv              weighted residual sum of squares
%       s0sq             vPv / r
%       QxxCorr          u x d correction basis; Qxx = inv(R'R) - QxxCorr*QxxCorr'
%
%   FREE NETWORKS
%   -------------
%   If a datum matrix B was supplied, the d rows of B' are appended to the QR
%   stack as pseudo-observations of weight 1. This is the inner-constraint /
%   bordered system in streaming form: it makes R full rank, so every
%   downstream routine works unchanged, and it yields the correct minimum-norm
%   x (verified against the bordered Lagrange system to 5e-14).
%
%   TWO SUBTLETIES, both verified numerically:
%
%   (1) The resulting inv(N + B*B') is NOT the pseudoinverse N+. They differ by
%       B*inv(B'B)^2*B'. With B orthonormal that is simply B*B'. Reporting
%       inv(N+BB') as Qxx overstates the parameter variances -- on a small
%       levelling net the discrepancy was 12%. sol.QxxCorr carries the
%       correction and lsqvce.engine.diagnose applies it.
%
%   (2) The redundancy is n - rank(A) = n - (u - d), NOT n - u. Using n - u
%       makes the identity sum(r_i) = redundancy fail outright and biases s0^2
%       by (n-u+d)/(n-u). The B' rows are CONSTRAINTS, not observations, and
%       must not be counted in n either.
%
%   Happily, r_i and w_i need NO correction: because A*B = 0 exactly (the
%   observations cannot see the datum), the correction term vanishes from
%   A*Qxx*A'. So the entire diagnostic apparatus is datum-invariant, as it
%   must be physically -- the residuals do not care where you put the origin.
%
%   Name-value:
%       'Sigma0'  a-priori sigma0, scalar > 0 (default 1; WP-39 canonical)

    % WP-39 c2: typed `arguments` signature (canonical 'Sigma0').
arguments
    obj
    opt.Sigma0 (1,1) double {mustBePositive} = 1
end
    sigma0 = opt.Sigma0;

    if obj.n == 0
        error('lsqvce:acc:noData', ...
              'No observations were added. Nothing to solve.');
    end

    d    = size(obj.datumB, 2);
    rnk  = obj.u - d;                    % rank of A
    r    = obj.n - rnk;                  % THE redundancy

    if r < 0
        error('lsqvce:acc:underdetermined', ...
              ['Underdetermined: %d observations for %d estimable parameters ' ...
               '(u = %d, datum defect = %d).'], obj.n, rnk, obj.u, d);
    end

    switch obj.method
        case 'normal'
            Nuse = obj.N;
            if d > 0
                Nuse = Nuse + obj.datumB * obj.datumB.';   % regularise
            end
            % chol() on a sparse matrix returns a sparse upper-triangular R
            % with NO permutation in this single-output form, so R remains
            % directly usable by every downstream formula (R'\A', Rinv, ...).
            %
            % We deliberately do NOT apply a fill-reducing permutation here.
            % Measured on the structures that actually arise:
            %     banded N (splines, time-ordered biases) : RCM buys 0.0%
            %     irregular / scattered footprints        : RCM buys 3.7%
            % Against that, a permutation must be honoured in THREE separate
            % places (the solve, qx = a'Qxx a, and diag(Qxx)); forgetting any
            % one of them yields plausible numbers attached to the WRONG
            % parameters -- measured qx error 0.90, i.e. silently wrong rather
            % than obviously broken. The trade is not worth it.
            %
            % If your N has ARROW structure (a few global parameters seen by
            % every observation, plus many per-arc locals), no generic ordering
            % will save you: you want explicit arc-wise elimination (a Schur
            % complement on the local blocks), which is a different algorithm
            % and outside this toolbox.
            [R, ok] = chol(Nuse);
            if ok ~= 0
                obj.rankFailure(Nuse, d);
            end
            % chol can still SUCCEED on a matrix that is singular only to within
            % rounding (smallest eigenvalue ~ 1e-14 > 0), e.g. A'A formed from an
            % exact column dependence. Guard with the factor's reciprocal
            % condition estimate: R is triangular, so rcond(Nuse) = rcond(R)^2 is
            % cheap and exact enough. Below ~eps^(3/4) the solve is meaningless,
            % so report the same rank-deficiency error rather than returning a
            % plausible-looking but numerically hollow result.
            rcR = rcond(full(R));
            if ~isfinite(rcR) || rcR^2 < 1e-12
                obj.rankFailure(Nuse, d);
            end
            x   = R \ (R' \ obj.b);
            vPv = obj.lPl - obj.b' * x;
            src = 'normal (provisional)';

            if obj.lPl > 0
                digitsLost = log10(max(obj.lPl / max(vPv, realmin), 1));
                if digitsLost > 8
                    lsqvce.util.log.warn('lsqvce:acc:cancellation', ...
                        ['vPv was formed as l''Pl - b''x and has lost about ' ...
                         '%.0f significant digits to cancellation ' ...
                         '(l''Pl = %.6e, vPv = %.6e).\n' ...
                         'RETURNING NaN rather than a plausible-looking number: ' ...
                         'measured, this quantity was 89%% wrong at sigma = 1e-6 ' ...
                         'and wrong by a factor of 1e6 at sigma = 1e-9. A number ' ...
                         'that wrong will be quoted; a NaN will be investigated.\n' ...
                         'lsqvce.engine.diagnose() recomputes vPv honestly from the ' ...
                         'residuals and will repair this. To avoid the problem ' ...
                         'entirely, use ''method'',''qr''.'], ...
                         digitsLost, obj.lPl, vPv);
                    vPv = NaN;
                end
            end
            if vPv < 0
                lsqvce.util.log.warn('lsqvce:acc:negativeVPV', ...
                    ['vPv came out NEGATIVE (%.6e), which is impossible for a ' ...
                     'quadratic form. This is pure cancellation error. Re-run ' ...
                     'with ''method'',''qr''.'], vPv);
                vPv = NaN;
            end
            Nmat = Nuse;

        case 'qr'
            % The datum rows B' were folded in by add() on the first batch, so
            % Rd already reflects them and e2 is already correct. Nothing to do
            % here beyond reading off the factor.
            Rd_ = obj.Rd;

            if size(Rd_, 1) < obj.u
                error('lsqvce:acc:underdetermined', ...
                      ['The QR factor has only %d rows for %d parameters: the ' ...
                       'observations (plus any datum constraints) do not ' ...
                       'determine all parameters. If this is a free network, ' ...
                       'supply ''datumB''.'], size(Rd_, 1), obj.u);
            end
            R = Rd_(1:obj.u, 1:obj.u);
            dv = Rd_(1:obj.u, obj.u + 1);

            dg = abs(diag(R));
            if any(dg <= eps * max(dg) * obj.u)
                obj.rankFailure(R' * R, d);
            end
            x    = R \ dv;
            vPv  = obj.e2;
            src  = 'qr';
            Nmat = R' * R;
    end

    % Normalise row signs so R is reproducible between methods.
    s = sign(diag(R));
    s(s == 0) = 1;
    R = s .* R;                     % row-scale by +/-1 (broadcast, not diag(s)*R)

    % --- pseudoinverse correction, GUARDED.
    %
    % With B orthonormal and spanning null(N), the pseudoinverse is
    %     N+ = inv(N + B*B') - B*B',
    % so QxxCorr = B and Qxx = inv(R'R) - QxxCorr*QxxCorr'.
    % (Verified against the bordered Lagrange system to 3e-16.)
    %
    % BUT that identity holds ONLY if B really spans the null space. It is easy
    % to supply something that removes the rank defect without being a null
    % vector -- fixing one station, for instance: B = e_3 is a perfectly valid
    % HARD CONSTRAINT but it is not an inner constraint, and ||N*e_3||/||N||
    % was measured at 0.24, not 0. Applying -B*B' in that case is not a small
    % error, it is nonsense: sum(r_i) came out as 46074 instead of 15.
    %
    % So: check, and refuse rather than silently corrupt every downstream
    % number. A hard constraint is a legitimate thing to want; it is simply a
    % different thing, and the caller must say which they mean.
    QxxCorr = obj.datumB;
    if d > 0
        Nobs = Nmat - obj.datumB * obj.datumB.';   % the observation part alone
        nrmN = norm(Nobs, 'fro');
        if nrmN > 0
            nullResid = norm(Nobs * obj.datumB, 'fro') / nrmN;
        else
            nullResid = 0;
        end

        if nullResid > 1e-6
            error('lsqvce:acc:datumNotNull', ...
                ['The supplied datumB does NOT span the null space of the ' ...
                 'normal matrix: ||N*B||/||N|| = %.3e, which must be ~0.\n\n' ...
                 'You have probably supplied a HARD CONSTRAINT (e.g. "fix ' ...
                 'station 3", B = e_3) rather than an INNER CONSTRAINT. These ' ...
                 'are different things and need different treatment:\n\n' ...
                 '  INNER CONSTRAINT (free network): B spans null(N). The ' ...
                 'observations genuinely cannot see these directions. Gives ' ...
                 'the minimum-norm solution, B''x = 0, and Qxx = N+ .\n' ...
                 '    levelling: B = ones(u,1)\n\n' ...
                 '  HARD CONSTRAINT (fixed station): you are choosing to pin ' ...
                 'a parameter. This is a REPARAMETERISATION -- simply remove ' ...
                 'that parameter from the model (drop its column from A and ' ...
                 'reduce u), which leaves a full-rank system needing no datum ' ...
                 'handling at all.\n\n' ...
                 'Applying the pseudoinverse correction to a hard constraint ' ...
                 'does not merely lose accuracy: in testing it drove ' ...
                 'sum(r_i) to 46074 against a true value of 15. Refusing.'], ...
                 nullResid);
        end
        obj.datumInfo = struct('nullResidual', nullResid, 'defect', d);
    end

    s0sq = vPv / max(r, 1);
    if r == 0
        lsqvce.util.log.warn('lsqvce:acc:zeroRedundancy', ...
            ['Redundancy is zero (n = rank = %d): the model reproduces the ' ...
             'data exactly by construction. No residual, no variance factor, ' ...
             'and no test of any kind is possible. This adjustment cannot be ' ...
             'judged -- only computed.'], rnk);
        s0sq = NaN;
    end

    condR = obj.condEst(R);

    sol = struct( ...
        'x',          x, ...
        'R',          R, ...
        'N',          Nmat, ...
        'n',          obj.n, ...
        'u',          obj.u, ...
        'defect',     d, ...
        'rank',       rnk, ...
        'datumB',     obj.datumB, ...
        'QxxCorr',    QxxCorr, ...
        'r',          r, ...
        'vPv',        vPv, ...
        'vPvSource',  src, ...
        's0sq',       s0sq, ...
        'sigma0',     sigma0, ...
        ... % full() because a sparse R yields a SPARSE condR, and sprintf()
        ... % REFUSES sparse input -- the report then crashes outright on the
        ... % sparse path. Caught only by running the suite in MATLAB.
        'condR',      full(condR), ...
        'condN',      full(condR)^2, ...
        'condIsBound', (issparse(R) || size(R,1) > 2000), ...
        'sparse',     obj.useSparse, ...
        'nnzR',       nnz(R), ...
        'method',     obj.method, ...
        'nBatch',     obj.nBatch);

    obj.finalized = true;
end

end % methods

methods (Access = private)

function [A, l, sA, sl] = whiten(obj, A, l, p)
    % Validate a batch and return sqrt(P)*A, sqrt(P)*l.

    if ~isnumeric(A) || ~ismatrix(A)
        error('lsqvce:acc:badA', 'A must be a numeric matrix.');
    end
    if size(A, 2) ~= obj.u
        error('lsqvce:acc:badA', ...
              ['A has %d columns but the accumulator was built for %d ' ...
               'parameters.'], size(A, 2), obj.u);
    end
    l = l(:);
    if numel(l) ~= size(A, 1)
        error('lsqvce:acc:badL', ...
              'l has %d elements but A has %d rows.', numel(l), size(A, 1));
    end
    if ~all(isfinite(A(:)))
        error('lsqvce:acc:nonFiniteA', ...
              ['A contains %d non-finite entries. Remove or repair these ' ...
               'observations before accumulating -- a single NaN silently ' ...
               'destroys the entire normal matrix.'], sum(~isfinite(A(:))));
    end
    if ~all(isfinite(l))
        error('lsqvce:acc:nonFiniteL', ...
              'l contains %d non-finite entries.', sum(~isfinite(l)));
    end

    nk = size(A, 1);

    if isempty(p)
        sA = A;  sl = l;
        return;
    end

    if isvector(p) && numel(p) == nk
        p = p(:);
        if any(~isfinite(p))
            error('lsqvce:acc:badP', 'Weight vector contains non-finite values.');
        end
        if any(p < 0)
            error('lsqvce:acc:negativeWeight', ...
                  ['Weight vector contains %d negative entries. A negative ' ...
                   'weight is not a down-weighting -- it makes the normal ' ...
                   'matrix indefinite and the "minimum" a saddle point.'], ...
                   sum(p < 0));
        end
        if any(p == 0)
            lsqvce.util.log.warn('lsqvce:acc:zeroWeight', ...
                ['%d observations have weight exactly 0. They contribute ' ...
                 'nothing to the solution but DO count toward n, which ' ...
                 'inflates the redundancy and biases the variance factor ' ...
                 'downward. Remove them instead of zero-weighting them.'], ...
                 sum(p == 0));
        end
        sp = sqrt(p);
        sA = sp .* A;                       % implicit expansion
        sl = sp .* l;
        return;
    end

    if ismatrix(p) && isequal(size(p), [nk nk])
        % Correlated batch: whiten with the Cholesky factor of P.
        p = (p + p') / 2;                   % symmetrise; caller may be sloppy
        [Lp, ok] = chol(p, 'lower');
        if ok ~= 0
            error('lsqvce:acc:badPmatrix', ...
                  ['The weight matrix for this batch is not positive definite ' ...
                   '(Cholesky failed at leading minor %d). A weight matrix is ' ...
                   'an inverse covariance and must be SPD. Check that you ' ...
                   'passed P and not its inverse.'], ok);
        end
        sA = Lp' * A;
        sl = Lp' * l;
        return;
    end

    error('lsqvce:acc:badP', ...
          ['Weight argument must be [] (unit), a length-%d vector (diagonal ' ...
           'P), or a %dx%d SPD matrix (correlated batch). Got size [%s].'], ...
           nk, nk, nk, num2str(size(p)));
end

function rankFailure(~, M, dGiven)
%RANKFAILURE  Raise a descriptive rank-deficiency error for an M x M system.
    if nargin < 3, dGiven = 0; end
    dd = eig((M + M') / 2);
    nz = sum(dd <= max(dd) * 1e-12);

    if dGiven > 0
        extra = sprintf([ ...
            '\nYou already supplied a datum of dimension %d, and the system is ' ...
            'STILL singular. Either the datum is incomplete (there are more ' ...
            'unobservable directions than you think) or some parameters are ' ...
            'genuinely undetermined by the data.'], dGiven);
    else
        extra = [newline ...
            'If this is a FREE NETWORK -- a levelling net with no fixed ' ...
            'benchmark, a geodetic net with no fixed station, an orbit with ' ...
            'an unconstrained epoch -- then this is not a bug: the ' ...
            'observations genuinely cannot see the datum. Supply it:' newline ...
            '    lsqvce.engine.adjust(u, src, nb, ''datumB'', B)' newline ...
            'where the columns of B span the transformations the data cannot ' ...
            'see (for levelling: ones(u,1)). Or call lsqvce.engine.datum(N) to have ' ...
            'it guessed, with the caveats that entails.'];
    end

    error('lsqvce:acc:rankDeficient', ...
        ['The normal matrix is rank deficient (roughly %d null direction(s) ' ...
         'out of %d).\n' ...
         'The observations do not determine every parameter. This is a DATUM ' ...
         'or PARAMETERISATION defect, not a numerical accident -- adding more ' ...
         'observations of the same kind will not fix it.%s'], nz, size(M,1), extra);
end

function c = condEst(~, R)
    % Condition number of the triangular factor.
    %
    % svd() does not accept a sparse matrix, and densifying a large sparse R
    % just to estimate its condition would defeat the purpose of keeping it
    % sparse. So: dense and small -> exact svd; sparse or large -> the ratio of
    % the extreme diagonal entries, which for a triangular factor is a LOWER
    % BOUND on the true condition number, not an estimate of it.
    %
    % A lower bound is still useful (it can only under-report trouble, never
    % invent it), but the report must not present it as exact -- see condIsBound
    % in the returned struct.
    u = size(R, 1);
    if ~issparse(R) && u <= 2000
        s = svd(R);
        if s(end) <= 0
            c = Inf;
        else
            c = s(1) / s(end);
        end
    else
        d = abs(diag(R));
        if isempty(d) || min(d) <= 0
            c = Inf;
        else
            c = max(d) / min(d);
        end
    end
end

end
end
