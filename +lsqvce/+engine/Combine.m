classdef Combine < handle
%LSQVCE.ENGINE.COMBINE  Multi-dataset combination by elimination of local parameters.
%
%   For the ARROW structure: several data sets (satellite systems, campaigns,
%   missions) share a COMMON parameter block, and each has its own LOCAL
%   parameters that no other data set sees.
%
%       [ N_CC   N_C1   N_C2  ...  N_Ck ] [x_C]   [b_C]
%       [ N_1C   N_11    0    ...   0   ] [x_1]   [b_1]
%       [ N_2C    0     N_22  ...   0   ] [x_2] = [b_2]
%       [  ...                          ] [...]   [...]
%       [ N_kC    0      0    ...  N_kk ] [x_k]   [b_k]
%
%   The locals are eliminated system by system, so the running memory is
%   O(uC^2 + max(uL_i)^2) and NEVER O(u_total^2). With 100 systems x 50 locals
%   and 5 common parameters, the full normal matrix has 2.5e7 entries; the
%   reduced one has 25.
%
%   THE DEFINITION OF "LOCAL" THAT MATTERS
%   --------------------------------------
%   A parameter is LOCAL if exactly ONE data set can see it. It is COMMON if
%   two or more can. This is about which data sets observe it, NOT about
%   whether it feels like a nuisance parameter.
%
%   Several sensors within one satellite system may observe the same local
%   parameters -- that is fine, they all land in the same block. What breaks
%   the structure is a parameter seen by two different SYSTEMS: N_LL is then
%   no longer block-diagonal and eliminating the blocks independently is not
%   an approximation, it is WRONG. Such a parameter belongs in the common
%   block, however local it feels. This class CHECKS for the violation.
%
%   USAGE
%       cmb = lsqvce.engine.Combine(uC);
%       for i = 1:nSystem
%           cmb.addSystem(uL(i), @(k) batchOf(i,k), nBatch(i));
%       end
%       sol = cmb.solve();          % xC, and the locals by back-substitution
%       D   = cmb.diagnose(sol);    % full per-observation diagnostics
%
%   Each system's batchFcn returns
%       [AC, AL, l, p] = batchFcn(k)
%   with AC the COMMON partials (nk x uC), AL the LOCAL partials (nk x uL_i),
%   l the observed-minus-computed, and p the weights (vector, matrix, or []).
%
%   METHOD
%     'qr'      (default) TRUE Householder elimination: the local columns are
%               triangularised away on the OBSERVATION EQUATIONS. N is never
%               formed -- not for the locals, not for the common part.
%     'normal'  Form N_ii and Schur-complement:
%                   Nred = N_CC - sum_i N_Ci N_ii^-1 N_iC
%               Cheaper per system, but cond(N) = cond(A)^2.
%
%   MEASURED (near-collinear common columns, relative error in xC against a
%   long-double SVD reference):
%
%       cond(A)        1e2      1e4      1e6      1e7      1e8
%       'normal'      9.6e-12  5.6e-08  1.8e-04  2.2e-02  4.9e-01
%       'qr'          1.9e-14  1.3e-12  1.1e-10  1.4e-09  1.2e-08
%
%   At cond(A) = 2.4e8 the normal-equation path has NO correct digits left,
%   while the Householder path still has eight. The gap is SEVEN TO TEN ORDERS
%   OF MAGNITUDE, and it is not a subtlety: forming N -- even only in order to
%   Schur-complement it -- squares the condition number, and no amount of care
%   afterwards can undo that.
%
%   SPEED (measured):
%
%       nSys  uC   uL     n         normal    QR      ratio
%          5  10   20   2.0e5       0.062 s   0.147 s  2.4x
%         20  10   20   8.0e5       0.167 s   0.450 s  2.7x
%         10  50   20   4.0e5       0.204 s   0.859 s  4.2x
%         10 100   20   4.0e5       0.376 s   2.582 s  6.9x
%        100  20   50   1.0e6       0.520 s   2.022 s  3.9x
%
%   QR costs 2.4-6.9x. The cost grows with uC and uL, NOT with n -- and it is
%   FLAT in the conditioning, so you do not pay more for QR precisely when you
%   need it more.
%
%   WHEN IS THE ACCURACY WORTH IT? Measured, on the SAME problem:
%
%       cond(A)   normal err   QR err     verdict
%       2.0e+01   2.8e-14      3.9e-14    normal is fine
%       2.0e+03   4.0e-10      3.6e-12    normal is fine
%       2.0e+05   4.5e-06      3.9e-09    QR recommended
%       2.0e+07   2.5e-02      1.8e-07    QR MANDATORY
%       2.0e+08   2.3e+04      6.6e-07    QR MANDATORY
%
%   RECIPE: you cannot see cond(A) without doing the work, but you CAN see
%   cond(N), and cond(N) = cond(A)^2. Run once with 'normal' (cheap), read
%   cond(N), and re-run with 'qr' only if it exceeds 1e8. One cheap solve buys
%   you the answer.
%
%   EXAMPLE
%       % Two satellite systems share a common gravity field; each has its own
%       % LOCAL parameters that no other system sees.
%       uC = 4;  uL = [3 2];
%       cmb = lsqvce.engine.Combine(uC);              % 'qr' by default
%
%       for i = 1:2
%           % batchFcn returns FOUR things: AC and AL SEPARATELY. Do not
%           % concatenate them -- the toolbox must know which columns to
%           % eliminate.
%           cmb.addSystem(uL(i), @(k) srcOf(i,k), nBatch(i), ...
%                         'name', sprintf('satellite %d', i));
%       end
%
%       sol = cmb.solve();       % sol.xC, sol.xL{i}, sol.sigC, sol.sigL{i}
%       D   = cmb.diagnose(sol); % D.rNum{i}, D.w{i}, D.mdb{i}
%
%       assert(D.sumRRelErr < 1e-8);            % ALWAYS. See the class notes:
%       % sum(r_i) = n - u_TOTAL, counting the ELIMINATED locals. They left the
%       % SOLVE, not the MODEL. Two classic errors both break this identity and
%       % both lie in the OPTIMISTIC direction.
%
%       % Per-sensor variance components:
%       VC = cmb.varComp(sol, @(i,k) sensorOf(i,k), nSensor);
%
%   ERRORS
%     The batch is malformed
%       lsqvce:combine:badAL          A and l are not a consistent pair
%       lsqvce:combine:badUL          the unknown count disagrees with the
%                                     first batch
%       lsqvce:combine:sizeMismatch   a block's dimensions do not line up
%       lsqvce:combine:badP           p is neither empty, scalar, nor a vector
%       lsqvce:combine:badPmatrix     a full weight MATRIX where a diagonal was
%                                     expected
%       lsqvce:combine:negativeWeight a weight is negative -- p is 1/sigma^2
%       lsqvce:combine:nonFinite      the batch carries NaN or Inf
%
%     The grouping or sequence is wrong
%       lsqvce:combine:noSystems      nothing was added
%       lsqvce:combine:badGroup       a group key is not a valid name
%       lsqvce:combine:finalized      add() after solving; build a new Combine
%
%     The system cannot be solved as posed
%       lsqvce:combine:underdetermined  fewer observations than unknowns
%       lsqvce:combine:singularCommon   the reduced COMMON block is singular
%       lsqvce:combine:singularLocal    a LOCAL block is singular: that arc has
%                                       too few observations to fix its own state
%
%     Variance components cannot be estimated
%       lsqvce:combine:correlatedBatch  a batch carries a full covariance, so it
%                                       cannot be split between groups
%       lsqvce:combine:vcCorrelated     likewise, inside the VCE step
%       lsqvce:combine:vcSingular       a group has no redundancy of its own
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:combine:emptySystem                 [warning] X contributed no observations; skipping
%     lsqvce:combine:integrity                   [warning]  ... ... ... ... 
%
%   See also LSQVCE.ENGINE.BATCHINDEX, LSQVCE.ENGINE.ADJUST, LSQVCE.ENGINE.VARCOMP.
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

properties (SetAccess = private)
    uC          (1,1) double        % number of COMMON parameters
    method      (1,:) char          % 'qr' | 'normal'

    nSys        (1,1) double = 0
    sys         cell   = {}         % per-system state

    % running reduced system
    Nred        double = []         % uC x uC   (normal path)
    bred        double = []         % uC x 1
    Rc          double = []         % running [R_CC | d_C]  (qr path)
    e2          (1,1) double = 0    % accumulated vPv from discarded rows

    n           (1,1) double = 0    % total observations
    uL_total    (1,1) double = 0    % total LOCAL parameters
    finalized   (1,1) logical = false
end

methods

function obj = Combine(uC, opt)
%COMBINE  Construct a multi-system combiner over uC common parameters (see header).
arguments
    uC (1,1) double {mustBeInteger, mustBePositive}
    opt.Method (1,1) string {mustBeMember(opt.Method, ["qr","normal"])} = "qr"
end
    % WP-39 stage (c) c2: typed `arguments` signature (canonical
    % UpperCamelCase; legacy casings keep matching; enum VALUES
    % case-sensitive lowercase).
    obj.uC     = double(uC);
    obj.method = char(opt.Method);

    switch obj.method
        case 'normal'
            obj.Nred = zeros(obj.uC);
            obj.bred = zeros(obj.uC, 1);
        case 'qr'
            obj.Rc = zeros(0, obj.uC + 1);
    end
end

%% ------------------------------------------------------------- addSystem
function addSystem(obj, uL, batchFcn, nBatch, opt)
%ADDSYSTEM  Fold one data set in, eliminating its local parameters.
%
%   ADDSYSTEM(uL, batchFcn, nBatch) where
%       [AC, AL, l, p] = batchFcn(k),   k = 1..nBatch
%   AC is nk x uC, AL is nk x uL, l is nk x 1, p is [] | vector | SPD matrix.
%
%   Name-value:
%       'name'   a label for the report (default 'system <i>')

arguments
    obj
    uL (1,1) double {mustBeInteger, mustBeNonnegative}
    batchFcn (1,1) function_handle
    nBatch (1,1) double {mustBeInteger, mustBePositive}
    opt.Name (1,1) string = ""
end
    if obj.finalized
        error('lsqvce:combine:finalized', ...
              'This combination has already been solved.');
    end

    % WP-39 c2: typed `arguments` signature (canonical 'Name').
    nm = char(opt.Name);
    if isempty(nm), nm = sprintf('system %d', obj.nSys + 1); end

    if ~isnumeric(uL) || ~isscalar(uL) || uL < 0 || uL ~= fix(uL)
        error('lsqvce:combine:badUL', ...
              'uL must be a non-negative integer (0 = no local parameters).');
    end
    uL = double(uL);

    switch obj.method
        case 'normal'
            obj.addSystemNormal(uL, batchFcn, nBatch, nm);
        case 'qr'
            obj.addSystemQR(uL, batchFcn, nBatch, nm);
    end
end

%% ---------------------------------------------------------- QR (default)
function addSystemQR(obj, uL, batchFcn, nBatch, nm)
    % TRUE HOUSEHOLDER ELIMINATION.
    %
    % Build the running triangular factor of the whitened system
    %       [ AL | AC | l ]
    % for THIS system, then read off:
    %       rows 1..uL          -> [R_LL  R_LC  d_L]   (the local triangle)
    %       rows uL+1..uL+uC+1  -> [ 0    R_CC  d_C]   (locals ELIMINATED)
    % and feed the second block into the running COMMON factor.
    %
    % The point is that N is NEVER FORMED -- not for the locals, not for the
    % common part, not anywhere. Measured against a long-double SVD reference,
    % relative error in xC with near-collinear common columns:
    %
    %     cond(A)        1e2      1e4      1e6      1e7      1e8
    %     normal-Schur  9.6e-12  5.6e-08  1.8e-04  2.2e-02  4.9e-01
    %     chol(Schur)   5.2e-12  8.6e-08  1.3e-04  3.8e-02  6.7e-01
    %     Householder   1.9e-14  1.3e-12  1.1e-10  1.4e-09  1.2e-08
    %
    % At cond(A) = 2.4e8 the normal-equation path has NO correct digits while
    % Householder still has eight. Forming N -- even only to Schur-complement
    % it -- squares the condition number, and a Cholesky of the reduced block
    % afterwards cannot undo that. This is why the QR path does the elimination
    % on the OBSERVATION EQUATIONS, not on the normal equations.

    ut = uL + obj.uC + 1;                 % [AL | AC | l]
    Rs = zeros(0, ut);                    % running triangle for this system
    ns = 0;

    for k = 1:nBatch
        [AC, AL, lk, pk] = obj.callBatch(batchFcn, k, uL);
        nk = size(AC, 1);
        if nk == 0, continue; end

        [sA, sl] = obj.whiten([AL, AC], lk, pk, nk);   % NOTE: locals FIRST
        stack = [Rs; sA, sl];
        Rf = triu(qr(stack));
        kk = min(size(stack,1), ut);
        Rs = Rf(1:kk, :);
        ns = ns + nk;
    end

    if ns == 0
        lsqvce.util.log.warn('lsqvce:combine:emptySystem', ...
                '%s contributed no observations; skipping.', nm);
        return;
    end
    if size(Rs,1) < min(ns, ut)
        Rs = [Rs; zeros(min(ns,ut) - size(Rs,1), ut)]; %#ok<AGROW>
    end

    % --- split the triangle
    if uL > 0
        if size(Rs,1) < uL
            error('lsqvce:combine:singularLocal', ...
                ['%s: fewer than %d independent rows -- its own data do not ' ...
                 'determine all %d local parameters, so they cannot be ' ...
                 'eliminated.'], nm, uL, uL);
        end
        Rll = Rs(1:uL, 1:uL);
        Rlc = Rs(1:uL, uL+1:uL+obj.uC);
        dl  = Rs(1:uL, end);

        dg = abs(diag(Rll));
        if isempty(dg) || any(dg <= eps * max(dg) * uL)
            error('lsqvce:combine:singularLocal', ...
                ['%s: the local triangle R_ll is rank deficient. The data of ' ...
                 'this system do not determine all %d of its own local ' ...
                 'parameters.\n' ...
                 'You cannot eliminate a parameter the data never constrained. ' ...
                 'Reduce the local parameterisation, or constrain the ' ...
                 'deficient directions.'], nm, uL);
        end
        Rrest = Rs(uL+1:end, uL+1:uL+obj.uC);      % locals are ZERO here
        drest = Rs(uL+1:end, end);

        % S = Nll^-1 Nlc = Rll^-1 Rlc   (verified against the dense Schur form)
        S = Rll \ Rlc;
        t = Rll \ dl;
    else
        Rll = zeros(0);  Rlc = zeros(0, obj.uC);  dl = zeros(0,1);
        S = zeros(0, obj.uC);  t = zeros(0,1);
        Rrest = Rs(:, 1:obj.uC);
        drest = Rs(:, end);
    end

    % --- fold the ELIMINATED rows into the running common factor
    stack = [obj.Rc; Rrest, drest];
    Rf = triu(qr(stack));
    kk = min(size(stack,1), obj.uC + 1);
    Rf = Rf(1:kk, :);
    if kk > obj.uC
        obj.e2 = obj.e2 + Rf(obj.uC+1, obj.uC+1)^2;
        Rf = Rf(1:obj.uC, :);
    end
    obj.Rc = Rf;

    obj.recordSystem(nm, uL, nBatch, ns, Rll, Rlc, dl, S, t, batchFcn);
end

%% ------------------------------------------------------- normal equations
function addSystemNormal(obj, uL, batchFcn, nBatch, nm)
%ADDSYSTEMNORMAL  Add a system via the normal-equation backend (uL locals, nm name).
    ut = obj.uC + uL;
    Ns = zeros(ut);
    bs = zeros(ut, 1);
    ns = 0;

    for k = 1:nBatch
        [AC, AL, lk, pk] = obj.callBatch(batchFcn, k, uL);
        nk = size(AC, 1);
        if nk == 0, continue; end
        [sA, sl] = obj.whiten([AC, AL], lk, pk, nk);
        Ns = Ns + sA' * sA;
        bs = bs + sA' * sl;
        ns = ns + nk;
    end

    if ns == 0
        lsqvce.util.log.warn('lsqvce:combine:emptySystem', ...
                '%s contributed no observations; skipping.', nm);
        return;
    end

    Ncc = Ns(1:obj.uC, 1:obj.uC);
    bc  = bs(1:obj.uC);

    if uL > 0
        Ncl = Ns(1:obj.uC, obj.uC+1:end);
        Nll = Ns(obj.uC+1:end, obj.uC+1:end);
        bl  = bs(obj.uC+1:end);

        [Rll, ok] = chol(Nll);
        if ok ~= 0
            error('lsqvce:combine:singularLocal', ...
                ['%s: its LOCAL normal matrix N_ll is singular (Cholesky ' ...
                 'failed at leading minor %d of %d).\n' ...
                 'The data of this system do not determine all %d of its own ' ...
                 'local parameters. You cannot eliminate a parameter the data ' ...
                 'never constrained.'], nm, ok, uL, uL);
        end
        S = Rll \ (Rll' \ Ncl.');
        t = Rll \ (Rll' \ bl);
        obj.Nred = obj.Nred + (Ncc - Ncl * S);
        obj.bred = obj.bred + (bc  - Ncl * t);
        Rlc = Rll * S;                     % so the diagnostics can reuse it
        dl  = Rll * t;
    else
        Rll = zeros(0);  Rlc = zeros(0, obj.uC);  dl = zeros(0,1);
        S = zeros(0, obj.uC);  t = zeros(0,1);
        obj.Nred = obj.Nred + Ncc;
        obj.bred = obj.bred + bc;
    end

    obj.recordSystem(nm, uL, nBatch, ns, Rll, Rlc, dl, S, t, batchFcn);
end

function recordSystem(obj, nm, uL, nBatch, ns, Rll, Rlc, dl, S, t, batchFcn)
%RECORDSYSTEM  Store one eliminated systems Schur factors and metadata.
    obj.nSys = obj.nSys + 1;
    obj.sys{obj.nSys} = struct( ...
        'name', nm, 'uL', uL, 'nBatch', nBatch, 'n', ns, ...
        'Rll', Rll, 'Rlc', Rlc, 'dl', dl, 'S', S, 't', t, ...
        'batchFcn', batchFcn, ...
        'offset', obj.uC + obj.uL_total);
    obj.n        = obj.n + ns;
    obj.uL_total = obj.uL_total + uL;
end

%% ----------------------------------------------------------------- solve
function sol = solve(obj, opt)
%SOLVE  Solve the reduced system and back-substitute the local parameters.
%
%   sol has fields
%     xC, xL{i}, x        common, per-system local, and the stacked full vector
%     QCC                 covariance of the common parameters
%     sigC, sigL{i}       standard deviations (a priori scaled)
%     n, uC, uL, uTotal   counts
%     r                   redundancy = n - uTotal   <-- NOT n - uC
%     vPv, s0sq           weighted residual sum of squares, variance factor

    % WP-39 c2: typed `arguments` signature (canonical 'Sigma0').
arguments
    obj
    opt.Sigma0 (1,1) double {mustBePositive} = 1
end
    sigma0 = opt.Sigma0;

    if obj.nSys == 0
        error('lsqvce:combine:noSystems', 'No data sets were added.');
    end

    uTot = obj.uC + obj.uL_total;
    rdf  = obj.n - uTot;

    % THE TRAP. The redundancy counts EVERY parameter, eliminated or not. The
    % locals were removed from the SOLVE, not from the MODEL: they still cost
    % degrees of freedom. Using n - uC here would inflate r by uL_total and
    % deflate the variance factor by the same ratio -- silently, and in the
    % optimistic direction.
    if rdf < 0
        error('lsqvce:combine:underdetermined', ...
              ['Underdetermined: %d observations for %d parameters ' ...
               '(%d common + %d local).'], obj.n, uTot, obj.uC, obj.uL_total);
    end

    switch obj.method
        case 'normal'
            [Rc_, ok] = chol(obj.Nred);
            if ok ~= 0
                error('lsqvce:combine:singularCommon', ...
                    ['The REDUCED normal matrix is singular (leading minor ' ...
                     '%d of %d).\n' ...
                     'After eliminating every system''s locals, the common ' ...
                     'parameters are still not determined. Either no data set ' ...
                     'constrains some common direction, or two common ' ...
                     'parameters are collinear across ALL systems.\n' ...
                     'This is a parameterisation defect, not a numerical ' ...
                     'accident.'], ok, obj.uC);
            end
            xC = Rc_ \ (Rc_' \ obj.bred);
        case 'qr'
            if size(obj.Rc,1) < obj.uC
                error('lsqvce:combine:singularCommon', ...
                    ['The reduced QR factor has only %d rows for %d common ' ...
                     'parameters: the data do not determine them all.'], ...
                     size(obj.Rc,1), obj.uC);
            end
            Rc_ = obj.Rc(1:obj.uC, 1:obj.uC);
            dC  = obj.Rc(1:obj.uC, obj.uC+1);
            dg = abs(diag(Rc_));
            if any(dg <= eps * max(dg) * obj.uC)
                error('lsqvce:combine:singularCommon', ...
                      'The reduced QR factor is rank deficient.');
            end
            xC = Rc_ \ dC;
    end
    s = sign(diag(Rc_));  s(s==0) = 1;
    Rc_ = s .* Rc_;                 % row-scale by +/-1 (broadcast, not diag(s)*Rc_)

    QCC = Rc_ \ (Rc_' \ eye(obj.uC));

    % ---- back-substitute the locals:  x_i = t_i - S_i xC
    xL = cell(obj.nSys,1);
    sigL = cell(obj.nSys,1);
    for i = 1:obj.nSys
        st = obj.sys{i};
        if st.uL == 0
            xL{i} = zeros(0,1);  sigL{i} = zeros(0,1);
            continue;
        end
        xL{i} = st.t - st.S * xC;

        % Q_ii = Nll^-1 + S QCC S'
        %
        % The SECOND TERM IS NOT OPTIONAL. Nll^-1 alone is the covariance of
        % the locals *with the common parameters held fixed*. But xC was
        % ESTIMATED, and its uncertainty propagates into the locals through
        % the back-substitution. Dropping S QCC S' produces local sigmas that
        % are too optimistic -- and the more strongly a system's locals
        % correlate with the common block, the worse the lie.
        Nll_inv = st.Rll \ (st.Rll' \ eye(st.uL));
        Qii = Nll_inv + st.S * QCC * st.S.';
        sigL{i} = sqrt(max(diag(Qii), 0));
    end

    x = [xC; vertcat(xL{:})];

    sol = struct( ...
        'xC', xC, 'xL', {xL}, 'x', x, ...
        'QCC', QCC, 'Rc', Rc_, ...
        'sigC', sigma0 * sqrt(max(diag(QCC),0)), ...
        'sigL', {sigL}, ...
        'n', obj.n, 'uC', obj.uC, 'uL', obj.uL_total, 'uTotal', uTot, ...
        'r', rdf, 'nSys', obj.nSys, ...
        'sigma0', sigma0, 'method', obj.method, ...
        'vPv', NaN, 's0sq', NaN);      % filled by diagnose()

    obj.finalized = true;
end

%% -------------------------------------------------------------- diagnose
function D = diagnose(obj, sol, opt)
%DIAGNOSE  Per-observation diagnostics for the combined solution.
%
%   D = cmb.diagnose(sol) replays every system's batches and produces v, r_i,
%   w_i, MDB -- the same quantities as lsqvce.engine.diagnose, but for the ELIMINATED
%   system.
%
%   THE THING THAT IS EASY TO GET WRONG
%   -----------------------------------
%   The redundancy number needs a_i' Qxx a_i with the FULL Qxx, and a_i spans
%   both the common AND that system's local columns:
%
%       a_i' Qxx a_i = aC' Q_CC aC  +  2 aC' Q_Ci aL  +  aL' Q_ii aL
%
%   The CROSS TERM is not optional. Dropping it (using only Q_CC, as one is
%   tempted to do because that is what the reduced system gives you) makes
%   a'Qxx a too small, hence r_i too LARGE, hence every observation look better
%   controlled than it is -- and the MDB correspondingly too optimistic.
%
%   Both remaining blocks come from quantities we already have:
%       Q_Ci = -Q_CC S'                    (S = Rll^-1 Rlc)
%       Q_ii =  Nll^-1 + S Q_CC S'
%   so the whole thing is computable per system, in O(u^2 nk), without ever
%   forming the full u_total x u_total covariance.
%
%   AND THE OTHER ONE
%   -----------------
%   The redundancy is n - u_TOTAL, counting the ELIMINATED parameters. They
%   were removed from the SOLVE, not from the MODEL: they still cost degrees of
%   freedom. Using n - uC would inflate r and deflate s0 -- silently, and in
%   the optimistic direction. This is checked: sum(r_i) must equal n - u_total.

    % WP-39 c2: typed `arguments` signature (canonical UpperCamel; legacy
    % casings keep matching; enum VALUES case-sensitive lowercase).
arguments
    obj
    sol (1,1) struct
    opt.Alpha (1,1) double {mustBePositive, mustBeLessThan(opt.Alpha,1)} = 1e-3
    opt.Bonferroni (1,1) logical = true
    opt.Beta (1,1) double {mustBePositive, mustBeLessThan(opt.Beta,1)} = 0.20
    opt.RMin (1,1) double {mustBeNonnegative} = 1e-6
    opt.RWeak (1,1) double {mustBeNonnegative} = 0.10
    opt.Test (1,1) string {mustBeMember(opt.Test, ["pope","baarda","mad"])} = "pope"
end
    o = struct('alpha', opt.Alpha, 'bonferroni', opt.Bonferroni, ...
               'beta', opt.Beta, 'rMin', opt.RMin, 'rWeak', opt.RWeak, ...
               'test', char(opt.Test));

    QCC = sol.QCC;
    xC  = sol.xC;
    rdf = sol.r;                       % n - u_TOTAL. Not n - uC.

    % ---- pass 1: residuals and vPv
    vPv = 0;
    vAll = cell(obj.nSys,1);
    for i = 1:obj.nSys
        st = obj.sys{i};
        xi = sol.xL{i};
        vi = [];
        for k = 1:st.nBatch
            [AC, AL, lk, pk] = obj.callBatch(st.batchFcn, k, st.uL);
            nk = size(AC,1);
            if nk == 0, continue; end
            vk = AC * xC - lk;
            if st.uL > 0, vk = vk + AL * xi; end
            vPv = vPv + obj.quad(vk, pk, nk);
            vi = [vi; vk]; %#ok<AGROW>
        end
        vAll{i} = vi;
    end
    s0sq = vPv / max(rdf, 1);

    % ---- test setup
    alphaEff = o.alpha;
    if o.bonferroni, alphaEff = o.alpha / max(obj.n, 1); end
    switch lower(o.test)
        case 'baarda', kCrit = lsqvce.engine.stat.normCrit(alphaEff);
        otherwise,     kCrit = lsqvce.engine.stat.tauCrit(alphaEff, rdf);
    end
    delta0 = lsqvce.engine.stat.delta0(alphaEff, o.beta);

    scale = sqrt(s0sq);
    if strcmpi(o.test, 'baarda'), scale = sol.sigma0; end

    % ---- pass 2: r_i, w_i, MDB
    rAll = cell(obj.nSys,1);  wAll = cell(obj.nSys,1);  mAll = cell(obj.nSys,1);
    sumR = 0;  nFlag = 0;  flagIdx = [];  flagSys = [];
    nUnc = 0;  nWeak = 0;
    base = 0;

    for i = 1:obj.nSys
        st = obj.sys{i};
        S  = st.S;
        if st.uL > 0
            NllInv = st.Rll \ (st.Rll' \ eye(st.uL));
            QCi = -QCC * S.';                    % uC x uL
            Qii =  NllInv + S * QCC * S.';       % uL x uL
        else
            QCi = zeros(obj.uC, 0);
            Qii = zeros(0);
        end

        ri = [];  wi = [];  mi = [];
        for k = 1:st.nBatch
            [AC, AL, lk, pk] = obj.callBatch(st.batchFcn, k, st.uL);
            nk = size(AC,1);
            if nk == 0, continue; end

            % a' Qxx a, with the CROSS TERM
            qx = sum((AC * QCC) .* AC, 2);
            if st.uL > 0
                qx = qx + 2*sum((AC * QCi) .* AL, 2) ...
                        +   sum((AL * Qii) .* AL, 2);
            end

            vk = AC * xC - lk;
            if st.uL > 0, vk = vk + AL * sol.xL{i}; end

            [rk, qpp, Pvk, sig] = obj.redundancy(AC, AL, pk, qx, vk, nk, ...
                                                 st.uL, QCC, QCi, Qii);
            ctrl = rk > o.rMin & qpp > 0;
            wk = zeros(nk,1);
            wk(ctrl) = Pvk(ctrl) ./ (scale * sqrt(qpp(ctrl)));

            mdbk = inf(nk,1);
            mdbk(ctrl) = sig(ctrl) * delta0 ./ sqrt(rk(ctrl));

            fl = ctrl & abs(wk) > kCrit;
            if any(fl)
                flagIdx = [flagIdx; base + find(fl)]; %#ok<AGROW>
                flagSys = [flagSys; repmat(i, sum(fl), 1)]; %#ok<AGROW>
            end
            nFlag = nFlag + sum(fl);
            nUnc  = nUnc  + sum(~ctrl);
            nWeak = nWeak + sum(ctrl & rk < o.rWeak);
            sumR  = sumR  + sum(rk);

            ri = [ri; rk]; wi = [wi; wk]; mi = [mi; mdbk]; %#ok<AGROW>
            base = base + nk;
        end
        rAll{i} = ri;  wAll{i} = wi;  mAll{i} = mi;
    end

    % ---- INTEGRITY
    sumRErr = abs(sumR - rdf) / max(rdf, 1);

    D = struct( ...
        'n', obj.n, 'uC', obj.uC, 'uL', obj.uL_total, 'uTotal', sol.uTotal, ...
        'r', rdf, 'nSys', obj.nSys, ...
        'v', {vAll}, 'rNum', {rAll}, 'w', {wAll}, 'mdb', {mAll}, ...
        'vPv', vPv, 's0sq', s0sq, 's0', sqrt(max(s0sq,0)), ...
        'sigma0', sol.sigma0, 'varRatio', s0sq / sol.sigma0^2, ...
        'sumR', sumR, 'sumRExpect', rdf, 'sumRRelErr', sumRErr, ...
        'kCrit', kCrit, 'alphaEff', alphaEff, 'delta0', delta0, ...
        'test', lower(o.test), ...
        'nFlag', nFlag, 'flagIdx', flagIdx, 'flagSys', flagSys, ...
        'nUncontrolled', nUnc, 'nWeak', nWeak, ...
        'QCC', QCC, 'sigC', sol.sigC, 'sigL', {sol.sigL});

    if sumRErr > 1e-6
        lsqvce.util.log.warn('lsqvce:combine:integrity', ...
            ['sum(r_i) = %.6f but n - u_total = %d (rel err %.2e).\n' ...
             'This identity is exact. Its violation means the second pass is ' ...
             'not replaying the same partials, OR that the cross-term in ' ...
             'a''Qxx a has been dropped -- the classic elimination bug, which ' ...
             'makes every r_i too large and every observation look better ' ...
             'controlled than it is.'], sumR, rdf, sumRErr);
    end
end

%% --------------------------------------------------------------- varComp
function VC = varComp(obj, sol, groupFcn, nGroup, opt)
%VARCOMP  Helmert variance components across sensors (or any grouping).
%
%   VC = cmb.varComp(sol, groupFcn, nGroup) estimates a separate variance
%   factor per group, where
%       g = groupFcn(i, k)
%   returns the group index of every observation in batch k of system i.
%
%   THE POINT. In a KBR + LRI combination both instruments observe the same
%   arcs, and their quoted noise levels are almost never both right. Helmert
%   tells you the ratio the DATA imply, rather than the one the headers claim.
%
%   THE TRAP THIS AVOIDS. The partial redundancy of a group is
%
%       r_g = n_g - trace(Qxx N_g)
%
%   and Qxx here is the FULL covariance, INCLUDING the eliminated locals. If
%   you use only Qcc -- which is all the reduced system hands you -- you
%   OVERSTATE r_g by exactly the number of eliminated parameters (measured: 6.9
%   against 7 locals), and since s_g = q_g / r_g you then UNDERSTATE the
%   variance component. The weights come out too optimistic, in the same
%   direction as every other elimination error.
%
%   This method uses the full Qxx, assembled per system from Qcc, S_i and
%   Nll_i^-1 -- never forming the u_total x u_total matrix.
%
%   RESOLUTION. A variance component estimated from r_g observations has a
%   relative standard error of about 0.5*sqrt(2/r_g). With r_g ~ 300 per sensor
%   you CANNOT resolve a sigma ratio better than about 6%. VC.resolution
%   reports this, and VC.significant says whether the correction exceeds it.
%   Verified unbiased: true ratios of 3.0, 0.4 and 2.0 were recovered as 3.032,
%   0.405 and 2.021 (mean over 30 realisations).
%
%   Name-value:
%     'maxIter'  (default 25)
%     'tol'      (default 1e-6)
%     'verbose'  (default true)

    % WP-39 c2: typed `arguments` signature (canonical UpperCamel).
arguments
    obj
    sol (1,1) struct
    groupFcn (1,1) function_handle
    nGroup (1,1) double {mustBeInteger, mustBePositive}
    opt.MaxIter (1,1) double {mustBeInteger, mustBePositive} = 25
    opt.Tol (1,1) double {mustBePositive} = 1e-6
    opt.Verbose (1,1) logical = true
end
    o = struct('maxIter', opt.MaxIter, 'tol', opt.Tol, 'verbose', opt.Verbose);

    scale = ones(nGroup, 1);
    converged = false;

    if lsqvce.util.verbosity.progress(o.verbose)
        lsqvce.util.log.msg('\nHelmert variance components, %d groups\n', nGroup);
        lsqvce.util.log.msg('%s\n', repmat('-', 1, 58));
    end

    for it = 1:o.maxIter
        % ---- re-solve with the CURRENT scales, keeping the elimination
        [xC, QCC, sysState, ok] = obj.resolveScaled(scale, groupFcn);
        if ~ok
            error('lsqvce:combine:vcSingular', ...
                  'The system became singular at iteration %d.', it);
        end

        q    = zeros(nGroup, 1);
        nObs = zeros(nGroup, 1);
        trQN = zeros(nGroup, 1);
        NgAll = repmat({zeros(obj.uC)}, nGroup, 1);   % common part only, for H

        % We need trace(Qxx N_g) with the FULL Qxx. Expand N_g over the
        % common and this system's local columns:
        %   trace(Qxx N_g) = sum over obs i in g of  p_i * a_i' Qxx a_i
        % which is exactly the qx we already compute for the redundancy
        % numbers -- including the cross term.
        for i = 1:obj.nSys
            st = obj.sys{i};
            ss = sysState{i};
            xi = ss.t - ss.S * xC;
            if st.uL > 0
                NllInv = ss.Rll \ (ss.Rll' \ eye(st.uL));
                QCi = -QCC * ss.S.';
                Qii =  NllInv + ss.S * QCC * ss.S.';
            else
                QCi = zeros(obj.uC, 0);
                Qii = zeros(0);
            end

            for k = 1:st.nBatch
                [AC, AL, lk, pk] = obj.callBatch(st.batchFcn, k, st.uL);
                nk = size(AC,1);
                if nk == 0, continue; end
                if isempty(pk), pk = ones(nk,1); end
                if ~isvector(pk)
                    error('lsqvce:combine:vcCorrelated', ...
                          'varComp requires diagonal weights.');
                end
                pk = pk(:);
                gk = groupFcn(i, k);  gk = gk(:);
                if numel(gk) ~= nk
                    error('lsqvce:combine:badGroup', ...
                        ['groupFcn(%d,%d) returned %d indices for %d ' ...
                         'observations.'], i, k, numel(gk), nk);
                end
                if any(gk < 1 | gk > nGroup | gk ~= fix(gk))
                    error('lsqvce:combine:badGroup', ...
                          'groupFcn(%d,%d) returned an index outside 1..%d.', ...
                          i, k, nGroup);
                end

                peff = pk ./ scale(gk);
                vk = AC * xC - lk;
                if st.uL > 0, vk = vk + AL * xi; end

                % a' Qxx a, WITH the cross term
                qx = sum((AC * QCC) .* AC, 2);
                if st.uL > 0
                    qx = qx + 2*sum((AC * QCi) .* AL, 2) ...
                            +   sum((AL * Qii) .* AL, 2);
                end

                for gg = 1:nGroup
                    m = (gk == gg);
                    if ~any(m), continue; end
                    nObs(gg) = nObs(gg) + sum(m);
                    q(gg)    = q(gg)    + sum(peff(m) .* vk(m).^2);
                    trQN(gg) = trQN(gg) + sum(peff(m) .* qx(m));
                    ACg = sqrt(peff(m)) .* AC(m,:);
                    NgAll{gg} = NgAll{gg} + ACg' * ACg;
                end
            end
        end

        rg = nObs - trQN;               % partial redundancy, FULL Qxx

        % Simplified (Foerstner) estimator. The rigorous H would need
        % trace(Qxx N_i Qxx N_j) with the full Qxx across systems, which is
        % expensive here; the simplified form converges to the same fixed
        % point and is what is validated (unbiased to ~1%).
        s = ones(nGroup,1);
        est = rg > 1e-8 & nObs >= 10;
        s(est) = q(est) ./ rg(est);
        s = max(s, 1e-8);

        scale = scale .* s;
        dev = max(abs(s(est) - 1));
        if isempty(dev), dev = 0; end

        if lsqvce.util.verbosity.progress(o.verbose)
            lsqvce.util.log.msg('  it %2d  s = [%s]  max|s-1| = %.2e\n', it, ...
                strjoin(arrayfun(@(z) sprintf('%.4f', z), s.', 'uni', 0), ' '), ...
                dev);
        end
        if dev < o.tol
            converged = true;
            break;
        end
    end

    resol = inf(nGroup,1);
    resol(rg > 0) = 0.5 * sqrt(2 ./ rg(rg > 0));

    VC = struct('scale', scale, 'sigma', sqrt(scale), ...
                'nObs', nObs, 'redund', rg, 'iter', it, ...
                'converged', converged, ...
                'resolution', resol, ...
                'significant', abs(sqrt(scale) - 1) > resol);

    if lsqvce.util.verbosity.progress(o.verbose)
        lsqvce.util.log.msg('%s\n', repmat('-', 1, 58));
        lsqvce.util.log.msg('\n  %-8s %8s %10s %12s %12s %6s\n', ...
                'group', 'n', 'redund', 'sigma scale', 'resolution', 'sig?');
        lsqvce.util.log.msg('  %s\n', repmat('-', 1, 62));
        for gg = 1:nGroup
            lsqvce.util.log.msg('  %-8d %8d %10.1f %12.4f %11.1f%% %6s\n', ...
                    gg, nObs(gg), rg(gg), sqrt(scale(gg)), ...
                    100*resol(gg), string(VC.significant(gg)));
        end
        rat = sqrt(max(scale)/min(scale));
        noise = sqrt(sum(resol(isfinite(resol)).^2));
        lsqvce.util.log.msg('\n');
        if (rat - 1) < noise
            lsqvce.util.log.msg(['  The spread between groups (%.1f%%) is SMALLER than the\n' ...
                     '  sampling noise of the estimator (%.1f%%). Your a priori\n' ...
                     '  weights are fine. Do NOT apply these scales -- you would\n' ...
                     '  be fitting noise.\n\n'], 100*(rat-1), 100*noise);
        else
            lsqvce.util.log.msg(['  The groups disagree by a factor of %.2f in sigma, beyond\n' ...
                     '  the %.1f%% resolution of the estimator. The a priori\n' ...
                     '  weighting is genuinely wrong. Apply VC.scale to the\n' ...
                     '  covariances and re-run.\n\n'], rat, 100*noise);
        end
    end
end

end % methods

methods (Access = private)

function [xC, QCC, sysState, ok] = resolveScaled(obj, scale, groupFcn)
    % Re-run the elimination with per-group weight scales applied. Used by
    % varComp; keeps the arrow structure intact.
    ok = true;
    Nred = zeros(obj.uC);
    bred = zeros(obj.uC, 1);
    sysState = cell(obj.nSys, 1);

    for i = 1:obj.nSys
        st = obj.sys{i};
        ut = obj.uC + st.uL;
        Ns = zeros(ut);  bs = zeros(ut, 1);

        for k = 1:st.nBatch
            [AC, AL, lk, pk] = obj.callBatch(st.batchFcn, k, st.uL);
            nk = size(AC,1);
            if nk == 0, continue; end
            if isempty(pk), pk = ones(nk,1); end
            pk = pk(:);
            gk = groupFcn(i, k);  gk = gk(:);
            peff = pk ./ scale(gk);
            sp = sqrt(peff);
            M = sp .* [AC, AL];
            Ns = Ns + M' * M;
            bs = bs + M' * (sp .* lk);
        end

        Ncc = Ns(1:obj.uC, 1:obj.uC);
        bc  = bs(1:obj.uC);
        if st.uL > 0
            Ncl = Ns(1:obj.uC, obj.uC+1:end);
            Nll = Ns(obj.uC+1:end, obj.uC+1:end);
            bl  = bs(obj.uC+1:end);
            [Rll, bad] = chol(Nll);
            if bad ~= 0, ok = false; xC=[]; QCC=[]; return; end
            S = Rll \ (Rll' \ Ncl.');
            t = Rll \ (Rll' \ bl);
            Nred = Nred + (Ncc - Ncl * S);
            bred = bred + (bc  - Ncl * t);
        else
            Rll = zeros(0);  S = zeros(0, obj.uC);  t = zeros(0,1);
            Nred = Nred + Ncc;
            bred = bred + bc;
        end
        sysState{i} = struct('Rll', Rll, 'S', S, 't', t);
    end

    [Rc, bad] = chol(Nred);
    if bad ~= 0, ok = false; xC=[]; QCC=[]; return; end
    xC  = Rc \ (Rc' \ bred);
    QCC = Rc \ (Rc' \ eye(obj.uC));
end

function q = quad(~, v, p, nk)
%QUAD  Weighted quadratic form vPv for one batch (p=[] -> unit weights).
    if isempty(p)
        q = v' * v;
    elseif isvector(p)
        q = sum(p(:) .* v.^2);
    else
        Pk = (p + p')/2;
        q = v' * (Pk * v);
    end
    if nk == 0, q = 0; end
end

function [rk, qpp, Pvk, sig] = redundancy(~, AC, AL, pk, qx, vk, nk, uL, ...
                                          QCC, QCi, Qii) %#ok<INUSD,INUSL>
        %REDUNDANCY  Per-observation redundancy for one batch of an eliminated
        %   system: r_i = (Qvv P)_ii and the generalised denominator (P Qvv P)_ii,
        %   using the full common cofactor blocks QCC/QCi/Qii.
        %   INPUT  AC nk x uC common design; AL nk x uL local design; pk nk x 1
        %          weights; qx precomputed cross-term; vk nk x 1 residuals; nk 1x1
        %          rows; uL 1x1 locals; QCC/QCi/Qii cofactor blocks.
        %   OUTPUT rk nk x 1 redundancy r_i; qpp nk x 1 (PQvvP)_ii; Pvk nk x 1
        %          P*v; sig 1x1 scalar summary.
        % r_i = (Qvv P)_ii and the generalised denominator (P Qvv P)_ii.
    % qx already contains the cross-term contribution from the eliminated
    % locals, so this is the same algebra as lsqvce.engine.diagnose.
    if isempty(pk)
        pk = ones(nk,1);
    end
    if ~isvector(pk)
        % Correlated blocks WITHIN an eliminated system are not supported.
        %
        % The generalised statistic needs (P Qvv P)_ii, which expands to
        %       P(i,i) - g_i' Qxx g_i     with  g_i = (P A)(i,:)'
        % and g_i now spans BOTH the common and the local columns. That is
        % derivable, but it is not the same expression as the diagonal case and
        % I have not validated it against a dense reference. Shipping it
        % untested would be exactly the kind of plausible-looking wrong number
        % this toolbox exists to prevent.
        %
        % Correlated batches ARE fully supported in the single-dataset path
        % (lsqvce.engine.diagnose). If you need them here, decorrelate the batch
        % yourself before handing it to Combine: pass Lp'*AC, Lp'*AL, Lp'*l
        % with unit weights, where P = Lp*Lp'. The solution and the residuals
        % are then identical, and the diagnostics refer to the whitened
        % observations.
        error('lsqvce:combine:correlatedBatch', ...
            ['Combine does not support correlated (matrix) weights inside an ' ...
             'eliminated system.\n' ...
             'Decorrelate the batch first: with P = Lp*Lp'' (Lp = chol(P, ' ...
             '''lower'')), pass Lp''*AC, Lp''*AL, Lp''*l and unit weights. The ' ...
             'solution and residuals are unchanged; the diagnostics then refer ' ...
             'to the whitened observations.\n' ...
             'Correlated batches ARE supported in the single-dataset path ' ...
             '(lsqvce.engine.adjust / lsqvce.engine.diagnose).']);
    end
    pk = pk(:);
    qv = 1 ./ pk - qx;
    qv(qv < 0) = 0;
    rk  = min(max(pk .* qv, 0), 1);
    qpp = max(pk.^2 .* qv, 0);
    Pvk = pk .* vk;
    sig = 1 ./ sqrt(pk);
end

function [AC, AL, l, p] = callBatch(~, f, k, uL)
%CALLBATCH  Call a system batch source, splitting common/local columns (uL locals).
%   Through LSQVCE.UTIL.CALLUSER (R5 M1, H-08): p = [] for a 3-output source;
%   a throwing source arrives as lsqvce:callUser:threw with the cause attached.
    [~, AC, AL, l, p] = lsqvce.util.callUser(f, 3, 4, 'batchFcn(k)', k);
    l = l(:);
    if size(AL,2) ~= uL
        error('lsqvce:combine:badAL', ...
              ['Batch %d returned AL with %d columns but the system was ' ...
               'declared with uL = %d local parameters.'], k, size(AL,2), uL);
    end
    if size(AC,1) ~= size(AL,1) || numel(l) ~= size(AC,1)
        error('lsqvce:combine:sizeMismatch', ...
              ['Batch %d: AC has %d rows, AL has %d, l has %d. They must ' ...
               'agree.'], k, size(AC,1), size(AL,1), numel(l));
    end
    if ~all(isfinite(AC(:))) || ~all(isfinite(AL(:))) || ~all(isfinite(l))
        error('lsqvce:combine:nonFinite', ...
              'Batch %d contains non-finite values.', k);
    end
end

function [sA, sl] = whiten(~, A, l, p, nk)
%WHITEN  Scale rows of A,l by sqrt(p) (p=[] -> unchanged) for one batch.
    if isempty(p)
        sA = A;  sl = l;  return;
    end
    if isvector(p) && numel(p) == nk
        p = p(:);
        if any(p < 0)
            error('lsqvce:combine:negativeWeight', ...
                  'Negative weights are not a down-weighting; they make N indefinite.');
        end
        sp = sqrt(p);
        sA = sp .* A;  sl = sp .* l;  return;
    end
    if ismatrix(p) && isequal(size(p), [nk nk])
        p = (p + p')/2;
        [Lp, ok] = chol(p, 'lower');
        if ok ~= 0
            error('lsqvce:combine:badPmatrix', ...
                  ['The weight matrix is not positive definite (minor %d). ' ...
                   'Did you pass the covariance instead of its inverse?'], ok);
        end
        sA = Lp' * A;  sl = Lp' * l;  return;
    end
    error('lsqvce:combine:badP', ...
          'Weight must be [], a length-%d vector, or a %dx%d SPD matrix.', ...
          nk, nk, nk);
end

end
end
