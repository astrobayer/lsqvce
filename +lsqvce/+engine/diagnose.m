function D = diagnose(sol, batchFcn, opt)
%LSQVCE.ENGINE.DIAGNOSE  Second-pass per-observation diagnostics for a batch LSQ.
%
%   D = LSQVCE.ENGINE.DIAGNOSE(sol, batchFcn) sweeps the observations a second time
%   and produces, for every observation, the quantities you actually need in
%   order to judge the adjustment:
%
%       v      residual
%       qv     cofactor of the residual, diag(Qvv)
%       r      redundancy number,  r_i = (Qvv P)_ii
%       w      normalised residual (Baarda) or studentised (Pope)
%       mdb    minimal detectable bias in the units of the observation
%       ext    external reliability: |dx| induced by an undetected MDB,
%              expressed as a multiple of the parameter's own sigma
%
%   Nothing of size n-by-n or n-by-u is ever formed. Per batch the cost is
%   O(u^2 * nk) and the memory is O(u^2 + nk).
%
%   INPUTS
%     sol       the struct returned by lsqvce.engine.Accumulator/solve
%     batchFcn  function handle, [A,l,p] = batchFcn(k), k = 1..nBatch.
%               Must reproduce EXACTLY the batches given to the accumulator,
%               in the same order. This is checked (see 'verify').
%
%   NAME-VALUE
%     'nBatch'      number of batches (default sol.nBatch)
%     'alpha'       per-test significance          (default 0.001)
%     'bonferroni'  correct alpha for n tests      (default true)
%     'beta'        MDB type-II error              (default 0.20)
%     'rMin'        redundancy below which an observation is declared
%                   UNCONTROLLED and its w is not reported (default 1e-6)
%     'rWeak'       redundancy below which an observation is flagged as
%                   weakly controlled                        (default 0.10)
%     'test'        'auto' | 'pope' | 'baarda'   (default 'auto')
%                   'baarda' assumes sigma0 is known a priori; 'pope'
%                   studentises with the estimated s0. 'auto' picks pope,
%                   because in practice sigma0 is essentially never known.
%     'store'       'all' | 'flagged' | 'none'   (default 'all')
%                   'flagged' keeps only observations that trip a test plus
%                   the worst offenders -- use this when n is enormous.
%     'nWorst'      how many worst-|w| observations to retain under
%                   'store','flagged'                        (default 1000)
%     'verify'      re-check sum(r) == n-u and vPv consistency (default true)
%
%   ON BONFERRONI
%   -------------
%   Data snooping runs one hypothesis test PER OBSERVATION. With n = 1e7 and an
%   uncorrected alpha = 0.001, roughly 10,000 perfectly good observations will
%   exceed the critical value by chance alone. Reporting those as "outliers" is
%   not conservatism, it is noise. The default therefore corrects alpha to
%   alpha/n. This raises the critical value from ~3.29 to ~5.7 at n = 1e7, and
%   is the difference between a usable outlier list and a meaningless one.
%
%   OUTPUT
%     D  1x1 struct  residuals, redundancy numbers r_i, leverage h_i, the
%                    standardised w_i / tau_i, MDB, flagIdx and the
%                    contamination gauge -- everything LSQVCE.REPORT reads
%
%   EXAMPLE
%       % The SAME batch source must be replayed -- this is checked.
%       [sol, D] = lsqvce.engine.adjust(u, src, 4, 'verbose', false);
%
%       % Re-diagnose with the ROBUST scale instead of Pope's tau. Do this
%       % whenever D.contamLevel is 'moderate' or worse: Pope goes BLIND at
%       % ~10%% contamination because the blunders inflate the s0 it divides by.
%       Dm = lsqvce.engine.diagnose(sol, src, 'nBatch', 4, 'test', 'mad');
%
%       % Dm.flagIdx    -> [40 150]
%       % Dm.quad       -> the (w, r) quadrants; read UNTESTED, not just w
%       % Dm.mdb        -> largest error that could hide in each observation
%
%   ERRORS
%     lsqvce:diagnose:nBatch          nBatch disagrees with what the solve used
%     lsqvce:diagnose:batchMismatch   batchFcn returned DIFFERENT data on this
%                                     second pass. It is called twice -- once to
%                                     accumulate, once to diagnose -- and must
%                                     be deterministic
%     lsqvce:diagnose:badPmatrix      a full weight matrix where a diagonal was
%                                     expected
%     lsqvce:diagnose:badTest         'test' is not 'auto', 'pope' or 'baarda'
%     lsqvce:diagnose:noRobustScale   a robust scale was requested but the
%                                     residuals cannot support one
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (b) b2 (F-07):
%   typed `arguments` signature, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   See also LSQVCE.ENGINE.ADJUST, LSQVCE.REPORT, LSQVCE.ENGINE.SNOOP.
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

% WP-39 stage (b) b2: typed `arguments` signature (canonical
% UpperCamelCase; legacy casings keep matching -- MATLAB NV
% semantics; enum VALUES case-sensitive lowercase).
arguments
    sol
    batchFcn
    opt.NBatch double {mustBeScalarOrEmpty} = []
    opt.Alpha (1,1) double {mustBePositive, mustBeLessThan(opt.Alpha,1)} = 1e-3
    opt.Bonferroni (1,1) logical = true
    opt.Beta (1,1) double {mustBePositive, mustBeLessThan(opt.Beta,1)} = 0.20
    opt.RMin (1,1) double {mustBeNonnegative} = 1e-6
    opt.RWeak (1,1) double {mustBeNonnegative, mustBeLessThanOrEqual(opt.RWeak,1)} = 0.10
    opt.Test (1,1) string {mustBeMember(opt.Test, ["auto","pope","baarda","mad"])} = "auto"
    opt.RobustCheck (1,1) logical = true
    opt.Store (1,1) string {mustBeMember(opt.Store, ["all","flagged","none"])} = "all"
    opt.NWorst (1,1) double {mustBeNonnegative} = 1000
    opt.Verify (1,1) logical = true
    opt.Progress (1,1) logical = false
end

o = struct('nBatch', opt.NBatch, 'alpha', opt.Alpha, 'bonferroni', opt.Bonferroni, 'beta', opt.Beta, 'rMin', opt.RMin, 'rWeak', opt.RWeak, 'test', char(opt.Test), 'robustCheck', opt.RobustCheck, 'store', char(opt.Store), 'nWorst', opt.NWorst, 'verify', opt.Verify, 'progress', opt.Progress);

nB = o.nBatch;
if isempty(nB)
    if isfield(sol, 'nBatch') && sol.nBatch > 0
        nB = sol.nBatch;
    else
        error('lsqvce:diagnose:nBatch', ...
              'Cannot infer the batch count; pass ''nBatch''.');
    end
end

R  = sol.R;
x  = sol.x;
n  = sol.n;
u  = sol.u;
rdf = sol.r;                          % redundancy = n - u

% -------------------------------------------------------------- test setup
alphaEff = o.alpha;
if o.bonferroni
    alphaEff = o.alpha / max(n, 1);
end

testName = lower(o.test);
if strcmp(testName, 'auto'), testName = 'pope'; end

switch testName
    case 'pope'
        kCrit = lsqvce.engine.stat.tauCrit(alphaEff, rdf);
    case 'mad'
        % The MAD statistic is still a STUDENTISED residual -- only the scale in
        % the denominator is robust. It remains tau-distributed, so the critical
        % value is Pope's.
        %
        % THIS CASE WAS MISSING. kCrit was simply never assigned on the 'mad'
        % path, and diagnose died with "Unrecognized function or variable
        % kCrit". The entire robust-scale cure -- the toolbox's headline fix for
        % Pope's breakdown -- was therefore DEAD ON ARRIVAL. The Python port
        % computed kCrit separately and so structurally could not see it.
        % Caught only by finally running the suite in MATLAB.
        kCrit = lsqvce.engine.stat.tauCrit(alphaEff, rdf);
    case 'baarda'
        kCrit = lsqvce.engine.stat.normCrit(alphaEff);
    otherwise
        error('lsqvce:diagnose:badTest', ...
              'Unknown test "%s". Use pope, baarda or mad.', testName);
end
delta0 = lsqvce.engine.stat.delta0(alphaEff, o.beta);

% Pope's tau is bounded by sqrt(rdf). Under a hard Bonferroni correction the
% requested critical value can exceed that bound, in which case NO residual can
% ever be flagged and the snooping is vacuous. Say so rather than silently
% returning an empty outlier list.
tauBound = sqrt(rdf);
snoopVacuous = strcmp(testName, 'pope') && kCrit >= tauBound;

% ---------------------------------------------------------- allocate output
storeAll = strcmpi(o.store, 'all');
storeNone = strcmpi(o.store, 'none');
if storeAll
    v_all  = zeros(n, 1);
    r_all  = zeros(n, 1);
    w_all  = zeros(n, 1);
    qv_all = zeros(n, 1);
    mdb_all = zeros(n, 1);
else
    v_all = []; r_all = []; w_all = []; qv_all = []; mdb_all = [];
end

% Running scalars -- these are the honest recomputations.
vPv_direct = 0;      % sum(p .* v.^2)   -- NOT l'Pl - b'x
sumR       = 0;      % must equal n - u
nSeen      = 0;

% Streaming collectors that never grow with n
flagIdx = [];  flagW = [];  flagR = [];  flagV = [];  flagMDB = [];
worstW  = -inf(min(o.nWorst, n), 1);
worstI  = zeros(numel(worstW), 1);

nUncontrolled = 0;   idxUncontrolled = [];

% (w, r) quadrant counters and the bounded MDB reservoir.
nQ_cleanOut = 0;  nQ_ambiguous = 0;  nQ_verified = 0;  nQ_untested = 0;
mdbCap   = 5000;          % reservoir size: quantiles to ~1% with 5000 samples
mdbPool  = zeros(0, 1);
nMdbSeen = 0;
nMdbInf  = 0;
nWeak         = 0;
rmin_seen     = inf;  rmin_idx = 0;
maxAbsW       = 0;    maxAbsW_idx = 0;

% External-reliability accumulator: the worst effect on any parameter, in
% units of that parameter's own sigma. Needs Qxx diagonal, which we can get
% once from R without ever forming Qxx fully at n-scale.
%
% COST NOTE. This is the only O(u^3) step in the toolbox and the only O(u^2)
% allocation. Everything else is O(u^2 * nk) time and O(u^2 + nk) memory, so
% for large u this -- not the observation count -- is the bottleneck:
%
%     u = 1e3   ->  Qxx is    8 MB   (trivial)
%     u = 1e4   ->  Qxx is  800 MB   (painful)
%     u = 3e4   ->  Qxx is  7.2 GB   (fatal)
%
% If R is SPARSE, inv(R) is generally DENSE, so forming Rinv here would throw
% away the entire benefit of the sparse path. We therefore warn loudly rather
% than quietly allocating a matrix that will not fit.
if u > 5000
    bytes = 8 * u^2;
    lsqvce.util.log.warn('lsqvce:diagnose:largeU', ...
        ['u = %d parameters. Forming Qxx requires a dense %dx%d matrix ' ...
         '(%.2f GB) and an O(u^3) inversion.\n' ...
         'This is the ONLY step in the toolbox that scales with u^3 -- the ' ...
         'observation handling is fully streaming. If this is too slow or too ' ...
         'large, you do not need full Qxx: r_i and w_i need only the ' ...
         'triangular solve R''\\A'', which is already streaming. Pass ' ...
         '''store'',''none'' to skip the correlation matrix, or extract ' ...
         'diag(Qxx) by selected inversion (Takahashi) if you have a sparse R.'], ...
         u, u, u, bytes / 2^30);
end
if issparse(R)
    % inv(sparse) is dense. Say so, then do it -- but only once, and only
    % because Qxx is genuinely needed for the correlation report.
    Rinv = full(R) \ eye(u);
else
    Rinv = R \ eye(u);
end

% --- FREE-NETWORK CORRECTION.
% On a free network R factorises N + B*B', so inv(R'R) is NOT the pseudoinverse
% N+: they differ by B*B' (with B orthonormal). Reporting inv(R'R) as Qxx
% overstates the parameter variances -- measured at 12% on a small levelling
% net. Subtract the correction.
%
% Crucially this affects ONLY Qxx. It does NOT affect r_i or w_i, because
% A*B = 0 exactly (the observations cannot see the datum), so the correction
% term drops out of A*Qxx*A'. Verified: sum(r_i) is identical with either Qxx.
% That is the formal statement that the residuals -- and therefore every
% diagnostic built on them -- are datum-invariant, as they must be physically.
if isfield(sol, 'QxxCorr') && ~isempty(sol.QxxCorr)
    Bc = sol.QxxCorr;
else
    Bc = zeros(u, 0);
end
dDefect = size(Bc, 2);

Qxx     = Rinv * Rinv.' - Bc * Bc.';
qxxDiag = diag(Qxx);

% The correction can push a diagonal entry a few eps negative for a parameter
% that lies almost entirely in the datum null space (i.e. one that is barely
% estimable at all). Clamp rather than emit a complex sigma.
negQ = qxxDiag < 0;
if any(negQ)
    qxxDiag(negQ) = 0;
end
sigX    = sqrt(max(qxxDiag, 0));
extMax  = 0;  extMax_idx = 0;  extMax_par = 0;

s0sq = sol.s0sq;
if ~isfinite(s0sq) || s0sq <= 0
    % Zero redundancy, or a poisoned vPv. Fall back to sigma0 so the studentised
    % residual is at least defined; flag it in the report.
    s0sq_use = sol.sigma0^2;
    s0Fallback = true;
else
    s0sq_use = s0sq;
    s0Fallback = false;
end

% =========================================================================
% ROBUST variance factor -- the fix for Pope's breakdown
% =========================================================================
% Pope's tau studentises the residual by s0, which is estimated FROM THE DATA.
% When several blunders are present they inflate s0, and since each blunder is
% then divided by that inflated s0, THEY MASK EACH OTHER. This is not a leverage
% effect and iterating does not cure it: if nothing crosses the threshold on the
% first pass, the rejection loop terminates immediately having removed nothing.
%
% Measured breakdown (n=200, u=8, blunders of 20 sigma, fraction of TRUE
% blunders detected):
%
%     contamination     Pope      MAD-scaled
%          2.5 %        100 %        100 %
%          5   %         23 %        100 %      <- Pope already collapsing
%         10   %          0 %        100 %      <- Pope totally blind
%         20   %          0 %         99 %
%         30   %          0 %         75 %
%
% At 10 % contamination a 100-sigma blunder is UNDETECTABLE by Pope: max|w|
% saturates around 3.3 while the critical value sits at 3.9. Tau is bounded by
% sqrt(r), so it cannot grow without limit however gross the error.
%
% The cure is a scale estimate with a high breakdown point. The MAD of the
% standardised residuals v_i/sqrt(qv_i) costs nothing on clean data (measured:
% zero extra false alarms over 200 clean realisations; mean scale 1.0036 vs
% s0's 0.9976) and is merely less EFFICIENT (sd 0.085 vs 0.051). That is a
% trade every adjustment should be willing to make.
%
% We therefore always compute BOTH, and cross-check them (see robustCheck).
% A large disagreement between s0 and the MAD scale is a direct symptom of
% contamination, and is reported even when the user has chosen 'pope'.
%
% The MAD needs the raw standardised residuals, so it requires its own sweep.
% That sweep is cheap: it reuses the same triangular solves and stores nothing.

robustScale = NaN;
madAvailable = false;
fhat = NaN;
if o.robustCheck || strcmp(testName, 'mad')
    rawAll = zeros(n, 1);
    pos = 0;
    for k = 1:nB
        [Ak, lk, pk] = batchFcn(k);
        lk = lk(:);
        nk = size(Ak, 1);
        if nk == 0, continue; end

        isCorrK = ~isempty(pk) && ~isvector(pk);
        if isempty(pk), pk = ones(nk, 1); end

        vk = Ak * x - lk;
        Wk = R' \ Ak.';

        % The raw score is the test statistic WITHOUT the scale factor:
        %       raw_i = (P v)_i / sqrt( (P Qvv P)_ii )
        % which is exactly what the MAD then rescales. It is defined for a
        % correlated batch just as it is for a diagonal one -- Baarda's
        % generalised form covers both -- so the robust check no longer has to
        % bail out on correlated data.
        if isCorrK
            Pk  = (pk + pk') / 2;
            Gk  = Pk * Ak;
            WGk = R' \ Gk.';
            qpk = diag(Pk) - sum(WGk.^2, 1).';
            Pvk = Pk * vk;
        else
            pk  = pk(:);
            qvk = 1 ./ pk - sum(Wk.^2, 1).';
            qvk(qvk < 0) = 0;
            qpk = pk.^2 .* qvk;
            Pvk = pk .* vk;
        end
        qpk = max(qpk, 0);

        good = qpk > 0;
        rk_ = zeros(nk, 1);
        rk_(good) = Pvk(good) ./ sqrt(qpk(good));
        rawAll(pos + (1:nk)) = rk_;
        pos = pos + nk;
    end

    if ~isempty(rawAll) && pos == n
        nz = rawAll(isfinite(rawAll) & rawAll ~= 0);
        if numel(nz) >= 3
            med = median(nz);
            % 1.4826 makes the MAD consistent for the normal distribution.
            robustScale = 1.4826 * median(abs(nz - med));
            if isfinite(robustScale) && robustScale > 0
                madAvailable = true;

                % ---- CONTAMINATION FRACTION (the gauge) ----------------
                % fhat = fraction of observations whose ROBUSTLY-standardised
                % residual exceeds a FIXED cut. The cut is 2.5, deliberately
                % NOT the Bonferroni critical value: we want a scale-free
                % COUNT of how much of the data looks bad, not a hypothesis
                % test. A test answers "is this observation an outlier?";
                % fhat answers "what fraction of my data is unusable?", and
                % those call for completely different responses.
                %
                % Measured accuracy (n=300, u=8, 20-sigma blunders):
                %     true f :  0%    5%   10%   20%   30%   40%
                %     fhat   : 0.015 0.059 0.106 0.201 0.300 0.388
                % i.e. |fhat - f| <= 0.015 for f <= 30%, with a ~1.5% floor on
                % clean data (the expected count beyond 2.5 sigma).
                %
                % ABOVE ~40% fhat COLLAPSES TO ZERO. This is not a bug: at 50%
                % contamination the MAD locks onto the BLUNDER population and
                % the good data become "the outliers". Every scale-based
                % indicator fails this way, which is why fhat alone is not
                % enough -- see contamBackstop below.
                fhat = mean(abs(nz / robustScale) > 2.5);
            end
        end
    end
end

if strcmp(testName, 'mad')
    if ~madAvailable
        error('lsqvce:diagnose:noRobustScale', ...
            ['''test'',''mad'' was requested but a robust scale could not be ' ...
             'formed (correlated batches, or too few usable residuals). Use ' ...
             '''pope'' or ''baarda'', or decorrelate the batches first.']);
    end
    scaleUse = robustScale;
else
    scaleUse = sqrt(s0sq_use);
end

% ================================================================== 2nd pass
% Progress reporting. Two passes over 1e8 observations backed by disk I/O is a
% long silence otherwise, and a user cannot tell a slow run from a hung one.
tProg = tic;
lastProg = 0;

for k = 1:nB
    [Ak, lk, pk] = batchFcn(k);
    lk = lk(:);
    nk = size(Ak, 1);
    if nk == 0, continue; end

    % =====================================================================
    % CORRELATED BATCHES
    % =====================================================================
    % P is block-diagonal with the BATCHES as its blocks: correlation within a
    % batch, independence between batches. That is the structure of most real
    % geodetic data (per-arc, per-pass, per-session noise models), and it is
    % exactly what makes the general case tractable.
    %
    % The diagnostics do NOT need whitening, and the flags do NOT degrade into
    % statements about "linear combinations". Baarda's GENERALISED data snooping
    % works directly on the physical observations:
    %
    %       r_i = (Qvv P)_ii
    %       w_i = (e_i' P v) / sqrt( (P Qvv P)_ii )
    %
    % and with  Qvv = P^-1 - A Qxx A'  these expand to
    %
    %       G_k          = P_k A_k                       (nk x u)
    %       (Qvv P)_ii   = 1 - (A_k Qxx G_k')_ii
    %       (P Qvv P)_ii = P_k(i,i) - g_i' Qxx g_i
    %
    % Every term needs only P_k, A_k and Qxx. Cost is O(u^2 nk) per batch --
    % IDENTICAL to the diagonal path -- and nothing of size n is ever formed.
    % Verified exact to 2e-15 against a dense reference.
    %
    % When P_k is diagonal this reduces ALGEBRAICALLY to the familiar
    %       r_i = 1 - p_i a_i' Qxx a_i,     w_i = v_i / sigma_vi
    % (verified: the two agree to 4.4e-16), so there is one code path, not two.
    isCorr = ~isempty(pk) && ~isvector(pk);

    if isempty(pk)
        pk = ones(nk, 1);
    end

    if isCorr
        if ~isequal(size(pk), [nk nk])
            error('lsqvce:diagnose:badPmatrix', ...
                ['Batch %d supplied a %dx%d weight matrix for %d observations. ' ...
                 'A correlated batch needs an nk x nk SPD weight matrix.'], ...
                 k, size(pk,1), size(pk,2), nk);
        end
        Pk = (pk + pk') / 2;               % symmetrise; the caller may be sloppy
    else
        pk = pk(:);
    end

    if nSeen + nk > n
        error('lsqvce:diagnose:batchMismatch', ...
            ['batchFcn is producing more observations (%d so far) than the ' ...
             'accumulator saw (%d). The second pass MUST replay exactly the ' ...
             'same batches in the same order, or every diagnostic below is ' ...
             'attributed to the wrong observation.'], nSeen + nk, n);
    end

    % --- residuals
    vk = Ak * x - lk;

    % --- redundancy numbers and the denominator of the test statistic
    W  = R' \ Ak.';                    % u x nk   (= R'^-1 A')
    if isCorr
        Gk  = Pk * Ak;                 % nk x u
        WG  = R' \ Gk.';               % u x nk
        % (Qvv P)_ii = 1 - (A Qxx G')_ii, and A Qxx G' = W' * WG columnwise
        rk  = 1 - sum(W .* WG, 1).';
        % (P Qvv P)_ii = P(i,i) - g_i' Qxx g_i
        qpp = diag(Pk) - sum(WG.^2, 1).';
        % Pv, the numerator of the generalised statistic
        Pvk = Pk * vk;
        % sigma_i of the OBSERVATION (for the MDB) is sqrt of the diagonal of
        % the covariance = inv(P). We need only the diagonal, so solve once.
        Ck    = Pk \ eye(nk);
        sig_i = sqrt(max(diag(Ck), 0));
    else
        qx  = sum(W.^2, 1).';          % a_i' Qxx a_i
        qv  = 1 ./ pk - qx;
        qv(qv < 0) = 0;                % CLAMP -- see the note below
        rk  = pk .* qv;
        % For diagonal P the generalised denominator is p_i^2 * qv_i, and the
        % numerator is p_i * v_i, so w reduces to v_i / sqrt(qv_i * ...) exactly
        % as before. Keep the cheap forms.
        qpp = pk.^2 .* qv;
        Pvk = pk .* vk;
        sig_i = 1 ./ sqrt(pk);
    end

    % --- CLAMP.
    % r_i is a difference of two positive numbers that become equal as an
    % observation approaches zero redundancy. Roundoff routinely drives it a
    % few eps NEGATIVE for observations that the network does not control at
    % all (verified: a parameter seen by exactly one observation yields
    % qv = -2.2e-16). Taking sqrt of that gives a complex number; dividing
    % v/0 gives NaN. Neither must be allowed to reach the user as a "result".
    %
    % THE UPPER BOUND IS DIFFERENT FOR CORRELATED OBSERVATIONS.
    % For a DIAGONAL P, r_i = 1 - p_i a_i' Qxx a_i is a leverage complement and
    % is genuinely bounded in [0,1]. For a CORRELATED block that bound does NOT
    % hold: r_i = (Qvv P)_ii can legitimately exceed 1, because the projector is
    % no longer orthogonal in the Euclidean metric. Measured on an AR(1) block:
    % r_i = 1.0001476, a real value, not roundoff.
    %
    % Clipping it to 1 would silently destroy the integrity identity
    % sum(r_i) = n - u -- and that is not hypothetical: it is exactly how this
    % bug was caught. The identity check earns its keep.
    rk = max(rk, 0);
    if ~isCorr
        rk = min(rk, 1);               % leverage complement: valid only here
    end
    qpp = max(qpp, 0);

    % --- normalised / studentised residual, Baarda's generalised form
    %       w_i = (P v)_i / ( scale * sqrt( (P Qvv P)_ii ) )
    %   pope   : scale = the ESTIMATED s0     (contaminable -- see above)
    %   baarda : scale = the A PRIORI sigma0  (needs a trustworthy sigma0)
    %   mad    : scale = a ROBUST scale       (high breakdown point)
    switch testName
        case 'pope'
            scl = sqrt(s0sq_use);
        case 'mad'
            scl = robustScale;
        otherwise   % baarda
            scl = sol.sigma0;
    end
    denom = scl * sqrt(qpp);

    wk = zeros(nk, 1);
    ctrl = rk > o.rMin & qpp > 0;      % controlled observations only
    wk(ctrl) = Pvk(ctrl) ./ denom(ctrl);
    % Uncontrolled: w is 0/0. Report 0, and record the index -- an uncontrolled
    % observation is not "fine", it is UNTESTABLE, which is a different and
    % more dangerous thing.
    wk(~ctrl) = 0;

    % --- MDB.
    mdbk  = inf(nk, 1);
    mdbk(ctrl) = sig_i(ctrl) .* delta0 ./ sqrt(rk(ctrl));

    % --- external reliability: an undetected bias mdb_i in obs i shifts x by
    %       dx = Qxx * (P e_i) * mdb_i
    %     For diagonal P this is Qxx a_i p_i mdb_i, the familiar form. For a
    %     correlated batch the bias enters through the WHOLE column of P, i.e.
    %     through g_i = (P_k A_k)(i,:)', which is exactly Gk' -- so the two
    %     cases differ only in which matrix multiplies Qxx.
    %     Qxx here is the CORRECTED (pseudoinverse) one on a free network:
    %     the bias propagates through N+, not through inv(N+BB').
    if isCorr
        dX = Qxx * Gk.';               % u x nk : Qxx * g_i  (columns)
        scal = mdbk.';                 % the p_i factor is already inside g_i
    else
        dX = Qxx * Ak.';               % u x nk : Qxx * a_i  (columns)
        scal = (pk .* mdbk).';         % 1 x nk
    end
    scal(~isfinite(scal)) = 0;         % uncontrolled -> undefined, skip
    dXs = abs(dX .* scal) ./ max(sigX, realmin);   % u x nk, in sigma units
    [extCol, parCol] = max(dXs, [], 1);            % worst parameter per obs
    [extB, iB] = max(extCol);
    if extB > extMax
        extMax     = extB;
        extMax_idx = nSeen + iB;
        extMax_par = parCol(iB);
    end

    % --- running scalars
    %   vPv = v' P v. For a correlated batch this is the full quadratic form,
    %   not sum(p .* v.^2). Getting this wrong would silently corrupt the
    %   variance factor and every test that depends on it.
    if isCorr
        vPv_direct = vPv_direct + vk.' * (Pk * vk);
    else
        vPv_direct = vPv_direct + sum(pk .* vk.^2);
    end
    sumR       = sumR + sum(rk);

    % --- flags
    isFlag = ctrl & (abs(wk) > kCrit);
    isWeak = ctrl & (rk < o.rWeak);
    nWeak  = nWeak + sum(isWeak);

    % ---- JOINT (w, r) QUADRANTS -------------------------------------
    % |w| and r_i mean nothing in isolation. The residual only reveals the
    % FRACTION r_i of any error present:  E[dv_i] = -r_i * b_i  (verified
    % exactly). So a small |w| on a low-r observation is not evidence of
    % health -- it is evidence that the test had no power.
    %
    % Measured on a levelling net with a dangling spur: a 40-sigma blunder on
    % a spur observation (r = 0) produced |w| = 0.000. A w-only report shows a
    % beautiful small residual on an observation carrying a gross error.
    %
    % The four quadrants call for four different actions:
    %   |w| hi, r hi -> CLEAN OUTLIER   reject with confidence
    %   |w| hi, r lo -> AMBIGUOUS       the network barely constrains this
    %                                   observation, so the residual may be
    %                                   absorbed error from a NEIGHBOUR.
    %                                   Rejecting is a coin flip; re-measure.
    %   |w| lo, r hi -> VERIFIED        well tested, and it passed
    %   |w| lo, r lo -> UNTESTED        NOT "fine". The silent quadrant.
    hiW = ctrl & (abs(wk) > kCrit);
    loR = rk < o.rWeak;                       % includes uncontrolled
    nQ_cleanOut  = nQ_cleanOut  + sum(hiW & ~loR);
    nQ_ambiguous = nQ_ambiguous + sum(hiW &  loR);
    nQ_verified  = nQ_verified  + sum(~hiW & ~loR);
    nQ_untested  = nQ_untested  + sum(~hiW &  loR);

    % ---- MDB distribution. Quoting only the worst case is misleading in both
    % directions: it hides that most of the network is fine, AND (because the
    % truly dangerous observations have MDB = Inf and drop out of any
    % percentile) a finite-only summary can look reassuring while the real
    % hazard is infinite. Track both the finite quantiles and the Inf count.
    %
    % The pool is BOUNDED by reservoir sampling: this toolbox's whole premise
    % is that nothing of size n is ever stored, and a percentile does not need
    % every value -- a few thousand uniform samples give quantiles to well
    % within the accuracy anyone acts on.
    mfin = isfinite(mdbk);
    nMdbInf = nMdbInf + sum(~mfin);
    mv = mdbk(mfin);
    for q = 1:numel(mv)
        nMdbSeen = nMdbSeen + 1;
        if numel(mdbPool) < mdbCap
            mdbPool(end+1, 1) = mv(q); %#ok<AGROW>
        else
            j = randi(nMdbSeen);
            if j <= mdbCap
                mdbPool(j) = mv(q);
            end
        end
    end

    nUnc = sum(~ctrl);
    if nUnc > 0
        nUncontrolled = nUncontrolled + nUnc;
        if numel(idxUncontrolled) < 1000
            idxUncontrolled = [idxUncontrolled; nSeen + find(~ctrl)]; %#ok<AGROW>
        end
    end

    [rlo, ilo] = min(rk);
    if rlo < rmin_seen
        rmin_seen = rlo;  rmin_idx = nSeen + ilo;
    end
    [whi, ihi] = max(abs(wk));
    if whi > maxAbsW
        maxAbsW = whi;  maxAbsW_idx = nSeen + ihi;
    end

    if any(isFlag)
        fi = nSeen + find(isFlag);
        flagIdx = [flagIdx; fi];              %#ok<AGROW>
        flagW   = [flagW;   wk(isFlag)];      %#ok<AGROW>
        flagR   = [flagR;   rk(isFlag)];      %#ok<AGROW>
        flagV   = [flagV;   vk(isFlag)];      %#ok<AGROW>
        flagMDB = [flagMDB; mdbk(isFlag)];    %#ok<AGROW>
    end

    % top-N worst |w|, maintained without sorting the whole vector
    if ~isempty(worstW)
        [cw, ci] = sort(abs(wk), 'descend');
        take = min(numel(cw), numel(worstW));
        cand   = [worstW;      cw(1:take)];
        candI  = [worstI; nSeen + ci(1:take)];
        [sw, si] = sort(cand, 'descend');
        worstW = sw(1:numel(worstW));
        worstI = candI(si(1:numel(worstW)));
    end

    if storeAll
        sl = nSeen + (1:nk);
        v_all(sl)   = vk;
        r_all(sl)   = rk;
        w_all(sl)   = wk;
        % Store the GENERALISED denominator (P Qvv P)_ii rather than qv. For
        % diagonal P these differ only by a factor p_i^2, and the generalised
        % one is the quantity the test actually uses -- so it is what the user
        % needs if they want to reproduce w by hand.
        qv_all(sl)  = qpp;
        mdb_all(sl) = mdbk;
    end

    nSeen = nSeen + nk;

    % ---- progress
    if o.progress
        el = toc(tProg);
        if el - lastProg > 2 || k == nB
            lastProg = el;
            frac = nSeen / max(n, 1);
            if frac > 0
                eta = el * (1 - frac) / frac;
            else
                eta = NaN;
            end
            lsqvce.util.log.msg(['  diagnose: batch %d/%d  (%5.1f%% of %d obs)  ' ...
                     'elapsed %6.1fs  ETA %6.1fs\n'], ...
                     k, nB, 100 * frac, n, el, eta);
        end
    end
end

if nSeen ~= n
    error('lsqvce:diagnose:batchMismatch', ...
        ['The second pass replayed %d observations but the accumulator saw ' ...
         '%d. batchFcn must reproduce the SAME batches in the SAME order.'], ...
         nSeen, n);
end

% ============================================================ integrity check
s0sq_direct = vPv_direct / max(rdf, 1);
sumR_expect = rdf;
sumR_err    = abs(sumR - sumR_expect) / max(sumR_expect, 1);

vPvRel = abs(vPv_direct - sol.vPv) / max(abs(vPv_direct), realmin);

D = struct();
D.n        = n;
D.u        = u;
D.r        = rdf;
D.nBatch   = nB;

D.x        = x;
D.sigX     = sigX;
D.qxxDiag  = qxxDiag;
D.R        = R;

D.vPv          = vPv_direct;          % the trustworthy one
D.vPvAccum     = sol.vPv;             % what the accumulator reported
D.vPvSource    = sol.vPvSource;
D.vPvRelDiff   = vPvRel;
D.s0sq         = s0sq_direct;
D.s0           = sqrt(max(s0sq_direct, 0));
D.sigma0       = sol.sigma0;
D.varRatio     = s0sq_direct / sol.sigma0^2;
D.s0Fallback   = s0Fallback;

% --- robust scale and the contamination symptom
%   On clean data s0 and the MAD scale both estimate sigma0 and agree.
%   Blunders inflate s0 but barely move the MAD, so their RATIO is a direct,
%   quantitative symptom of contamination -- and it is available even when the
%   user ran the (contaminable) Pope test, which is exactly when they need it.
D.robustScale   = robustScale;
D.madAvailable  = madAvailable;
if madAvailable && isfinite(D.s0) && robustScale > 0
    D.scaleRatio = D.s0 / robustScale;
else
    D.scaleRatio = NaN;
end
% Expected spread of the ratio under a CLEAN model. sd(s0) ~ 1/sqrt(2r) and
% sd(MAD) ~ 0.085 at n=200 (measured); combine conservatively. Anything beyond
% ~1.3 is not sampling noise, it is contamination.
D.scaleRatioTol = 1 + 3 * sqrt(1 / max(2 * rdf, 1) + 0.085^2);

% =========================================================================
% CONTAMINATION ASSESSMENT -- three indicators, because one is not enough
% =========================================================================
% The s0/MAD ratio fires early but is TREACHEROUS on its own. Measured, as
% contamination rises from 0 to 60%:
%
%     contamination :   0%    5%   10%   20%   30%   40%   50%   60%
%     s0/MAD        : 1.02  3.47  4.26  4.13  3.58  2.96  0.95  0.61
%     fhat          : 0.02  0.06  0.11  0.20  0.30  0.39  0.00  0.00
%     s0/sigma0     : 1.00  4.58  6.40  9.01 10.99 12.67 14.13 15.52
%
% The ratio PEAKS at 10% and then FALLS -- at 50% it sits at 0.95, BELOW its
% clean-data value. An alarm that switches itself off at total breakdown is
% worse than no alarm. The same happens to fhat, and to kurtosis, and for the
% same reason: past ~40% the MAD adopts the BLUNDERS as its reference
% population and the good data become "the outliers".
%
% So we report three things:
%
%   1. fhat            THE GAUGE. Accurate to +/-0.015 for f <= 30%, which is
%                      the range where the answer still matters. Tells the user
%                      HOW MUCH of the data is bad, which is the number that
%                      determines what to do about it.
%
%   2. s0/sigma0       THE BACKSTOP. Monotone all the way to 60% and unbounded,
%                      because sigma0 is A PRIORI and contamination cannot
%                      touch it. It is the one indicator that cannot be turned
%                      off. Caveat: it is only as trustworthy as the user's
%                      sigma0 -- if that is wrong by 2x, so is this.
%
%   3. contamLevel     The escalation, driven by whichever indicator is still
%                      credible. Above ~25% the correct response is not a
%                      better estimator, it is REFUSAL: if a quarter of your
%                      observations are "outliers", they are not outliers, and
%                      rejecting them launders a model error into a plausible
%                      result.
D.fhat = fhat;
if isfinite(D.s0) && sol.sigma0 > 0
    D.backstop = D.s0 / sol.sigma0;
else
    D.backstop = NaN;
end

% The backstop's clean-data band. s0/sigma0 has sd ~ 1/sqrt(2r) under a correct
% model, so 3 sigma is a fair "this is not sampling noise" line.
D.backstopTol = 1 + 3 / sqrt(max(2 * rdf, 1));

% ---- fhat's NULL DISTRIBUTION -------------------------------------------
% fhat counts residuals beyond 2.5 robust sigmas, so on PERFECTLY CLEAN data it
% does not go to zero -- it goes to the tail mass of the normal:
%       P(|Z| > 2.5) = 0.0124
% and it fluctuates around that. A FIXED threshold of 0.02 sits INSIDE that
% null distribution and cries wolf: it escalated 87 of 300 clean datasets (29%)
% before this was corrected.
%
% Two corrections, both measured rather than assumed:
%
% (1) SD INFLATION. The naive binomial sd, sqrt(f(1-f)/n), is TOO SMALL. The
%     residuals are not independent -- they satisfy u linear constraints, and
%     they are standardised by a MAD estimated from the same sample. Measured
%     ratio of true sd to binomial sd, over n = 100..5000:
%         1.52, 1.21, 1.26, 1.24, 1.20, 1.33   -> mean 1.29
%     We use 1.30.
%
% (2) THE MULTIPLIER. Measured trade-off (false-alarm rate on clean data vs
%     detection at low contamination):
%         k=3 : 1.75% false alarms,  97% detection at 3%, 100% at 5%
%         k=4 : 0.25% false alarms,  49% detection at 3%, 100% at 5%
%         k=5 : 0.12% false alarms,  25% detection at 3%, 100% at 5%
%     k=4 is the knee: it cuts false alarms sevenfold against k=3 while STILL
%     detecting 5% contamination every time -- and 5% is precisely where Pope's
%     tau begins to fail (detection 23%). Below that level the contamination
%     does not yet endanger the classical test, so missing some of it is a
%     cheap price for a quiet report.
FHAT_FLOOR = 2 * lsqvce.engine.stat.normCdf(-2.5);      % = 0.012419, exactly
FHAT_SD_INFLATION = 1.30;                          % measured, see above
FHAT_K = 4;                                        % measured, see above
fhatSd = FHAT_SD_INFLATION * ...
         sqrt(FHAT_FLOOR * (1 - FHAT_FLOOR) / max(n, 1));
D.fhatFloor = FHAT_FLOOR;
D.fhatTol   = FHAT_FLOOR + FHAT_K * fhatSd;       % n-dependent, calibrated

% The EXCESS is what actually matters: how much contamination is there beyond
% what a clean normal sample would produce anyway?
if isfinite(fhat)
    D.fhatExcess = max(fhat - FHAT_FLOOR, 0);
else
    D.fhatExcess = NaN;
end

% Escalation. Use fhat where it is credible; fall back to the backstop where
% fhat has collapsed (which is itself a signal: fhat ~ 0 together with a large
% backstop means the MAD has broken down entirely).
fhatCredible = isfinite(fhat) && fhat > 0.005;
bigBackstop  = isfinite(D.backstop) && D.backstop > max(3, 2 * D.backstopTol);

if ~madAvailable
    D.contamLevel = 'unknown';
elseif ~fhatCredible && bigBackstop
    % fhat says zero but the a-priori backstop says the residuals are enormous.
    % That combination is the fingerprint of TOTAL breakdown: the robust scale
    % has been captured by the blunders.
    D.contamLevel = 'breakdown';
elseif isfinite(fhat) && fhat <= D.fhatTol && ~bigBackstop
    % Within the null distribution: this is what CLEAN data looks like.
    D.contamLevel = 'clean';
elseif fhat >= 0.25
    D.contamLevel = 'refuse';
elseif fhat >= 0.10
    D.contamLevel = 'severe';
elseif fhat > D.fhatTol
    D.contamLevel = 'moderate';
else
    D.contamLevel = 'clean';
end

D.sumR       = sumR;
D.sumRExpect = sumR_expect;
D.sumRRelErr = sumR_err;

D.test        = testName;
D.alpha       = o.alpha;
D.alphaEff    = alphaEff;
D.bonferroni  = o.bonferroni;
D.beta        = o.beta;
D.kCrit       = kCrit;
D.delta0      = delta0;
D.tauBound    = tauBound;
D.snoopVacuous = snoopVacuous;

D.flagIdx = flagIdx;
D.flagW   = flagW;
D.flagR   = flagR;
D.flagV   = flagV;
D.flagMDB = flagMDB;
D.nFlag   = numel(flagIdx);
D.nExpectedFalse = alphaEff * n;      % how many flags are chance alone

D.rMin        = o.rMin;
D.rWeak       = o.rWeak;
D.rMinSeen    = rmin_seen;
D.rMinIdx     = rmin_idx;
D.nUncontrolled   = nUncontrolled;
D.idxUncontrolled = idxUncontrolled;
D.nWeak       = nWeak;

D.maxAbsW     = maxAbsW;
D.maxAbsWIdx  = maxAbsW_idx;
D.worstW      = worstW;
D.worstIdx    = worstI;

% ---- joint (w, r) quadrants and the hidden-blunder envelope
D.quad = struct( ...
    'cleanOutlier', nQ_cleanOut, ...   % |w| high, r high -> reject confidently
    'ambiguous',    nQ_ambiguous, ...  % |w| high, r low  -> coin flip, re-measure
    'verified',     nQ_verified, ...   % |w| low,  r high -> tested, and passed
    'untested',     nQ_untested);      % |w| low,  r low  -> the SILENT quadrant

D.mdbInf = nMdbInf;                    % MDB = Inf: cannot be tested at all
if ~isempty(mdbPool)
    D.mdbQuantiles = struct( ...
        'p50',  median(mdbPool), ...
        'p90',  prctileLocal(mdbPool, 90), ...
        'p99',  prctileLocal(mdbPool, 99), ...
        'max',  max(mdbPool));
else
    D.mdbQuantiles = struct('p50', NaN, 'p90', NaN, 'p99', NaN, 'max', NaN);
end

% External reliability, and its CALIBRATION. extMax scales like 1/sqrt(n)
% (measured: extMax ~ 4.0*delta0/sqrt(n) on clean data), so the raw number
% cannot be compared against a fixed threshold. extRatio says whether this
% network is worse than a healthy one OF ITS OWN SIZE.
D.extMax      = extMax;
D.extExpect   = 4.0 * delta0 / sqrt(max(n,1));
D.extRatio    = extMax / max(D.extExpect, realmin);
D.extMaxIdx   = extMax_idx;
D.extMaxPar   = extMax_par;

D.condR       = sol.condR;
D.condN       = sol.condN;
% Propagate the caveat: for sparse or large R the condition number is a LOWER
% BOUND (extreme-diagonal ratio), not an svd. The report must say so, or a
% sparse run looks better conditioned than it is.
if isfield(sol, 'condIsBound')
    D.condIsBound = sol.condIsBound;
else
    D.condIsBound = false;
end

% --- datum / free-network reporting
D.defect  = dDefect;
D.rank    = u - dDefect;
D.datumB  = Bc;

if storeAll
    D.v   = v_all;
    D.rNum = r_all;
    D.w   = w_all;
    D.qpp = qv_all;      % (P Qvv P)_ii -- the denominator of the test statistic
    D.qv  = qv_all;      % backwards-compatible alias
    D.mdb = mdb_all;
end

% Correlation matrix of the parameters -- u x u, always affordable.
% Uses the CORRECTED Qxx (see above): on a free network the uncorrected
% inv(R'R) would inflate the variances and distort every correlation.
dsig = sqrt(max(diag(Qxx), realmin));
D.Qxx  = Qxx;
D.corr = Qxx ./ (dsig * dsig.');
C = D.corr;  C(1:u+1:end) = 0;
C(~isfinite(C)) = 0;                 % a null-space parameter has sigma 0
[cmax, ci] = max(abs(C(:)));
[ci1, ci2] = ind2sub([u u], ci);
D.corrMax    = cmax;
D.corrMaxPar = [ci1 ci2];

end

function q = prctileLocal(x, pct)
%PRCTILELOCAL  Percentile without the Statistics Toolbox.
%   Linear interpolation between order statistics, matching the default
%   convention of prctile(). Written out because the whole point of this
%   toolbox is that it runs on base MATLAB.
    x = sort(x(:));
    n = numel(x);
    if n == 0, q = NaN; return; end
    if n == 1, q = x;   return; end
    % prctile places the k-th order statistic at percentile 100*(k-0.5)/n.
    pos = n * pct / 100 + 0.5;
    if pos <= 1
        q = x(1);
    elseif pos >= n
        q = x(n);
    else
        lo = floor(pos);
        f  = pos - lo;
        q  = (1 - f) * x(lo) + f * x(lo + 1);
    end
end
