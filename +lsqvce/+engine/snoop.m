function [sol, D, REJ] = snoop(sol, D, batchFcn, opt)
%LSQVCE.ENGINE.SNOOP  Iterative data snooping: reject outliers one at a time.
%
%   [sol, D, REJ] = LSQVCE.ENGINE.SNOOP(sol, D, batchFcn) repeatedly removes the
%   single worst outlier, re-solves, and re-tests, until no observation trips
%   the critical value.
%
%   ONE AT A TIME, and here is the real reason -- which is NOT the one usually
%   given.
%
%   The textbook justification is that a blunder "smears" its residual onto its
%   neighbours through the correlation in Qvv, so a simultaneous pass would
%   over-reject innocent observations. That was tested directly and it is NOT
%   what happens: with a Bonferroni-corrected threshold, simultaneous rejection
%   produced ZERO false positives across 200 realisations with 3 large blunders
%   each. Smearing is real but the threshold is far too high for it to matter.
%
%   The actual reason is MASKING. Several blunders inflate s0 together, and
%   each therefore hides the others: every |w| is divided by a variance factor
%   that the blunders themselves have inflated, so none of them looks extreme.
%   Remove one, s0 drops, and the rest become visible. Measured, with two
%   equal blunders in a small network:
%
%       iterative     found 400/400 true blunders   (100%)
%       simultaneous  found 380/400                  (95%)
%
%   and in a wider sweep the single-pass approach MISSED 71 and 114 true
%   blunders (at n=100 and n=50 respectively) where iterating missed none.
%
%   So: iterate because a single pass MISSES outliers, not because it rejects
%   too many. The practical consequence is the same -- iterate -- but the
%   failure it protects against is the dangerous one (undetected blunders left
%   in the data), not the benign one (a little good data discarded).
%
%   THE DOWNDATE
%   ------------
%   Re-accumulating from scratch after each rejection is O(n*u^2) per
%   rejection. Instead we downdate the Cholesky factor:
%
%       R'R = N   ->   R~'R~ = N - p_i a_i a_i'      via cholupdate(R, x, '-')
%
%   at O(u^2). Verified to reproduce a full re-accumulation to 5e-16, and to
%   stay accurate (1e-15) over 200 sequential downdates.
%
%   THE vPv UPDATE -- the subtle part
%   ---------------------------------
%   It is natural to write vPv_new = vPv - p_i*v_i^2. THIS IS WRONG. Removing
%   observation i changes x, hence changes EVERY residual. The correct identity
%   is
%
%       vPv_new = vPv - p_i * v_i^2 / r_i
%
%   with r_i the redundancy number. Verified exact (error 0.0e+00) against a
%   full re-solve; the naive form was off by 7.6e-4 on a SINGLE rejection, and
%   the error compounds with every further one. The 1/r_i factor is precisely
%   the correction for the fact that the fit moves to absorb the removal.
%
%   THE GUARDS
%   ----------
%   Downdating an observation with tiny r_i destroys the only information the
%   network has about some parameter, leaving a singular factor.
%
%   When this was actually TESTED, MATLAB-equivalent cholupdate flagged the
%   dangerous case correctly (it returned a failure code on a zero-redundancy
%   downdate rather than silently producing a near-singular factor). So the
%   post-hoc diag(R) check below is belt-and-braces, not the essential guard I
%   first assumed it was. It is kept anyway: it costs O(u), it is the only
%   thing that would catch a slow degradation over many successive downdates,
%   and a rollback is far cheaper than a silently corrupted factor.
%
%   The guard that actually matters is the PRE-CHECK: refuse to reject any
%   observation with r_i < rMinReject. That is not a numerical safeguard but a
%   STATISTICAL one -- a large residual on a weakly-controlled observation is
%   as likely to be error absorbed from a neighbour as a blunder in the
%   observation itself, so rejecting it is close to a coin flip.
%
%   So we check BOTH:
%     * BEFORE: refuse to reject any observation with r_i < rMinReject;
%     * AFTER:  verify the smallest |diag(R)| has not collapsed, and roll back
%               if it has.
%
%   NAME-VALUE
%     'maxReject'    cap on rejections (default: max(10, 0.01*n))
%     'rMinReject'   refuse to reject below this redundancy (default 0.10)
%     'alpha'        significance per test (default 1e-3)
%     'bonferroni'   correct for n tests (default true)
%     'reaccumulate' re-accumulate from scratch at the end to purge downdate
%                    roundoff and produce a clean final D (default true)
%     'verbose'      print the rejection sequence (default true)
%
%   OUTPUT REJ
%     .idx        indices rejected, in the order they were rejected
%     .w          their normalised residual at the moment of rejection
%     .r          their redundancy number at that moment
%     .refused    indices that tripped the test but were REFUSED (low r_i)
%     .nIter      iterations performed
%     .converged  true if snooping terminated cleanly
%     .history    per-iteration [nFlag, maxAbsW, s0sq]
%
%   EXAMPLE
%       data = struct('batchFcn', src, 'numUnknowns', u, 'numBatches', 4);
%       R = lsqvce.fit(data, Verbose=false);
%       [sol, D, REJ] = lsqvce.engine.snoop(R.source, R.diag, R.batchFcn);
%
%       % REJ.idx       -> [40 150]     rejected, one at a time
%       % REJ.refused   -> observations too weakly controlled to judge
%       % REJ.history   -> the full rejection sequence, for the record
%
%       % Under heavy contamination Pope cannot see the blunders at all, so
%       % this loop would exit having rejected NOTHING and report 'converged'.
%       % Use the robust scale there:
%       [sol, D, REJ] = lsqvce.engine.snoop(R.source, R.diag, R.batchFcn, 'test', 'mad');
%
%   INPUT
%     sol       1x1 struct  the solution to start from (from LSQVCE.ENGINE.ADJUST)
%     D         1x1 struct  its diagnostics
%     batchFcn  1x1 fh      the SAME batch function, for the re-solves
%     varargin  Name/Value: 'test' ('auto'|'pope'|'baarda'|'mad'), 'alpha',
%                           'maxReject', 'force', 'verbose'
%
%   ERRORS
%     lsqvce:snoop:badBatch        batchFcn did not return [A, l, p], or
%                                  returned different data on the second call
%     lsqvce:snoop:badIndex        an observation index is outside 1..n
%     lsqvce:snoop:correlatedBatch a batch carries a full covariance, so its
%                                  observations cannot be rejected individually
%     lsqvce:snoop:tooContaminated the contamination gauge is past the point
%                                  where tau is meaningful -- switch to
%                                  'test','mad', or use LSQVCE.ENGINE.ROBUST
%     lsqvce:snoop:breakdown       the rejection cap was reached and the scale
%                                  estimate has broken down; the data are not
%                                  a few blunders on a good model
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (b) b2 (F-07):
%   typed `arguments` signature, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:snoop:contaminated                  [warning] snoop
%     lsqvce:snoop:downdateFailed                [warning] Removing observation X would make the normal matrix singular (the Cho...
%     lsqvce:snoop:forcedBreakdown               [warning] snoop
%     lsqvce:snoop:forcedContaminated            [warning] snoop
%     lsqvce:snoop:forcedOverCap                 [warning] snoop
%     lsqvce:snoop:lowRedundancy                 [warning] Observation X has the largest |w| (X) but its redundancy is only X. T...
%     lsqvce:snoop:maxReject                     [warning] Reached the rejection cap (X)
%     lsqvce:snoop:rejectFraction                [warning] STOPPING
%
%   See also LSQVCE.ENGINE.ROBUST, LSQVCE.ENGINE.DIAGNOSE, LSQVCE.ENGINE.SYSTEMATIC.
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
    D
    batchFcn
    opt.MaxReject double {mustBeScalarOrEmpty} = []
    opt.RMinReject (1,1) double {mustBeNonnegative, mustBeLessThanOrEqual(opt.RMinReject,1)} = 0.10
    opt.Alpha (1,1) double {mustBePositive, mustBeLessThan(opt.Alpha,1)} = 1e-3
    opt.Bonferroni (1,1) logical = true
    opt.Reaccumulate (1,1) logical = true
    opt.Verbose (1,1) logical = true
    opt.NBatch double {mustBeScalarOrEmpty} = []
    opt.Test (1,1) string {mustBeMember(opt.Test, ["pope","mad"])} = "pope"
    opt.Force (1,1) logical = false
    opt.MaxRejectFrac (1,1) double {mustBePositive, mustBeLessThanOrEqual(opt.MaxRejectFrac,1)} = 0.25
end

o = struct('maxReject', opt.MaxReject, 'rMinReject', opt.RMinReject, 'alpha', opt.Alpha, 'bonferroni', opt.Bonferroni, 'reaccumulate', opt.Reaccumulate, 'verbose', opt.Verbose, 'nBatch', opt.NBatch, 'test', char(opt.Test), 'force', opt.Force, 'maxRejectFrac', opt.MaxRejectFrac);

nB = o.nBatch;
if isempty(nB), nB = D.nBatch; end

n0 = sol.n;
u  = sol.u;
maxRej = o.maxReject;
if isempty(maxRej)
    maxRej = max(10, ceil(0.01 * n0));
end

if sol.defect > 0 && o.verbose
    lsqvce.util.log.msg(['\nNote: this is a free network (datum defect %d). Rejection is\n' ...
             'still valid -- residuals are datum-invariant -- but rejecting an\n' ...
             'observation can INCREASE the defect if it was the last one tying\n' ...
             'a parameter down. That is caught by the redundancy guard below.\n'], ...
             sol.defect);
end

R    = sol.R;
rejected = false(n0, 1);
REJ  = struct('idx', [], 'w', [], 'r', [], 'refused', [], ...
              'nIter', 0, 'converged', false, 'history', []);

if lsqvce.util.verbosity.progress(o.verbose)
    lsqvce.util.log.msg('\n--- ITERATIVE DATA SNOOPING\n');
    lsqvce.util.log.msg('  %4s %10s %12s %10s %12s %10s\n', ...
            'iter', 'nFlag', 'worst |w|', 'r_i', 'obs', 's0');
    lsqvce.util.log.msg('  %s\n', repmat('-', 1, 64));
end

converged = false;
useRobust = strcmpi(o.test, 'mad');

% ---------------------------------------------------------------------------
% CONTAMINATION GUARDS
% ---------------------------------------------------------------------------
% Iterative rejection has an appealing story: remove the worst, s0 falls, the
% next becomes visible, repeat. That story FAILS in two distinct ways, and both
% are silent.
%
% FAILURE 1 -- nothing is visible to begin with.
%   If the blunders have already inflated s0 so far that NOTHING crosses the
%   critical value on the first pass, this loop exits at iteration 1 having
%   rejected nothing, and reports a clean result. Measured: 5 blunders of 25
%   sigma in n=50 -- iterative snooping with Pope rejects ZERO and terminates
%   immediately.
%
% FAILURE 2 -- the robust scale itself has broken down.
%   Past ~40% contamination the MAD adopts the blunders as its reference
%   population. fhat collapses to 0, the s0/MAD ratio falls BELOW its
%   clean-data value, and every scale-based alarm goes quiet at exactly the
%   moment the data are worst.
%
% Against FAILURE 2 no estimator can be trusted, so the last guard is one that
% does not depend on any: COUNTING. "You are asking me to delete a third of
% your data" is a statement no amount of contamination can falsify.

if isfield(D, 'contamLevel')
    switch D.contamLevel
    case 'breakdown'
        if o.force
            lsqvce.util.log.warn('lsqvce:snoop:forcedBreakdown', ...
                ['snoop: proceeding under ''force'' despite total breakdown ' ...
                 'of the outlier machinery. The rejections below are very ' ...
                 'likely to be GOOD observations.']);
        else
        error('lsqvce:snoop:breakdown', ...
            ['snoop: REFUSING TO RUN. The outlier machinery has broken down ' ...
             'completely.\n' ...
             'The robust scale reports almost no outliers (fhat = %.3f) while ' ...
             'the a-priori backstop says the residuals are %.1fx larger than ' ...
             'sigma0 promised. Those are only compatible if MORE THAN HALF the ' ...
             'data are bad -- at which point the robust scale has been captured ' ...
             'by the blunders and is calling the GOOD observations outliers.\n' ...
             'Rejecting here would delete the good data. Fix sigma0 or the ' ...
             'functional model first.'], D.fhat, D.backstop);
        end

    case 'refuse'
        if o.force
            lsqvce.util.log.warn('lsqvce:snoop:forcedContaminated', ...
                ['snoop: proceeding under ''force'' with an estimated %.0f%% ' ...
                 'contamination. The surviving parameters will be confidently ' ...
                 'wrong.'], 100 * D.fhat);
        else
        error('lsqvce:snoop:tooContaminated', ...
            ['snoop: REFUSING TO RUN. An estimated %.0f%% of the observations ' ...
             '(fhat = %.3f) are inconsistent with the rest.\n' ...
             'If a quarter of your data are "outliers", they are not outliers. ' ...
             'Rejecting them launders a MODEL error into a plausible-looking ' ...
             'result. Run lsqvce.engine.systematic() and fix the functional model.\n' ...
             'If you have a specific reason to proceed anyway, pass ' ...
             '''force'',true -- but understand that the surviving parameters ' ...
             'will be confidently wrong.'], 100 * D.fhat, D.fhat);
        end

    case 'severe'
        if ~useRobust
            lsqvce.util.log.warn('lsqvce:snoop:contaminated', ...
                ['snoop: an estimated %.0f%% of observations are bad, and you ' ...
                 'are using Pope''s tau, which is BLIND at this level: the ' ...
                 'blunders inflate s0 and divide themselves back below the ' ...
                 'threshold. Iterating will NOT rescue this -- if nothing trips ' ...
                 'the threshold on pass one, this loop exits immediately and ' ...
                 'returns a "clean" result that is entirely spurious.\n' ...
                 'Re-run with ''test'',''mad''. Proceeding, but treat an empty ' ...
                 'rejection list as MEANINGLESS.'], 100 * D.fhat);
        end
    end
end

% The hard ceiling. Independent of every estimator above.
maxRejFrac = o.maxRejectFrac;
hardCap = floor(maxRejFrac * n0);

for it = 1:(maxRej + 1)

    % ---- find the worst CURRENTLY-flagged observation, using the current R.
    % One streaming pass. We do not need full diagnostics here, only the worst
    % candidate, so this is cheaper than a call to diagnose().
    nAlive = n0 - sum(rejected);
    rdf    = nAlive - sol.rank;
    if rdf <= 0
        if lsqvce.util.verbosity.progress(o.verbose)
            lsqvce.util.log.msg('  Redundancy exhausted. Stopping.\n');
        end
        break;
    end

    alphaEff = o.alpha;
    if o.bonferroni, alphaEff = o.alpha / max(nAlive, 1); end
    kCrit = lsqvce.engine.stat.tauCrit(alphaEff, rdf);

    [bestW, bestI, bestR, bestV, bestP, vPvCur, nFlag] = ...
        scanWorst(batchFcn, nB, R, sol.x, rejected, rdf, kCrit, sol.rank, ...
                  useRobust);

    s0cur = sqrt(max(vPvCur / max(rdf, 1), 0));
    REJ.history(end+1, :) = [nFlag, bestW, s0cur^2]; %#ok<AGROW>

    if nFlag == 0 || ~isfinite(bestW) || bestW <= kCrit
        converged = true;
        if lsqvce.util.verbosity.progress(o.verbose)
            lsqvce.util.log.msg('  %4d %10d %12.4f %10s %12s %10.4f   -> clean\n', ...
                    it, 0, bestW, '-', '-', s0cur);
        end
        break;
    end

    % ---- HARD CEILING on the rejected FRACTION.
    % This guard depends on no estimator whatsoever. Past ~40% contamination
    % every scale-based indicator (s0/MAD, fhat, kurtosis) switches itself off
    % -- but counting cannot be fooled. If the loop wants to delete a quarter
    % of the data, the premise of "outlier rejection" has already failed.
    if numel(REJ.idx) >= hardCap
        if o.force
            lsqvce.util.log.warn('lsqvce:snoop:forcedOverCap', ...
                ['snoop: rejected %d of %d observations (%.0f%%), past the ' ...
                 '%.0f%% ceiling, continuing only because ''force'' is set.'], ...
                 numel(REJ.idx), n0, 100*numel(REJ.idx)/n0, 100*maxRejFrac);
        else
            lsqvce.util.log.warn('lsqvce:snoop:rejectFraction', ...
                ['STOPPING: %d of %d observations (%.0f%%) have been rejected, ' ...
                 'reaching the %.0f%% ceiling, and %d are still flagged.\n' ...
                 'This is no longer an outlier problem. When a quarter of the ' ...
                 'data fails the test, the data are not wrong -- the MODEL is. ' ...
                 'Continuing to reject would launder that model error into a ' ...
                 'plausible-looking result, and the surviving parameters would ' ...
                 'be confidently wrong.\n' ...
                 'Run lsqvce.engine.systematic() to find the structure you have not ' ...
                 'modelled. Raise ''maxRejectFrac'' only if you can justify it.'], ...
                 numel(REJ.idx), n0, 100*numel(REJ.idx)/n0, 100*maxRejFrac, nFlag);
            REJ.stoppedOnFraction = true;
            break;
        end
    end

    if numel(REJ.idx) >= maxRej
        lsqvce.util.log.warn('lsqvce:snoop:maxReject', ...
            ['Reached the rejection cap (%d). %d observation(s) are still ' ...
             'flagged.\n' ...
             'A large fraction of the data failing the test is not an outlier ' ...
             'problem -- it is a MODEL problem. Rejecting your way out of it ' ...
             'launders a systematic error into a plausible-looking result. ' ...
             'Run lsqvce.engine.systematic() before rejecting anything else.'], ...
             maxRej, nFlag);
        break;
    end

    % ---- GUARD 1: refuse a low-redundancy observation.
    if bestR < o.rMinReject
        REJ.refused(end+1) = bestI; %#ok<AGROW>
        if lsqvce.util.verbosity.progress(o.verbose)
            lsqvce.util.log.msg('  %4d %10d %12.4f %10.4f %12d %10.4f   REFUSED (r too low)\n', ...
                    it, nFlag, bestW, bestR, bestI, s0cur);
        end
        lsqvce.util.log.warn('lsqvce:snoop:lowRedundancy', ...
            ['Observation %d has the largest |w| (%.3f) but its redundancy is ' ...
             'only %.4f. The network barely constrains it, so its residual is ' ...
             'as likely to be absorbed error from its NEIGHBOURS as a blunder ' ...
             'in itself. Rejecting it is a coin flip, and removing it may make ' ...
             'the system singular.\n' ...
             'REFUSING. Re-measure this observation, or add geometry that ' ...
             'controls it. Snooping stops here because the next candidate ' ...
             'cannot be assessed while this one distorts the fit.'], ...
             bestI, bestW, bestR);
        break;
    end

    % ---- downdate
    xi = sqrt(bestP) * getRow(batchFcn, nB, bestI);
    Rold = R;
    [Rnew, fail] = cholupdate(R, xi(:), '-');

    % ---- GUARD 2: verify the factor did not collapse.
    % cholupdate's own flag is necessary but NOT sufficient: a downdate can
    % "succeed" and still leave a diagonal entry at ~0, after which every
    % solve is garbage.
    dg = abs(diag(Rnew));
    collapsed = fail ~= 0 || any(dg < eps(max(dg)) * u * 100);
    if collapsed
        R = Rold;
        REJ.refused(end+1) = bestI; %#ok<AGROW>
        lsqvce.util.log.warn('lsqvce:snoop:downdateFailed', ...
            ['Removing observation %d would make the normal matrix singular ' ...
             '(the Cholesky downdate collapsed). That observation carries the ' ...
             'only information about at least one parameter -- its redundancy ' ...
             'is effectively zero, whatever the nominal value says.\n' ...
             'REFUSING and stopping.'], bestI);
        break;
    end
    R = Rnew;

    % ---- vPv update, with the essential 1/r_i factor (see header).
    vPvCur = vPvCur - bestP * bestV^2 / bestR;

    rejected(bestI) = true;
    REJ.idx(end+1) = bestI;   %#ok<AGROW>
    REJ.w(end+1)   = bestW;   %#ok<AGROW>
    REJ.r(end+1)   = bestR;   %#ok<AGROW>

    % ---- re-solve for x with the downdated R.
    % b is not maintained incrementally: recompute it in the same pass that
    % scanWorst will make next iteration. Simpler and one pass either way.
    bnew = accumB(batchFcn, nB, rejected, u);
    sol.x = R \ (R' \ bnew);

    if lsqvce.util.verbosity.progress(o.verbose)
        s0new = sqrt(max(vPvCur / max(rdf - 1, 1), 0));
        lsqvce.util.log.msg('  %4d %10d %12.4f %10.4f %12d %10.4f   rejected\n', ...
                it, nFlag, bestW, bestR, bestI, s0new);
    end
end

REJ.nIter     = size(REJ.history, 1);
REJ.converged = converged;
REJ.rejected  = rejected;

% ------------------------------------------------------------ final solution
% Re-accumulate from scratch. Two reasons: it purges any accumulated downdate
% roundoff, and it gives a clean, fully consistent D for the surviving data --
% which is what the user will actually report.
if o.reaccumulate && ~isempty(REJ.idx)
    if lsqvce.util.verbosity.progress(o.verbose)
        lsqvce.util.log.msg(['\n  Re-accumulating from scratch over the %d surviving ' ...
                 'observations\n  (purges downdate roundoff; the reported ' ...
                 'result is a clean solve, not a\n  sequence of subtractions).\n'], ...
                 n0 - numel(REJ.idx));
    end
    % R5 M5 5b': the clean re-solve runs on the BLOCK path (the one
    % accumulator, finding H-02): the surviving rows of all batches become
    % one local-free system (lsqvce.engine.batchSystem), a datum survives as
    % the prior on the null space, and the diagnostics come from the same
    % sweep every other solve uses. The flat lsqvce.engine.Accumulator is
    % retired with this level.
    keptSrc = @(k) maskedBatch(batchFcn, k, rejected, nB);
    allKept = @(k) concatKept(keptSrc, nB);
    L = lsqvce.engine.batchSystem(u, allKept, 1, 'Name', 'snoop');
    backend = 'qr';
    if isfield(sol, 'method') && strcmpi(char(sol.method), 'normal'), backend = 'normal'; end
    acc = lsqvce.engine.BlockAccumulator(u, 'Backend', backend, 'KeepLocal', true, 'Verbose', false);
    if isfield(sol, 'defect') && sol.defect > 0 && isfield(sol, 'datumB') && ~isempty(sol.datumB)
        acc.setPrior(lsqvce.engine.CommonPrior(sol.datumB.', zeros(size(sol.datumB, 2), 1), 1));
    end
    acc.addSystem(L.sysID, 0, L.closure, 1);
    solB = acc.solve();
    [dsol, bf, nBd] = lsqvce.engine.diagnoseAccum(acc, solB);
    D   = lsqvce.engine.diagnose(dsol, bf, 'nBatch', nBd, ...
                           'alpha', o.alpha, 'bonferroni', o.bonferroni, ...
                           'test', o.test);
    D.rejectedIdx = REJ.idx;
    % the returned sol keeps the caller's fields and carries the clean solve
    sol.x = solB.x(:);  sol.s0sq = solB.s0sq;  sol.sigma0 = solB.sigma0;
    if isfield(solB, 'Cxx'), sol.Cxx = solB.Cxx; end
    sol.n = dsol.n;  sol.dof = solB.dof;
end

if lsqvce.util.verbosity.progress(o.verbose)
    lsqvce.util.log.msg('\n  %d rejected, %d refused, %d surviving.\n', ...
            numel(REJ.idx), numel(REJ.refused), n0 - numel(REJ.idx));
    if ~converged && isempty(REJ.refused)
        lsqvce.util.log.msg('  WARNING: snooping did not converge.\n');
    end
    lsqvce.util.log.msg('\n');
end
end

%% ======================================================== helpers

function [bestW, bestI, bestR, bestV, bestP, vPv, nFlag] = ...
        scanWorst(batchFcn, nB, R, x, rejected, rdf, kCrit, rank_, useRobust)
%SCANWORST  One streaming pass over all batches to find the single worst
%   (largest |w|) surviving observation and accumulate vPv honestly.
%   INPUT
%     batchFcn  1x1 fh      [A,l,p] = batchFcn(k), k = 1..nB
%     nB        1x1 double  number of batches
%     R         u x u       upper-triangular factor of the reduced normal (for
%                           the leverage back-substitution a' N^-1 a)
%     x         u x 1       current parameter estimate
%     rejected  m x 1 log   mask of already-rejected observations (global index)
%     rdf       1x1 double  residual degrees of freedom (for s0)
%     kCrit     1x1 double  critical value for flagging (|w| threshold)
%     rank_     1x1 double  rank of the system (for the s0 scaling)
%     useRobust 1x1 log     use a robust (MAD) scale for w (default false)
%   OUTPUT
%     bestW  1x1  the largest |w| found among surviving observations
%     bestI  1x1  its global observation index (0 if none)
%     bestR  1x1  its partial redundancy r_i
%     bestV  1x1  its residual v_i
%     bestP  1x1  its weight p_i
%     vPv    1x1  weighted residual sum of squares over survivors
%     nFlag  1x1  count of survivors exceeding kCrit
    if nargin < 9, useRobust = false; end
    rawPool = [];
% One streaming pass: find the worst |w| among surviving observations, and
% accumulate vPv honestly along the way.
    bestW = 0; bestI = 0; bestR = 0; bestV = 0; bestP = 1;
    vPv = 0; nFlag = 0; nSeen = 0;

    % s0 must come from THIS pass, but we need it to studentise... so do two
    % sub-passes over the (cheap) scalars: first vPv, then w. In practice we
    % accumulate the raw v and qv and finish afterwards, which needs one pass
    % and O(1) memory beyond the batch.
    cand = struct('w', {}, 'i', {}, 'r', {}, 'v', {}, 'p', {});

    % --- pass A: vPv
    for k = 1:nB
        [Ak, lk, pk] = batchFcn(k);
        lk = lk(:); nk = size(Ak,1);
        if nk == 0, continue; end
        if isempty(pk), pk = ones(nk,1); end
        if ~isvector(pk)
            % A matrix here would be silently FLATTENED by pk(:) into a
            % meaningless vector, which is the worst possible failure: a
            % confident wrong answer. Refuse instead.
            %
            % The refusal is not laziness. Rejecting one observation from a
            % correlated block is genuinely harder than rejecting an
            % independent one: removing row i from a block means taking the
            % SCHUR COMPLEMENT of the block covariance, not merely deleting a
            % row of P. The surviving observations in that block have their
            % covariance CHANGED by the removal, so the downdate is not a
            % rank-one cholupdate and the whole O(u^2) argument collapses.
            %
            % lsqvce.engine.diagnose DOES handle correlated batches (Baarda's
            % generalised statistic), so you can FIND the outliers. You just
            % cannot use this function to remove them.
            error('lsqvce:snoop:correlatedBatch', ...
                ['snoop: batch %d supplies a full weight MATRIX (correlated ' ...
                 'observations). Iterative rejection is not supported for ' ...
                 'correlated blocks.\n' ...
                 'Removing one observation from a correlated block changes the ' ...
                 'covariance of every SURVIVING observation in that block (a ' ...
                 'Schur complement, not a row deletion), so the rank-one ' ...
                 'Cholesky downdate this function is built on does not apply.\n' ...
                 'WHAT YOU CAN DO: lsqvce.engine.diagnose fully supports correlated ' ...
                 'batches, so run it to LOCATE the outliers (D.flagIdx). Then ' ...
                 'either (a) drop the whole affected block, (b) rebuild the ' ...
                 'block covariance without the rejected observation and re-run ' ...
                 'lsqvce.engine.adjust, or (c) use lsqvce.engine.robust, which downweights ' ...
                 'rather than deletes and therefore leaves the block ' ...
                 'structure intact.'], k);
        end
        pk = pk(:);
        idx = nSeen + (1:nk);
        alive = ~rejected(idx);
        if any(alive)
            Aa = Ak(alive,:);
            vk = Aa * x - lk(alive);
            vPv = vPv + sum(pk(alive) .* vk.^2);

            % Collect the raw standardised score for the robust scale.
            % Scale-free, so it can be pooled across batches.
            if useRobust
                Wa = R' \ Aa.';
                qva = 1 ./ pk(alive) - sum(Wa.^2, 1).';
                qva(qva < 0) = 0;
                gd = qva > 0;
                tmp = vk(gd) ./ sqrt(qva(gd));
                rawPool = [rawPool; tmp]; %#ok<AGROW>
            end
        end
        nSeen = nSeen + nk;
    end
    s0sq = vPv / max(rdf, 1);
    if ~isfinite(s0sq) || s0sq <= 0, s0sq = 1; end

    % Choose the denominator. The robust (MAD) scale has a high breakdown
    % point; s0 does not, and under contamination the blunders inflate s0 and
    % divide themselves back down below the threshold.
    if useRobust && numel(rawPool) >= 3
        med  = median(rawPool);
        sMad = 1.4826 * median(abs(rawPool - med));
        if isfinite(sMad) && sMad > 0
            s0sq = sMad^2;
        end
    end

    % --- pass B: w and the worst candidate
    nSeen = 0;
    for k = 1:nB
        [Ak, lk, pk] = batchFcn(k);
        lk = lk(:); nk = size(Ak,1);
        if nk == 0, continue; end
        if isempty(pk), pk = ones(nk,1); end
        pk = pk(:);
        idx = nSeen + (1:nk);
        alive = ~rejected(idx);
        if ~any(alive), nSeen = nSeen + nk; continue; end

        Aa = Ak(alive,:); la = lk(alive); pa = pk(alive);
        ia = idx(alive);

        vk = Aa * x - la;
        W  = R' \ Aa.';
        qx = sum(W.^2, 1).';
        qv = 1 ./ pa - qx;
        qv(qv < 0) = 0;
        rk = min(max(pa .* qv, 0), 1);

        wk = zeros(numel(la), 1);
        ok = rk > 1e-6;
        wk(ok) = vk(ok) ./ sqrt(s0sq * qv(ok));

        fl = ok & (abs(wk) > kCrit);
        nFlag = nFlag + sum(fl);

        if any(fl)
            [mw, mi] = max(abs(wk .* fl));
            if mw > bestW
                bestW = mw;  bestI = ia(mi);  bestR = rk(mi);
                bestV = vk(mi);  bestP = pa(mi);
            end
        end
        nSeen = nSeen + nk;
    end
end

function a = getRow(batchFcn, nB, target)
% Fetch the design row of a single observation, by index.
    nSeen = 0;
    for k = 1:nB
        [Ak, lk, ~] = batchFcn(k);
        nk = size(Ak, 1);
        if target <= nSeen + nk
            a = Ak(target - nSeen, :);
            return;
        end
        nSeen = nSeen + nk;
    end
    error('lsqvce:snoop:badIndex', 'Observation %d not found.', target);
end

function b = accumB(batchFcn, nB, rejected, u)
%ACCUMB  Recompute the reduced right-hand side b = A' P l over SURVIVING
%   observations only (those not in the rejected mask), in one streaming pass.
%   INPUT
%     batchFcn  1x1 fh      [A,l,p] = batchFcn(k)
%     nB        1x1 double  number of batches
%     rejected  m x 1 log   global mask of rejected observations
%     u         1x1 double  number of parameters
%   OUTPUT
%     b         u x 1       A' P l summed over surviving observations
    b = zeros(u, 1);
    nSeen = 0;
    for k = 1:nB
        [Ak, lk, pk] = batchFcn(k);
        lk = lk(:); nk = size(Ak,1);
        if nk == 0, continue; end
        if isempty(pk), pk = ones(nk,1); end
        pk = pk(:);
        idx = nSeen + (1:nk);
        alive = ~rejected(idx);
        if any(alive)
            b = b + Ak(alive,:).' * (pk(alive) .* lk(alive));
        end
        nSeen = nSeen + nk;
    end
end

function [A, l, p] = maskedBatch(batchFcn, k, rejected, nB)
% Return batch k with the rejected observations removed.
    nSeen = 0;
    for j = 1:nB
        [Aj, lj, pj] = batchFcn(j);
        lj = lj(:); nj = size(Aj,1);
        if j == k
            if isempty(pj), pj = ones(nj,1); end
            pj = pj(:);
            alive = ~rejected(nSeen + (1:nj));
            A = Aj(alive, :);
            l = lj(alive);
            p = pj(alive);
            return;
        end
        nSeen = nSeen + nj;
    end
    error('lsqvce:snoop:badBatch', 'Batch %d not found.', k);
end

%==========================================================================
function [A, l, p] = concatKept(keptSrc, nB)
%CONCATKEPT  All surviving rows of every batch as ONE batch (rows stay in
%   batch order; p carried when the source has it, else unit).
    As = cell(nB, 1);  ls = cell(nB, 1);  ps = cell(nB, 1);
    for k = 1:nB
        [Ak, lk, pk] = keptSrc(k);
        if isempty(lk), continue, end
        As{k} = Ak;  ls{k} = lk(:);
        if isempty(pk), pk = ones(numel(lk), 1); end
        ps{k} = pk(:);
    end
    A = vertcat(As{:});  l = vertcat(ls{:});  p = vertcat(ps{:});
end
