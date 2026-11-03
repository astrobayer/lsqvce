classdef verdict
%LSQVCE.ENGINE.VERDICT  The verdict of an adjustment: diagnostics, snooping, systematics, tests.
%
%   The ONE engine unit that judges a solved adjustment (R5 M4 WP-86, the
%   "verdict from the sweep" row; part 4). It absorbs the seven units that
%   did it before -- diagnoseAccum, diagnose, snoop, systematic, ftest,
%   ellipse and the statistics half of stat (the distributions moved to
%   LSQVCE.KERNEL.STAT) -- as static methods of one class. Users reach it
%   through lsqvce.Results: R.diag is its diagnose, R.systematic / R.snoop
%   / R.reliability / R.ellipse / R.compare call it, and LSQVCE.REPORT prints
%   it. The report never names an engine function.
%
%   METHODS (static)
%     [D, dsol, src, nB] = fromAccumulator(acc, sol, Name, Value, ...)
%         the five-level verdict of a solved block accumulator in one call:
%         source + diagnose (what every driver runs; Name-Value forwarded)
%     sol = supersede(sol, D, acc)
%         a refused (cancelled) v'Pv replaced by D's residual pass: vtv,
%         s0sq, sigma0 and Cxx (what lsqvce.Results stores; WP-149)
%     [dsol, src, nB] = source(acc, sol)
%         the full-parameter solution (its cofactor as blocks,
%         LSQVCE.ENGINE.BLOCKCOFACTOR) and the replayable batch source
%     D = diagnose(dsol, src, Name, Value, ...)
%         per-observation residuals, redundancy, w / tau, MDB, external
%         reliability, the contamination gauge, conditioning, correlation
%     [sol, D, REJ] = snoop(dsol, D, src, Name, Value, ...)
%         iterative data snooping (rank-one downdates of the block cofactor)
%     S = systematic(D, Name, Value, ...)
%         Ljung-Box / runs / covariate tests for unmodelled signal
%     F = ftest(solFull, solRed, Name, Value, ...)
%         nested-model F test
%     E = ellipse(D, [i j], Name, Value, ...)
%         confidence ellipse of a parameter pair
%   Each method's own help: help lsqvce.engine.verdict.<method>.
%
%   EXAMPLE
%       R = lsqvce.fit(task, Verbose=false);          % R.diag is the verdict
%       [D, dsol, src, nB] = lsqvce.engine.verdict.fromAccumulator(R.acc, R.sol, 'Test', "mad");
%       S = lsqvce.engine.verdict.systematic(D);      % what R.systematic() runs
%
%   See also LSQVCE.RESULTS, LSQVCE.REPORT, LSQVCE.ENGINE.BLOCKCOFACTOR, LSQVCE.KERNEL.STAT
%
%   Updated 2026-11-03 (CET, UTC+01:00), WP-162: the examples fit a
%   lsqvce.TaskDefinition (lsqvce.fit takes no struct since 7.0.0). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-17 (CEST, UTC+02:00), WP-149: supersede, the refused v'Pv
%   replaced by the verdict's residual pass. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-138 (code review of WP-135):
%   every method states its ACCURACY; the contamination gauge -- fhat's null
%   law and the ladder -- is the private contaminationOf (the class stays one,
%   WP-86 part 4); one redundancy formula (redundancyNumbers) for diagnose's
%   passes, its datum pre-pass and snoop's scan; one list of the prior batches
%   (priorBatchesOf); snoop reads D's thresholds. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-14 (CEST, UTC+02:00), WP-135 (review 3 findings 5, 19 and
%   the diagnostics lows): local prior groups get the datum rule and the no-
%   rejection rule; the contamination ladder is calibrated for every n;
%   snoopVacuous is a practical criterion; systematic leaves the priors' rows
%   out; snoop keeps the caller's options. CI round: P(K = 0) of the ladder's
%   null law is the closed form at every n (the normal limit read an empty
%   table). With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Created 2026-10-11 (CEST, UTC+02:00), R7 WP-86 part 4 (user decision:
%   the full absorption): diagnoseAccum, diagnose, snoop, systematic, ftest
%   and ellipse are this class's static methods, their bodies unchanged but
%   for the renamed references; stat moved to lsqvce.kernel.stat. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.

methods (Static)

function [D, diagSol, batchFcn, nBatch] = fromAccumulator(acc, sol, varargin)
%LSQVCE.ENGINE.VERDICT.FROMACCUMULATOR  The verdict of a solved block accumulator, in one call.
%
%   [D, diagSol, batchFcn, nBatch] = LSQVCE.ENGINE.VERDICT.FROMACCUMULATOR(acc, sol, ...)
%   is source (the full-parameter solution and the replayable batch source)
%   followed by diagnose on them -- the pair every driver ran (fit,
%   adjustBlock, adjustRobust, adjustNonlinear, snoop's clean re-solve).
%
%   INPUT
%     acc       1x1 lsqvce.engine.BlockAccumulator  solved, 'KeepLocal' true
%     sol       1x1 struct  the solution it returned
%     varargin  Name-Value pairs forwarded to LSQVCE.ENGINE.VERDICT.DIAGNOSE
%               (Alpha, Bonferroni, Beta, RMin, RWeak, Test, RobustCheck,
%               Store, NWorst, Progress; NBatch is set here)
%   OUTPUT
%     D         1x1 struct  the diagnostics (R.diag)
%     diagSol   1x1 struct  the full-parameter solution (R.source)
%     batchFcn  1x1 fh      the batch source (R.batchFcn)
%     nBatch    1x1 double  its batch count (R.nBatch)
%
%   A REFUSED v'Pv (WP-150, F6): when the backend refused a cancelled v'Pv
%   (sol.s0sq NaN), diagnose studentises with the fallback sigma0 and flags
%   against the wrong scale -- reproduced: 28 rows flagged on Backend "normal"
%   against 1 (the blunder) on "qr", R.sigma0 3.27 on both. The verdict then
%   runs diagnose once more with the s0^2 of its own residual pass (D.s0sq,
%   the number LSQVCE.ENGINE.VERDICT.SUPERSEDE gives R.sigma0), and diagSol
%   carries it. Only in that case: the flags are decided inside the data
%   pass, before its v'Pv exists.
%
%   ACCURACY
%     As LSQVCE.ENGINE.VERDICT.SOURCE and LSQVCE.ENGINE.VERDICT.DIAGNOSE: it
%     adds no arithmetic of its own (WP-138).
%
%   EXAMPLE
%       [D, dsol, src, nB] = lsqvce.engine.verdict.fromAccumulator(R.acc, R.sol, 'Test', "mad");
%
%   See also LSQVCE.ENGINE.VERDICT.SOURCE, LSQVCE.ENGINE.VERDICT.DIAGNOSE
%
%   Updated 2026-10-17 (CEST, UTC+02:00), WP-150 (F6): a refused v'Pv is
%   diagnosed again with the pass's own s0^2. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-138 (code review of WP-135):
%   ACCURACY; it adds no arithmetic of its own. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
    [diagSol, batchFcn, nBatch] = lsqvce.engine.verdict.source(acc, sol);
    D = lsqvce.engine.verdict.diagnose(diagSol, batchFcn, 'nBatch', nBatch, varargin{:});
    if ~isfinite(diagSol.s0sq) && D.s0Fallback && isfinite(D.s0sq) && D.s0sq > 0
        diagSol.s0sq = D.s0sq;                       % WP-150 (F6): the pass's own s0^2
        D = lsqvce.engine.verdict.diagnose(diagSol, batchFcn, 'nBatch', nBatch, varargin{:});
    end
end

function sol = supersede(sol, D, acc)
%LSQVCE.ENGINE.VERDICT.SUPERSEDE  A refused v'Pv replaced by the verdict's residual pass.
%
%   sol = LSQVCE.ENGINE.VERDICT.SUPERSEDE(sol, D, acc) returns sol unchanged
%   unless the backend refused a cancellation-destroyed v'Pv (sol.s0sq NaN,
%   and with it sol.sigma0 and every entry of sol.Cxx; lsqvce:accum:cancellation,
%   lsqvce:normalEq:cancellation). Then vtv, s0sq and sigma0 are the verdict's
%   own -- D.vPv, the DIRECT sum of the whitened residuals of its second pass,
%   which is what those warnings promise supersedes the refused number
%   (WP-49) -- and Cxx is the solve's cofactor scaled by it,
%   D.s0sq * acc.commonCofactor(), the rule lsqvce.engine.residuals applies to
%   its local covariances (WP-132). lsqvce.Results calls it, so every producer
%   of a Results (fit, adjustBlock, adjustRobust, adjustNonlinear) returns
%   the honest numbers rather than the NaN; the pass is the one the verdict
%   ran anyway, no data is read again.
%
%   INPUT
%     sol   1x1 struct  the solution (.s0sq .sigma0 .vtv .Cxx)
%     D     1x1 struct  its verdict (LSQVCE.ENGINE.VERDICT.FROMACCUMULATOR's D:
%                       .vPv .s0sq .s0); anything without them leaves sol as is
%     acc   lsqvce.engine.BlockAccumulator | []  the accumulator that solved it
%   OUTPUT
%     sol   1x1 struct  as the input, or with .vtv .s0sq .sigma0 .Cxx replaced
%
%   ACCURACY
%     No arithmetic beyond the one scaling. On a tight fit refused on every
%     solve, sigma0 agrees with lsqvce.engine.residuals' to 1e-12 relative
%     (measured 2.3e-13) and sigma0 and the formal errors with the QR
%     backend's on the same data to 1e-9 (measured <= 2.7e-11, robust runs
%     included; test_vcePolicy, WP-149). D.s0sq carries a robust hard
%     solution's truncation divisor, as sol.s0sq would have.
%
%   EXAMPLE
%       [D, dsol] = lsqvce.engine.verdict.fromAccumulator(acc, sol);
%       sol = lsqvce.engine.verdict.supersede(sol, D, acc);   % what Results does
%
%   See also LSQVCE.RESULTS, LSQVCE.ENGINE.RESIDUALS, LSQVCE.ENGINE.BLOCKACCUMULATOR
%
%   Created 2026-10-17 (CEST, UTC+02:00), WP-149: the fit's own Results carried
%   R.sigma0 = NaN and NaN formal errors after a refused v'Pv, although the
%   warning said the residual pass supersedes it. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
    if ~isa(acc, 'lsqvce.engine.BlockAccumulator') || isfinite(lsqvce.util.getfielddef(sol, 's0sq', 1)) ...
            || ~all(isfield(D, {'vPv', 's0sq', 's0'})) || ~isfinite(D.s0sq)
        return
    end
    sol.vtv    = D.vPv;
    sol.s0sq   = D.s0sq;
    sol.sigma0 = D.s0;
    sol.Cxx    = D.s0sq * acc.commonCofactor();
end

function [diagSol, batchFcn, nBatch] = source(acc, sol)
%LSQVCE.ENGINE.VERDICT.SOURCE  Bridge a block-arrow lsqvce.engine.BlockAccumulator adjustment into LSQVCE diagnostics.
%
%   [diagSol, batchFcn, nBatch] = LSQVCE.ENGINE.VERDICT.SOURCE(acc, sol) returns exactly what
%   lsqvce.engine.verdict.diagnose / lsqvce.engine.verdict.snoop expect:
%       diagSol   struct with fields {x, Q, n, u, r, rank, defect, sigma0, s0sq,
%                 vPv, vPvSource, condR, condN, condIsBound, nBatch, nExcluded,
%                 datumRows, prior, priorBatch} (n: every replayed row, r: the
%                 live rows' redundancy -- n - nExcluded - rank plus a ridge's
%                 defect sum tr(Q_ll Lambda), WP-133; nExcluded: the rows of
%                 weight 0, which batchFcn replays with p = 0;
%                 datumRows: the prior's datum rows, WP-127;
%                 prior / priorBatch: the CommonPrior and its batch index, 0
%                 without one, WP-128; localPriorBatches: the batches of
%                 the local prior groups' pseudo-arcs, accOpt: the solve's
%                 Backend / Reg options for snoop's re-solve, truncation: a
%                 hard robust solution's tau (1 otherwise), WP-135); Q is the
%                 LSQVCE.ENGINE.BLOCKCOFACTOR of the full parameter set
%       batchFcn  @(k) -> [A, l, p]   the k-th WHITENED batch (one arc), unit
%                 weight, p = 0 on the arc's excluded rows (WP-133)
%       nBatch    number of batches (the arcs, then the prior's pseudo-batch)
%
%   so you can run the full five-level lsqvce verdict on a gravity-recovery
%   adjustment built with the streaming accumulator:
%
%       R   = lsqvce.fit(task, Estimate="noise", Pool="obsType");
%       acc = R.acc;   sol = R.sol;               % fit did the next three lines
%       [dsol, src, nB] = lsqvce.engine.verdict.source(acc, sol);   % already: R.source / R.batchFcn / R.nBatch
%       D = lsqvce.engine.verdict.diagnose(dsol, src, 'nBatch', nB);
%       V = lsqvce.report(D);              % the five levels
%       S = lsqvce.engine.verdict.systematic(D);          % Ljung-Box / runs / covariate
%
%   WHY THE FULL PARAMETER SET, NOT THE REDUCED ONE
%   -----------------------------------------------
%   lsqAccum pre-eliminates the local parameters per system by Schur complement
%   and keeps only reduced common normals. lsqvce per-observation redundancy
%   r_i = 1 - a_i' Qxx a_i is only correct if a_i and Qxx span the SAME
%   parameter set the observation actually constrains -- i.e. the FULL set
%   [common | all locals]. Diagnosing against the reduced common block alone
%   would miscount the degrees of freedom the local parameters absorb and
%   inflate every r_i. So this bridge works on the full block-arrow system:
%   columns are [x_common | y_sys1 | y_sys2 | ...], and each batch places its
%   whitened A into the common columns and its whitened B into that system's
%   local columns. Its cofactor Qxx is never formed: BlockCofactor answers
%   a'*Qxx*a from the reduced commons factor and the per-system Schur blocks
%   the solve used (WP-86 part 2), exactly the block-arrow inverse.
%
%   The whitening replays the loader's own Wfcn/Cfcn, batch by batch, with the
%   same overlap discard as accumulation -- so the batches are bit-identical to
%   what the accumulator saw.
%
%   COST
%   ----
%   Nothing u_full x u_full (u_full = nc + sum(nl)) is formed (WP-86 part 2):
%   the cofactor is the nc x nc commons factor plus, per system, its nl x nl
%   local factor and the nl x nc coupling G -- memory O(nc^2 + sum nl*nc),
%   and lsqvce.engine.verdict.diagnose costs O(nc^2 + nl*nc) per row with it. Before,
%   a dense QR of every row at full width built a u_full x u_full R (1.5 GB at
%   a GRAVFIRE month, ~u_full^2/2 flops per row and pass, and diagnose held
%   four such matrices). The batches are NOT stored (WP-86 part 1, R7-24):
%   batchFcn re-reads an arc through its loader and re-whitens it on demand,
%   caching one arc, so the memory is one arc at full width. The row count
%   comes from the accumulator and the index from the loaders' meta answers,
%   so the bridge itself reads no data. The price is a re-read of the data
%   per diagnostic pass, and it RELIES ON THE LOADER'S CONTRACT that an arc is
%   returned identically on every call (a changed count is still caught:
%   lsqvce:diagnose:batchMismatch).
%
%   REQUIRES  acc constructed with 'KeepLocal', true (the default).
%
%   INPUT
%     acc       1x1 lsqvce.engine.BlockAccumulator  the accumulator that solved
%     sol       1x1 struct  the solution it returned
%
%   OUTPUT
%     diagSol   1x1 struct   the FULL-parameter solution (fields above): the
%                            common redundancies are CORRECTED for the
%                            eliminated locals -- reading the reduced common
%                            system alone inflates every reliability number
%     batchFcn  1x1 fh       the replay source, one arc per batch
%     nBatch    1x1 double   number of batches
%
%   ACCURACY
%     The batches are no computation: each arc is re-read and re-whitened by
%     the one iterator, so the rows are bit-identical to accumulation (a
%     replay that differs is refused by diagnose, batchMismatch). The
%     cofactor is LSQVCE.ENGINE.BLOCKCOFACTOR's: its blocks against the dense
%     inverse of the same system to 1e-12 relative (test_lsqAccum). WP-138.
%
%   ERRORS
%     lsqvce:diagnoseAccum:noSystems     the accumulator carries no systems
%     lsqvce:diagnoseAccum:needKeepLocal the local blocks were not retained, so
%                                   the common redundancies cannot be corrected;
%                                   re-run keeping them
%
%   See also LSQVCE.ENGINE.VERDICT.DIAGNOSE, LSQVCE.ENGINE.VERDICT.SNOOP, LSQVCE.REPORT, LSQVCE.ENGINE.REDUNDANCY,
%   LSQVCE.ENGINE.BLOCKCOFACTOR
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-138 (code review of WP-135):
%   ACCURACY. With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-14 (CEST, UTC+02:00), WP-135 (review 3 finding 5, lane 1):
%   diagSol.localPriorBatches (a local prior group's rows: a datum in
%   diagnose, never a snoop candidate) and diagSol.accOpt (snoop re-solves
%   with the solve's Backend and Reg). With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-133: the rows of weight 0
%   (arc.excluded) are replayed with p = 0 and leave r (diagSol.nExcluded);
%   r adds a regularised local block's defect. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-130: r and rank count a
%   rank-deficient local block by its rank. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-128: diagSol.prior and
%   .priorBatch -- snoop re-solves with the prior instead of as data rows.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-127: diagSol.datumRows -- the
%   global indices of the prior's datum rows (BlockAccumulator.datumRows),
%   which diagnose does not test or count. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-86 part 2 (R7-24): the dense
%   full-width factor replaced by the block cofactor; one batch per arc,
%   indexed from the meta answers, n from the accumulator -- the bridge reads
%   no data. With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-86 part 1 (R7-24): the
%   batches are rebuilt on demand from their address, one arc cached, instead
%   of being stored at full width (n x u_full). With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-118 part 1 (R7-34): the
%   deleted engine.adjust reference replaced. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    % The bridge needs the per-system Schur blocks (the cofactor's local
    % factors, the back-substitution, the arc replay). They are only retained
    % when the accumulator was built with 'KeepLocal', true; refuse loudly
    % rather than return a silently degenerate (zero-batch) diagSol.
    if ~lsqvce.util.getfielddef(acc.opt, 'KeepLocal', true)
        error('lsqvce:diagnoseAccum:needKeepLocal', ...
            ['diagnoseAccum needs the per-system Schur blocks to rebuild the ' ...
             'full block-arrow system, but the accumulator was built with ' ...
             '''KeepLocal'', false. Reconstruct it with ''KeepLocal'', true ' ...
             '(the default) before bridging into diagnostics.']);
    end

    sysIds = acc.systems();
    nSys   = numel(sysIds);
    if nSys == 0
        error('lsqvce:diagnoseAccum:noSystems', ...
            ['The accumulator has no stored systems to diagnose. Add at ' ...
             'least one system and solve before calling diagnoseAccum.']);
    end
    nc     = numel(sol.x);

    % ---- global column layout: [common | locals(sys1) | locals(sys2) | ...]
    localOffset = zeros(1, nSys);
    localCount  = zeros(1, nSys);
    col         = nc;
    stores      = cell(1, nSys);
    n           = 0;
    defL        = 0;                                % WP-130: locals the solve could not estimate
    nExcl       = 0;                                % WP-133: rows of weight 0 (replayed, not observations)
    for iSys = 1:nSys
        st = acc.systemStore(sysIds(iSys));
        stores{iSys}     = st;
        localOffset(iSys) = col;
        localCount(iSys)  = st.nl;
        kept = acc.getStore(sysIds(iSys));
        mEx = lsqvce.util.getfielddef(kept, 'mExcl', 0);
        n   = n + kept.mObs + mEx;                     % rows the replay serves
        nExcl = nExcl + mEx;
        defL = defL + st.nl - lsqvce.engine.BlockAccumulator.localRankOf(kept);
        col = col + st.nl;
    end
    uFull = col;                                    % nc + sum(nl)

    % ---- full parameter vector: x_common then each system's back-substituted y
    xFull = zeros(uFull, 1);
    xFull(1:nc) = sol.x;
    for iSys = 1:nSys
        y = acc.backSubstitute(sysIds(iSys), sol.x);
        xFull(localOffset(iSys) + (1:localCount(iSys))) = y;
    end

    % ---- the cofactor, as blocks (WP-86 part 2) -------------------------
    % Qxx of the full parameter set from the reduced commons factor and the
    % per-system blocks: nothing u_full x u_full is factored or held (it was a
    % dense QR of every whitened row at full width, 1.5 GB per copy and
    % ~u_full^2/2 flops per row at a GRAVFIRE month).
    Q = lsqvce.engine.BlockCofactor.fromAccumulator(acc);   % same system order, same layout

    % ---- the batch index: one batch = one ARC, from the meta calls alone --
    % batchFcn rebuilds an arc's whitened rows from its address on demand
    % (WP-86 part 1, R7-24); the row count comes from the accumulator, so
    % building the index reads no data (part 2).
    index = enumerateArcs(stores);
    % braces around the cells: a cell VALUE, not a struct array (the struct()-cell pitfall)
    src = struct('stores', {stores}, 'localOffset', localOffset, 'localCount', localCount, ...
                 'nc', nc, 'uFull', uFull, 'index', {index}, 'prior', [], ...
                 'cache', containers.Map('KeyType', 'char', 'ValueType', 'any'));
    % R5 M5 5b: the common-block PRIOR is one more batch of pseudo-observations
    % (rows K/sigma on the common columns, observations x0/sigma), so the hat
    % matrix, the redundancies and the residual sum include it exactly as the
    % solve did (Koch: the prior's redundancy is what it absorbs).
    nBatch = size(index, 1);
    if ~isempty(acc.prior)
        K  = full(acc.prior.K) / acc.prior.sigma;
        Ap = [K, zeros(size(K, 1), uFull - nc)];
        lp = acc.prior.x0(:) / acc.prior.sigma;
        src.prior = struct('A', Ap, 'l', lp);       % nk x uFull: the one stored batch
        nBatch = nBatch + 1;
        n     = n + size(Ap, 1);
    end

    % ---- diagSol in lsqvce contract -------------------------------------
    % cond: a LOWER bound from the blocks (an SVD of the u_full factor is the
    % O(u_full^3) step part 2 removes); report prints it as one.
    kN = Q.condLowerBound();
    % R5 M5 5b': sigma0 below is the A-PRIORI sigma0 of the law-whitened rows,
    % 1 by construction, so diagnose's varRatio and its a-priori backstop are
    % s0sq/1. The accumulator's sol.sigma0 is the a-posteriori one -- passing
    % it here made both identically 1 (the breakdown backstop could never fire).
    diagSol = struct( ...
        'x',      xFull, ...
        'Q',      Q, ...
        'n',      n, ...
        'u',      uFull, ...
        'r',      n - nExcl - uFull + defL + acc.regularisationDefect(), ...
        'rank',   uFull - defL, ...
        'defect', 0, ...
        'sigma0', 1, ...
        's0sq',   sol.s0sq, ...
        'vPv',    sol.vtv, ...
        'vPvSource', 'accumulator (v''v)', ...
        'condR',  sqrt(kN), ...
        'condN',  kN, ...
        'condIsBound', true, ...
        'nBatch', nBatch, ...
        'nExcluded', nExcl, ...
        'datumRows', datumRowsOf(acc, n), ...
        'prior', {acc.prior}, 'priorBatch', nBatch * ~isempty(acc.prior), ...   % WP-128: snoop re-solves with it
        'localPriorBatches', find(lsqvce.util.priorStream.is(index(:, 3))), ...   % WP-135: the local groups' arcs
        'truncation', lsqvce.util.getfielddef(sol, 'truncation', 1), ...   % WP-135: a hard robust solution's tau
        'accOpt', struct('Backend', acc.backend, 'Reg', acc.opt.Reg, 'RegLambda', acc.opt.RegLambda, ...
                         'RegMatrix', acc.opt.RegMatrix, 'TruncTol', acc.opt.TruncTol));   % WP-135: snoop's re-solve

    % ---- batchFcn: replay the k-th whitened batch --------------------------
    % A local function (not deal) so a caller may take 1, 2, or 3 outputs:
    % lsqvce.engine.verdict.diagnose asks for [A,l,p], but a redundancy-only sweep wants just
    % A. deal() would throw on nargout < 3; batchOf never does.
    batchFcn = @(k) batchOf(src, k);
end

function D = diagnose(sol, batchFcn, opt)
%LSQVCE.ENGINE.VERDICT.DIAGNOSE  Second-pass per-observation diagnostics for a batch LSQ.
%
%   D = LSQVCE.ENGINE.VERDICT.DIAGNOSE(sol, batchFcn) sweeps the observations a second time
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
%   Nothing of size n-by-n, n-by-u or u-by-u is formed: the cofactor of the
%   FULL parameter set (commons + every system's locals) is sol.Q, an
%   LSQVCE.ENGINE.BLOCKCOFACTOR, so per row the cost is O(nc^2 + nl*nc) and the
%   memory O(nc^2 + sum nl*nc) (WP-86 part 2; a dense u x u factor before). The
%   batch source every fit uses, LSQVCE.ENGINE.VERDICT.SOURCE, rebuilds the
%   batches on demand (WP-86 part 1, R7-24).
%
%   INPUTS
%     sol       the diagnosed solution (R.source; LSQVCE.ENGINE.VERDICT.SOURCE builds
%               it): x, Q (the block cofactor), n, u, r, sigma0, s0sq, vPv, ...
%     batchFcn  function handle, [A,l,p] = batchFcn(k), k = 1..nBatch.
%               Must reproduce EXACTLY the batches given to the accumulator,
%               in the same order. This is checked: a replay that sees a
%               different observation count is lsqvce:diagnose:batchMismatch,
%               and D.sumRRelErr / D.vPvRelDiff report the integrity numbers.
%
%   NAME-VALUE
%     'nBatch'      number of batches (default sol.nBatch)
%     'alpha'       per-test significance          (default 0.001)
%     'bonferroni'  correct alpha for n tests      (default true)
%     'beta'        MDB type-II error              (default 0.20)
%     'rMin'        redundancy below which an observation is declared
%                   UNCONTROLLED and its w is not reported (default 0.01,
%                   the same as Options.Diagnostics.RMin -- R7-35)
%     'rWeak'       redundancy below which an observation is flagged as
%                   weakly controlled                        (default 0.10)
%     'test'        'auto' | 'pope' | 'baarda' | 'mad'   (default 'auto')
%                   'baarda' assumes sigma0 is known a priori; 'pope'
%                   studentises with the estimated s0. 'auto' picks pope,
%                   because in practice sigma0 is essentially never known.
%     'store'       'all' | 'flagged' | 'none'   (default 'all')
%                   'flagged' keeps only observations that trip a test plus
%                   the worst offenders -- use this when n is enormous.
%     'nWorst'      how many worst-|w| observations to retain under
%                   'store','flagged'                        (default 1000)
%     'robustCheck' also form the MAD scale and the contamination gauge
%                   (in the main pass; 'test','mad' needs the scale first and
%                   sweeps once more; always on for 'mad') (default true)
%     'progress'    print a progress line every ~2 s     (default false)
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
%                    contamination gauge -- everything LSQVCE.REPORT reads.
%                    The parameter cofactor is D.Q (D.Q.block(idx) is
%                    Qxx(idx, idx)); D.sigX / D.qxxDiag its diagonal, and
%                    D.corrMax / D.corrMaxPar the strongest correlation.
%                    A row of weight 0 (p = 0 in its batch: a robust
%                    rejection, a zero data weight -- WP-133) keeps its index
%                    but is no observation: r = 0, not tested, not in sum(r)
%                    or D.nTests (sol.nExcluded), D.quad.rejected counts it
%                    and D.excludedIdx names it
%
%   EXAMPLE
%       % The SAME batch source must be replayed -- this is checked. On the
%       % block path (v5) R.source is the diagnosed solution and R.batchFcn /
%       % R.nBatch the replayable source LSQVCE.ENGINE.VERDICT.SOURCE built.
%       data = struct('batchFcn', src, 'numUnknowns', u, 'numBatches', 4);
%       R = lsqvce.fit(task, Verbose=false);           % R.diag is this D
%
%       % Re-diagnose with the ROBUST scale instead of Pope's tau. Do this
%       % whenever D.contamLevel is 'moderate' or worse: Pope goes BLIND at
%       % ~10%% contamination because the blunders inflate the s0 it divides by.
%       Dm = lsqvce.engine.verdict.diagnose(R.source, R.batchFcn, 'nBatch', R.nBatch, 'test', 'mad');
%
%       % Dm.flagIdx    -> [40 150]
%       % Dm.quad       -> the (w, r) quadrants; read UNTESTED, not just w
%       % Dm.mdb        -> largest error that could hide in each observation
%
%   ACCURACY
%     Per observation, from the block cofactor batch by batch (no n x n
%     matrix): the partition identity sum(r_i) = n - u holds to 1e-10
%     relative, the sparse and the dense path agree to 1e-10 in r and 1e-8
%     in w, a diagonal weight MATRIX reproduces the vector form to 1e-10
%     (test_lsqvce, test_lsqvce_ext, test_lsqvce_opts), and the external
%     reliability matches a brute-force injection of the MDB to 1e-7
%     relative (test_lsqvce, test_external_reliability_is_exact). The
%     critical values and delta0 are the Baarda / Pope textbook values at
%     the quoted decimals and scipy.stats to 1e-9 relative down to alpha =
%     1e-17 (HANDOVER 6.1.4, the "Baarda / Pope" and "scipy.stats" rows).
%     The contamination gauge's null band (fhatTol) keeps false alarms <=
%     0.25 % for n = 20..2000 (the third review's mirrors row; the private
%     contaminationOf). WP-138.
%
%   ERRORS
%     lsqvce:diagnose:nBatch          nBatch not given and sol carries none
%     lsqvce:diagnose:batchMismatch   batchFcn returned DIFFERENT data on this
%                                     second pass. It is called twice -- once to
%                                     accumulate, once to diagnose -- and must
%                                     be deterministic
%     lsqvce:diagnose:badPmatrix      a full weight matrix where a diagonal was
%                                     expected
%     lsqvce:diagnose:noRobustScale   a robust scale was requested but the
%                                     residuals cannot support one
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (b) b2 (F-07):
%   typed `arguments` signature, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   See also LSQVCE.ENGINE.VERDICT.SOURCE, LSQVCE.REPORT, LSQVCE.ENGINE.VERDICT.SNOOP,
%   LSQVCE.ENGINE.BLOCKCOFACTOR.
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-138 (code review of WP-135):
%   ACCURACY; a local prior batch with a full weight matrix is judged by the
%   generalised redundancy (it was skipped without a word, its datum row
%   counted uncontrolled); one redundancy formula for every pass. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-14 (CEST, UTC+02:00), WP-135 (review 3 findings 5, 19, lane
%   1): a local prior group's row with r_i <= 1e-6 is a datum row (D.datumRows,
%   not tested or counted; it made the verdict SUSPECT); D.priorRows; the
%   contamination ladder from fhat's measured null law at every n (D.fhatTol,
%   D.collapseDetectable, the level 'indeterminate' -- n = 30 clean data raised
%   severe in 3 %, a sigma0 3.5x too small read 'breakdown' in 50 %);
%   D.snoopVacuous from kCrit^2 / r >= 0.9 (D.vacuityShare; the old test could
%   not fire); a mad replay that differs is batchMismatch; the dead badTest
%   branch and the unused inv(P) diagonal removed. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-133 (review 3 finding 14): a row
%   of weight 0 is "rejected" -- not tested, not counted, its own quadrant
%   (D.quad.rejected, D.excludedIdx); a robust fit's rejections were
%   "verified" and in the dof. With assistance from Claude (Anthropic), model
%   Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-131: the MDB of a correlated
%   batch is delta0 / sqrt((P Qvv P)_ii), the generalised test's own (it was
%   the diagonal-P formula). With assistance from Claude (Anthropic), model
%   Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-130: Test 'mad' leaves the datum
%   rows out of its scale (as the fused path does since WP-127); D.rank is the
%   solve's rank. With assistance from Claude (Anthropic), model Claude Opus
%   5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-128: Test 'mad' uses
%   lsqvce.kernel.stat.madCrit (Student t, nu = r/2), not Pope's tau.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-127: a prior's datum rows
%   (sol.datumRows; r = 0 by construction) stay in sum(r) and vPv but are not
%   tested or counted -- uncontrolled, MDB, quadrants, the smallest r and the
%   Bonferroni n (D.nTests, D.datumRows); every free network was SUSPECT.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-12 (CEST, UTC+02:00), R7 WP-121: the robust check's raw
%   scores are formed in the main pass -- one read of the data and one
%   halfRows per row instead of two; only 'test','mad' sweeps first. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-86 part 2 (R7-24): the
%   cofactor comes as blocks (sol.Q); the dense R, Qxx and correlation matrix
%   (D.R, D.Qxx, D.corr) and the largeU warning about them are retired. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-118 part 3 (R7-24): the help
%   separates this function's own working set from its batch source's --
%   diagnoseAccum keeps all whitened rows at full width (WP-86). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-118 part 1 (R7-34..36): the
%   dead 'verify' option retired, RobustCheck / Progress documented, rMin
%   defaults to Options.Diagnostics.RMin (0.01), the MDB reservoir draws
%   from its own stream. With assistance from Claude (Anthropic), model
%   Claude Opus 5.5.
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
    opt.RMin (1,1) double {mustBeNonnegative} = 0.01   % R7-35: Options.Diagnostics.RMin's value
    opt.RWeak (1,1) double {mustBeNonnegative, mustBeLessThanOrEqual(opt.RWeak,1)} = 0.10
    opt.Test (1,1) string {mustBeMember(opt.Test, ["auto","pope","baarda","mad"])} = "auto"
    opt.RobustCheck (1,1) logical = true
    opt.Store (1,1) string {mustBeMember(opt.Store, ["all","flagged","none"])} = "all"
    opt.NWorst (1,1) double {mustBeNonnegative} = 1000
    opt.Progress (1,1) logical = false
end

o = struct('nBatch', opt.NBatch, 'alpha', opt.Alpha, 'bonferroni', opt.Bonferroni, 'beta', opt.Beta, 'rMin', opt.RMin, 'rWeak', opt.RWeak, 'test', char(opt.Test), 'robustCheck', opt.RobustCheck, 'store', char(opt.Store), 'nWorst', opt.NWorst, 'progress', opt.Progress);

nB = o.nBatch;
if isempty(nB)
    if isfield(sol, 'nBatch') && sol.nBatch > 0
        nB = sol.nBatch;
    else
        error('lsqvce:diagnose:nBatch', ...
              'Cannot infer the batch count; pass ''nBatch''.');
    end
end

Q  = sol.Q;
x  = sol.x;
n  = sol.n;
u  = sol.u;
rdf = sol.r;                          % redundancy = n - u

% -------------------------------------------------------------- test setup
% WP-127: the datum rows of a prior (r = 0 by construction) are constraints,
% not observations: kept in sum(r) and vPv, left out of the tests and counts
% WP-133 (review 3 finding 14): rows of weight 0 (a batch's p = 0: a robust
% rejection, a zero data weight) are not observations either -- not tested,
% not in sum(r) or nTests; sol.r already leaves them out (sol.nExcluded)
exIdx = lsqvce.util.getfielddef(sol, 'datumRows', zeros(0, 1));
% WP-135 (review 3 finding 5): a LOCAL prior group's rows ride as pseudo-arcs
% of their system (sol.localPriorBatches). A row the data cannot see is a
% datum like the common prior's (WP-127): r_i = 0 by construction, and it was
% counted as an uncontrolled observation (REL_UNCONTROLLED, verdict SUSPECT).
% Found before the pass -- the prior batches are cheap -- so Bonferroni's n
% is right; tol: r_i <= 1e-6 -- a true datum's r is 0 to rounding (~1e-11 on
% the factor, ~1e-7 on squared normals), and a prior a million times
% stronger than the data along its row is a constraint in practice.
[locDatum, nLocDatum] = localDatumOf(batchFcn, Q, lsqvce.util.getfielddef(sol, 'localPriorBatches', []));
nTest = n - numel(exIdx) - nLocDatum - lsqvce.util.getfielddef(sol, 'nExcluded', 0);
alphaEff = o.alpha;
if o.bonferroni
    alphaEff = o.alpha / max(nTest, 1);
end

testName = lower(o.test);
if strcmp(testName, 'auto'), testName = 'pope'; end

switch testName
    case 'pope'
        kCrit = lsqvce.kernel.stat.tauCrit(alphaEff, rdf);
    case 'mad'
        % WP-128: NOT Pope's tau. tau holds only when the scale contains the
        % tested residual; a MAD scale does not, and its own noise makes the
        % statistic Student t with nu = r/2 (lsqvce.kernel.stat.madCrit, calibrated):
        % the tau threshold gave 63x the nominal false-alarm rate at n = 30.
        % (Previously: "The MAD statistic is still a STUDENTISED residual --
        % only the scale in the denominator is robust. It remains
        % tau-distributed, so the critical value is Pope's.")
        %
        % THIS CASE WAS MISSING. kCrit was simply never assigned on the 'mad'
        % path, and diagnose died with "Unrecognized function or variable
        % kCrit". The entire robust-scale cure -- the toolbox's headline fix for
        % Pope's breakdown -- was therefore DEAD ON ARRIVAL. The Python port
        % computed kCrit separately and so structurally could not see it.
        % Caught only by finally running the suite in MATLAB.
        kCrit = lsqvce.kernel.stat.madCrit(alphaEff, rdf);
    otherwise   % 'baarda' -- the arguments block admits nothing else (WP-135:
                % the unknown-test error branch here was unreachable)
        kCrit = lsqvce.kernel.stat.normCrit(alphaEff);
end
delta0 = lsqvce.kernel.stat.delta0(alphaEff, o.beta);

% Pope's tau is bounded by sqrt(rdf): under a hard Bonferroni correction on a
% weak network the critical value creeps up to that bound, and the snooping
% becomes vacuous. Say so rather than silently returning an empty outlier
% list. WP-135 (review 3, lane 1): tauCrit < sqrt(r) for every r > 1, so the
% old test kCrit >= sqrt(r) could never fire. The practical measure: a
% residual of full redundancy (r_i = 1) is flagged only if it carries
% kCrit^2 / r of v'Pv; from 90 % on only a blunder that dominates the whole
% adjustment can be seen (r = 5, n = 20: 98 %; r = 10, n = 40: 87 %).
tauBound = sqrt(rdf);
vacuityShare = kCrit^2 / max(rdf, 1);
snoopVacuous = strcmp(testName, 'pope') && (rdf <= 1 || vacuityShare >= 0.9);

% ---------------------------------------------------------- allocate output
storeAll = strcmpi(o.store, 'all');   % 'none' stores what 'flagged' does
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
flagIdx = [];  flagW = [];  flagR = [];  flagV = [];  flagMDB = [];  exclIdx = zeros(0, 1);
locIdx = zeros(0, 1);                          % WP-135: the local prior groups' datum rows
priIdx = zeros(0, 1);                          % WP-135: every prior row (systematic leaves them out)
priBatches = priorBatchesOf(sol);              % WP-138: one list with snoop's
worstW  = -inf(min(o.nWorst, n), 1);
worstI  = zeros(numel(worstW), 1);

nUncontrolled = 0;   idxUncontrolled = [];

% (w, r) quadrant counters and the bounded MDB reservoir.
nQ_cleanOut = 0;  nQ_ambiguous = 0;  nQ_verified = 0;  nQ_untested = 0;
mdbCap   = 5000;          % reservoir size: quantiles to ~1% with 5000 samples
% R7-36: the reservoir draws from its OWN fixed-seed stream. randi on the
% global stream advanced the user's RNG on every fit with n > 5000 (against
% WP-36: the global RNG is never touched) and made the reported MDB
% quantiles differ run to run on identical input.
mdbRng   = RandStream('mt19937ar', 'Seed', 0);
mdbPool  = zeros(0, 1);
nMdbSeen = 0;
nMdbInf  = 0;
nWeak         = 0;
rmin_seen     = inf;  rmin_idx = 0;
maxAbsW       = 0;    maxAbsW_idx = 0;

% External-reliability accumulator: the worst effect on any parameter, in
% units of that parameter's own sigma. Needs diag(Qxx), which the block
% cofactor gives without forming Qxx (WP-86 part 2: the dense u x u Qxx, its
% inverse factor and the correlation matrix were the toolbox's O(u^3) step and
% its only O(u^2) allocations -- 1.5 GB each at u = 13 729).
%
% A free network enters as a prior on its null space (data.priors), folded
% into the cofactor like every other prior, so Qxx is the cofactor of the
% solve actually made and needs no correction here.
qxxDiag = Q.diagonal();

% A parameter that lies almost entirely in the datum null space can come out
% a few eps negative. Clamp rather than emit a complex sigma.
qxxDiag(qxxDiag < 0) = 0;
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
% The MAD needs the raw standardised residuals of EVERY observation before its
% scale exists. For 'pope' / 'baarda' they are formed in the main pass (WP-121);
% only 'mad', whose test needs that scale, sweeps the data once more first.

robustScale = NaN;
madAvailable = false;
fhat = NaN;
nScored = 0;
% WP-121 (the residue §12.48 recorded): for every test but 'mad' the scale the
% classification needs is known before the sweep (s0, or the a-priori sigma0),
% so the raw scores are formed IN the main pass from its own (P Qvv P)_ii and
% P v -- one read of the data and one halfRows per row instead of two, values
% bitwise those of the separate sweep. Only 'mad' needs the robust scale before
% it can classify, and keeps its own first sweep.
fuseRobust = o.robustCheck && ~strcmp(testName, 'mad');
if fuseRobust, rawAll = zeros(n, 1); end
if strcmp(testName, 'mad')
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

        % The raw score is the test statistic WITHOUT the scale factor:
        %       raw_i = (P v)_i / sqrt( (P Qvv P)_ii )
        % which is exactly what the MAD then rescales. It is defined for a
        % correlated batch just as it is for a diagonal one -- Baarda's
        % generalised form covers both -- so the robust check no longer has to
        % bail out on correlated data.
        % (redundancyNumbers forms the cofactor products, WP-138)
        if isCorrK
            Pk  = (pk + pk') / 2;
            [~, ~, WGk] = redundancyNumbers(Q, Ak, Pk, true);
            qpk = diag(Pk) - sum(WGk.^2, 1).';
            Pvk = Pk * vk;
        else
            pk  = pk(:);
            [~, ~, ~, ~, qvk] = redundancyNumbers(Q, Ak, pk, false);
            qvk(qvk < 0 | pk == 0) = 0;              % WP-133: a weight-0 row is no score
            qpk = pk.^2 .* qvk;
            Pvk = pk .* vk;
        end
        qpk = max(qpk, 0);

        % WP-130: a datum row (r = 0, qpk ~ +-eps) is no score -- as in the
        % fused path (WP-127), or its near-zero raw value pulls the MAD down
        good = qpk > 0 & ~ismember(pos + (1:nk).', exIdx);
        if isKey(locDatum, k), good = good & ~locDatum(k); end   % WP-135
        rk_ = zeros(nk, 1);
        rk_(good) = Pvk(good) ./ sqrt(qpk(good));
        rawAll(pos + (1:nk)) = rk_;
        pos = pos + nk;
    end

    if pos ~= n          % WP-135: a replay that differs is batchMismatch, not noRobustScale
        error('lsqvce:diagnose:batchMismatch', ...
            ['The second pass replayed %d observations but the accumulator saw ' ...
             '%d. batchFcn must reproduce the SAME batches in the SAME order.'], pos, n);
    end
    [robustScale, madAvailable, fhat, nScored] = robustScaleOf(rawAll, true);
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
    % Every term needs only P_k, A_k and Qxx. Cost is O((nc^2 + nl*nc) nk) per batch --
    % IDENTICAL to the diagonal path -- and nothing of size n is ever formed.
    % Verified exact to 2e-15 against a dense reference.
    %
    % When P_k is diagonal this reduces ALGEBRAICALLY to the familiar
    %       r_i = 1 - p_i a_i' Qxx a_i,     w_i = v_i / sigma_vi
    % (verified: the two agree to 4.4e-16), so there is one code path, not two.
    isCorr = ~isempty(pk) && ~isvector(pk);

    if isempty(pk), pk = ones(nk, 1); end

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
    ex = ismember(nSeen + (1:nk).', exIdx);   % the prior's datum rows (WP-127)
    if ismember(k, priBatches), priIdx = [priIdx; nSeen + (1:nk).']; end   %#ok<AGROW> prior rows
    if isKey(locDatum, k)                      % a local prior group's datum rows (WP-135)
        ex = ex | locDatum(k);
        locIdx = [locIdx; nSeen + find(locDatum(k))];   %#ok<AGROW> a handful of rows
    end
    rej = false(nk, 1);  if ~isCorr, rej = pk == 0; end   % WP-133: the rows of weight 0
    exclIdx = [exclIdx; nSeen + find(rej)];  ex = ex | rej;   %#ok<AGROW> neither is tested or counted

    % --- redundancy numbers and the denominator of the test statistic
    %     (r_i: redundancyNumbers, the one formula -- WP-138)
    if isCorr
        [rk, ~, WG, Gk] = redundancyNumbers(Q, Ak, Pk, true);   % Gk = P A
        % (P Qvv P)_ii = P(i,i) - g_i' Qxx g_i
        qpp = diag(Pk) - sum(WG.^2, 1).';
        % Pv, the numerator of the generalised statistic
        Pvk = Pk * vk;
        % (WP-135: the inv(P) diagonal formed here -- O(nk^3) -- went unused
        % since WP-131's MDB delta0 / sqrt(qpp); removed)
    else
        [rk, ~, ~, ~, qv] = redundancyNumbers(Q, Ak, pk, false);   % qv = 1/p - a'Qxx a
        qv(qv < 0 | pk == 0) = 0;      % CLAMP -- see the note below; weight 0: no r (WP-133)
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
    if fuseRobust                      % the robust check's raw score (WP-121)
        rawK = zeros(nk, 1);  okq = qpp > 0 & ~ex;
        rawK(okq) = Pvk(okq) ./ sqrt(qpp(okq));
        rawAll(nSeen + (1:nk)) = rawK;
    end

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
    ctrl = rk > o.rMin & qpp > 0 & ~ex;  % controlled observations only
    wk(ctrl) = Pvk(ctrl) ./ denom(ctrl);
    % Uncontrolled: w is 0/0. Report 0, and record the index -- an uncontrolled
    % observation is not "fine", it is UNTESTABLE, which is a different and
    % more dangerous thing.
    wk(~ctrl) = 0;

    % --- MDB.
    mdbk  = inf(nk, 1);
    % WP-131: the w-test's noncentrality under a bias nabla_i is nabla_i sqrt((P Qvv P)_ii),
    % so MDB_i = delta0 / sqrt(qpp_i); sigma_i delta0 / sqrt(r_i) is that only for diagonal P
    if isCorr, mdbk(ctrl) = delta0 ./ sqrt(qpp(ctrl)); else, mdbk(ctrl) = sig_i(ctrl) .* delta0 ./ sqrt(rk(ctrl)); end

    % --- external reliability: an undetected bias mdb_i in obs i shifts x by
    %       dx = Qxx * (P e_i) * mdb_i
    %     For diagonal P this is Qxx a_i p_i mdb_i, the familiar form. For a
    %     correlated batch the bias enters through the WHOLE column of P, i.e.
    %     through g_i = (P_k A_k)(i,:)', which is exactly Gk' -- so the two
    %     cases differ only in which matrix multiplies Qxx.
    if isCorr
        dX = Q.qxxTimes(Gk);           % u x nk : Qxx * g_i  (columns)
        scal = mdbk.';                 % the p_i factor is already inside g_i
    else
        dX = Q.qxxTimes(Ak);           % u x nk : Qxx * a_i  (columns)
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
    nQ_verified  = nQ_verified  + sum(~hiW & ~loR & ~ex);
    nQ_untested  = nQ_untested  + sum(~hiW &  loR & ~ex);

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
    nMdbInf = nMdbInf + sum(~mfin & ~ex);
    mv = mdbk(mfin);
    for q = 1:numel(mv)
        nMdbSeen = nMdbSeen + 1;
        if numel(mdbPool) < mdbCap
            mdbPool(end+1, 1) = mv(q); %#ok<AGROW>
        else
            j = randi(mdbRng, nMdbSeen);
            if j <= mdbCap
                mdbPool(j) = mv(q);
            end
        end
    end

    nUnc = sum(~ctrl & ~ex);
    if nUnc > 0
        nUncontrolled = nUncontrolled + nUnc;
        if numel(idxUncontrolled) < 1000
            idxUncontrolled = [idxUncontrolled; nSeen + find(~ctrl & ~ex)]; %#ok<AGROW>
        end
    end

    rkT = rk;  rkT(ex) = inf;  [rlo, ilo] = min(rkT);   % a datum row's r = 0 is not a finding
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
if fuseRobust
    [robustScale, madAvailable, fhat, nScored] = robustScaleOf(rawAll, true);
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
D.Q        = Q;                       % the cofactor, as blocks (WP-86 part 2)

D.vPv          = vPv_direct;          % the trustworthy one
D.vPvAccum     = sol.vPv;             % what the accumulator reported
D.vPvSource    = sol.vPvSource;
D.vPvRelDiff   = vPvRel;
% WP-135 (review 3, lane 4): a joint robust run's hard solution corrects its
% s0^2 for the truncation of the rejected epochs (sol.truncation = tau,
% adjustVCEjoint>hardSolution; WP-133 gave it to residuals) -- the report's
% s0, its chi^2 test and the scaled PAR_TABLE column read R.sigma0 / sqrt(tau)
D.truncation   = lsqvce.util.getfielddef(sol, 'truncation', 1);
D.s0sq         = s0sq_direct / D.truncation;
D.s0           = sqrt(max(D.s0sq, 0));
D.sigma0       = sol.sigma0;
D.varRatio     = D.s0sq / sol.sigma0^2;
D.s0Fallback   = s0Fallback;

% --- robust scale and the contamination symptom
%   On clean data s0 and the MAD scale both estimate sigma0 and agree.
%   Blunders inflate s0 but barely move the MAD, so their RATIO is a direct,
%   quantitative symptom of contamination -- and it is available even when the
%   user ran the (contaminable) Pope test, which is exactly when they need it.
D.robustScale   = robustScale;
D.madAvailable  = madAvailable;
D.scaleRatio = NaN;
if madAvailable && isfinite(D.s0) && robustScale > 0, D.scaleRatio = D.s0 / robustScale; end
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

% ---- fhat's null law and the escalation (WP-135; private contaminationOf, its
%      own file since WP-138)
D = contaminationOf(D, fhat, nScored, madAvailable);

D.sumR       = sumR;
D.sumRExpect = sumR_expect;
D.sumRRelErr = sumR_err;

D.test        = testName;  D.store = o.store;   % store: what R.reliability can read
D.nWorst      = o.nWorst;  D.robustCheck = o.robustCheck;   % WP-135: snoop's re-diagnosis keeps them
D.alpha       = o.alpha;
D.alphaEff    = alphaEff;
D.bonferroni  = o.bonferroni;
D.beta        = o.beta;
D.kCrit       = kCrit;
D.delta0      = delta0;
D.tauBound    = tauBound;
D.vacuityShare = vacuityShare;     % WP-135: the share of v'Pv a residual must carry
D.snoopVacuous = snoopVacuous;

D.flagIdx = flagIdx;
D.flagW   = flagW;
D.flagR   = flagR;
D.flagV   = flagV;
D.flagMDB = flagMDB;
D.nFlag   = numel(flagIdx);
D.nExpectedFalse = alphaEff * nTest;  % how many flags are chance alone
% WP-127: n less the prior's datum rows, which are not tested; WP-133: the
% rows of weight 0 (snoop's own D.rejectedIdx are its deletions)
D.nTests = nTest;  D.datumRows = sort([exIdx(:); locIdx]);  D.excludedIdx = exclIdx;
D.priorRows = priIdx;                          % WP-135: pseudo-observations, not data

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
    'untested',     nQ_untested, ...   % |w| low,  r low  -> the SILENT quadrant
    'rejected',     numel(exclIdx));   % weight 0: not an observation (WP-133)

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

% --- datum / free-network reporting: a datum is a prior inside Q (above)
D.defect  = 0;
D.rank    = lsqvce.util.getfielddef(sol, 'rank', u);   % WP-130: < u when a local block is rank-deficient
D.datumB  = zeros(u, 0);

if storeAll
    D.v   = v_all;
    D.rNum = r_all;
    D.w   = w_all;
    D.qpp = qv_all;      % (P Qvv P)_ii -- the denominator of the test statistic
    D.qv  = qv_all;      % backwards-compatible alias
    D.mdb = mdb_all;
end

% The strongest parameter correlation, from column panels of Qxx: the u x u
% correlation matrix is never held (it was 1.5 GB at u = 13 729).
[cmax, cpair] = Q.maxCorrelation();
D.corrMax    = cmax;
D.corrMaxPar = cpair;

end

function [sol, D, REJ] = snoop(sol, D, batchFcn, opt)
%LSQVCE.ENGINE.VERDICT.SNOOP  Iterative data snooping: reject outliers one at a time.
%
%   [sol, D, REJ] = LSQVCE.ENGINE.VERDICT.SNOOP(sol, D, batchFcn) repeatedly removes the
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
%   rejection. Instead we downdate the cofactor, N -> N - p_i a_i a_i', on its
%   BLOCKS (LSQVCE.ENGINE.BLOCKCOFACTOR, WP-86 part 2): the Schur complement of
%   a rank-one downdate is a rank-one downdate, so the commons factor takes
%   one cholupdate(Rc, x, '-') at O(nc^2) and the observation's own system its
%   O(nl^2 + nl*nc) local update -- no u x u factor exists. The estimate
%   follows in closed form, x <- x + p_i v_i Qxx~ a_i (Qxx~ after the
%   downdate), so no pass over the data is needed to re-solve. Verified
%   against a dense re-solve to 1e-15 (mirror wp86/blockcofactor.py).
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
%     * AFTER:  verify the smallest |diag(Rc)| has not collapsed and the
%               observation is not its local block's only information, and
%               roll back if either fails.
%
%   NAME-VALUE
%     'maxReject'    cap on rejections (default: max(10, 0.01*n))
%     'rMinReject'   refuse to reject below this redundancy (default 0.10)
%     'alpha'        significance per test (default 1e-3)
%     'bonferroni'   correct for n tests (default true)
%     'reaccumulate' re-accumulate from scratch at the end to purge downdate
%                    roundoff and produce a clean final D (default true): the
%                    original block-arrow structure is rebuilt -- each system
%                    from its surviving batches, with its own locals
%     'verbose'      print the rejection sequence (default true)
%
%   OUTPUT REJ
%     .idx        indices rejected, in the order they were rejected
%     .w          their normalised residual at the moment of rejection
%     .r          their redundancy number at that moment
%     .refused    indices that tripped the test but were REFUSED (low r_i)
%     .nIter      iterations performed
%     .converged  true if snooping terminated cleanly
%     .history    per-iteration [nFlag, maxAbsW, s0sq, kCrit] (kCrit: the
%                 critical value tested against, WP-130)
%
%   EXAMPLE
%       data = struct('batchFcn', src, 'numUnknowns', u, 'numBatches', 4);
%       R = lsqvce.fit(task, Verbose=false);
%       [sol, D, REJ] = lsqvce.engine.verdict.snoop(R.source, R.diag, R.batchFcn);
%
%       % REJ.idx       -> [40 150]     rejected, one at a time
%       % REJ.refused   -> observations too weakly controlled to judge
%       % REJ.history   -> the full rejection sequence, for the record
%
%       % Under heavy contamination Pope cannot see the blunders at all, so
%       % this loop would exit having rejected NOTHING and report 'converged'.
%       % Use the robust scale there:
%       [sol, D, REJ] = lsqvce.engine.verdict.snoop(R.source, R.diag, R.batchFcn, 'test', 'mad');
%
%   INPUT
%     sol       1x1 struct  the solution to start from (R.source): x, Q (its
%                           LSQVCE.ENGINE.BLOCKCOFACTOR), n, u, rank, defect
%     D         1x1 struct  its diagnostics
%     batchFcn  1x1 fh      the SAME batch function, for the re-solves
%     Name/Value: 'Test' ('pope'|'mad'), 'Alpha', 'Bonferroni', 'MaxReject',
%                 'MaxRejectFrac' (the breakdown cap, default 0.25), 'RMinReject',
%                 'Reaccumulate', 'NBatch' (default D.nBatch), 'Force', 'Verbose'
%
%   ACCURACY
%     The rank-one downdates of the block cofactor and the closed-form x
%     against a fresh QR solve on the survivors: 1e-8 relative after three
%     rejections (test_vceFunction, Reaccumulate false against true); the
%     clean re-solve keeps a regularised local block's Reg, to 1e-10 against
%     the regularised solve with the rejected rows at weight 0
%     (test_reliability). Its critical values are diagnose's. WP-138.
%
%   ERRORS
%     lsqvce:snoop:correlatedBatch a batch carries a full covariance, so its
%                                  observations cannot be rejected individually
%     lsqvce:snoop:tooContaminated the contamination gauge is past the point
%                                  where tau is meaningful -- switch to
%                                  'test','mad', or use LSQVCE.ENGINE.ADJUSTROBUST
%     lsqvce:snoop:breakdown       the contamination level is 'breakdown': fhat
%                                  collapsed below its clean null while the
%                                  backstop is large -- the robust scale has
%                                  adopted the blunders ('Force' warns instead)
%   WARNINGS
%     lsqvce:snoop:indeterminate   the level is 'indeterminate' (n too small to
%                                  tell a captured robust scale from a sigma0
%                                  declared too small): snoop proceeds under its
%                                  hard cap (WP-137, user decision)
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-137: the level 'indeterminate'
%   (WP-135) has its case -- a warning, lsqvce:snoop:indeterminate, and snoop
%   proceeds (it ran on in silence); ERRORS names breakdown's real trigger.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-133: a row of weight 0 (p = 0,
%   sol.nExcluded) is no rejection candidate and is out of the Bonferroni n,
%   the caps and the MAD pool (its score was a 0 there); rdf is sol.r less the
%   rejections, so a ridge's defect stays in. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-130: the Bonferroni n and the caps
%   leave the datum rows out (diagnose's nTests); 'mad' forms no datum score;
%   sol (n, r, vPv, s0sq, datumRows) is the solution after rejection on both
%   paths, and Reaccumulate=false re-diagnoses D on the downdated cofactor over
%   the survivors (.rejectedIdx); REJ.history carries kCrit. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-128: the clean re-solve keeps a
%   prior as a prior (it went in as data rows: a datum's cofactor grew with its
%   sigma); a prior row is never a rejection candidate; the returned D counts
%   in the ORIGINAL observation indices (D.obsIdx maps its per-observation
%   vectors); Test 'mad' uses madCrit. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
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
%     lsqvce:snoop:indeterminate                 [warning] snoop
%     lsqvce:snoop:lowRedundancy                 [warning] Observation X has the largest |w| (X) but its redundancy is only X. T...
%     lsqvce:snoop:maxReject                     [warning] Reached the rejection cap (X)
%     lsqvce:snoop:rejectFraction                [warning] STOPPING
%
%   See also LSQVCE.ENGINE.ADJUSTROBUST, LSQVCE.ENGINE.VERDICT.DIAGNOSE, LSQVCE.ENGINE.VERDICT.SYSTEMATIC.
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-138 (code review of WP-135):
%   ACCURACY; the re-diagnosis and the candidate floor read D's thresholds
%   (they restated diagnose's defaults). With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-14 (CEST, UTC+02:00), WP-135 (review 3 finding 5, lane 1): a
%   local prior group's row is never a candidate (the common prior's rule,
%   WP-128), its datum rows count like the common prior's; candidates are the
%   rows diagnose calls controlled (r > D.rMin, not 1e-6); the re-diagnosis
%   keeps D's Store / NWorst / RobustCheck and the re-solve the solve's Backend
%   and Reg. With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-86 part 3 (R7-24): the clean
%   re-solve rebuilds the block-arrow system (each system from its surviving
%   batches, with its locals) instead of re-solving every parameter as a
%   common one. With assistance from Claude (Anthropic), model Claude Opus
%   5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-86 part 2 (R7-24): the
%   downdate runs on the block cofactor (sol.Q) and the estimate follows in
%   closed form; the row search and the b re-accumulation passes, and the
%   badBatch / badIndex refusals of the retired helpers, are gone; the clean
%   re-solve reads one batch at a time. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-118 part 1 (R7-34, R7-35):
%   the runtime remedy names live units; the Name-Value list matches the
%   arguments block; the re-diagnosis inherits D's thresholds. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
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
% WP-130: the datum rows; WP-133: the rows of weight 0, replayed (p = 0) but no
% observations -- not candidates (scanWorst), not in the ceiling or Bonferroni's n
dRows = lsqvce.util.getfielddef(sol, 'datumRows', zeros(0, 1));
% WP-135 (review 3 finding 5): D.datumRows adds the local prior groups' datum
% rows (diagnose finds them; sol carries the common prior's only); no prior
% row, common or local, is a candidate; candidates are rows above D.rMin
[dRowsAll, priorBatches, rMinCand] = snoopScope(sol, D, dRows);
nExcl = lsqvce.util.getfielddef(sol, 'nExcluded', 0);  nObs0 = n0 - numel(dRowsAll) - nExcl;
maxRej = o.maxReject;  if isempty(maxRej), maxRej = max(10, ceil(0.01 * nObs0)); end

if sol.defect > 0 && o.verbose
    lsqvce.util.log.msg(['\nNote: this is a free network (datum defect %d). Rejection is\n' ...
             'still valid -- residuals are datum-invariant -- but rejecting an\n' ...
             'observation can INCREASE the defect if it was the last one tying\n' ...
             'a parameter down. That is caught by the redundancy guard below.\n'], ...
             sol.defect);
end

Q    = sol.Q;                      % the cofactor, as blocks (WP-86 part 2)
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

guardContamination(D, o, useRobust, nObs0);   % WP-137: lifted out of snoop

% The hard ceiling. Independent of every estimator above.
maxRejFrac = o.maxRejectFrac;
hardCap = floor(maxRejFrac * nObs0);

for it = 1:(maxRej + 1)

    % ---- find the worst CURRENTLY-flagged observation, using the current R.
    % One streaming pass. We do not need full diagnostics here, only the worst
    % candidate, so this is cheaper than a call to diagnose().
    nAlive = n0 - nExcl - sum(rejected);  rdf = sol.r - sum(rejected);   % WP-133: sol.r has a ridge's defect
    if rdf <= 0
        if lsqvce.util.verbosity.progress(o.verbose)
            lsqvce.util.log.msg('  Redundancy exhausted. Stopping.\n');
        end
        break;
    end

    alphaEff = o.alpha;
    if o.bonferroni, alphaEff = o.alpha / max(nAlive - numel(dRowsAll), 1); end   % WP-130: diagnose's nTests
    kCrit = snoopCrit(useRobust, alphaEff, rdf);   % WP-128: madCrit for 'mad', tau otherwise

    [bestW, bestI, bestR, bestV, bestP, vPvCur, nFlag, bestA, sizes, sysOf] = ...
        scanWorst(batchFcn, nB, Q, sol.x, rejected, rdf, kCrit, useRobust, priorBatches, dRowsAll, rMinCand);

    s0cur = sqrt(max(vPvCur / max(rdf, 1), 0));
    REJ.history(end+1, :) = [nFlag, bestW, s0cur^2, kCrit]; %#ok<AGROW> (kCrit: WP-130)

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
                 'Run R.systematic() to find the structure you have not ' ...
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
             'Run R.systematic() before rejecting anything else.'], ...
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

    % ---- downdate the cofactor's blocks (WP-86 part 2).
    % GUARD 2: Q.downdate refuses (ok = false, Q unchanged) when the removal
    % would make the system singular -- cholupdate's own flag, a commons
    % factor whose smallest diagonal collapsed (a downdate can "succeed" and
    % still leave ~0 there, after which every solve is garbage), or an
    % observation that is its local block's only information.
    [Qnew, ok] = Q.downdate(bestA, bestP);
    if ~ok
        REJ.refused(end+1) = bestI; %#ok<AGROW>
        lsqvce.util.log.warn('lsqvce:snoop:downdateFailed', ...
            ['Removing observation %d would make the normal matrix singular ' ...
             '(the Cholesky downdate collapsed). That observation carries the ' ...
             'only information about at least one parameter -- its redundancy ' ...
             'is effectively zero, whatever the nominal value says.\n' ...
             'REFUSING and stopping.'], bestI);
        break;
    end
    Q = Qnew;

    % ---- vPv update, with the essential 1/r_i factor (see header).
    vPvCur = vPvCur - bestP * bestV^2 / bestR;

    rejected(bestI) = true;
    REJ.idx(end+1) = bestI;   %#ok<AGROW>
    REJ.w(end+1)   = bestW;   %#ok<AGROW>
    REJ.r(end+1)   = bestR;   %#ok<AGROW>

    % ---- the estimate without observation i, in closed form: with v_i = a_i'x - l_i,
    % x~ = x + p_i v_i Qxx~ a_i (Sherman-Morrison; Qxx~ is the downdated cofactor).
    sol.x = sol.x + bestP * bestV * Q.qxxTimes(bestA);

    if lsqvce.util.verbosity.progress(o.verbose)
        s0new = sqrt(max(vPvCur / max(rdf - 1, 1), 0));
        lsqvce.util.log.msg('  %4d %10d %12.4f %10.4f %12d %10.4f   rejected\n', ...
                it, nFlag, bestW, bestR, bestI, s0new);
    end
end

REJ.nIter     = size(REJ.history, 1);
REJ.converged = converged;
REJ.rejected  = rejected;
sol.Q         = Q;               % the downdated cofactor, consistent with sol.x
if ~isempty(REJ.idx), [sol, D] = afterRejection(sol, D, REJ.idx, rejected, vPvCur, dRows, batchFcn, nB, sizes, o); end

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
    [acc, obsIdx, order] = blockResolveAccumulator(sol, Q, batchFcn, nB, rejected, sizes, sysOf);
    solB = acc.solve();
    [dsol, bf, nBd] = lsqvce.engine.verdict.source(acc, solB);
    % WP-135: the re-solve's batches are the kept originals in a new order; the
    % local prior groups' batches keep their role there
    dsol.localPriorBatches = find(ismember(order, lsqvce.util.getfielddef(sol, 'localPriorBatches', [])));
    % R7-35: the re-diagnosis classifies with the SAME thresholds as the D it
    % was handed (rMin / rWeak / beta), not diagnose's defaults -- before, R.diag
    % and snoop's D could call one observation controlled and uncontrolled.
    % WP-135: and with its Store / NWorst / RobustCheck (Store "flagged" came
    % back "all": five n-sized vectors at the scale it was chosen to avoid)
    nvR = rediagnoseNV(D, o);
    D   = lsqvce.engine.verdict.diagnose(dsol, bf, 'nBatch', nBd, nvR{:});
    D = toOriginalIndex(D, obsIdx);  D.rejectedIdx = REJ.idx;   % WP-128: the original indices
    % the returned sol keeps the caller's fields and carries the clean solve
    sol.x = dsol.x(:);  sol.s0sq = solB.s0sq;  sol.sigma0 = solB.sigma0;   % commons + locals
    if isfield(solB, 'Cxx'), sol.Cxx = solB.Cxx; end
    sol.n = dsol.n;  sol.dof = solB.dof;  sol.Q = dsol.Q;  sol.r = dsol.r;  sol.vPv = dsol.vPv;
    sol.datumRows = dsol.datumRows;                    % WP-130: r, vPv, datum rows follow too
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

function S = systematic(D, opt)
%LSQVCE.ENGINE.VERDICT.SYSTEMATIC  Detect unmodelled systematic errors in the residuals.
%
%   S = LSQVCE.ENGINE.VERDICT.SYSTEMATIC(D) runs three complementary tests on the normalised
%   residuals from lsqvce.engine.verdict.diagnose. This closes the single largest blind spot
%   in the rest of the toolbox.
%
%   WHY THIS EXISTS
%   ---------------
%   A smooth unmodelled signal -- a thermal drift, an unmodelled bias that
%   varies with elevation, a once-per-revolution term you forgot -- produces NO
%   single large residual. Data snooping therefore reports a perfectly clean
%   result while the parameters are quietly corrupted. Meanwhile the global
%   chi-square test rejects, because the residuals ARE too large in aggregate.
%
%   That combination -- inflated variance factor, zero outliers -- is the
%   signature of a systematic error, and until now the toolbox could only tell
%   you to go and plot something. Demonstrated on a sinusoid of amplitude 0.8
%   sigma: s0^2 = 1.33 (inflated), Bonferroni snooping flagged 0 observations,
%   and the three tests below returned p = 5e-5, p < 1e-16 and p = 3e-9.
%
%   THE THREE TESTS
%   ---------------
%   (a) RUNS TEST on sign(w). Distribution-free. Counts sign changes: a drift
%       produces too few. Weak on smooth signals, strong on steps.
%
%   (b) LJUNG-BOX on w, lags 1..h. Detects autocorrelation of any shape.
%       Powerful, but see the caveat below.
%
%   (c) BINNED MEANS of w against each supplied covariate. In bin b the mean
%       of w has standard deviation 1/sqrt(n_b) under H0, so z_b = mean_b *
%       sqrt(n_b). Bonferroni over the bins. This is effectively a score test
%       for an omitted parameter, and it is the WORKHORSE: in testing it
%       reached 84-100% power where the runs test managed 1-9%.
%
%   (a) and (b) require the residuals to be ORDERED along the systematic --
%   they see the sequence you hand them. (c) does not: it works on any
%   covariate directly. Pass covariates and let (c) do the work.
%
%   CAVEAT ON LJUNG-BOX -- and it is an honest one.
%   Least-squares residuals are NOT independent, even when the observation
%   errors are: cov(v) = sigma0^2 * Qvv, and Qvv is not diagonal. Normalising
%   to w removes the variance but not the correlation. Ljung-Box therefore
%   OVER-REJECTS, and the effect grows with u/n. Measured false-alarm rate at a
%   nominal 1%:
%       u/n = 0.060  ->  1.6%
%       u/n = 0.012  ->  1.3%
%       u/n = 0.003  ->  0.6%   (within noise of nominal)
%   Geodetic problems live at u/n << 0.01, where it is harmless. The toolbox
%   warns when u/n > 0.05, and it leads with the binned test, which is exact.
%
%   NAME-VALUE
%     'covariates'  struct: each field is an n x 1 vector to bin against
%                   (or one entry per DATA row: the priors' pseudo-
%                   observations, D.priorRows, are no residual sequence and
%                   are left out of every test -- WP-135).
%                   e.g. struct('time', t, 'elev', e, 'azim', az)
%                   THIS IS THE IMPORTANT ONE. Supply everything you have.
%     'nbins'       bins per covariate (default 10)
%     'lags'        Ljung-Box lags h (default 10)
%     'alpha'       significance (default 0.01)
%     'verbose'     print (default true)
%
%   OUTPUT S
%     .runs      struct(z, p, nRuns, expected)
%     .ljungBox  struct(Q, p, h, reliable)
%     .bins      struct array, one per covariate:
%                (name, z, p, worstBin, worstZ, edges, counts, means)
%     .detected  true if any test rejects at alpha
%     .worst     name of the most significant finding
%
%   INPUT
%     D         1x1 struct  diagnostics from LSQVCE.ENGINE.VERDICT.DIAGNOSE
%     varargin  Name/Value: 'covariate' (a vector to test against), 'alpha', ...
%
%   OUTPUT
%     S  1x1 struct  .runs, .durbinWatson, .groupOffsets, .verdict -- the
%                    structure tests that distinguish a MODEL defect from
%                    scattered blunders
%
%   EXAMPLE
%       % The failure that snooping CANNOT see. A smooth drift produces no
%       % single large residual: the variance factor inflates, the flag list
%       % stays clean, and THAT COMBINATION is the signature.
%       data = struct('batchFcn', src, 'numUnknowns', u, 'numBatches', 4);
%       R = lsqvce.fit(task, Verbose=false);
%
%       t = (1:n).';
%       S = lsqvce.engine.verdict.systematic(R.diag, 'covariates', ...
%                              struct('time', t, 'elev', elev));
%
%       % S.detected     did anything fire?
%       % S.ljungBox.p   autocorrelation in the residual sequence
%       % S.runs.p       distribution-free sign test
%       % S.bins         one entry per covariate: .name, .p, .worstBin
%
%       % The COVARIATE ranking lives in S.bins -- S.worst ranks every test,
%       % and a drift in time necessarily also trips Ljung-Box.
%       [~, iw] = min([S.bins.p]);
%       lsqvce.util.log.msg('worst covariate: %s\n', S.bins(iw).name);
%
%   ACCURACY
%     The sizes are claims of the statistics, measured by simulation: the
%     binned covariate test is exact under H0 (z_b ~ N(0,1)); Ljung-Box
%     over-rejects with u/n (1.6 % at a nominal 1 % for u/n = 0.06, 0.6 % at
%     0.003 -- the CAVEAT above); the pooled Ljung-Box is divided by the
%     pooled pair counts (R7-38: 0.048 at a nominal 0.05 on 64 arcs of 30,
%     it was 0.001; dev/reviews/R7_2026-10-09/mirrors/wp118/ljungbox_pooled.py).
%     The tails are LSQVCE.KERNEL.STAT's (the "scipy.stats" row of HANDOVER
%     6.1.4). WP-138.
%
%   ERRORS
%     lsqvce:systematic:noResiduals   the diagnostics carry no residuals
%     lsqvce:systematic:badCovariate  the covariate is not a vector matching the
%                                     residual length
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (b) b2 (F-07):
%   typed `arguments` signature, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:systematic:lbUnreliable             [warning] Ljung-Box is anti-conservative here (u/n = X, n = X, h = X)
%     lsqvce:systematic:tooFew                   [warning] Only X controlled observations
%
%   See also LSQVCE.ENGINE.VERDICT.DIAGNOSE, LSQVCE.ENGINE.VERDICT.FTEST.
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-138 (code review of WP-135):
%   ACCURACY. With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-14 (CEST, UTC+02:00), WP-135 (review 3, lane 1): the priors'
%   pseudo-observations (D.priorRows) leave every sequence test, and a covariate
%   may be given over the data rows only. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

% WP-39 stage (b) b2: typed `arguments` signature (canonical
% UpperCamelCase; legacy casings keep matching -- MATLAB NV
% semantics; enum VALUES case-sensitive lowercase).
arguments
    D
    opt.Covariates (1,1) struct = struct()
    opt.NBins (1,1) double {mustBeGreaterThanOrEqual(opt.NBins,2)} = 10
    opt.Lags (1,1) double {mustBeGreaterThanOrEqual(opt.Lags,1)} = 10
    opt.Alpha (1,1) double {mustBePositive, mustBeLessThan(opt.Alpha,1)} = 0.01
    opt.Verbose (1,1) logical = true
end

o = struct('covariates', opt.Covariates, 'nbins', opt.NBins, 'lags', opt.Lags, 'alpha', opt.Alpha, 'verbose', opt.Verbose);

if ~isfield(D, 'w')
    error('lsqvce:systematic:noResiduals', ...
        ['The diagnostics struct has no residual vector. Re-run ' ...
         'the fit with Diagnostics.Store = "all" -- the systematic tests ' ...
         'need every residual, not just the flagged ones.']);
end

w = D.w(:);
n = numel(w);
u = D.u;

% Only observations the network actually controls carry information about a
% systematic. An uncontrolled observation has w = 0 by construction, and
% including those zeros would dilute every test toward "no signal".
if isfield(D, 'rNum')
    ctrl = D.rNum(:) > D.rMin;
else
    ctrl = true(n, 1);
end
% WP-135 (review 3, lane 1): the priors' pseudo-observations (the common
% prior's rows last, a local group's inside its system) are no residual
% sequence -- they entered the runs and Ljung-Box statistics
pr = lsqvce.util.getfielddef(D, 'priorRows', zeros(0, 1));
if isfield(D, 'obsIdx') && ~isempty(D.obsIdx), pr = find(ismember(D.obsIdx, pr)); end   % a snooped D
isData = true(n, 1);  isData(pr) = false;
ctrl = ctrl & isData;
wc = w(ctrl);
nc = numel(wc);

S = struct();
S.n = n; S.nControlled = nc; S.alpha = o.alpha;

if nc < 20
    lsqvce.util.log.warn('lsqvce:systematic:tooFew', ...
        ['Only %d controlled observations. The systematic-error tests are ' ...
         'asymptotic and are meaningless at this size.'], nc);
    S.detected = false;  S.worst = 'n/a';
    S.runs = struct('z', NaN, 'p', NaN);
    S.ljungBox = struct('Q', NaN, 'p', NaN, 'reliable', false);
    S.bins = struct([]);
    return;
end

%% ------------------------------------------------------------- (a) runs test
sg = sign(wc);
sg = sg(sg ~= 0);
ns = numel(sg);
npos = sum(sg > 0);
nneg = ns - npos;

if npos == 0 || nneg == 0 || ns < 3
    S.runs = struct('z', NaN, 'p', NaN, 'nRuns', NaN, 'expected', NaN);
else
    nRuns = 1 + sum(sg(2:end) ~= sg(1:end-1));
    muR   = 2 * npos * nneg / ns + 1;
    varR  = 2*npos*nneg * (2*npos*nneg - ns) / (ns^2 * (ns - 1));
    if varR <= 0
        S.runs = struct('z', NaN, 'p', NaN, 'nRuns', nRuns, 'expected', muR);
    else
        z = (nRuns - muR) / sqrt(varR);
        p = 2 * (1 - lsqvce.kernel.stat.normCdf(abs(z)));
        S.runs = struct('z', z, 'p', p, 'nRuns', nRuns, 'expected', muR);
    end
end

%% -------------------------------------------------------------- (b) Ljung-Box
h = min(o.lags, floor(nc / 5));
if h < 1
    S.ljungBox = struct('Q', NaN, 'p', NaN, 'h', 0, 'reliable', false);
else
    wz = wc - mean(wc);
    c0 = (wz' * wz) / nc;
    if c0 <= 0
        S.ljungBox = struct('Q', NaN, 'p', NaN, 'h', h, 'reliable', false);
    else
        Q = 0;
        for k = 1:h
            ck = (wz(1:end-k)' * wz(k+1:end)) / nc;
            rk = ck / c0;
            Q  = Q + rk^2 / (nc - k);
        end
        Q = nc * (nc + 2) * Q;
        p = 1 - lsqvce.kernel.stat.chi2Cdf(Q, h);

        % Reliability flags, both established by Monte Carlo:
        %   * residual correlation makes LB anti-conservative as u/n grows;
        %   * the asymptotic chi2(h) needs n >> h (at n=100, h=20 the variance
        %     of Q was 52 against a nominal 40).
        ratio    = u / nc;
        reliable = (ratio <= 0.05) && (nc >= 20 * h);

        S.ljungBox = struct('Q', Q, 'p', p, 'h', h, 'reliable', reliable, ...
                            'uOverN', ratio);
        if ~reliable
            lsqvce.util.log.warn('lsqvce:systematic:lbUnreliable', ...
                ['Ljung-Box is anti-conservative here (u/n = %.3f, n = %d, ' ...
                 'h = %d). Least-squares residuals are correlated by ' ...
                 'construction -- cov(v) = s0^2*Qvv is not diagonal -- and the ' ...
                 'over-rejection grows with u/n. Treat its p-value as a screen, ' ...
                 'not a verdict, and lead with the binned test.'], ...
                 ratio, nc, h);
        end
    end
end

%% ----------------------------------------------------------- (c) binned means
names = fieldnames(o.covariates);
bins = struct('name', {}, 'z', {}, 'p', {}, 'worstBin', {}, 'worstZ', {}, ...
              'edges', {}, 'counts', {}, 'means', {});

for ii = 1:numel(names)
    cv = o.covariates.(names{ii});
    cv = cv(:);
    if numel(cv) == sum(isData) && sum(isData) < n   % WP-135: over the DATA rows
        cvAll = nan(n, 1);  cvAll(isData) = cv;  cv = cvAll;
    end
    if numel(cv) ~= n
        error('lsqvce:systematic:badCovariate', ...
            ['Covariate ''%s'' has %d elements but there are %d ' ...
             'observations (%d without the priors'' rows). It must be aligned ' ...
             'with the observation order the batch function produced.'], ...
             names{ii}, numel(cv), n, sum(isData));
    end
    cvc = cv(ctrl);
    good = isfinite(cvc);
    if sum(good) < 2 * o.nbins
        continue;
    end

    % Quantile edges: equal-count bins, so every bin has power. Equal-WIDTH
    % bins would leave the tails nearly empty and waste the Bonferroni budget
    % on bins that cannot reject anything.
    q  = linspace(0, 1, o.nbins + 1);
    ed = quantile(cvc(good), q);
    ed = unique(ed);
    nb = numel(ed) - 1;
    if nb < 2, continue; end
    ed(end) = ed(end) + eps(ed(end)) * 10;

    zb = zeros(nb, 1);  cnt = zeros(nb, 1);  mu = zeros(nb, 1);
    for b = 1:nb
        if b < nb
            m = good & cvc >= ed(b) & cvc < ed(b+1);
        else
            m = good & cvc >= ed(b) & cvc <= ed(b+1);
        end
        cnt(b) = sum(m);
        if cnt(b) < 2
            zb(b) = 0;  mu(b) = 0;
            continue;
        end
        mu(b) = mean(wc(m));
        % Under H0, w ~ N(0,1) so the bin mean has sd 1/sqrt(n_b).
        zb(b) = mu(b) * sqrt(cnt(b));
    end

    [~, ib] = max(abs(zb));
    praw = 2 * (1 - lsqvce.kernel.stat.normCdf(abs(zb(ib))));
    padj = min(1, praw * nb);          % Bonferroni over the bins

    bins(end+1) = struct('name', names{ii}, 'z', zb, 'p', padj, ...
                         'worstBin', ib, 'worstZ', zb(ib), 'edges', ed, ...
                         'counts', cnt, 'means', mu); %#ok<AGROW>
end
S.bins = bins;

%% ------------------------------------------------------------------ verdict
ps = [];  lbl = {};
if isfinite(S.runs.p),     ps(end+1) = S.runs.p;     lbl{end+1} = 'runs'; end
if isfinite(S.ljungBox.p), ps(end+1) = S.ljungBox.p; lbl{end+1} = 'Ljung-Box'; end
for ii = 1:numel(bins)
    ps(end+1) = bins(ii).p; %#ok<AGROW>
    lbl{end+1} = sprintf('bins:%s', bins(ii).name); %#ok<AGROW>
end

if isempty(ps)
    S.detected = false;  S.worst = 'n/a';  S.pMin = NaN;
else
    [pmin, im] = min(ps);
    S.pMin     = pmin;
    S.worst    = lbl{im};
    S.detected = pmin < o.alpha;
end

%% ------------------------------------------------------------------- print
if ~o.verbose, return; end

lsqvce.util.log.msg('\n--- SYSTEMATIC ERROR TESTS  (alpha = %.3g)\n', o.alpha);

if isfinite(S.runs.p)
    lsqvce.util.log.msg('  runs test  : z = %+7.3f   p = %.3e   (%d runs, expected %.1f)%s\n', ...
        S.runs.z, S.runs.p, S.runs.nRuns, S.runs.expected, mark(S.runs.p, o.alpha));
end
if isfinite(S.ljungBox.p)
    rel = '';
    if ~S.ljungBox.reliable
        rel = '   [anti-conservative here: treat as a screen]';
    end
    lsqvce.util.log.msg('  Ljung-Box  : Q = %8.2f   p = %.3e   (h = %d)%s%s\n', ...
        S.ljungBox.Q, S.ljungBox.p, S.ljungBox.h, ...
        mark(S.ljungBox.p, o.alpha), rel);
end

if isempty(bins)
    lsqvce.util.log.msg(['  binned     : no covariates supplied.\n' ...
             '               This is the most powerful test and you have not\n' ...
             '               given it anything to work with. Pass every ordering\n' ...
             '               variable you have:\n' ...
             '                 ''covariates'', struct(''time'',t,''elev'',e,''azim'',az)\n']);
else
    for ii = 1:numel(bins)
        b = bins(ii);
        lsqvce.util.log.msg('  bins:%-8s: worst bin %d/%d, mean w = %+.4f (n=%d), z = %+.3f, p = %.3e%s\n', ...
            b.name, b.worstBin, numel(b.z), b.means(b.worstBin), ...
            b.counts(b.worstBin), b.worstZ, b.p, mark(b.p, o.alpha));
    end
end

lsqvce.util.log.msg('\n');
if S.detected
    lsqvce.util.log.msg(['  ==> SYSTEMATIC ERROR DETECTED (%s, p = %.2e).\n' ...
             '      The residuals carry structure that the model does not.\n' ...
             '      This is NOT an outlier problem and rejecting observations\n' ...
             '      will not fix it -- it will only hide it. Find the missing\n' ...
             '      parameter.\n'], S.worst, S.pMin);
    if ~isempty(bins)
        [~, iw] = min([bins.p]);
        b = bins(iw);
        lsqvce.util.log.msg(['      The signal is strongest against ''%s''. Plot w against\n' ...
                 '      it and look at the shape -- that shape IS the missing\n' ...
                 '      term.\n'], b.name);
    end
else
    lsqvce.util.log.msg(['  ==> No systematic structure detected at alpha = %.3g.\n' ...
             '      Note this is only as good as the covariates you supplied:\n' ...
             '      a systematic in a variable you did not pass is invisible\n' ...
             '      to every test here.\n'], o.alpha);
end
lsqvce.util.log.msg('\n');
end

function F = ftest(solFull, solRed, opt)
%LSQVCE.ENGINE.VERDICT.FTEST  Is the extra parameter earning its place?
%
%   F = LSQVCE.ENGINE.VERDICT.FTEST(solFull, solRed) tests whether the additional parameters
%   in the FULL model give a statistically significant improvement in fit over
%   the REDUCED (nested) model.
%
%       F = [ (vPv_red - vPv_full) / (u_full - u_red) ] / s0sq_full
%           ~ F(u_full - u_red,  r_full)   under H0 (the extra params are zero)
%
%   The reduced model must be strictly NESTED inside the full one -- same
%   observations, same weights, a subset of the parameters. This is checked as
%   far as it can be: n must match, and vPv_red >= vPv_full (adding freedom can
%   only improve the fit, so a reduced model that fits BETTER means the two
%   solutions are not what you think they are).
%
%   WHY YOU NEED THIS
%   -----------------
%   The report tells you that a LOW variance factor may mean over-
%   parameterisation -- the model has enough freedom to absorb the noise, which
%   produces a beautiful fit and meaningless parameters. That warning is useless
%   without a way to act on it. This is the way to act on it.
%
%   The temptation is to look at the t-statistic x_j / sigma_j for each
%   parameter and drop the insignificant ones. That is fine for ONE parameter,
%   but it is wrong for several at once when they are correlated: two parameters
%   can each be individually insignificant while their JOINT contribution is
%   highly significant (and vice versa). The F-test asks the right question --
%   does this GROUP of parameters, taken together, earn its keep.
%
%   INPUT
%     solFull  1x1 struct  solution of the FULL model (from lsqvce.engine.adjustBlock /
%                          Accumulator.solve); uses .vtv (or .vPv), .m/.n,
%                          .u/.rank, .r (the redundancy; n - rank without
%                          it) and .sigma0/.s0 as available.
%     solRed   1x1 struct  solution of the nested REDUCED model, same fields.
%
%   NAME-VALUE
%     'alpha'    significance (default 0.05 -- model selection is not blunder
%                detection, and 0.001 here would retain almost nothing)
%     'verbose'  print (default true)
%
%   OUTPUT F
%     .F, .p, .df1, .df2, .crit
%     .significant   true if the extra parameters are justified
%     .dvPv          the improvement in vPv
%
%   EXAMPLE
%       % Is the extra parameter justified, or is it fitting noise?
%       Rfull = lsqvce.fit(lsqvce.TaskDefinition.fromMatrices([A extra], l, Weights=p), Verbose=false);
%       Rred  = lsqvce.fit(lsqvce.TaskDefinition.fromMatrices(A, l, Weights=p), Verbose=false);
%
%       F = lsqvce.engine.verdict.ftest(Rfull.sol, Rred.sol);   % = Rfull.compare(Rred)
%
%       % F.significant   false for a pure-noise column (measured p = 0.548)
%       % F.p             the p-value
%
%       % The test assumes the REDUCED model is correct. If a systematic error
%       % is present, both models are wrong and the p-value means nothing --
%       % pass 'checkD' and it will refuse:
%       F = lsqvce.engine.verdict.ftest(sF, sR, 'checkD', S);
%
%   ACCURACY
%     F = ((vPv_red - vPv_full) / df1) / (vPv_full / df2) is rounding only;
%     the p-value is betainc's UPPER tail (the 1 - lower form lost every p
%     below ~1e-16, WP-135) and the critical value LSQVCE.KERNEL.STAT.FINV
%     (betaincinv). Against the closed form of F(2, r), P(F > f) =
%     (1 + 2f/r)^(-r/2) and F_crit = (r/2)(alpha^(-2/r) - 1): p to 1e-12,
%     crit to 1e-10 relative (test_errorDoctrine,
%     testFtestMatchesTheClosedFormOfF2). WP-138.
%
%   ERRORS
%     lsqvce:ftest:notNested       the two models are not nested, so the F ratio
%                                  has no meaning
%     lsqvce:ftest:differentData    they were not fitted to the same observations
%     lsqvce:ftest:noRedundancy     the fuller model has no degrees of freedom left
%     lsqvce:ftest:zeroVariance     the residual sum of squares is zero
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (b) b2 (F-07):
%   typed `arguments` signature, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:ftest:invalidPremise                [warning]  ... ... ... ... ... 
%
%   See also LSQVCE.ENGINE.VERDICT.SYSTEMATIC, LSQVCE.ENGINE.ADJUSTBLOCK.
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-138 (code review of WP-135):
%   ACCURACY, the p-value and the critical value pinned against the closed
%   form of F(2, r). With assistance from Claude (Anthropic), model Claude
%   Opus 5.5.
%   Updated 2026-10-14 (CEST, UTC+02:00), WP-135 (review 3, lane 1): the p-value
%   is betainc's upper tail (1 - lower lost every p below ~1e-16). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-118 part 1 (R7-34): the
%   deleted engine.adjust reference replaced. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

% WP-39 stage (b) b2: typed `arguments` signature (canonical
% UpperCamelCase; legacy casings keep matching -- MATLAB NV
% semantics; enum VALUES case-sensitive lowercase).
arguments
    solFull
    solRed
    opt.Alpha (1,1) double {mustBePositive, mustBeLessThan(opt.Alpha,1)} = 0.05
    opt.Verbose (1,1) logical = true
    opt.CheckD = []
end

o = struct('alpha', opt.Alpha, 'verbose', opt.Verbose, 'checkD', opt.CheckD);

uF = solFull.rank;   if isempty(uF), uF = solFull.u; end
uR = solRed.rank;    if isempty(uR), uR = solRed.u;  end

if solFull.n ~= solRed.n
    error('lsqvce:ftest:differentData', ...
        ['The two models were fitted to different numbers of observations ' ...
         '(%d and %d). An F-test compares two models on the SAME data; ' ...
         'otherwise the difference in vPv reflects the data, not the model.'], ...
         solFull.n, solRed.n);
end

df1 = uF - uR;
% WP-135: the full model's redundancy where it carries one (sol.r: the rows
% of weight 0 and a ridge's defect counted, WP-133) -- n - rank otherwise
df2 = lsqvce.util.getfielddef(solFull, 'r', solFull.n - uF);

if df1 <= 0
    error('lsqvce:ftest:notNested', ...
        ['The "full" model has %d estimable parameters and the "reduced" one ' ...
         'has %d. The full model must have MORE.'], uF, uR);
end
if df2 <= 0
    error('lsqvce:ftest:noRedundancy', ...
        'The full model has no redundancy (n = %d, rank = %d).', solFull.n, uF);
end

dvPv = solRed.vPv - solFull.vPv;

if dvPv < 0
    if abs(dvPv) < 1e-8 * max(solFull.vPv, 1)
        dvPv = 0;      % roundoff
    else
        error('lsqvce:ftest:notNested', ...
            ['The REDUCED model fits BETTER than the full one ' ...
             '(vPv %.6e vs %.6e). That is impossible for genuinely nested ' ...
             'models -- adding parameters cannot worsen a least-squares fit.\n' ...
             'The two models are not nested, or they were not fitted to the ' ...
             'same data with the same weights, or one of them has not ' ...
             'converged.'], solRed.vPv, solFull.vPv);
    end
end

s0sqF = solFull.vPv / df2;
if s0sqF <= 0
    error('lsqvce:ftest:zeroVariance', 'The full model has zero residual.');
end

Fstat = (dvPv / df1) / s0sqF;
crit  = lsqvce.kernel.stat.fInv(1 - o.alpha, df1, df2);
% p-value from the F CDF, via the incomplete beta
xb = (df1 * Fstat) / (df1 * Fstat + df2);
pval = betainc(xb, df1/2, df2/2, 'upper');   % WP-135: 1 - lower lost every p below ~1e-16

F = struct('F', Fstat, 'p', pval, 'df1', df1, 'df2', df2, ...
           'crit', crit, 'alpha', o.alpha, 'dvPv', dvPv, ...
           'significant', Fstat > crit, ...
           'vPvFull', solFull.vPv, 'vPvRed', solRed.vPv);

% ---------------------------------------------------------------------------
% VALIDITY GUARD
% ---------------------------------------------------------------------------
% The F-test asks "does the extra parameter improve the fit MORE THAN CHANCE?"
% and answers it by comparing dvPv against the full model's residual variance.
% That comparison is only meaningful if the residual variance is NOISE.
%
% If a systematic error is present, or the data are contaminated, then s0sqF
% measures the systematic error rather than the noise, and the F-statistic is
% the ratio of two things that are both wrong. It will happily declare a
% parameter "significant" because it soaks up part of an unmodelled signal --
% which is precisely how a spurious parameter gets into a model and stays there.
F.valid = true;
F.validityWarning = '';
if ~isempty(o.checkD) && isstruct(o.checkD)
    bad = '';
    if isfield(o.checkD, 'contamLevel') && ...
       any(strcmp(o.checkD.contamLevel, {'moderate','severe','refuse','breakdown'}))
        bad = sprintf('contamination (level "%s")', o.checkD.contamLevel);
    elseif isfield(o.checkD, 'detected') && o.checkD.detected
        bad = 'an unmodelled systematic error';
    end
    if ~isempty(bad)
        F.valid = false;
        F.validityWarning = bad;
        lsqvce.util.log.warn('lsqvce:ftest:invalidPremise', ...
            ['ftest: the residuals show %s, so this F-test is NOT valid.\n' ...
             'The F-statistic divides the improvement in fit by the full ' ...
             'model''s residual variance, which is only meaningful when that ' ...
             'variance is NOISE. Here it is measuring the model error itself, ' ...
             'so the test compares two wrong things.\n' ...
             'A parameter can come out "significant" simply because it soaks up ' ...
             'part of an unmodelled signal -- which is exactly how a spurious ' ...
             'parameter enters a model and stays there. Fix the model first.'], ...
             bad);
    end
end

if lsqvce.util.verbosity.progress(o.verbose)
    lsqvce.util.log.msg('\n--- F-TEST: are the %d extra parameter(s) justified?\n', df1);
    lsqvce.util.log.msg('  vPv (reduced, %d par) = %.6e\n', uR, solRed.vPv);
    lsqvce.util.log.msg('  vPv (full,    %d par) = %.6e\n', uF, solFull.vPv);
    lsqvce.util.log.msg('  improvement          = %.6e\n', dvPv);
    lsqvce.util.log.msg('  F = %.4f   ~ F(%d, %d)   critical = %.4f   p = %.4e\n', ...
            Fstat, df1, df2, crit, pval);
    lsqvce.util.log.msg('\n');
    if F.significant
        lsqvce.util.log.msg(['  ==> SIGNIFICANT. The extra parameter(s) reduce the residual\n' ...
                 '      by more than chance would. Keep them.\n']);
    else
        lsqvce.util.log.msg(['  ==> NOT significant. The extra parameter(s) buy an\n' ...
                 '      improvement no bigger than noise. They are absorbing\n' ...
                 '      random error, not signal -- which makes the fit look\n' ...
                 '      better while making the parameters MEAN less. Drop them.\n']);
    end
    lsqvce.util.log.msg('\n');
end
end

function E = ellipse(D, pair, opt)
%LSQVCE.ENGINE.VERDICT.ELLIPSE  Confidence ellipse for a pair of correlated parameters.
%
%   E = LSQVCE.ENGINE.VERDICT.ELLIPSE(D, [i j]) returns the confidence ellipse of parameters
%   i and j from the 2x2 submatrix of Qxx.
%
%   WHEN THE REPORT SAYS "look at the error ellipse instead", THIS IS IT.
%
%   Two parameters correlated at 0.99 may each have a small standard deviation
%   while being almost completely inseparable. The individual sigmas describe
%   the two SIDES of a bounding box; the ellipse inside it can be a needle. The
%   sigmas are then true but deeply misleading: they suggest you know each
%   parameter to that precision, when in fact you only know one combination of
%   them well and the orthogonal combination hardly at all.
%
%   The ellipse says what you actually know. Its semi-major axis is the
%   direction you know WORST, and its length is how badly.
%
%   INPUT
%     D     1x1 struct  diagnostics from lsqvce.engine.verdict.diagnose (R.diag); uses
%                       the parameter estimate and the cofactor D.Q plus sigma0/s0.
%     pair  1x2 double  the two parameter indices [i j] to plot the ellipse for.
%
%   NAME-VALUE
%     'conf'     confidence level (default 0.95)
%     'scale'    'apriori' (use sigma0) | 'aposteriori' (use s0)  [default
%                'aposteriori', matching what the report tabulates]. With
%                the estimated s0 the region is 2 F(2, r) (r = D.r, the
%                redundancy), with sigma0 chi2(2) -- R7-39
%     'npts'     points on the returned outline (default 200)
%
%   OUTPUT E
%     .semiMajor, .semiMinor    axis half-lengths at the requested confidence
%     .k2                       the squared scale used: 2 F(2, r) or chi2(2)
%     .angle                    orientation of the major axis, radians CCW from
%                               the +i axis
%     .ratio                    semiMajor/semiMinor -- the elongation. Above
%                               ~10 the parameters are effectively one.
%     .xy                       2 x npts outline, centred on the estimate
%     .centre                   [x_i; x_j]
%     .corr                     the correlation between them
%     .sigInd                   [sigma_i; sigma_j] -- the misleading numbers
%     .worstDir                 unit vector of the semi-major axis: the
%                               combination of the two parameters you know least
%
%   EXAMPLE
%       data = struct('batchFcn', src, 'numUnknowns', u, 'numBatches', 4);
%       R = lsqvce.fit(task, Verbose=false);
%       E = lsqvce.engine.verdict.ellipse(R.diag, [1 2]);       % parameters 1 and 2
%
%       plot(E.xy(:,1), E.xy(:,2), 'LineWidth', 1.5);  hold on
%       plot(E.centre(1), E.centre(2), 'r+');
%       axis equal; grid on
%
%       % E.xy         npts x 2, ready for plot()
%       % E.semiMajor  / E.semiMinor / E.angle
%       % E.ratio      elongation -- large means the pair is poorly separated
%
%       % Draw this whenever |rho| > 0.95: the individual sigmas are then
%       % meaningless in isolation, however small they look.
%
%   ACCURACY
%     k2 = 2 F(2, r) (a posteriori) or chi2(2) (a priori), both against
%     their closed forms for two dimensions -- r ((1 - p)^(-2/r) - 1) and
%     -2 log(1 - p) -- to 1e-10 relative; the semi-axes sqrt(k2 s0^2) of a
%     unit cofactor to 1e-12 (test_lsqvce_opts,
%     test_ellipse_uses_2F_with_the_estimated_variance_factor). WP-138.
%
%   ERRORS
%     lsqvce:ellipse:badPair        the parameter pair is not two valid indices
%     lsqvce:ellipse:samePair       both indices are the same parameter
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (b) b2 (F-07):
%   typed `arguments` signature, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   See also LSQVCE.ENGINE.VERDICT.DIAGNOSE, LSQVCE.ENGINE.BLOCKCOFACTOR.
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-138 (code review of WP-135):
%   ACCURACY. With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-86 part 2 (R7-24): the 2 x 2
%   cofactor comes from D.Q.block (the dense D.Qxx is retired). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-118 part 1 (R7-34, R7-39): 2
%   F(2, r) with the estimated s0^2 (chi2(2) under-covered at small r);
%   output .k2. With assistance from Claude (Anthropic), model Claude Opus
%   5.5.
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

% WP-39 stage (b) b2: typed `arguments` signature (canonical
% UpperCamelCase; legacy casings keep matching -- MATLAB NV
% semantics; enum VALUES case-sensitive lowercase).
arguments
    D
    pair
    opt.Conf (1,1) double {mustBePositive, mustBeLessThan(opt.Conf,1)} = 0.95
    opt.Scale (1,1) string {mustBeMember(opt.Scale, ["apriori","aposteriori"])} = "aposteriori"
    opt.NPts (1,1) double {mustBeGreaterThanOrEqual(opt.NPts,8)} = 200
end

o = struct('conf', opt.Conf, 'scale', char(opt.Scale), 'npts', opt.NPts);

if numel(pair) ~= 2 || any(pair < 1) || any(pair > D.u) || any(pair ~= fix(pair))
    error('lsqvce:ellipse:badPair', ...
          'pair must be two integer parameter indices in 1..%d.', D.u);
end
i = pair(1);  j = pair(2);
if i == j
    error('lsqvce:ellipse:samePair', 'The two indices must differ.');
end

Q = D.Q.block([i j]);          % Qxx([i j], [i j]) from the block cofactor (WP-86 part 2)
Q = (Q + Q') / 2;

estimated = false;
if strcmpi(o.scale, 'aposteriori')
    fac = D.s0sq;
    estimated = isfinite(fac) && fac > 0;
    if ~estimated, fac = D.sigma0^2; end
else
    fac = D.sigma0^2;
end
C = fac * Q;

[V, Lam] = eig(C);
lam = diag(Lam);
[lam, ix] = sort(lam, 'descend');
V = V(:, ix);

if any(lam < 0)
    lam = max(lam, 0);
end

% Chi-square with 2 dof gives the scale factor for a joint confidence region
% when the variance factor is KNOWN. (Using 1 dof would give the marginal
% interval on each axis, which is a DIFFERENT and smaller region -- the
% classic factor-of-2.45-vs-1.96 confusion.) R7-39: with the ESTIMATED s0^2
% the exact region is 2 F(2, r): chi2(2) under-covers at small redundancy
% (a "95 %" ellipse covered 86 % at r = 5, 90.5 % at r = 10; 20 000 draws,
% R7 mirror m7_ellipse.py) and converges to it as r grows.
r = lsqvce.util.getfielddef(D, 'r', Inf);
if estimated && isfinite(r) && r >= 1
    k2 = 2 * lsqvce.kernel.stat.fInv(o.conf, 2, r);
else
    k2 = lsqvce.kernel.stat.chi2Inv(o.conf, 2);
end

a = sqrt(k2 * lam(1));      % semi-major
b = sqrt(k2 * lam(2));      % semi-minor
ang = atan2(V(2,1), V(1,1));

th = linspace(0, 2*pi, o.npts);
% npts x 2, NOT 2 x npts: this is the orientation plot() expects, so the caller
% can write plot(E.xy(:,1), E.xy(:,2)) directly.
%
% It WAS 2 x npts, which silently broke demo_workflow: E.xy(:,1) plotted a
% single POINT instead of the ellipse, and no error was raised. Caught only by
% running in MATLAB.
xy = (V * [a*cos(th); b*sin(th)]).';

cx = [D.x(i); D.x(j)];
si = sqrt(fac * [Q(1,1); Q(2,2)]);
rho = Q(1,2) / sqrt(max(Q(1,1)*Q(2,2), realmin));

E = struct('semiMajor', a, 'semiMinor', b, 'angle', ang, ...
           'ratio', a / max(b, realmin), ...
           'xy', xy + cx.', 'centre', cx, 'corr', rho, ...
           'sigInd', si, 'conf', o.conf, 'pair', [i j], ...
           'worstDir', V(:,1), 'bestDir', V(:,2), ...
           'lambda', lam, 'scale', o.scale, 'k2', k2);
end

end   % methods (Static)
end   % classdef

% ---------------------------------------------- source (formerly lsqvce.engine.verdict.source)

%==========================================================================
function k = snoopCrit(useRobust, alphaEff, rdf)
%SNOOPCRIT  The snooping threshold: Student t on the MAD's efficiency for a
%   robust scale (lsqvce.kernel.stat.madCrit, WP-128), Pope's tau otherwise.
    if useRobust
        k = lsqvce.kernel.stat.madCrit(alphaEff, rdf);
    else
        k = lsqvce.kernel.stat.tauCrit(alphaEff, rdf);
    end
end

function D = toOriginalIndex(D, obsIdx)
%TOORIGINALINDEX  Map the re-diagnosis' row indices back to the original
%   observations (WP-128). The clean re-solve numbers its rows in its own
%   order (survivors only, systems in order, local-free last, prior last),
%   while D.rejectedIdx and every caller count the original rows: after one
%   rejection "observation 399" was observation 400. The index fields are
%   mapped; the per-observation vectors (v, rNum, w, mdb, ...) stay in the
%   re-solve's order and D.obsIdx names their original index.
%   INPUT  D the re-diagnosis; obsIdx the original index of each of its rows
    D.obsIdx = obsIdx(:);
    for f = {'flagIdx', 'idxUncontrolled', 'worstIdx', 'datumRows', 'excludedIdx', ...
             'rMinIdx', 'maxAbsWIdx', 'extMaxIdx', 'priorRows'}   % priorRows: WP-135
        if ~isfield(D, f{1}) || isempty(D.(f{1})), continue, end
        v = D.(f{1});  ok = v >= 1 & v <= numel(obsIdx);   % 0 = "none"
        v(ok) = obsIdx(v(ok));  D.(f{1}) = v;
    end
end

function idx = datumRowsOf(acc, n)
%DATUMROWSOF  Global row indices of the prior's DATUM rows (WP-127), [] if none.
%   The prior is the LAST batch (rows n-nk+1..n); a datum row has r = 0 by
%   construction -- a constraint on the null space, not an observation -- so
%   diagnose keeps it in sum(r) and vPv but not in the per-observation
%   verdict (it was counted as an uncontrolled observation: every free
%   network came back SUSPECT).
%   INPUT  acc the accumulator; n the total row count (data + prior)
%   OUTPUT idx k x 1 double
    idx = zeros(0, 1);
    if isempty(acc.prior), return, end
    isD = acc.datumRows();
    idx = (n - numel(isD)) + find(isD(:));
end

function [mask, nD] = localDatumOf(batchFcn, Q, kList)
%LOCALDATUMOF  The datum rows of the local prior groups' batches (WP-135):
%   r_i = (Qvv P)_ii <= 1e-6 (the threshold is diagnose's), by diagnose's
%   formula -- redundancyNumbers. WP-138: a batch with a full weight matrix
%   was skipped without a word; it is judged by the generalised form. Only a
%   malformed one is passed over here: the main pass refuses it
%   (lsqvce:diagnose:badPmatrix), so nothing is skipped in silence.
%   INPUT  batchFcn, Q the cofactor, kList the batches (sol.localPriorBatches)
%   OUTPUT mask containers.Map batch -> nk x 1 logical (only batches with a
%          datum row); nD 1x1 their count
    mask = containers.Map('KeyType', 'double', 'ValueType', 'any');  nD = 0;
    for k = kList(:).'
        [Ak, ~, pk] = batchFcn(k);
        nk = size(Ak, 1);
        if nk == 0, continue, end
        isCorr = ~isempty(pk) && ~isvector(pk);
        if isCorr && ~isequal(size(pk), [nk nk]), continue, end   % badPmatrix, main pass
        if isempty(pk), pk = ones(nk, 1); end
        if isCorr, pk = (pk + pk') / 2; else, pk = pk(:); end
        rk = redundancyNumbers(Q, Ak, pk, isCorr);
        dm = rk <= 1e-6;                    % tol: see diagnose (WP-135)
        if ~isCorr, dm = dm & pk > 0; end   % a row of weight 0 is no observation
        if any(dm), mask(k) = dm;  nD = nD + sum(dm); end
    end
end

function [rk, W, WG, G, qv] = redundancyNumbers(Q, Ak, P, isCorr)
%REDUNDANCYNUMBERS  r_i = (Qvv P)_ii of one batch, UNCLAMPED -- the one
%   formula of diagnose's main pass, its datum pre-pass (localDatumOf) and
%   snoop's scan (WP-138: three copies of the diagonal form, and the pre-pass
%   skipped a full weight matrix).
%   Full P (nk x nk, symmetric): r_i = 1 - (A Qxx G')_ii with G = P A, and
%   A Qxx G' = W' * WG columnwise. An r_i may exceed 1 there (the projector
%   is not orthogonal in the Euclidean metric; diagnose's CLAMP note).
%   Diagonal p (nk x 1): r_i = p_i (1/p_i - a_i' Qxx a_i), 0 for a row of
%   weight 0 (WP-133). The two forms agree to 4.4e-16 when P is diagonal.
%   INPUT  Q BlockCofactor; Ak nk x u; P nk x nk | nk x 1; isCorr 1x1 logical
%   OUTPUT rk nk x 1; W u x nk = Q.halfRows(A); WG u x nk = Q.halfRows(G) and
%          G nk x u = P A (full P; [] otherwise); qv nk x 1 = 1/p - a'Qxx a
%          (diagonal p; [] otherwise), for the caller's clamp
    W = Q.halfRows(Ak);
    if isCorr
        G  = P * Ak;
        WG = Q.halfRows(G);
        rk = 1 - sum(W .* WG, 1).';
        qv = [];
    else
        G = [];  WG = [];
        qv = 1 ./ P - sum(W.^2, 1).';
        rk = P .* qv;
        rk(P == 0) = 0;
    end
end

function [A, l, p] = batchOf(src, k)
%BATCHOF  The k-th whitened batch (one arc), unit weight. Tolerates nargout 1..3.
%   Rebuilt from its address (WP-86 part 1): the arc is re-read through the
%   system's loader and re-whitened by the ONE iterator, so the rows are
%   bit-identical to accumulation. The last arc served is cached (src.cache
%   is a handle), so asking for the same batch twice reads it once.
    if k > size(src.index, 1)                      % the prior pseudo-batch, last
        A = src.prior.A;  l = src.prior.l;  p = ones(size(A, 1), 1);
        return
    end
    e = src.index(k, :);                           % {iSys, day, obsType, arcId}
    key = sprintf('%d|%s|%s|%s', e{1}, num2str(e{2}, 17), e{3}, num2str(e{4}, 17));
    c = src.cache;
    if ~(isKey(c, 'key') && strcmp(c('key'), key))
        [Ak, lk, pk] = arcRows(src, e{1}, e{2}, e{3}, e{4});
        c('key') = key;  c('A') = Ak;  c('l') = lk;  c('p') = pk;
    end
    A = c('A');  l = c('l');  p = c('p');
end

function [A, l, p] = arcRows(src, iSys, day, obsType, arcId)
%ARCROWS  One arc's whitened rows placed into the columns [common | this
%   system's locals], every kept row of every iterator batch, in order; p is
%   1 per row and 0 on the rows of weight 0 (arc.excluded, WP-133) -- the
%   iterator keeps every row once and in order, so the arc's mask is the
%   replay's (a whitener is causal: whitened row i is observation i's).
    st = src.stores{iSys};
    uo = st.uo;  uo.obsType = obsType;  uo.day = day;  uo.arc = arcId;
    arc = st.dataFcn(day, obsType, arcId, uo);
    A = zeros(0, src.uFull);  l = zeros(0, 1);  p = zeros(0, 1);
    if isempty(arc) || ~isfield(arc, 'l') || isempty(arc.l)
        return
    end
    colLocal = src.localOffset(iSys) + (1:src.localCount(iSys));
    % the iterator's batches, gathered in a handle (a local function of a
    % classdef file takes no nested function): key = batch number, in order
    got = containers.Map('KeyType', 'double', 'ValueType', 'any');
    arcBatchWhiten(arc, uo, @(~, keepRows, Ablk, Bblk, lblk) ...
        keepBatch(got, keepRows, Ablk, Bblk, lblk, colLocal, src.nc, src.uFull));
    parts = values(got);                           % ascending keys = iterator order
    As = cellfun(@(c) c{1}, parts, 'UniformOutput', false);
    ls = cellfun(@(c) c{2}, parts, 'UniformOutput', false);
    A = vertcat(A, As{:});  l = vertcat(l, ls{:});
    p = ones(size(A, 1), 1);
    ex = lsqvce.util.getfielddef(arc, 'excluded', []);
    if any(ex), p(logical(ex(:))) = 0; end
end

function keepBatch(got, keepRows, Ablk, Bblk, lblk, colLocal, nc, uFull)
%KEEPBATCH  One iterator batch, its kept rows placed full width, appended to
%   the handle GOT under the next batch number.
    Ab = zeros(numel(keepRows), uFull);
    Ab(:, 1:nc) = Ablk(keepRows, :);
    Ab(:, colLocal) = Bblk(keepRows, :);
    got(got.Count + 1) = {Ab, lblk(keepRows)};
end

function index = enumerateArcs(stores)
%ENUMERATEARCS  The address {iSys, day, obsType, arcId} of every arc, in the
%   accumulator's order (system, day, obsType, arc), from the loaders' meta
%   answers only -- no arc is read (the loader memoises its discovery).
    index = cell(0, 4);
    for iSys = 1:numel(stores)
        st = stores{iSys};
        for day = st.dayList(:).'
            meta = st.dataFcn('meta', day, [], st.uo);
            for iType = 1:numel(meta.obsTypes)
                for arcId = meta.arcs{iType}(:).'
                    index(end+1, :) = {iSys, day, meta.obsTypes{iType}, arcId};   %#ok<AGROW> addresses only
                end
            end
        end
    end
end

% ---------------------------------------------- diagnose (formerly lsqvce.engine.verdict.diagnose)

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

% ---------------------------------------------- snoop (formerly lsqvce.engine.verdict.snoop)

%% ======================================================== helpers

function guardContamination(D, o, useRobust, nObs0)
%GUARDCONTAMINATION  snoop's guard on diagnose's contamination level (WP-137:
%   lifted out of snoop, which gained the 'indeterminate' case): breakdown and
%   refuse stop snoop unless o.force (then a warning), indeterminate warns and
%   proceeds (user decision), severe warns under Pope's tau. Raises snoop's
%   identifiers unchanged.
%   INPUT  D (contamLevel, fhat, backstop), o (snoop's options: force,
%          maxRejectFrac), useRobust 1x1 logical (Test 'mad'), nObs0 the
%          observations snoop may reject from
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

        case 'indeterminate'
            % WP-137 (user decision: warn and proceed). fhat lies in its clean band
            % while the backstop is large, and at this n a robust scale captured by
            % the blunders cannot be told from a sigma0 declared too small
            % (contaminationOf: P(fhat = 0) > 1 % on clean data below n ~ 450).
            % Refusing would stop every clean run with a mis-declared sigma0; the
            % hard cap below bounds the damage of the other cause. It ran on in
            % silence (WP-135 added the level, not its case here).
            lsqvce.util.log.warn('lsqvce:snoop:indeterminate', ...
                ['snoop: the residuals are %.1fx larger than sigma0 promised, yet the ' ...
                 'robust scale sees no excess of outliers (fhat = %.3f). With %d ' ...
                 'observations the two causes cannot be told apart: sigma0 declared too ' ...
                 'small (the rejections below are sound) or more than half the data bad, ' ...
                 'the robust scale captured by them (the rejections are GOOD ' ...
                 'observations). Proceeding, at most %.0f%% of the data; check sigma0 ' ...
                 'first.'], D.backstop, D.fhat, nObs0, 100 * o.maxRejectFrac);

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
                 'result. Run R.systematic() and fix the functional model.\n' ...
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
end

function [dRowsAll, priorBatches, rMin] = snoopScope(sol, D, dRows)
%SNOOPSCOPE  What snoop leaves out (WP-135, review 3 finding 5 and lane 1):
%   the datum rows of the common prior (sol) and of the local groups (D);
%   the batches of every prior -- a prior row is a declaration, never a
%   candidate (snoop rejected an inconsistent soft local prior row, |w| 7.3);
%   and the rows diagnose calls uncontrolled (r <= D.rMin, not 1e-6).
%   WP-138: D is diagnose's output (snoop reads D.nBatch too), so D.rMin is
%   read, not defaulted -- the default was a second copy of diagnose's.
    dRowsAll = unique([dRows(:); D.datumRows(:)]);
    priorBatches = priorBatchesOf(sol);
    rMin = D.rMin;
end

function pb = priorBatchesOf(sol)
%PRIORBATCHESOF  The batches of every prior: the common prior's
%   (sol.priorBatch, 0 without one) and the local groups'
%   (sol.localPriorBatches). ONE list (WP-138) for diagnose, which leaves
%   their rows out of systematic, and snoop, which never rejects them --
%   each assembled it (WP-135).
    c = lsqvce.util.getfielddef(sol, 'priorBatch', 0);
    pb = [reshape(c(c > 0), [], 1); ...
          reshape(lsqvce.util.getfielddef(sol, 'localPriorBatches', []), [], 1)];
end

function nv = rediagnoseNV(D, o)
%REDIAGNOSENV  The Name-Value list of snoop's re-diagnosis: snoop's alpha /
%   Bonferroni / test, and the thresholds and storage of the D it was handed
%   (R7-35: rMin / rWeak / beta; WP-135: store / nWorst / robustCheck).
%   WP-138: read from D, which diagnose always fills -- the fallbacks here
%   restated diagnose's defaults, a second authority for six values.
    nv = {'alpha', o.alpha, 'bonferroni', o.bonferroni, 'test', o.test, ...
          'rMin', D.rMin, 'rWeak', D.rWeak, 'beta', D.beta, 'store', D.store, ...
          'nWorst', D.nWorst, 'robustCheck', D.robustCheck};
end

function [sol, D] = afterRejection(sol, D, rejIdx, rejected, vPvCur, dRows, batchFcn, nB, sizes, o)
%AFTERREJECTION  WP-130: the returned sol is the solution AFTER rejection on
%   either path -- vPvCur (the 1/r_i identity) was computed and dropped, so
%   with Reaccumulate=false s0sq, n, r and vPv still held the blunders. Without
%   the re-solve, D is re-diagnosed on the downdated cofactor over the
%   SURVIVING rows (masked batches, no re-accumulation), in the original
%   indices; before, D was the input D -- rejected rows still flagged, no
%   .rejectedIdx. (The re-solve overwrites both.)
    obsIdx = find(~rejected);
    sol.r = sol.r - numel(rejIdx);   % WP-133: keeps the weight-0 rows and a regularisation out
    sol.n = numel(obsIdx);  sol.vPv = vPvCur;
    sol.s0sq = vPvCur / max(sol.r, 1);
    sol.datumRows = find(ismember(obsIdx, dRows));   % in the survivors' numbering
    if o.reaccumulate, return, end
    edges = [0; cumsum(sizes(:))];
    nvR = rediagnoseNV(D, o);
    D = lsqvce.engine.verdict.diagnose(sol, @(k) maskedBatch(batchFcn, k, rejected, edges), ...
                           'nBatch', nB, nvR{:});
    D = toOriginalIndex(D, obsIdx);  D.rejectedIdx = rejIdx;   % in rejection order, as REJ.idx
end

function [bestW, bestI, bestR, bestV, bestP, vPv, nFlag, bestA, sizes, sysOf] = ...
        scanWorst(batchFcn, nB, Q, x, rejected, rdf, kCrit, useRobust, priorBatches, dRows, rMin)
%SCANWORST  One streaming pass over all batches to find the single worst
%   (largest |w|) surviving observation and accumulate vPv honestly.
%   INPUT
%     batchFcn  1x1 fh      [A,l,p] = batchFcn(k), k = 1..nB
%     nB        1x1 double  number of batches
%     Q         1x1         lsqvce.engine.BlockCofactor of the current system
%                           (the leverage a' Qxx a, WP-86 part 2)
%     x         u x 1       current parameter estimate
%     rejected  m x 1 log   mask of already-rejected observations (global index)
%     rdf       1x1 double  residual degrees of freedom (for s0)
%     kCrit     1x1 double  critical value for flagging (|w| threshold)
%     useRobust 1x1 log     use a robust (MAD) scale for w
%     priorBatches p x 1    the batches of the common and the local priors
%                           (never a candidate; WP-128, WP-135)
%     dRows     k x 1       the datum rows (global index; no robust score, WP-130)
%     rMin      1x1         diagnose's controlled threshold: a row at r <= rMin
%                           is no candidate (WP-135: r > 1e-6 made an
%                           uncontrolled row the worst one, refused, and the
%                           snooping stopped)
%   OUTPUT
%     bestW  1x1  the largest |w| found among surviving observations
%     bestI  1x1  its global observation index (0 if none)
%     bestR  1x1  its partial redundancy r_i
%     bestV  1x1  its residual v_i
%     bestP  1x1  its weight p_i
%     vPv    1x1  weighted residual sum of squares over survivors
%     nFlag  1x1  count of survivors exceeding kCrit
%     bestA  1xu  its design row (the downdate needs it; no second search pass)
%     sizes  nBx1 the row count of every batch
%     sysOf  nBx1 the system whose locals each batch touches (0: none) --
%                 the clean re-solve rebuilds the systems from it
    rawPool = [];
    sizes = zeros(nB, 1);
    sysOf = zeros(nB, 1);
    bestA = [];
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
        sizes(k) = nk;
        if nk == 0, continue; end
        sysOf(k) = Q.systemOf(Ak);
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
            % lsqvce.engine.verdict.diagnose DOES handle correlated batches (Baarda's
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
                 'WHAT YOU CAN DO: the diagnostics (R.diag) fully support correlated ' ...
                 'batches, so run it to LOCATE the outliers (D.flagIdx). Then ' ...
                 'either (a) drop the whole affected block, (b) rebuild the ' ...
                 'block covariance without the rejected observation and re-run ' ...
                 'lsqvce.fit, or (c) use lsqvce.engine.adjustRobust, which downweights ' ...
                 'rather than deletes and therefore leaves the block ' ...
                 'structure intact.'], k);
        end
        pk = pk(:);
        idx = nSeen + (1:nk);
        alive = ~rejected(idx) & pk > 0;            % WP-133: a weight-0 row is no observation
        if any(alive)
            Aa = Ak(alive,:);
            vk = Aa * x - lk(alive);
            vPv = vPv + sum(pk(alive) .* vk.^2);

            % Collect the raw standardised score for the robust scale.
            % Scale-free, so it can be pooled across batches.
            if useRobust
                Wa = Q.halfRows(Aa);
                qva = 1 ./ pk(alive) - sum(Wa.^2, 1).';
                qva(qva < 0) = 0;
                gd = qva > 0 & ~ismember(idx(alive).', dRows);   % WP-130: no datum score
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
        alive = ~rejected(idx) & pk > 0;            % WP-133: a weight-0 row is no candidate
        if ~any(alive), nSeen = nSeen + nk; continue; end

        if ismember(k, priorBatches), nSeen = nSeen + nk; continue; end   % WP-128 / WP-135: a prior
                                                                           % row is a declaration
        Aa = Ak(alive,:); la = lk(alive); pa = pk(alive);
        ia = idx(alive);

        vk = Aa * x - la;
        [rk, ~, ~, ~, qv] = redundancyNumbers(Q, Aa, pa, false);   % WP-138: diagnose's formula
        qv(qv < 0) = 0;
        rk = min(max(rk, 0), 1);

        wk = zeros(numel(la), 1);
        ok = rk > rMin;
        wk(ok) = vk(ok) ./ sqrt(s0sq * qv(ok));

        fl = ok & (abs(wk) > kCrit);
        nFlag = nFlag + sum(fl);

        if any(fl)
            [mw, mi] = max(abs(wk .* fl));
            if mw > bestW
                bestW = mw;  bestI = ia(mi);  bestR = rk(mi);
                bestV = vk(mi);  bestP = pa(mi);  bestA = Aa(mi, :);
            end
        end
        nSeen = nSeen + nk;
    end
end

%==========================================================================
function [robustScale, madAvailable, fhat, nScored] = robustScaleOf(rawAll, complete)
%ROBUSTSCALEOF  The MAD scale of the raw scores (P v)_i / sqrt((P Qvv P)_ii) and
%   the contamination gauge fhat, from either sweep of diagnose (WP-121: the
%   main pass for every test but 'mad', the separate first sweep for 'mad').
%   INPUT  rawAll   n x 1 double  the raw scores (0 where uncontrolled)
%          complete 1x1 logical   every observation was swept
%   OUTPUT robustScale (NaN if not formed), madAvailable, fhat (NaN if not
%          formed), nScored (the scores fhat counts over; WP-135: the n of its
%          null law)
robustScale = NaN;
madAvailable = false;
fhat = NaN;
nScored = 0;
if ~isempty(rawAll) && complete
    nz = rawAll(isfinite(rawAll) & rawAll ~= 0);
    nScored = numel(nz);
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
            % enough -- see contamBackstop in diagnose.
            fhat = mean(abs(nz / robustScale) > 2.5);
        end
    end
end
end

%==========================================================================
function [acc, obsIdx, order] = blockResolveAccumulator(sol, Q, batchFcn, nB, rejected, sizes, sysOf)
%BLOCKRESOLVEACCUMULATOR  The surviving rows, accumulated on the block-arrow
%   structure the cofactor Q describes (WP-86 part 3): one system per
%   original system (its surviving batches, its own local columns), the
%   batches that touch no locals as one local-free system, ready to solve.
%   INPUT  sol (method / defect / datumB), Q, batchFcn, nB, the rejected
%          mask, the batch sizes and each batch's system (scanWorst)
%   OUTPUT acc  1x1 lsqvce.engine.BlockAccumulator (KeepLocal), not solved
%          obsIdx  the ORIGINAL index of every row of acc, in acc's order
%                  (systems in order, local-free last, the prior's rows at
%                  the end) -- WP-128: the re-diagnosis numbers its rows in
%                  that order, the caller maps them back
%          order   the ORIGINAL batch of every batch of acc (WP-135)
%   WP-135: the accumulator carries the solve's Backend and Reg options
%   (sol.accOpt) -- the re-solve of a regularised fit was unregularised, and
%   a singular local block took the pinv fallback.
%   WP-128: the prior is not a batch of data here. It was re-accumulated
%   as rows of the local-free system (sol.defect is 0 from SOURCE, so the
%   setPrior branch never ran): the clean cofactor became inv(N + K'K/s^2)
%   instead of the datum's minimum-constraint one, growing with sigma.
    % The clean re-solve rebuilds the BLOCK-ARROW system the cofactor
    % describes (WP-86 part 3): each original system becomes one system
    % again -- its surviving batches, its own local columns -- and the
    % batches that touch no locals (a kind B/C source, the prior's rows)
    % one local-free system. So the clean solve has the original structure
    % and the original column layout; before, every parameter was re-solved
    % as a common one (a dense u_full factor). One batch per batch that
    % keeps a row, read through its row offset: the surviving rows are never
    % concatenated (n x u at full width -- R7-24 once more).
    edges = [0; cumsum(sizes)];
    alive = arrayfun(@(k) any(~rejected(edges(k)+1:edges(k+1))), (1:nB).');
    pb = lsqvce.util.getfielddef(sol, 'priorBatch', 0);
    if pb > 0, alive(pb) = false; end                   % the prior goes in by setPrior
    obsIdx = zeros(0, 1);  order = zeros(0, 1);
    backend = 'qr';
    if isfield(sol, 'method') && strcmpi(char(sol.method), 'normal'), backend = 'normal'; end
    ao = lsqvce.util.getfielddef(sol, 'accOpt', struct());
    regNV = {};
    if isfield(ao, 'Backend'), backend = ao.Backend; end
    if isfield(ao, 'Reg')
        regNV = {'Reg', ao.Reg, 'RegLambda', ao.RegLambda, 'RegMatrix', ao.RegMatrix, 'TruncTol', ao.TruncTol};
    end
    acc = lsqvce.engine.BlockAccumulator(Q.nc, 'Backend', backend, 'KeepLocal', true, ...
                                         'Verbose', false, regNV{:});
    if isfield(sol, 'defect') && sol.defect > 0 && isfield(sol, 'datumB') && ~isempty(sol.datumB)
        acc.setPrior(lsqvce.engine.CommonPrior(sol.datumB.', zeros(size(sol.datumB, 2), 1), 1));
    end
    nS = numel(Q.sys);
    for s = [1:nS, 0]                                   % the systems in order, local-free last
        kept = find(alive & sysOf(:) == s);
        if isempty(kept), continue, end
        cols = zeros(1, 0);  sid = nS + 1;
        if s > 0, cols = Q.sys(s).cols;  sid = s; end
        L = lsqvce.engine.batchSystem(Q.nc, @(j) maskedBatch(batchFcn, kept(j), rejected, edges), ...
                numel(kept), 'LocalColumns', cols, 'SystemId', sid, 'Name', 'snoop');
        acc.addSystem(L.sysID, numel(cols), L.closure, 1:numel(kept));   % one day per kept batch
        for k = kept(:).'
            obsIdx = [obsIdx; edges(k) + find(~rejected(edges(k)+1:edges(k+1)))]; %#ok<AGROW>
        end
        order = [order; kept(:)];                      %#ok<AGROW> one batch per kept batch
    end
    if pb > 0
        acc.setPrior(sol.prior);
        obsIdx = [obsIdx; (edges(pb)+1 : edges(pb+1)).'];
    end
end

function [A, l, p] = maskedBatch(batchFcn, k, rejected, edges)
%MASKEDBATCH  Batch k with the rejected observations removed; edges(k) is the
%   row offset of batch k, so only batch k is read.
    [A, l, p] = batchFcn(k);
    l = l(:);
    if isempty(p), p = ones(numel(l), 1); end
    alive = ~rejected(edges(k) + (1:numel(l)));
    A = A(alive, :);  l = l(alive);  p = p(alive);
    p = p(:);
end

% ---------------------------------------------- systematic (formerly lsqvce.engine.verdict.systematic)

function s = mark(p, a)
%MARK  Return a '<<< REJECT' tag when a p-value falls below the level a, else ''.
%   INPUT   p 1x1 double  a test p-value;  a 1x1 double  the significance level
%   OUTPUT  s 1xk char     the marker string ('' if not significant)
    if p < a, s = '   <<< REJECT'; else, s = ''; end
end
