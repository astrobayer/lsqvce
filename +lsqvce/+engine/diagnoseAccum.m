function [diagSol, batchFcn, nBatch] = diagnoseAccum(acc, sol)
%LSQVCE.ENGINE.DIAGNOSEACCUM  Bridge a block-arrow lsqvce.engine.BlockAccumulator adjustment into LSQVCE diagnostics.
%
%   [diagSol, batchFcn, nBatch] = LSQVCE.ENGINE.DIAGNOSEACCUM(acc, sol) returns exactly what
%   lsqvce.engine.diagnose / lsqvce.engine.snoop expect:
%       diagSol   struct with fields {x, Q, n, u, r, rank, defect, sigma0, s0sq,
%                 vPv, vPvSource, condR, condN, condIsBound, nBatch}; Q is the
%                 LSQVCE.ENGINE.BLOCKCOFACTOR of the full parameter set
%       batchFcn  @(k) -> [A, l, p]   the k-th WHITENED batch (one arc), unit weight
%       nBatch    number of batches (the arcs, then the prior's pseudo-batch)
%
%   so you can run the full five-level lsqvce verdict on a gravity-recovery
%   adjustment built with the streaming accumulator:
%
%       R   = lsqvce.fit(data, Estimate="noise", Pool="obsType");
%       acc = R.acc;   sol = R.sol;               % fit did the next three lines
%       [dsol, src, nB] = lsqvce.engine.diagnoseAccum(acc, sol);   % already: R.source / R.batchFcn / R.nBatch
%       D = lsqvce.engine.diagnose(dsol, src, 'nBatch', nB);
%       V = lsqvce.report(D);              % the five levels
%       S = lsqvce.engine.systematic(D);          % Ljung-Box / runs / covariate
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
%   and lsqvce.engine.diagnose costs O(nc^2 + nl*nc) per row with it. Before,
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
%   ERRORS
%     lsqvce:diagnoseAccum:noSystems     the accumulator carries no systems
%     lsqvce:diagnoseAccum:needKeepLocal the local blocks were not retained, so
%                                   the common redundancies cannot be corrected;
%                                   re-run keeping them
%
%   See also LSQVCE.ENGINE.DIAGNOSE, LSQVCE.ENGINE.SNOOP, LSQVCE.REPORT, LSQVCE.ENGINE.REDUNDANCY,
%   LSQVCE.ENGINE.BLOCKCOFACTOR
%
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
    blocks      = struct('cols', cell(1, nSys), 'Nll', [], 'Nlc', [], 'truncTol', []);
    n           = 0;
    for iSys = 1:nSys
        st = acc.systemStore(sysIds(iSys));
        stores{iSys}     = st;
        localOffset(iSys) = col;
        localCount(iSys)  = st.nl;
        kept = acc.getStore(sysIds(iSys));             % the Schur blocks the solve used
        blocks(iSys) = struct('cols', col + (1:st.nl), 'Nll', kept.Nll_r, 'Nlc', kept.Nlc, ...
                              'truncTol', acc.opt.TruncTol);   % the solve's pinv cut
        n   = n + kept.mObs;                           % rows the accumulator folded
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
    Q = lsqvce.engine.BlockCofactor(acc.commonFactor(), blocks);

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
        'r',      n - uFull, ...
        'rank',   uFull, ...
        'defect', 0, ...
        'sigma0', 1, ...
        's0sq',   sol.s0sq, ...
        'vPv',    sol.vtv, ...
        'vPvSource', 'accumulator (v''v)', ...
        'condR',  sqrt(kN), ...
        'condN',  kN, ...
        'condIsBound', true, ...
        'nBatch', nBatch);

    % ---- batchFcn: replay the k-th whitened batch --------------------------
    % A local function (not deal) so a caller may take 1, 2, or 3 outputs:
    % lsqvce.engine.diagnose asks for [A,l,p], but a redundancy-only sweep wants just
    % A. deal() would throw on nargout < 3; batchOf never does.
    batchFcn = @(k) batchOf(src, k);
end

%==========================================================================
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
        [Ak, lk] = arcRows(src, e{1}, e{2}, e{3}, e{4});
        c('key') = key;  c('A') = Ak;  c('l') = lk;
    end
    A = c('A');  l = c('l');  p = ones(size(A, 1), 1);
end

function [A, l] = arcRows(src, iSys, day, obsType, arcId)
%ARCROWS  One arc's whitened rows placed into the columns [common | this
%   system's locals], every kept row of every iterator batch, in order.
    st = src.stores{iSys};
    uo = st.uo;  uo.obsType = obsType;  uo.day = day;  uo.arc = arcId;
    arc = st.dataFcn(day, obsType, arcId, uo);
    A = zeros(0, src.uFull);  l = zeros(0, 1);
    if isempty(arc) || ~isfield(arc, 'l') || isempty(arc.l)
        return
    end
    colLocal = src.localOffset(iSys) + (1:src.localCount(iSys));
    As = {};  ls = {};
    arcBatchWhiten(arc, uo, @collect);
    A = vertcat(A, As{:});  l = vertcat(l, ls{:});
    function collect(~, keepRows, Ablk, Bblk, lblk)
        Ab = zeros(numel(keepRows), src.uFull);
        Ab(:, 1:src.nc) = Ablk(keepRows, :);
        Ab(:, colLocal) = Bblk(keepRows, :);
        As{end+1} = Ab;  ls{end+1} = lblk(keepRows);
    end
end

%==========================================================================
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
