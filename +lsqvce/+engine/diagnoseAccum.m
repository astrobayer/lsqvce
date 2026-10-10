function [diagSol, batchFcn, nBatch] = diagnoseAccum(acc, sol)
%LSQVCE.ENGINE.DIAGNOSEACCUM  Bridge a block-arrow lsqvce.engine.BlockAccumulator adjustment into LSQVCE diagnostics.
%
%   [diagSol, batchFcn, nBatch] = LSQVCE.ENGINE.DIAGNOSEACCUM(acc, sol) returns exactly what
%   lsqvce.engine.diagnose / lsqvce.engine.snoop expect:
%       diagSol   struct with fields {x, R, n, u, r, sigma0, s0sq, vPv, vPvSource, nBatch}
%       batchFcn  @(k) -> [A, l, p]   the k-th WHITENED batch, unit weight
%       nBatch    number of batches
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
%   inflate every r_i. So this bridge rebuilds the full block-arrow system:
%   columns are [x_common | y_sys1 | y_sys2 | ...], and each batch places its
%   whitened A into the common columns and its whitened B into that system's
%   local columns. The redundancy numbers this yields match lsqAccum's exact
%   block-arrow leverage to 1e-16 (the two formulations are identical).
%
%   The whitening replays the loader's own Wfcn/Cfcn, batch by batch, with the
%   same overlap discard as accumulation -- so the batches are bit-identical to
%   what the accumulator saw. This is the 'verify' contract lsqvce checks.
%
%   COST
%   ----
%   diagSol carries the full R (u_full x u_full, u_full = nc + sum(nl)), and
%   lsqvce.engine.diagnose costs O(u_full^2) per row with it. For GRACE-FO with
%   two satellites and a handful of local parameters each this is trivial; with
%   per-arc locals in the thousands both are large (WP-86 part 2: the
%   block-native verdict). The batches themselves are NOT stored (WP-86 part 1,
%   R7-24): batchFcn re-reads an arc through its loader and re-whitens it on
%   demand, caching one arc, so the memory is one arc's full-width batches.
%   Before, every whitened row was kept at full width -- n x u_full doubles in
%   the Results, ~228 GB for one GRAVFIRE-scale month. The price is a re-read
%   of the data per diagnostic pass, and it RELIES ON THE LOADER'S CONTRACT that
%   an arc is returned identically on every call (a changed count is still
%   caught: lsqvce:diagnose:batchMismatch).
%
%   REQUIRES  acc constructed with 'KeepLocal', true (the default).
%
%   INPUT
%     acc       1x1 lsqvce.engine.BlockAccumulator  the accumulator that solved
%     sol       1x1 struct  the solution it returned
%     varargin  Name/Value forwarded to LSQVCE.ENGINE.DIAGNOSE
%
%   OUTPUT
%     D  1x1 struct  the same diagnostics LSQVCE.ENGINE.DIAGNOSE produces, but with the
%                    common-parameter redundancies CORRECTED for the eliminated
%                    arc-local parameters. Reading the reduced common system
%                    alone inflates every reliability number.
%
%   ERRORS
%     lsqvce:diagnoseAccum:noSystems     the accumulator carries no systems
%     lsqvce:diagnoseAccum:needKeepLocal the local blocks were not retained, so
%                                   the common redundancies cannot be corrected;
%                                   re-run keeping them
%
%   See also LSQVCE.ENGINE.DIAGNOSE, LSQVCE.ENGINE.SNOOP, LSQVCE.REPORT, LSQVCE.ENGINE.REDUNDANCY
%
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

    % The bridge rebuilds the FULL block-arrow system from the per-system
    % Schur blocks. Those are only retained when the accumulator was built with
    % 'KeepLocal', true. Without them we could neither back-substitute the
    % locals nor re-enumerate the arcs, so refuse loudly rather than return a
    % silently degenerate (zero-batch) diagSol.
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
    for iSys = 1:nSys
        st = acc.systemStore(sysIds(iSys));
        stores{iSys}     = st;
        localOffset(iSys) = col;
        localCount(iSys)  = st.nl;
        col               = col + st.nl;
    end
    uFull = col;                                    % nc + sum(nl)

    % ---- full parameter vector: x_common then each system's back-substituted y
    xFull = zeros(uFull, 1);
    xFull(1:nc) = sol.x;
    for iSys = 1:nSys
        y = acc.backSubstitute(sysIds(iSys), sol.x);
        xFull(localOffset(iSys) + (1:localCount(iSys))) = y;
    end

    % ---- enumerate the WHITENED batches once: fold R, keep an INDEX ---------
    % One batch = one (system, arc, row-batch). The full-width root factor R
    % is folded batch by batch DURING the enumeration, and only the batch's
    % ADDRESS is kept (system, day, obsType, arc, sub-batch). batchFcn rebuilds
    % a batch from that address on demand (WP-86 part 1, R7-24): before, every
    % full-width batch was stored -- n x uFull doubles inside the Results,
    % ~228 GB for one GRAVFIRE-scale month.
    [Rfull, n, index] = enumerateBatches(stores, sysIds, localOffset, localCount, nc, uFull);
    src = struct('stores', {stores}, 'localOffset', localOffset, 'localCount', localCount, ...
                 'nc', nc, 'uFull', uFull, 'index', index, 'prior', [], ...
                 'cache', containers.Map('KeyType', 'char', 'ValueType', 'any'));
    % R5 M5 5b: the common-block PRIOR is one more batch of pseudo-observations
    % (rows K/sigma on the common columns, observations x0/sigma), so the hat
    % matrix, the redundancies and the residual sum include it exactly as the
    % solve did (Koch: the prior's redundancy is what it absorbs). Without it
    % the diagnostics rebuilt a rank-deficient N from the data alone.
    nBatch = size(index, 1);
    if ~isempty(acc.prior)
        K  = acc.prior.K / acc.prior.sigma;
        Ap = [K, zeros(size(K, 1), uFull - nc)];
        lp = acc.prior.x0(:) / acc.prior.sigma;
        src.prior = struct('A', Ap, 'l', lp);       % nk x uFull: the one stored batch
        nBatch = nBatch + 1;
        Rfull = foldR(Rfull, Ap, uFull);
        n     = n + size(Ap, 1);
    end

    % ---- diagSol in lsqvce contract -------------------------------------
    kR = cond(full(Rfull));            % reduced-system conditioning; N = R'R
    % R5 M5 5b': sigma0 below is the A-PRIORI sigma0 of the law-whitened rows,
    % 1 by construction, so diagnose's varRatio and its a-priori backstop are
    % s0sq/1. The accumulator's sol.sigma0 is the a-posteriori one -- passing
    % it here made both identically 1 (the breakdown backstop could never fire).
    diagSol = struct( ...
        'x',      xFull, ...
        'R',      Rfull, ...
        'n',      n, ...
        'u',      uFull, ...
        'r',      n - uFull, ...
        'sigma0', 1, ...
        's0sq',   sol.s0sq, ...
        'vPv',    sol.vtv, ...
        'vPvSource', 'accumulator (v''v)', ...
        'condR',  kR, ...
        'condN',  kR^2, ...
        'nBatch', nBatch);

    % ---- batchFcn: replay the k-th whitened batch --------------------------
    % A local function (not deal) so a caller may take 1, 2, or 3 outputs:
    % lsqvce.engine.diagnose asks for [A,l,p], but a redundancy-only sweep wants just
    % A. deal() would throw on nargout < 3; batchOf never does.
    batchFcn = @(k) batchOf(src, k);
end

%==========================================================================
function [A, l, p] = batchOf(src, k)
%BATCHOF  The k-th whitened batch, unit weight. Tolerates nargout 1..3.
%   Rebuilt from its address (WP-86 part 1): the arc is re-read through the
%   system's loader and re-whitened by the ONE iterator, so the rows are
%   bit-identical to the enumeration pass. One arc's batches are cached
%   (src.cache is a handle), so a sequential pass -- diagnose, the residual
%   replay -- reads every arc once; random access (snoop's getRow) may re-read.
    if k > size(src.index, 1)                      % the prior pseudo-batch, last
        A = src.prior.A;  l = src.prior.l;  p = ones(size(A, 1), 1);
        return
    end
    e = src.index(k, :);                           % {iSys, day, obsType, arcId, sub}
    key = sprintf('%d|%s|%s|%s', e{1}, num2str(e{2}, 17), e{3}, num2str(e{4}, 17));
    c = src.cache;
    if ~(isKey(c, 'key') && strcmp(c('key'), key))
        c('key') = key;
        c('batches') = arcBatches(src, e{1}, e{2}, e{3}, e{4});
    end
    b = c('batches');
    A = b{e{5}}.A;  l = b{e{5}}.l;  p = ones(size(A, 1), 1);
end

function batches = arcBatches(src, iSys, day, obsType, arcId)
%ARCBATCHES  One arc's whitened full-width batches, in iterator order.
    st = src.stores{iSys};
    uo = st.uo;  uo.obsType = obsType;  uo.day = day;  uo.arc = arcId;
    arc = st.dataFcn(day, obsType, arcId, uo);
    colLocal = src.localOffset(iSys) + (1:src.localCount(iSys));
    batches = {};
    arcBatchWhiten(arc, uo, @rebuild);
    function rebuild(~, keepRows, Ablk, Bblk, lblk)
        batches{end+1} = fullWidth(keepRows, Ablk, Bblk, lblk, colLocal, src.nc, src.uFull);
    end
end

function b = fullWidth(keepRows, Ablk, Bblk, lblk, colLocal, nc, uFull)
%FULLWIDTH  One whitened batch placed into the columns [common | this system's locals].
    fullA = zeros(numel(keepRows), uFull);
    fullA(:, 1:nc)     = Ablk(keepRows, :);
    fullA(:, colLocal) = Bblk(keepRows, :);
    b = struct('A', fullA, 'l', lblk(keepRows));
end

function R = foldR(R, A, uFull)
%FOLDR  Fold one batch into the full-width root factor (streamed QR).
    R = triu(qr([R; A], 0));
    R = R(1:min(size(R, 1), uFull), :);
end

%==========================================================================
function [Rfull, n, index] = enumerateBatches(stores, sysIds, localOffset, localCount, nc, uFull)
%ENUMERATEBATCHES  Whiten every arc of every system, batch by batch, place the
%   rows into the full-width columns, FOLD each into the root factor, and keep
%   only the batch's address (WP-86 part 1). Mirrors the accumulator's arc
%   loop and overlap-discard, so the batches are identical to accumulation.
%   OUTPUT  Rfull uFull x uFull (or fewer rows); n rows folded; index K x 5
%           cell {iSys, day, obsType, arcId, sub-batch}.
    Rfull = zeros(0, uFull);  n = 0;
    index = cell(0, 5);
    for iSys = 1:numel(stores)
        st = stores{iSys};
        uo = st.uo;
        colLocal = localOffset(iSys) + (1:localCount(iSys));

        for day = st.dayList(:).'
            meta = st.dataFcn('meta', day, [], uo);
            for iType = 1:numel(meta.obsTypes)
                uo.obsType = meta.obsTypes{iType};
                for arcId = meta.arcs{iType}(:).'
                    uo.day = day;  uo.arc = arcId;
                    arc = st.dataFcn(day, uo.obsType, arcId, uo);
                    if isempty(arc) || ~isfield(arc,'l') || isempty(arc.l)
                        continue
                    end
                    sub = 0;
                    arcBatchWhiten(arc, uo, @foldAndIndex);
                end
            end
        end
    end

    function foldAndIndex(~, keepRows, Ablk, Bblk, lblk)
        b = fullWidth(keepRows, Ablk, Bblk, lblk, colLocal, nc, uFull);
        Rfull = foldR(Rfull, b.A, uFull);
        n = n + size(b.A, 1);
        sub = sub + 1;
        index(end+1, :) = {iSys, day, uo.obsType, arcId, sub};   %#ok<AGROW> addresses only
    end
end

