function arcBatchWhiten(a, uo, cb)
%ARCBATCHWHITEN  Iterate one arc in whitened overlap batches (the ONE loop).
%   arcBatchWhiten(a, uo, cb)
%
%   The ONE batch-and-whiten iterator (WP-38, fixes F-08; rule 3). Since
%   R6 WP-91 (finding E-09) it is behind EVERY loop that reads an arc in
%   whitened batches: LSQVCE.ENGINE.BLOCKACCUMULATOR's accumArc,
%   LSQVCE.ENGINE.RESIDUALS' arcResiduals, LSQVCE.ENGINE.VERDICT.SOURCE's
%   whitenArcIntoBatches, LSQVCE.ENGINE.REDUNDANCY's arcLeverageAndResidual
%   and LSQVCE.ENGINE.STREAMCONTRIBUTION's arcCommonProbe -- before, three of
%   them were verbatim copies WITHOUT the accumulator's whole-epoch rounding,
%   so a multi-axis arc longer than one batch failed in the sweeps. The
%   iterator owns ONLY the skeleton: batch bounds, overlap, the rounding of
%   both to whole epochs for a multi-axis arc, and the Wfcn / Cfcn / dense-W
%   whitening of the FULL overlap block. The callback owns everything numeric
%   downstream -- keep-slicing, payload, accumulation -- so each caller's
%   summation order is untouched and the fixture numbers stay byte-identical.
%
%   INPUT
%     a    1x1 struct   the loaded arc: .A (m x nc), .B (m x nl), .l (m x 1),
%                       optional .Wfcn(idx,Ab,Bb,lb,uo) | .Cfcn(idx,uo) | .W
%                       (a dense m x m whitening matrix), optional
%                       .filterMemory (default overlap), optional .nAxis
%                       (rows epoch-major; batches hold whole epochs)
%                       optional .wholeArc (true: ONE batch, no overlap --
%                       a nugget arc, whose whitener is not causal, or a
%                       dense law composed with Cep). An arc whitened by a
%                       dense covariance (Cfcn or W, no Wfcn) is always ONE
%                       batch: a dense factor cannot be joined by an overlap
%     uo   1x1 struct   user options; honours .BatchSize (default 20000)
%                       and .BatchOverlap (default a.filterMemory, else 0)
%     cb   1x1 function_handle   cb(rows, keep, Ab, Bb, lb) called once per
%                       batch: rows = i0:i1 (the batch''s own rows),
%                       keep = the indices of those rows WITHIN the
%                       whitened overlap block, Ab/Bb/lb = the whitened
%                       FULL overlap block (idx = j0:i1). NOTE the Cfcn
%                       branch whitens lb as well (the leverage path''s
%                       behaviour); a callback that does not need lb
%                       simply ignores it -- outputs are unchanged.
%
%   OUTPUT
%     (none -- the callback accumulates)
%
%   ACCURACY
%     A dense law is factorised over the whole arc (exact). For a filter law,
%     with the overlap >= the filter memory a batched pass equals the
%     whole-arc pass up to summation order (test_lsqAccum: BatchSize
%     invariance at 1e-6 across backends, and 1e-12 on x / the whitened
%     residuals / the diagnostic rows of a 3-axis arc at BatchSize 100 vs one
%     batch, R6 WP-91).
%
%   ERRORS
%     lsqvce:accum:notWholeEpochs   a multi-axis arc whose row count is not a
%                                   multiple of its nAxis (rows must be
%                                   epoch-major); Wfcn / Cfcn / decorrFactor
%                                   errors propagate
%
%   EXAMPLE
%     h = zeros(size(a.A,1),1);
%     arcBatchWhiten(a, uo, @(rows,keep,Ab,Bb,lb) accum(rows,keep,Ab,Bb));
%
%   See also LSQVCE.ENGINE.BLOCKACCUMULATOR, LSQVCE.ENGINE.RESIDUALS,
%            LSQVCE.ENGINE.VERDICT.SOURCE, LSQVCE.ENGINE.REDUNDANCY,
%            LSQVCE.ENGINE.STREAMCONTRIBUTION, LSQVCE.KERNEL.DECORRFACTOR
%
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-125: a dense law, or an overlap the
%   batch cannot hold, is whitened as ONE batch. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-115 part 3b: an arc that
%   carries wholeArc (a nugget arc) is whitened in ONE batch. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-86 part 4: the diagnostics
%   replay is lsqvce.engine.verdict.source. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-09-07 (CEST, UTC+02:00), R6 WP-91 (E-09): the whole-epoch
%   rounding and the dense-W branch moved in from the accumulator; the
%   residual sweep and the diagnostic re-read became callbacks. With
%   assistance from Claude (Anthropic), model Claude Fable 5.1.
%   Created 2026-08-19 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.

    m   = size(a.A, 1);
    bs  = lsqvce.util.getfielddef(uo, 'BatchSize', 20000);
    nOv = lsqvce.util.getfielddef(uo, 'BatchOverlap', []);
    if isempty(nOv)
        nOv = lsqvce.util.getfielddef(a, 'filterMemory', 0);
    end
    % Multi-axis: rows are epoch-major, so a batch must never split an
    % epoch. Round BatchSize and the overlap DOWN to whole epochs (R6 WP-91:
    % this rounding lived in the accumulator only -- finding E-09).
    nAx = lsqvce.util.getfielddef(a, 'nAxis', 1);
    if nAx > 1
        if mod(m, nAx) ~= 0
            error('lsqvce:accum:notWholeEpochs', ...
                ['The arc has %d rows, not a multiple of nAxis = %d. ' ...
                 'Rows must be epoch-major.'], m, nAx);
        end
        bs  = max(nAx, floor(bs /nAx)*nAx);
        nOv =           floor(nOv/nAx)*nAx;
    end
    hasW = isfield(a, 'Wfcn') && ~isempty(a.Wfcn);
    hasC = isfield(a, 'Cfcn') && ~isempty(a.Cfcn);
    hasM = isfield(a, 'W')    && ~isempty(a.W);
    if lsqvce.util.getfielddef(a, 'wholeArc', false) || (~hasW && (hasC || hasM))
        bs = m;  nOv = 0;                            % a nugget arc, or a dense law (WP-125)
    end
    bs  = max(1, min(round(bs),  m));
    if round(nOv) >= bs                              % WP-125: an overlap the batch cannot
        bs = m;  nOv = 0;                            % hold would be clipped inexactly --
    end                                              % the arc is whitened whole instead
    nOv = max(0, min(round(nOv), bs-1));
    if nAx > 1, nOv = floor(nOv/nAx)*nAx; end        % the clip to bs-1 may split an epoch

    for i0 = 1:bs:m
        i1   = min(i0+bs-1, m);
        j0   = max(1, i0-nOv);
        idx  = j0:i1;
        keep = (i0-j0+1):numel(idx);
        rows = i0:i1;

        Ab = a.A(idx, :);
        Bb = a.B(idx, :);
        lb = a.l(idx);  lb = lb(:);

        if hasW
            [Ab, Bb, lb] = a.Wfcn(idx, Ab, Bb, lb, uo);
        elseif hasC
            W  = lsqvce.kernel.decorrFactor(a.Cfcn(idx, uo));
            Ab = W*Ab;  Bb = W*Bb;  lb = W*lb;
        elseif hasM
            W  = a.W(idx, idx);
            Ab = W*Ab;  Bb = W*Bb;  lb = W*lb;
        end
        % else: unit weights
        cb(rows, keep, Ab, Bb, lb);
    end
end
