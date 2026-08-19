function arcBatchWhiten(a, uo, cb)
%ARCBATCHWHITEN  Iterate one arc in whitened overlap batches (the ONE loop).
%   arcBatchWhiten(a, uo, cb)
%
%   The ONE batch-and-whiten iterator (WP-38, fixes F-08; rule 3) behind
%   LSQVCE.REDUNDANCY's arcLeverageAndResidual and
%   LSQVCE.STREAMCONTRIBUTION's arcCommonProbe -- previously two verbatim
%   copies of the batch/overlap/whitening skeleton that had to stay in
%   lockstep by hand ("whitening applied EXACTLY as the exact leverage
%   path"). The iterator owns ONLY the skeleton: batch bounds, overlap,
%   and the Wfcn/Cfcn whitening of the FULL overlap block. The callback
%   owns everything numeric downstream -- keep-slicing, payload,
%   accumulation -- so each caller''s summation order is untouched and the
%   fixture numbers stay byte-identical.
%
%   INPUT
%     a    1x1 struct   the loaded arc: .A (m x nc), .B (m x nl), .l (m x 1),
%                       optional .Wfcn(idx,Ab,Bb,lb,uo) | .Cfcn(idx,uo),
%                       optional .filterMemory (default overlap)
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
%   ERRORS
%     (none of its own; Wfcn/Cfcn/decorrFactor errors propagate)
%
%   EXAMPLE
%     h = zeros(size(a.A,1),1);
%     arcBatchWhiten(a, uo, @(rows,keep,Ab,Bb,lb) accum(rows,keep,Ab,Bb));
%
%   See also LSQVCE.REDUNDANCY, LSQVCE.STREAMCONTRIBUTION,
%            LSQVCE.DECORRFACTOR
%
%   Created 2026-08-19 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.

    m   = size(a.A, 1);
    bs  = lsqvce.util.getfielddef(uo, 'BatchSize', 20000);
    nOv = lsqvce.util.getfielddef(uo, 'BatchOverlap', []);
    if isempty(nOv)
        nOv = lsqvce.util.getfielddef(a, 'filterMemory', 0);
    end
    bs  = max(1, min(round(bs),  m));
    nOv = max(0, min(round(nOv), bs-1));

    hasW = isfield(a, 'Wfcn') && ~isempty(a.Wfcn);
    hasC = isfield(a, 'Cfcn') && ~isempty(a.Cfcn);

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
            W  = lsqvce.decorrFactor(a.Cfcn(idx, uo));
            Ab = W*Ab;  Bb = W*Bb;  lb = W*lb;
        end
        cb(rows, keep, Ab, Bb, lb);
    end
end
