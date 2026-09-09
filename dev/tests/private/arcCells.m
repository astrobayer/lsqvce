function D = arcCells(D, schedule)
%ARCCELLS  Cut a v5-shaped day (matrices + .arcIdx) into the v6 cell contract.
%   D = arcCells(D)             cells in .arcIdx order; ids are .arcId or 1:n
%   D = arcCells(D, schedule)   cells in SCHEDULE order (a 1xS id vector): an
%                               id the day does not carry becomes an EMPTY
%                               cell -- the v6 data gap
%
%   The ONE instrument of the test tree that speaks both day contracts (R6
%   WP-102): the fixtures keep generating a day as whole-day matrices with
%   row-index arcs (the cheapest way to plant a signal and its noise across
%   arcs), and this cut turns that into what the v6 loader accepts -- one
%   cell per arc for A, B, l, t and Cep. Every builder of the test tree calls
%   it last, so the mapping is written once and the parity test
%   (test_fit>testCellContractParity) pins it against the v5 numbers.
%
%   INPUT
%     D         1x1 struct   .A (mxnc) .B (mxnl) .l (mx1) .t (mx1),
%                            .arcIdx {1xn} row-index vectors, optional .arcId
%                            (1xn ids), .Cep (nEpoch x 1 | x d | x d x d, per
%                            DAY epoch), .nAxis
%     schedule  1xS double   optional; the ids the cells must follow
%   OUTPUT
%     D         1x1 struct   .A .B .l .t {1xS} cells (S = n without a
%                            schedule), .Cep {1xS} when the day had one,
%                            .nAxis kept; .arcIdx / .arcId removed
%
%   ERRORS  none of its own (a malformed fixture fails in the loader's
%           validateDay, which is what the contract tests exercise).
%
%   See also test_fit, test_adjustBlock, lsqvce.engine.makeLoader
%
%   Created 2026-09-09 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R6 WP-102.
    if ~isfield(D, 'arcIdx') || isempty(D.arcIdx)
        D.arcIdx = {1:size(D.A, 1)};
    end
    n = numel(D.arcIdx);
    if isfield(D, 'arcId') && ~isempty(D.arcId)
        ids = D.arcId(:).';
    else
        ids = 1:n;
    end
    if nargin < 2 || isempty(schedule), schedule = ids; end
    hasCep = isfield(D, 'Cep') && ~isempty(D.Cep);
    if hasCep
        d = size(D.Cep, 2);
    elseif isfield(D, 'nAxis') && ~isempty(D.nAxis)
        d = D.nAxis;
    else
        d = 1;
    end
    S = numel(schedule);
    A = cell(1, S);  B = cell(1, S);  l = cell(1, S);  t = cell(1, S);  Cep = cell(1, S);
    nc = size(D.A, 2);  nl = size(D.B, 2);
    for k = 1:S
        pos = find(ids == schedule(k), 1);
        if isempty(pos)                              % the gap: 0 rows, right widths
            A{k} = zeros(0, nc);  B{k} = zeros(0, nl);  l{k} = zeros(0, 1);  t{k} = zeros(0, 1);
            continue
        end
        r = D.arcIdx{pos};
        A{k} = D.A(r, :);  B{k} = D.B(r, :);  l{k} = D.l(r);  t{k} = D.t(r);
        if hasCep
            ep = unique(ceil(r(:) / d), 'stable');   % the arc's day-epochs (epoch-major rows)
            if ismatrix(D.Cep), Cep{k} = D.Cep(ep, :); else, Cep{k} = D.Cep(ep, :, :); end
        end
    end
    out = struct('A', {A}, 'B', {B}, 'l', {l}, 't', {t});
    if hasCep, out.Cep = Cep; end
    if isfield(D, 'nAxis') && ~isempty(D.nAxis), out.nAxis = D.nAxis; end
    D = out;
end
