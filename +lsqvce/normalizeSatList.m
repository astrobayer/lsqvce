function satListCell = normalizeSatList(satList)
%LSQVCE.NORMALIZESATLIST  Normalise a satellite list to a cell of scalar ids.
%   satListCell = lsqvce.normalizeSatList(satList)
%
%   Accepts the several natural ways of writing a satellite list and returns a
%   uniform 1xN cell of scalar identifiers, so downstream code can iterate with
%   satListCell{i} regardless of how the caller spelled it.
%
%   INPUT
%     satList   1xN | Nx1   numeric | char | string | cell   The list. Required.
%               Mappings:
%                 [1 2]        (numeric)      -> {1, 2}
%                 'CD'         (char row)     -> {'C', 'D'}
%                 'C'          (1x1 char)     -> {'C'}
%                 ["C" "D"]    (string array) -> {'C', 'D'}
%                 {'C','D'}    (cell)         -> {'C', 'D'}
%
%   OUTPUT
%     satListCell   1xN   cell   Row cell of scalar ids (char or numeric).
%
%   EXAMPLE
%     lsqvce.normalizeSatList([1 2])       % -> {[1] [2]}
%     lsqvce.normalizeSatList('CD')        % -> {'C' 'D'}
%
%   ERRORS
%     lsqvce:normalizeSatList:empty  the satellite list is empty
%     lsqvce:normalizeSatList:type   an entry is neither numeric nor char
%
%   See also LSQVCE.RESOLVESATSPECS, LSQVCE.ITERATESEPARABLE
%
%   Created 2026-07-15 13:40 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.
    narginchk(1, 1);
    if isempty(satList)
        error('lsqvce:normalizeSatList:empty', 'satList must be non-empty.');
    end
    if iscell(satList)
        satListCell = satList(:).';
    elseif ischar(satList)
        if isrow(satList)
            satListCell = num2cell(satList);              % 'CD' -> {'C','D'}
        else
            satListCell = num2cell(satList, 2).';         % char column
        end
    elseif isstring(satList)
        satListCell = cellstr(satList(:).');
    elseif isnumeric(satList)
        satListCell = num2cell(satList(:).');             % [1 2] -> {1,2}
    else
        error('lsqvce:normalizeSatList:type', ...
              'satList must be numeric, char, string or cell.');
    end
end
