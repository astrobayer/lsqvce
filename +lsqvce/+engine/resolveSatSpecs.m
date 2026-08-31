function [satSpecs, specFor] = resolveSatSpecs(satListCell, userSpec, defaults)
%LSQVCE.ENGINE.RESOLVESATSPECS  Per-satellite configuration: overrides over defaults.
%   [satSpecs, specFor] = lsqvce.engine.resolveSatSpecs(satListCell, userSpec, defaults)
%
%   Produces one resolved configuration struct per satellite by layering
%   per-satellite overrides (userSpec, matched by .satId) over shared defaults.
%   Field-agnostic: it never hard-codes a knob name, so a caller adds a knob
%   simply by adding a field to DEFAULTS.
%
%   INPUT
%     satListCell   1xN   cell     Satellite ids (see LSQVCE.ENGINE.NORMALIZESATLIST).
%                                   Required.
%     userSpec      1xM   struct|[] Struct array of overrides. Each element MUST
%                                   have a .satId matching an entry of
%                                   satListCell, plus any subset of the knobs in
%                                   DEFAULTS. [] means "no overrides". Required
%                                   (may be []).
%     defaults      1x1   struct   Fallback values; its field names ARE the
%                                   configurable knobs (e.g. obsSpec, nl,
%                                   arSelect, maxLag). Required.
%
%   OUTPUT
%     satSpecs   1xN   struct   One element per satellite, in satListCell order.
%                               Fields: .satId, .key (= satKeyOf(satId)), and
%                               every field of DEFAULTS (override or fallback).
%     specFor    1x1   function_handle   specFor(satId) -> the element of
%                               satSpecs whose .satId matches (errors if none).
%
%   ERRORS
%     lsqvce:resolveSatSpecs:badDefaults   defaults is not a struct
%     lsqvce:resolveSatSpecs:badSatSpec    userSpec lacks a .satId field
%     lsqvce:resolveSatSpecs:noSpec        specFor called with an unknown id
%
%   EXAMPLE
%     def = struct('nl', 12, 'maxLag', 128);
%     us(1).satId = 'D';  us(1).nl = 9;                 % GRACE-D: fewer locals
%     [S, specFor] = lsqvce.engine.resolveSatSpecs({'C','D'}, us, def);
%     specFor('D').nl    % -> 9 ;   specFor('C').nl -> 12 (default)
%
%   See also LSQVCE.ENGINE.NORMALIZESATLIST, LSQVCE.ITERATESEPARABLE
%
%   Created 2026-07-15 13:40 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.
    narginchk(3, 3);
    if nargin < 2, userSpec = []; end
    if ~isstruct(defaults)
        error('lsqvce:resolveSatSpecs:badDefaults', 'defaults must be a struct.');
    end
    if ~isempty(userSpec) && (~isstruct(userSpec) || ~isfield(userSpec, 'satId'))
        error('lsqvce:resolveSatSpecs:badSatSpec', ...
              'userSpec must be a struct array with a .satId field, or [].');
    end

    flds = fieldnames(defaults);
    n    = numel(satListCell);
    satSpecs = repmat(emptyTemplate(flds), 1, n);
    for i = 1:n
        sid = satListCell{i};
        s   = matchById(userSpec, sid);
        satSpecs(i).satId = sid;
        satSpecs(i).key   = lsqvce.engine.satKeyOf(sid);  % public package helper
        for k = 1:numel(flds)
            satSpecs(i).(flds{k}) = override(s, flds{k}, defaults.(flds{k}));
        end
    end
    specFor = @(satId) lookup(satSpecs, satId);
end

%==========================================================================
function t = emptyTemplate(flds)
%EMPTYTEMPLATE  An empty spec struct with fields {satId,key,flds{:}} all [].
%   INPUT   flds  1xF cell   extra field names to include
%   OUTPUT  t     1x1 struct  scalar template with every field empty
    names = [{'satId'; 'key'}; flds(:)];
    t = cell2struct(repmat({[]}, numel(names), 1), names);
end

function s = matchById(userSpec, satId)
%MATCHBYID  First entry of userSpec whose .satId equals satId ([] if none).
%   INPUT   userSpec 1xN struct array;  satId 1x1 the id to match
%   OUTPUT  s  1x1 struct | []   the matching spec or empty
    s = [];
    for j = 1:numel(userSpec)
        if isequal(userSpec(j).satId, satId), s = userSpec(j); return, end
    end
end

function v = override(s, field, fallback)
%OVERRIDE  s.(field) if s is a struct with a non-empty field, else fallback.
%   INPUT   s 1x1 struct | [];  field 1xk char;  fallback any
%   OUTPUT  v  the resolved value
    if ~isempty(s) && isstruct(s) && isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = fallback;
    end
end

function spec = lookup(satSpecs, satId)
%LOOKUP  The resolved spec for satId (errors if absent).
%   INPUT   satSpecs 1xN struct array;  satId 1x1 the id
%   OUTPUT  spec 1x1 struct  the matching entry
    for i = 1:numel(satSpecs)
        if isequal(satSpecs(i).satId, satId), spec = satSpecs(i); return, end
    end
    error('lsqvce:resolveSatSpecs:noSpec', 'No spec for satellite %s.', string(satId));
end
