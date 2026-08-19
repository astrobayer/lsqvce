function pool = poolGroups(vc, varargin)
%LSQVCE.POOLGROUPS  Suggest (or build) a coarser VCE grouping for groups whose
%   redundancy is too low to estimate a reliable variance component.
%
%   A variance component estimated from r_g observations has a relative standard
%   error of about 0.5*sqrt(2/r_g), so below a redundancy FLOOR the estimated
%   sigma_g is mostly noise. The usual remedy is to POOL such groups with a
%   neighbour (or with each other) so their combined redundancy clears the floor.
%   This helper reads a finished VCE result, identifies the under-redundant
%   groups, and returns a remapping from the old group keys to pooled keys, plus
%   a ready-to-use grouping function you can feed back into the next VCE pass.
%
%   IT DOES NOT RE-SOLVE. Pooling changes the grouping, which changes the whole
%   VCE fixed point, so the honest workflow is: pool -> re-run VCE with the new
%   grouping -> re-check. This function gives you the new grouping; you re-run.
%
%   pool = LSQVCE.POOLGROUPS(vc) uses the default floor and the default strategy
%   ('mergeLow': merge every below-floor group into ONE pooled group).
%
%   pool = LSQVCE.POOLGROUPS(vc, Name, Value, ...) sets the options below.
%
%   INPUT
%     vc        1x1 struct   a result from LSQVCE.VCE or LSQVCE.VARCOMP, carrying
%                            at least:
%                              .groups  1xG cell/│ the group keys (labels)
%                              .r       Gx1 double  redundancy per group (r_g)
%                              .sigma   Gx1 double  variance components (optional,
%                                       used only for the printed report)
%
%   NAME-VALUE
%     'floor'    1x1 double  redundancy floor r_min below which a group is judged
%                            under-redundant (default 50). A group with r_g <
%                            floor is a pooling candidate.
%     'strategy' char        how to pool the candidates:
%                  'mergeLow'  (default) merge ALL below-floor groups into one
%                              pooled group keyed 'pooledLow'. Simple and safe.
%                  'nearest'   merge each below-floor group into the ABOVE-floor
%                              group whose sigma is closest (needs vc.sigma). Keeps
%                              components physically similar together.
%     'poolKey'  char        the label for the merged group under 'mergeLow'
%                            (default 'pooledLow').
%     'verbose'  1x1 logical print a short report (default true).
%
%   OUTPUT  pool  (1x1 struct)
%     .map        containers.Map  oldKey (char) -> newKey (char). Every original
%                 group key maps to itself unless it was pooled.
%     .newGroups  1xGnew cell     the distinct new group keys after pooling.
%     .remap      function_handle  @(oldKey) -> newKey, a thin wrapper over .map
%                 you can compose with your existing groupBy to produce the new
%                 grouping. If your groupBy returns numeric indices, use
%                 .indexRemap instead.
%     .indexRemap Gx1 double      newIndex(oldIndex): maps an old 1..G group index
%                 to a new 1..Gnew index, for numeric-index groupings.
%     .pooled     1xP cell        the old keys that were pooled (r_g < floor).
%     .floor      1x1 double      the floor used.
%     .changed    1x1 logical     true if any pooling happened.
%
%   EXAMPLE
%     vc  = lsqvce.vce(accFcn, 'groupBy','obsType');   % some groups thin
%     pl  = lsqvce.poolGroups(vc, 'floor', 80);
%     if pl.changed
%         % re-run VCE with the coarser grouping (string-key example):
%         newGroupBy = @(varargin) pl.remap(oldGroupBy(varargin{:}));
%         vc2 = lsqvce.vce(accFcn, 'groupBy', newGroupBy);
%     end
%
%   ERRORS
%     lsqvce:poolGroups:badInput    the group specification is malformed
%     lsqvce:poolGroups:needSigma   pooling needs per-group variance factors
%     lsqvce:poolGroups:sizeMismatch  the group and sigma lists differ in length
%
%   See also LSQVCE.VCE, LSQVCE.VARCOMP, LSQVCE.REDUNDANCY.
%
%   Created 2026-07-17 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8. Implements HISTORY_v1.md roadmap item 4 (VCE group-pooling
%   helper): it SUGGESTS/builds a coarser grouping; the caller re-runs VCE.

    ip = inputParser;  ip.FunctionName = 'lsqvce.poolGroups';
    ip.addRequired('vc', @isstruct);
    ip.addParameter('floor',    50, @(v) isscalar(v) && v > 0);
    ip.addParameter('strategy', 'mergeLow', ...
                    @(s) any(strcmpi(s, {'mergeLow','nearest'})));
    ip.addParameter('poolKey',  'pooledLow', @(s) ischar(s) || isstring(s));
    ip.addParameter('verbose',  true, @(v) (islogical(v) || isnumeric(v)) && isscalar(v));
    ip.parse(vc, varargin{:});
    o = ip.Results;

    if ~isfield(vc, 'groups') || ~isfield(vc, 'r')
        error('lsqvce:poolGroups:badInput', ...
              'vc must carry .groups and .r (from lsqvce.vce / lsqvce.varComp).');
    end
    keys = cellfun(@charKey, ensureCell(vc.groups), 'UniformOutput', false);
    r    = vc.r(:);
    G    = numel(keys);
    if numel(r) ~= G
        error('lsqvce:poolGroups:sizeMismatch', ...
              'numel(vc.groups) = %d but numel(vc.r) = %d.', G, numel(r));
    end
    poolKey = char(o.poolKey);

    low = r < o.floor;                    % under-redundant groups
    map = containers.Map('KeyType','char','ValueType','char');
    for g = 1:G, map(keys{g}) = keys{g}; end   % identity by default

    if ~any(low)
        pool = finish(map, keys, {}, o.floor, false, G);
        if lsqvce.util.verbosity.progress(o.verbose)
            lsqvce.util.log.msg(['%sall %d groups clear the r >= %.0f ' ...
                'floor; no pooling needed.\n'], lsqvce.util.verbosity.banner('Group pooling'), G, o.floor);
        end
        return
    end

    switch lower(o.strategy)
        case 'mergelow'
            for g = find(low(:).')
                map(keys{g}) = poolKey;
            end

        case 'nearest'
            if ~isfield(vc, 'sigma') || numel(vc.sigma) ~= G
                error('lsqvce:poolGroups:needSigma', ...
                      '''nearest'' needs vc.sigma aligned with vc.groups.');
            end
            sig  = vc.sigma(:);
            hi   = find(~low);
            if isempty(hi)
                % nothing above the floor to attach to -> fall back to mergeLow
                for g = find(low(:).'), map(keys{g}) = poolKey; end
            else
                for g = find(low(:).')
                    [~, j] = min(abs(sig(hi) - sig(g)));
                    map(keys{g}) = keys{hi(j)};
                end
            end
    end

    pooledKeys = keys(low);
    pool = finish(map, keys, pooledKeys, o.floor, true, G);

    if lsqvce.util.verbosity.progress(o.verbose)
        lsqvce.util.log.msg(['%s%d of %d groups below r = %.0f: %s\n' ...
            '           strategy ''%s'' -> %d group(s) after pooling.\n'], lsqvce.util.verbosity.banner('Group pooling'), ...
            sum(low), G, o.floor, strjoin(pooledKeys, ', '), ...
            lower(o.strategy), numel(pool.newGroups));
    end
end

% ------------------------------------------------------------------------
function pool = finish(map, keys, pooledKeys, fl, changed, G)
%FINISH  Assemble the output struct, including the numeric index remap.
    % distinct new keys, in first-appearance order of the old keys
    newForOld = cellfun(@(k) map(k), keys, 'UniformOutput', false);
    newGroups = uniqueStable(newForOld);
    % numeric index remap: old 1..G -> new 1..Gnew
    idx = zeros(G, 1);
    for g = 1:G
        idx(g) = find(strcmp(newGroups, newForOld{g}), 1);
    end
    pool = struct( ...
        'map',        map, ...
        'newGroups',  {newGroups}, ...
        'remap',      @(k) map(charKey(k)), ...
        'indexRemap', idx, ...
        'pooled',     {pooledKeys}, ...
        'floor',      fl, ...
        'changed',    changed);
end

function c = ensureCell(g)
%ENSURECELL  Return g as a cell array: pass a cell through, wrap a numeric/other
%   array element-wise via num2cell.
%   INPUT  g   1xG cell | numeric array   group keys
%   OUTPUT c   1xG cell                    the same keys as a cell array
    if iscell(g), c = g; else, c = num2cell(g); end
end

function s = charKey(k)
%CHARKEY  Normalise a group key (char/string/numeric) to a char label.
    if ischar(k),         s = k;
    elseif isstring(k),   s = char(k);
    elseif isnumeric(k),  s = num2str(k);
    else,                 s = char(string(k));
    end
end

function u = uniqueStable(c)
%UNIQUESTABLE  Distinct cellstr in first-appearance order (no sort).
    u = {};  seen = containers.Map('KeyType','char','ValueType','logical');
    for i = 1:numel(c)
        if ~isKey(seen, c{i}), seen(c{i}) = true; u{end+1} = c{i}; end %#ok<AGROW>
    end
end
