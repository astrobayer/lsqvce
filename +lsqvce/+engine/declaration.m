function [spec, systems] = declaration(data, cacheDays)
%LSQVCE.ENGINE.DECLARATION  Validate a kind-A data struct and build its loaders.
%   [spec, systems] = lsqvce.engine.declaration(data)
%   [spec, systems] = lsqvce.engine.declaration(data, cacheDays)
%
%   THE one place the data contract (the Manual, chapter 3) is checked and turned
%   into loaders: field names and types by name (lsqvce:fit:* identifiers --
%   this unit keeps fit's identifiers because it IS fit's declaration step,
%   lifted out in R5-D8 round 4 so that Results.toTable can reuse the very
%   same expansion plans for the write-back), the scope tokens, the calendar,
%   the common plan (with validity, R5-D8 round 3) and one lsqvce.engine.makeLoader
%   per system (with the local plan, priors, processing options).
%
%   INPUT   data       1x1 struct   the kind-A declaration (the Manual, chapter 3)
%           cacheDays  1x1 double   resident built days per system (default 1)
%   OUTPUT  spec       1x1 struct   the validated declaration: .commonParams .days
%                                   .systems .ids .commonPlan .validityColumns .nc ...
%           systems    1xK cell     Loader structs (isLoader true; .plan, .nc, .nl, ...)
%   ERRORS  lsqvce:fit:noSystems / noCommon / unknownField / badScope / badId /
%           badParamColumns / ...  (the fit help lists them)
%
%   See also LSQVCE.FIT, LSQVCE.ENGINE.MAKELOADER, LSQVCE.RESULTS
%
%   Created 2026-09-03 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5-D8 round 4 (lifted from fit.m,
%   where these functions lived since R4 Level L2).
    if nargin < 2, cacheDays = 1; end
    spec = validateData(data);
    systems = buildSystems(spec, cacheDays);
end

function spec = validateData(data)
%VALIDATEDATA  Enforce Manual ch. 3 (§2) (v3) on the raw struct; fail fast.
%   OUTPUT spec 1x1 struct  .commonParams (RxC cell), .days (double),
%                           .systems (1xM), .ids (1xM cellstr, resolved)
    topAllowed = {'commonParams', 'days', 'systems', 'priors', 'paramColumns'};   % priors: R5 M5 5b; paramColumns: R5-D8 round 3
    rejectUnknown(fieldnames(data), topAllowed, 'data');

    if ~isfield(data, 'systems') || ~isstruct(data.systems) || isempty(data.systems)
        error('lsqvce:fit:noSystems', ...
            'data.systems must be a non-empty 1xM struct array (one entry per observing system).');
    end
    sys = data.systems(:).';

    % ---- the common unknowns: the SCOPED table is the declaration (v4).
    % Column 1 = scope token 'common' | 'day'; rows are TEMPLATES ('day' rows
    % expand over the global day calendar); nc is DERIVED by the common plan.
    if ~isfield(data, 'commonParams') || isempty(data.commonParams) || ...
            ~iscell(data.commonParams) || ~ismatrix(data.commonParams)
        error('lsqvce:fit:noCommon', ...
            ['data.commonParams is required: an R x C cell, one ROW per common ' ...
             'template, column 1 the scope token ''common'' | ''day''. ' ...
             'Minimal table: [repmat({''common''}, nc, 1), num2cell((1:nc).'')].']);
    end
    for r = 1:size(data.commonParams, 1)
        v = data.commonParams{r, 1};
        if ~(ischar(v) || (isstring(v) && isscalar(v))) || ...
                ~any(strcmp(char(v), {'common', 'day'}))
            error('lsqvce:fit:badScope', ...
                ['data.commonParams row %d: column 1 must be the scope token ' ...
                 '''common'' | ''day'' (''arc'' is system-specific and lives in ' ...
                 'localParams).'], r);
        end
    end

    % ---- days
    days = [];
    if isfield(data, 'days') && ~isempty(data.days)
        days = double(data.days(:).');
        if ~all(isfinite(days) & days == round(days))
            error('lsqvce:fit:noDays', 'data.days must be numeric integers (the loader cache is keyed on them).');
        end
    end

    % ---- systems
    sysRequired = {'buildDay', 'covariance'};
    sysAllowed  = [sysRequired, {'localParams', 'arcs', 'days', ...
                   'id', 'name', 'processingOptions', 'epochCovariance', 'priors'}];   % priors: R5-D8 round 5
    rejectUnknown(fieldnames(sys), sysAllowed, 'data.systems');
    ids = cell(1, numel(sys));
    for k = 1:numel(sys)
        s = sys(k);
        for f = sysRequired
            if ~isfield(s, f{1}) || isempty(s.(f{1}))
                error('lsqvce:fit:missingField', ...
                    'data.systems(%d) needs the field ''%s''.', k, f{1});
            end
        end
        if ~isa(s.buildDay, 'function_handle')
            error('lsqvce:fit:builderKind', ...
                'data.systems(%d): buildDay must be a function handle.', k);
        end
        if ~isstruct(s.covariance) || ~isscalar(s.covariance) || isempty(fieldnames(s.covariance))
            error('lsqvce:fit:missingField', ...
                ['data.systems(%d).covariance must be a scalar struct with one field per ' ...
                 'observation type (the fieldnames ARE the observation types).'], k);
        end
        if isfield(s, 'arcs') && (ischar(s.arcs) || isstring(s.arcs))
            error('lsqvce:fit:badArcs', ...
                ['data.systems(%d).arcs = ''%s'': the string form (''fromDay'') is ' ...
                 'ENGINE-tier only -- it reads data before the run. Declare a count, ' ...
                 'an id vector, or a resolver @(day,obsType).'], k, char(s.arcs));
        end
        if isfield(s, 'id') && ~isempty(s.id)
            id = s.id;
            if ~((isnumeric(id) && isscalar(id)) || (ischar(id) && isscalar(id)) || ...
                 (isstring(id) && isscalar(id) && strlength(id) == 1))
                error('lsqvce:fit:badId', ...
                    'data.systems(%d).id must be a numeric scalar or a single character.', k);
            end
            ids{k} = char(string(id));
        else
            ids{k} = char(string(k));               % the index default
        end
        if isempty(days) && ~(isfield(s, 'days') && ~isempty(s.days))
            error('lsqvce:fit:noDays', ...
                'no day list: set data.days or data.systems(%d).days.', k);
        end
    end
    if numel(unique(ids)) < numel(ids)
        error('lsqvce:fit:duplicateSystem', ...
            'system ids must be unique across data.systems (explicit ids and index defaults share one namespace).');
    end
    spec = struct('commonParams', {data.commonParams}, 'days', days, ...
                  'systems', sys, 'ids', {ids});
    spec.validityColumns = [];
    spec.nameColumn = 2;                                   % R5-D8 round 6: where the names are
    if isfield(data, 'paramColumns') && ~isempty(data.paramColumns)
        pc = data.paramColumns;
        if ~isstruct(pc)
            error('lsqvce:fit:badParamColumns', ...
                  'data.paramColumns must be a struct (fields start, end, name) naming columns of your tables.');
        end
        bad = setdiff(fieldnames(pc), {'start', 'end', 'name'});
        if ~isempty(bad)
            error('lsqvce:fit:badParamColumns', 'data.paramColumns has an unknown field ''%s''.', bad{1});
        end
        if isfield(pc, 'start') || isfield(pc, 'end')
            if ~(isfield(pc, 'start') && isfield(pc, 'end'))
                error('lsqvce:fit:badParamColumns', 'data.paramColumns needs BOTH start and end for validity.');
            end
            spec.validityColumns = [double(pc.start), double(pc.end)];
        end
        if isfield(pc, 'name') && ~isempty(pc.name), spec.nameColumn = double(pc.name); end
    end
    spec.commonPlan = buildCommonPlan(spec);
end

function plan = buildCommonPlan(spec)
%BUILDCOMMONPLAN  Expand the SCOPED commonParams over the GLOBAL day calendar
%   (L11 WP-75, v4). Layout, binding: the 'common' rows first (table order),
%   then one 'day' block per calendar day (sorted union of every declared day
%   list) whose SCHEDULE carries at least one arc in at least one system --
%   schedules are resolvable without reading data by contract, so the plan is
%   built before any read. nc is DERIVED.
%   OUTPUT plan 1x1 struct  .nc .commonCount .days (calendar) .dayOffset
%          .dayCount (1xD) .rows (.common / .day template row indices)
    cp  = spec.commonParams;
    tok = cellfun(@(v) char(string(v)), cp(:, 1), 'UniformOutput', false);
    iCommon = find(strcmp(tok, 'common')).';
    iDay    = find(strcmp(tok, 'day')).';
    nCom = numel(iCommon);  nDayRows = numel(iDay);
    cal = spec.days;
    for k = 1:numel(spec.systems)
        s = spec.systems(k);
        if isfield(s, 'days') && ~isempty(s.days)
            cal = union(cal, double(s.days(:).'));
        end
    end
    cal = sort(cal(:).');
    D = numel(cal);
    dayOffset = zeros(1, D);  dayCount = zeros(1, D);  dayRowsOf = cell(1, D);
    valid = commonValidity(cp, spec.validityColumns);        % R5-D8 round 3
    off = nCom;
    if nDayRows > 0
        for i = 1:D
            rows = iDay(arrayfun(@(r) valid(r, cal(i)), iDay));
            dayRowsOf{i} = rows;
            if calendarDayHasArcs(spec.systems, cal(i)) && ~isempty(rows)
                dayOffset(i) = off;  dayCount(i) = numel(rows);
                off = off + numel(rows);
            end
        end
    end
    plan = struct('nc', off, 'commonCount', nCom, 'days', cal, ...
                  'dayOffset', dayOffset, 'dayCount', dayCount, ...
                  'rows', struct('common', iCommon, 'day', iDay), 'dayRowsOf', {dayRowsOf});
end

function tf = calendarDayHasArcs(sys, day)
%CALENDARDAYHASARCS  Does at least one system's SCHEDULE carry an arc on DAY?
%   Mirrors the engine's arc resolution on the declared forms (default 1;
%   count; id vector; per-type struct; @(day,obsType) resolver) -- resolvable
%   pre-run by contract.
    tf = false;
    for k = 1:numel(sys)
        s = sys(k);
        if isfield(s, 'days') && ~isempty(s.days) && ~any(double(s.days) == day)
            continue                                   % this system skips the day
        end
        if ~isfield(s, 'arcs') || isempty(s.arcs)
            tf = true;  return                         % default: one arc per day
        end
        a = s.arcs;
        types = fieldnames(s.covariance);
        for j = 1:numel(types)
            aj = a;
            if isstruct(a)
                if ~isfield(a, types{j}), continue, end
                aj = a.(types{j});
            end
            if isa(aj, 'function_handle')
                if ~isempty(aj(day, types{j})), tf = true;  return, end
            elseif isnumeric(aj)
                if (isscalar(aj) && aj >= 1) || (~isscalar(aj) && ~isempty(aj))
                    tf = true;  return
                end
            end
        end
    end
end

function rejectUnknown(names, allowed, where)
%REJECTUNKNOWN  Unknown field names are errors (typo protection, rule 5).
    bad = setdiff(names, allowed);
    if ~isempty(bad)
        error('lsqvce:fit:unknownField', ...
            '%s has unknown field ''%s''. Allowed: %s.', where, bad{1}, strjoin(allowed, ', '));
    end
end

%==========================================================================
function systems = buildSystems(spec, cacheDays)
    if nargin < 2, cacheDays = 1; end
%BUILDSYSTEMS  One Loader per system through LSQVCE.ENGINE.MAKELOADER.
%   Derivations (v4): ObservationTypes = fieldnames(covariance); nc rides in
%   as the COMMON PLAN (fit expands the scoped commonParams over the global
%   calendar); nl is derived by makeLoader's scope plan from the scoped
%   LocalParams table, or NumLocal = 0 when none is declared.
    sys  = spec.systems;
    days = spec.days;
    % ---- the ONE mapping table: lowerCamelCase field -> makeLoader NV name
    map = {'covariance', 'Covariance'; 'arcs', 'Arcs'; ...
           'localParams', 'LocalParams'; ...
           'days', 'Days'; 'name', 'Name'; ...
           'processingOptions', 'ProcessingOptions'; ...
           'epochCovariance', 'EpochCovariance'; ...
           'priors', 'Priors'};                           % R5-D8 round 5
    systems = cell(1, numel(sys));
    for k = 1:numel(sys)
        s  = sys(k);
        nv = {};
        for r = 1:size(map, 1)
            f = map{r, 1};
            if isfield(s, f) && ~isempty(s.(f))
                nv = [nv, {map{r, 2}, s.(f)}];                       %#ok<AGROW>
            end
        end
        nv = [nv, {'ObservationTypes', fieldnames(s.covariance).'}]; %#ok<AGROW>
        nv = [nv, {'SystemId', idValue(spec.ids{k})}];               %#ok<AGROW>
        nv = [nv, {'CommonParams', spec.commonParams}];              %#ok<AGROW>
        nv = [nv, {'Scoped', true, 'CommonPlan', spec.commonPlan}];  %#ok<AGROW>
        nv = [nv, {'CacheDays', cacheDays}];                         %#ok<AGROW>  R5 M4 4d
        if ~isempty(spec.validityColumns)
            nv = [nv, {'ValidityColumns', spec.validityColumns}];    %#ok<AGROW>  R5-D8 round 3
        end
        nv = [nv, {'NameColumn', spec.nameColumn}];                   %#ok<AGROW>  R5-D8 round 6
        if ~(isfield(s, 'days') && ~isempty(s.days)) && ~isempty(days)
            nv = [nv, {'Days', days}];                               %#ok<AGROW>
        end
        if ~(isfield(s, 'localParams') && ~isempty(s.localParams))
            nv = [nv, {'NumLocal', 0}];                              %#ok<AGROW>
        end
        try
            systems{k} = lsqvce.engine.makeLoader(s.buildDay, nv{:});
        catch e
            if startsWith(e.identifier, 'lsqvce:makeLoader:')
                e2 = MException(strrep(e.identifier, 'lsqvce:makeLoader:', 'lsqvce:fit:'), ...
                                'data.systems(%d): %s', k, e.message);
                throw(addCause(e2, e));
            end
            rethrow(e);
        end
    end
end

function v = idValue(idChar)
%IDVALUE  A resolved id back to the type makeLoader stores: numeric where the
%   char is a number ('1' -> 1, index defaults and explicit numerics), the
%   single char otherwise ('C' -> 'C').
    n = str2double(idChar);
    if isscalar(idChar) && isnan(n)
        v = idChar;
    else
        v = n;
    end
end

%========================================================== kinds B / C
function f = commonValidity(tbl, valCols)
%COMMONVALIDITY  @(row, day) validity of a commonParams row on a calendar day
%   (R5-D8 round 3; the same rule as the loader's: [MJD sec] or fractional
%   MJD in the named columns, [] = always, half-open day overlap).
    if isempty(valCols), f = @(row, day) true;  return, end
    if numel(valCols) ~= 2 || any(valCols < 1) || any(valCols > size(tbl, 2))
        error('lsqvce:fit:badParamColumns', ...
              'data.paramColumns names columns %s; commonParams has %d columns.', mat2str(valCols), size(tbl, 2));
    end
    n = size(tbl, 1);  t0 = -inf(n, 1);  t1 = inf(n, 1);
    for r = 1:n
        a = tbl{r, valCols(1)};  b = tbl{r, valCols(2)};
        if ~isempty(a), t0(r) = mjdOf(a); end
        if ~isempty(b), t1(r) = mjdOf(b); end
    end
    f = @(row, day) t0(row) < day + 1 && t1(row) > day;
end

function m = mjdOf(v)
    if isscalar(v), m = double(v); else, m = double(v(1)) + double(v(2)) / 86400; end
end
