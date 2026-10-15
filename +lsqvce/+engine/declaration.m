function [spec, systems] = declaration(data, cacheDays)
%LSQVCE.ENGINE.DECLARATION  Validate a kind-A data struct and build its loaders.
%   [spec, systems] = lsqvce.engine.declaration(data)
%   [spec, systems] = lsqvce.engine.declaration(data, cacheDays)
%
%   THE one place the data contract (the API reference) is checked and turned
%   into loaders: field names and types by name (lsqvce:fit:* identifiers --
%   this unit keeps fit's identifiers because it IS fit's declaration step,
%   lifted out in R5-D8 round 4 so that Results.toTable can reuse the very
%   same expansion plans for the write-back), the scope tokens, the calendar,
%   the common plan (with validity, R5-D8 round 3) and one lsqvce.engine.makeLoader
%   per system (with the local plan, priors, processing options).
%
%   INPUT   data       1x1 struct   the kind-A declaration (the API reference)
%           cacheDays  1x1 double   resident built days per system (default 1)
%   OUTPUT  spec       1x1 struct   the validated declaration: .commonParams .days
%                                   .systems .ids .commonPlan .validityColumns .nc ...
%           systems    1xK cell     Loader structs (isLoader true; .plan, .nc, .nl, ...)
%   ACCURACY
%           Not applicable: validation and plan construction; the one number it
%           forms, a validity entry's MJD, is the private validityOf's.
%   ERRORS  lsqvce:fit:noSystems / noCommon / unknownField / badScope / badId /
%           badParamColumns / ...  (lsqvce.fit's ERRORS lists them)
%
%   See also LSQVCE.FIT, LSQVCE.ENGINE.MAKELOADER, LSQVCE.RESULTS
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-141: data.systems(k).prepareDay
%   is a known field, refused by builderKind when it is not a function handle,
%   and handed to makeLoader as PrepareDay. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-138 (code review of WP-135): the
%   validity rule is the private validityOf (mjdOf retired: one rule, the
%   loader's), the id and day rules the private idDayRules (the accumulator's
%   too); a non-integer day list is lsqvce:fit:badDays (it was noDays), ids of
%   two kinds lsqvce:fit:mixedIdTypes (it was badId); ACCURACY. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-14 (CEST, UTC+02:00), WP-135 (review 3 findings 11, 12 and
%   the declaration lows): day lists are unique integers
%   (lsqvce:fit:duplicateDay); the system ids are of one kind, 'i' / 'j' stay
%   characters and a numeric id is finite (lsqvce:fit:badId); a system without
%   its own days claims only data.days for the common 'day' blocks; a common
%   validity entry of another type, or a start after its end, is refused
%   (lsqvce:fit:badValidity). With assistance from Claude (Anthropic), model
%   Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-131: a one-argument arc
%   resolver @(day) works with day-scoped commonParams (it threw "too many
%   inputs"). With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-118 part 1 (R7-34): points
%   at lsqvce.fit's ERRORS, which now lists badParamColumns. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-09-09 (CEST, UTC+02:00), R6 WP-102 (v6): the arc schedule
%   is optional -- a system without `arcs` claims every day of its calendar
%   for the common 'day' rows; the string form 'fromDay' is refused at the
%   declaration (lsqvce:fit:badArcs). With assistance from Claude
%   (Anthropic), model Claude Fable 5.1.
%   Created 2026-09-03 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5-D8 round 4 (lifted from fit.m,
%   where these functions lived since R4 Level L2).
%
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-124: references to the data
%   contract point at the API reference. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
    if nargin < 2, cacheDays = 1; end
    spec = validateData(data);
    systems = buildSystems(spec, cacheDays);
end

function spec = validateData(data)
%VALIDATEDATA  Enforce the API reference (The data struct) on the raw struct; fail fast.
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
        days = dayListOf(data.days, 'data.days');
    end

    % ---- systems
    sysRequired = {'buildDay', 'covariance'};
    sysAllowed  = [sysRequired, {'localParams', 'arcs', 'days', ...
                   'id', 'name', 'processingOptions', 'priors', ...  % priors: R5-D8 round 5;
                   'prepareDay'}];                   % WP-141; epochCovariance retired R6 WP-90
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
        if isfield(s, 'prepareDay') && ~isempty(s.prepareDay) && ~isa(s.prepareDay, 'function_handle')
            error('lsqvce:fit:builderKind', ...                          % WP-141
                'data.systems(%d): prepareDay must be a function handle.', k);
        end
        if ~isstruct(s.covariance) || ~isscalar(s.covariance) || isempty(fieldnames(s.covariance))
            error('lsqvce:fit:missingField', ...
                ['data.systems(%d).covariance must be a scalar struct with one field per ' ...
                 'observation type (the fieldnames ARE the observation types).'], k);
        end
        if isfield(s, 'arcs') && (ischar(s.arcs) || isstring(s.arcs))
            error('lsqvce:fit:badArcs', ...
                ['data.systems(%d).arcs = ''%s'': the string form is retired (v6). Omit ' ...
                 'arcs -- the cells buildDay returns are the arcs of the day -- or declare ' ...
                 'a count, an id vector, a per-stream struct or a resolver @(day,obsType).'], ...
                k, char(s.arcs));
        end
        if isfield(s, 'id') && ~isempty(s.id)
            id = s.id;
            if ~((isnumeric(id) && isscalar(id) && isreal(id) && isfinite(id)) || ...
                 (ischar(id) && isscalar(id)) || ...
                 (isstring(id) && isscalar(id) && strlength(id) == 1))
                error('lsqvce:fit:badId', ...
                    'data.systems(%d).id must be a finite numeric scalar or a single character.', k);
            end
            ids{k} = char(string(id));
        else
            ids{k} = char(string(k));               % the index default
        end
        if isfield(s, 'days') && ~isempty(s.days)
            dayListOf(s.days, sprintf('data.systems(%d).days', k));
        elseif isempty(days)
            error('lsqvce:fit:noDays', ...
                'no day list: set data.days or data.systems(%d).days.', k);
        end
    end
    if numel(unique(ids)) < numel(ids)
        error('lsqvce:fit:duplicateSystem', ...
            'system ids must be unique across data.systems (explicit ids and index defaults share one namespace).');
    end
    oneIdKind(ids, arrayfun(@(t) isfield(t, 'id') && ~isempty(t.id), sys));
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
            if calendarDayHasArcs(spec.systems, cal(i), spec.days) && ~isempty(rows)
                dayOffset(i) = off;  dayCount(i) = numel(rows);
                off = off + numel(rows);
            end
        end
    end
    plan = struct('nc', off, 'commonCount', nCom, 'days', cal, ...
                  'dayOffset', dayOffset, 'dayCount', dayCount, ...
                  'rows', struct('common', iCommon, 'day', iDay), 'dayRowsOf', {dayRowsOf});
end

function tf = calendarDayHasArcs(sys, day, globalDays)
%CALENDARDAYHASARCS  Does at least one system's SCHEDULE carry an arc on DAY?
%   Mirrors the engine's arc resolution on the declared forms (count; id
%   vector; per-type struct; @(day,obsType) resolver) -- resolvable pre-run
%   by contract. A system WITHOUT a schedule (v6, R6 WP-102: the cells its
%   buildDay returns are its arcs) claims every day of its calendar -- the
%   calendar is the declaration, and a day that then returns no data is
%   refused at its 'meta' probe (lsqvce:makeLoader:emptyDay). A system
%   without its own days reads data.days only (buildSystems), so it claims
%   no other day (WP-135: it claimed every calendar day, and a 'day' block
%   nobody observed made the solve singular).
    tf = false;
    for k = 1:numel(sys)
        s = sys(k);
        if isfield(s, 'days') && ~isempty(s.days)
            if ~any(double(s.days) == day), continue, end   % this system skips the day
        elseif ~any(globalDays == day)
            continue                                   % its loader reads data.days only
        end
        if ~isfield(s, 'arcs') || isempty(s.arcs)
            tf = true;  return                         % no schedule: the calendar is
        end                                            % the declaration (v6)
        a = s.arcs;
        types = fieldnames(s.covariance);
        for j = 1:numel(types)
            aj = a;
            if isstruct(a)
                if ~isfield(a, types{j}), continue, end
                aj = a.(types{j});
            end
            if isa(aj, 'function_handle')
                % WP-131: a resolver may take the day alone, as makeLoader's
                % resolveArcSpec accepts -- @(d) ... threw "too many inputs" here
                if nargin(aj) == 1, arcs = aj(day); else, arcs = aj(day, types{j}); end
                if ~isempty(arcs), tf = true;  return, end
            elseif isnumeric(aj)
                if (isscalar(aj) && aj >= 1) || (~isscalar(aj) && ~isempty(aj))
                    tf = true;  return
                end
            end
        end
    end
end

function days = dayListOf(v, where)
%DAYLISTOF  A declared day list: numeric integers (lsqvce:fit:badDays --
%   WP-138: it was noDays, the reason for a list that is absent), each day
%   ONCE (lsqvce:fit:duplicateDay, WP-135: a repeated day was accumulated
%   twice). The repeat rule is the accumulator's too: private idDayRules.
%   INPUT  v 1xD numeric; where 1xN char (the field, for the message)
%   OUTPUT days 1xD double
    days = NaN;                                  % a cell / char is no day list
    if isnumeric(v), days = double(v(:).'); end
    if ~all(isfinite(days) & days == round(days))
        error('lsqvce:fit:badDays', '%s must be numeric integers (the loader cache is keyed on them).', where);
    end
    [~, twice] = idDayRules({}, days);
    if ~isempty(twice)
        error('lsqvce:fit:duplicateDay', ...
            '%s lists day %d more than once; every listing would be accumulated.', where, twice);
    end
end

function oneIdKind(ids, explicit)
%ONEIDKIND  The system ids of a declaration are of ONE kind (WP-135, finding
%   12; the rule is the accumulator's too: private idDayRules). WP-138: its
%   own identifier, lsqvce:fit:mixedIdTypes (the accumulator's reason name;
%   it was badId, which keeps an id that is no number or letter at all).
%   INPUT  ids 1xM cellstr (resolved); explicit 1xM logical (id declared)
    mixed = idDayRules(cellfun(@idValue, ids, 'UniformOutput', false), []);
    if ~isempty(mixed)
        kc = mixed(1);  kn = mixed(2);
        how = '';
        if ~explicit(kn), how = ' (its index default: it declares no id)'; end
        error('lsqvce:fit:mixedIdTypes', ...
            ['system ids must be all numbers or all single letters: data.systems(%d) has ' ...
             'the id ''%s'', data.systems(%d) the number %s%s. Give every system an id ' ...
             'of one kind.'], kc, ids{kc}, kn, ids{kn}, how);
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
           'priors', 'Priors'; ...                        % R5-D8 round 5
           'prepareDay', 'PrepareDay'};                   % WP-141
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
    if isscalar(idChar) && (isnan(n) || ~isreal(n))    % WP-135: 'i' / 'j' parse as 1i
        v = idChar;
    else
        v = n;
    end
end

%========================================================== kinds B / C
function f = commonValidity(tbl, valCols)
%COMMONVALIDITY  @(row, day) validity of a commonParams row on a calendar day
%   (R5-D8 round 3): the loader's rule on the named columns -- the private
%   validityOf, ONE definition since WP-138 (it was re-implemented here as
%   mjdOf). Its refusal is relabelled lsqvce:fit:badValidity, the rule's
%   error attached as the cause; the columns themselves are fit's to judge.
    if isempty(valCols), f = @(row, day) true;  return, end
    if numel(valCols) ~= 2 || any(valCols < 1) || any(valCols > size(tbl, 2)) || any(valCols ~= round(valCols))
        error('lsqvce:fit:badParamColumns', ...
              'data.paramColumns names columns %s; commonParams has %d columns.', mat2str(valCols), size(tbl, 2));
    end
    try
        f = validityOf(tbl, valCols);
    catch e
        if ~strcmp(e.identifier, 'lsqvce:makeLoader:badValidity'), rethrow(e), end
        e2 = MException('lsqvce:fit:badValidity', 'data.commonParams %s', e.message);
        throw(addCause(e2, e));
    end
end
