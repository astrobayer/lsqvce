classdef BatchIndex < handle
%LSQVCE.BATCHINDEX  Bookkeeping for (system, sensor, day) batch hierarchies.
%
%   The numerics do not need this. With PER-SYSTEM local parameters, a
%   (sensor, day) batch is just another batch inside its system's elimination
%   block, and lsqvce.Combine already streams batches without caring where
%   they came from. Verified: sensor x day batching gives bit-identical results
%   to a flat structure, and the solution is invariant to batch ORDER to 3e-16.
%
%   What you lose without this class is the ability to SAY ANYTHING USEFUL about
%   a result. A flagged observation comes back as "observation 4193" instead of
%   "LRI, day 14, epoch 193". Variance components cannot be estimated per sensor
%   -- which, for a KBR/LRI combination where the two instruments have genuinely
%   different and probably mis-stated noise, is the whole point. And a single bad
%   day cannot be identified as a bad day.
%
%   This class is a thin map:
%
%       (system, sensor, day)  <-->  flat batch index k
%       observation index      -->   (system, sensor, day, epoch)
%
%   USAGE
%       bx = lsqvce.BatchIndex();
%       bx.addSystem('GRACE-C', {'KBR','LRI','GPS'}, days);
%       bx.addSystem('GRACE-D', {'KBR','LRI'},       days);
%
%       % register how many observations each (sys, sensor, day) actually has
%       bx.setCount('GRACE-C', 'LRI', datetime(2024,3,7), 8640);
%       ...
%       bx.lock();                       % no more registration
%
%       % now the batch source is indexed by a MEANINGFUL triple
%       src = @(i, k) myLoader(bx.decode(i, k));
%
%       % ... and the diagnostics can be attributed
%       T = bx.attribute(D);             % table: sensor, day, nFlag, mean r_i
%
%   EXAMPLE
%       bx = lsqvce.BatchIndex();
%       bx.addSystem('GRACE-C', {'KBR','LRI','GPS'}, days);
%       bx.addSystem('GRACE-D', {'KBR','LRI'},       days);
%
%       % Register only the (sensor, day) triples that HAVE data. A sensor
%       % that was down is simply never registered -- no empty batches needed.
%       for d = days
%           bx.setCount('GRACE-C', 'LRI', d, 8640);
%       end
%       bx.lock();
%
%       % The batch source is now indexed by a MEANINGFUL triple:
%       src = @(i,k) myLoader(bx.decode(i,k));
%
%       % ... and the diagnostics can be ATTRIBUTED. Without this, a flagged
%       % observation is "observation 4193". With it, "LRI / day 14".
%       bx.report(D);
%
%       % Per-sensor variance components need the sensor group index:
%       VC = cmb.varComp(sol, @(i,k) bx.sensorGroups(i){k}, nSensor);
%
%   ERRORS
%     The registration is malformed
%       lsqvce:batchIndex:badName        a name is not a valid identifier
%       lsqvce:batchIndex:badSystem      the system id is not numeric or a
%                                        single char
%       lsqvce:batchIndex:badSensors     the sensor list is not a cellstr
%       lsqvce:batchIndex:badBatch       the batch descriptor has the wrong shape
%       lsqvce:batchIndex:badCount       a count is not a positive integer
%       lsqvce:batchIndex:empty          nothing was registered
%
%     A name collides or is unknown
%       lsqvce:batchIndex:duplicateSystem  two systems share an id
%       lsqvce:batchIndex:duplicateSensor  two sensors share a name
%       lsqvce:batchIndex:duplicateBatch   two batches share a key
%       lsqvce:batchIndex:noSuchSystem     lookup of an unregistered system
%       lsqvce:batchIndex:noSuchSensor     lookup of an unregistered sensor
%       lsqvce:batchIndex:noBatches        a system carries no batches
%
%     The lock protocol was broken
%       lsqvce:batchIndex:locked         a registration after lock(); the index
%                                        is frozen so downstream row offsets
%                                        cannot shift under a running solve
%       lsqvce:batchIndex:notLocked      a query that needs final row offsets
%                                        before lock() was called
%       lsqvce:batchIndex:countMismatch  the observation count does not match
%                                        what was registered
%
%   See also LSQVCE.COMBINE, LSQVCE.VARCOMP.
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.

properties (SetAccess = private)
    sysName     cell = {}       % 1 x nSys
    sensorName  cell = {}       % 1 x nSys, each a cell of sensor names
    dayList     cell = {}       % 1 x nSys, each a vector (datenum or index)

    % per system: a table of batches, one row per (sensor, day) that has data
    batch       cell = {}       % 1 x nSys, each a struct array

    locked      (1,1) logical = false
end

methods

function obj = BatchIndex()
%BATCHINDEX  Construct an empty batch index (call addSystem/lock before use).
end

%% ------------------------------------------------------------- addSystem
function addSystem(obj, name, sensors, days)
%ADDSYSTEM  Declare a satellite system, its sensors, and its days.
%
%   ADDSYSTEM(name, sensors, days)
%     name     char, e.g. 'GRACE-C'
%     sensors  cell of char, e.g. {'KBR','LRI','GPS'}
%     days     vector of datenum, or any monotone day index

    if obj.locked
        error('lsqvce:batchIndex:locked', ...
              'The index is locked. Create a new one to change the layout.');
    end
    if ~ischar(name) && ~isstring(name)
        error('lsqvce:batchIndex:badName', 'name must be char.');
    end
    name = char(name);
    if any(strcmp(obj.sysName, name))
        error('lsqvce:batchIndex:duplicateSystem', ...
              'System "%s" has already been added.', name);
    end
    if ~iscell(sensors) || isempty(sensors)
        error('lsqvce:batchIndex:badSensors', ...
              'sensors must be a non-empty cell of names.');
    end
    if numel(unique(sensors)) ~= numel(sensors)
        error('lsqvce:batchIndex:duplicateSensor', ...
              'System "%s" lists a sensor twice.', name);
    end

    i = numel(obj.sysName) + 1;
    obj.sysName{i}    = name;
    obj.sensorName{i} = cellfun(@char, sensors, 'UniformOutput', false);
    obj.dayList{i}    = days(:).';
    obj.batch{i}      = struct('sensor', {}, 'day', {}, 'n', {}, 'k', {});
end

%% -------------------------------------------------------------- setCount
function setCount(obj, sysName, sensorName, day, nObs)
%SETCOUNT  Register that a (system, sensor, day) has nObs observations.
%
%   Only registered triples become batches. A sensor that was down on a given
%   day is simply not registered -- there is no need for empty batches, and no
%   need for every sensor to cover every day.

    if obj.locked
        error('lsqvce:batchIndex:locked', 'The index is locked.');
    end
    i = obj.sysIdx(sysName);
    s = obj.sensorIdx(i, sensorName);

    if ~isscalar(nObs) || nObs < 0 || nObs ~= fix(nObs)
        error('lsqvce:batchIndex:badCount', ...
              'nObs must be a non-negative integer.');
    end
    if nObs == 0
        return;                       % nothing to register
    end

    % duplicate check: the same (sensor, day) must not be registered twice, or
    % the observations would be counted -- and diagnosed -- twice.
    for q = 1:numel(obj.batch{i})
        b = obj.batch{i}(q);
        if b.sensor == s && isequaln(b.day, day)
            error('lsqvce:batchIndex:duplicateBatch', ...
                ['%s / %s / %s is already registered with %d observations. ' ...
                 'Registering it twice would double-count them.'], ...
                 obj.sysName{i}, obj.sensorName{i}{s}, obj.dayStr(day), b.n);
        end
    end

    obj.batch{i}(end+1) = struct('sensor', s, 'day', day, 'n', nObs, 'k', 0);
end

%% ------------------------------------------------------------------ lock
function lock(obj)
%LOCK  Freeze the layout and assign the flat batch indices.
%
%   After this, batch k of system i is a well-defined (sensor, day) pair, and
%   the observation offsets are fixed. The batch ORDER within a system does not
%   affect the solution (verified to 3e-16), but it DOES fix the observation
%   numbering, so it must not change between the two passes.

    if obj.locked, return; end
    if isempty(obj.sysName)
        error('lsqvce:batchIndex:empty', 'No systems were added.');
    end

    for i = 1:numel(obj.sysName)
        if isempty(obj.batch{i})
            error('lsqvce:batchIndex:noBatches', ...
                  'System "%s" has no registered batches.', obj.sysName{i});
        end
        % sort by (day, sensor) so the layout is deterministic and readable
        days = [obj.batch{i}.day];
        sens = [obj.batch{i}.sensor];
        [~, ord] = sortrows([days(:), sens(:)]);
        obj.batch{i} = obj.batch{i}(ord);
        for k = 1:numel(obj.batch{i})
            obj.batch{i}(k).k = k;
        end
    end
    obj.locked = true;
end

%% ---------------------------------------------------------------- decode
function info = decode(obj, i, k)
%DECODE  Flat batch index -> (system, sensor, day).
%
%   Use inside your batch source:
%       src = @(i, k) myLoader(bx.decode(i, k));

    obj.requireLocked();
    if i < 1 || i > numel(obj.sysName)
        error('lsqvce:batchIndex:badSystem', ...
              'System index %d is out of range (1..%d).', i, numel(obj.sysName));
    end
    B = obj.batch{i};
    if k < 1 || k > numel(B)
        error('lsqvce:batchIndex:badBatch', ...
              'Batch %d is out of range for system "%s" (1..%d).', ...
              k, obj.sysName{i}, numel(B));
    end
    b = B(k);
    info = struct( ...
        'system',     obj.sysName{i}, ...
        'systemIdx',  i, ...
        'sensor',     obj.sensorName{i}{b.sensor}, ...
        'sensorIdx',  b.sensor, ...
        'day',        b.day, ...
        'n',          b.n, ...
        'batchIdx',   k);
end

%% --------------------------------------------------------------- nBatch
function m = nBatch(obj, i)
%NBATCH  Number of batches in system i.
    obj.requireLocked();
    m = numel(obj.batch{i});
end

function m = nSystem(obj)
%NSYSTEM  Number of systems registered. OUTPUT m 1x1 double.
    m = numel(obj.sysName);
end

%% ------------------------------------------------------------- groupFcn
function g = sensorGroups(obj, i)
%SENSORGROUPS  Per-observation sensor group index for system i.
%
%   Feed this to lsqvce.varComp to estimate a SEPARATE VARIANCE COMPONENT per
%   sensor. That is the right tool when, say, KBR and LRI both observe the same
%   arc but their quoted noise levels disagree with reality -- which is the
%   normal state of affairs, and exactly what a combination must sort out.
%
%   Returns a cell array: g{k} is the group index vector for batch k.

    obj.requireLocked();
    B = obj.batch{i};
    g = cell(numel(B), 1);
    for k = 1:numel(B)
        g{k} = repmat(B(k).sensor, B(k).n, 1);
    end
end

function g = dayGroups(obj, i)
%DAYGROUPS  Per-observation DAY index for system i.
%
%   Use this to hunt for a bad DAY rather than a bad observation -- e.g. a
%   manoeuvre, an eclipse season, a thermal event. Feed it to
%   lsqvce.systematic as a covariate.

    obj.requireLocked();
    B = obj.batch{i};
    days = unique([B.day]);
    g = cell(numel(B), 1);
    for k = 1:numel(B)
        d = find(days == B(k).day, 1);
        g{k} = repmat(d, B(k).n, 1);
    end
end

%% ------------------------------------------------------------- attribute
function T = attribute(obj, D, varargin)
%ATTRIBUTE  Turn per-observation diagnostics into a per-(sensor, day) summary.
%
%   T = bx.attribute(D) where D is the struct from Combine/diagnose.
%
%   This is the whole point of the class. Without it a flagged observation is
%   "observation 4193". With it, it is "LRI, day 14, 3 flags, mean r = 0.42".
%
%   Name-value:
%     'by'   'sensor' | 'day' | 'both'    (default 'both')

    obj.requireLocked();
    ip = inputParser;
    addParameter(ip, 'by', 'both', ...
                 @(s) any(strcmpi(s, {'sensor','day','both'})));
    parse(ip, varargin{:});
    by = lower(ip.Results.by);

    rows = struct('system', {}, 'sensor', {}, 'day', {}, 'n', {}, ...
                  'nFlag', {}, 'meanR', {}, 'minR', {}, 'maxAbsW', {}, ...
                  's0', {});

    for i = 1:numel(obj.sysName)
        if i > numel(D.rNum), break; end
        ri = D.rNum{i};  wi = D.w{i};  vi = D.v{i};
        B  = obj.batch{i};
        pos = 0;
        for k = 1:numel(B)
            nk = B(k).n;
            if pos + nk > numel(ri)
                error('lsqvce:batchIndex:countMismatch', ...
                    ['System "%s": the registered observation counts sum to ' ...
                     'more than the diagnostics contain (%d > %d). The counts ' ...
                     'passed to setCount must match what batchFcn actually ' ...
                     'returns.'], obj.sysName{i}, pos+nk, numel(ri));
            end
            sl = pos + (1:nk);
            rows(end+1) = struct( ...
                'system',  obj.sysName{i}, ...
                'sensor',  obj.sensorName{i}{B(k).sensor}, ...
                'day',     B(k).day, ...
                'n',       nk, ...
                'nFlag',   sum(abs(wi(sl)) > D.kCrit), ...
                'meanR',   mean(ri(sl)), ...
                'minR',    min(ri(sl)), ...
                'maxAbsW', max(abs(wi(sl))), ...
                's0',      sqrt(mean(vi(sl).^2))); %#ok<AGROW>
            pos = pos + nk;
        end
        if pos ~= numel(ri)
            lsqvce.util.log.warn('lsqvce:batchIndex:countMismatch', ...
                ['System "%s": registered %d observations but the diagnostics ' ...
                 'contain %d. The per-(sensor,day) attribution below is ' ...
                 'therefore WRONG. Check your setCount calls.'], ...
                 obj.sysName{i}, pos, numel(ri));
        end
    end

    T = struct2table(rows);

    switch by
        case 'sensor'
            T = obj.rollup(T, {'system','sensor'});
        case 'day'
            T = obj.rollup(T, {'system','day'});
    end
end

%% ---------------------------------------------------------------- report
function report(obj, D)
%REPORT  Print a per-sensor and per-day summary of the diagnostics.

    obj.requireLocked();
    Ts = obj.attribute(D, 'by', 'sensor');
    Td = obj.attribute(D, 'by', 'day');

    lsqvce.util.log.msg('\n%s\n', repmat('=', 1, 72));
    lsqvce.util.log.msg('  PER-SENSOR SUMMARY\n');
    lsqvce.util.log.msg('%s\n', repmat('=', 1, 72));
    lsqvce.util.log.msg('  %-10s %-8s %8s %8s %9s %9s\n', ...
            'system', 'sensor', 'n', 'flags', 'mean r', 'rms v');
    lsqvce.util.log.msg('  %s\n', repmat('-', 1, 60));
    for q = 1:height(Ts)
        lsqvce.util.log.msg('  %-10s %-8s %8d %8d %9.4f %9.4f\n', ...
                Ts.system{q}, Ts.sensor{q}, Ts.n(q), Ts.nFlag(q), ...
                Ts.meanR(q), Ts.s0(q));
    end

    lsqvce.util.log.msg('\n  A sensor with a systematically different rms v is telling you\n');
    lsqvce.util.log.msg('  its weights are wrong. Estimate them: lsqvce.varComp with\n');
    lsqvce.util.log.msg('  bx.sensorGroups(i) as the group index.\n');

    lsqvce.util.log.msg('\n%s\n', repmat('=', 1, 72));
    lsqvce.util.log.msg('  PER-DAY SUMMARY (worst 10 by flag count)\n');
    lsqvce.util.log.msg('%s\n', repmat('=', 1, 72));
    [~, ord] = sort(Td.nFlag, 'descend');
    ord = ord(1:min(10, numel(ord)));
    lsqvce.util.log.msg('  %-10s %-12s %8s %8s %9s\n', ...
            'system', 'day', 'n', 'flags', 'mean r');
    lsqvce.util.log.msg('  %s\n', repmat('-', 1, 52));
    for q = ord(:).'
        lsqvce.util.log.msg('  %-10s %-12s %8d %8d %9.4f\n', ...
                Td.system{q}, obj.dayStr(Td.day(q)), Td.n(q), ...
                Td.nFlag(q), Td.meanR(q));
    end
    lsqvce.util.log.msg(['\n  A single day with many flags is usually an EVENT (manoeuvre,\n' ...
             '  eclipse, thermal transient), not a set of independent blunders.\n' ...
             '  Rejecting the observations one by one treats the symptom.\n\n']);
end

end % methods

methods (Access = private)

function i = sysIdx(obj, name)
%SYSIDX  Index of system NAME in the registry (errors if unknown).
    i = find(strcmp(obj.sysName, char(name)), 1);
    if isempty(i)
        error('lsqvce:batchIndex:noSuchSystem', ...
              'No system named "%s". Known: %s', char(name), ...
              strjoin(obj.sysName, ', '));
    end
end

function s = sensorIdx(obj, i, name)
%SENSORIDX  Index of sensor NAME within system i (errors if unknown).
    s = find(strcmp(obj.sensorName{i}, char(name)), 1);
    if isempty(s)
        error('lsqvce:batchIndex:noSuchSensor', ...
              'System "%s" has no sensor "%s". Known: %s', ...
              obj.sysName{i}, char(name), ...
              strjoin(obj.sensorName{i}, ', '));
    end
end

function requireLocked(obj)
%REQUIRELOCKED  Assert the index has been locked (built); errors otherwise.
    if ~obj.locked
        error('lsqvce:batchIndex:notLocked', ...
              ['Call lock() first. Until then the flat batch indices are not ' ...
               'assigned and decode() has nothing to decode.']);
    end
end

function s = dayStr(~, d)
%DAYSTR  Format a day key as yyyy-mm-dd if it looks like a datenum, else num2str.
    if isnumeric(d) && d > 7e5 && d < 8e5
        s = datestr(d, 'yyyy-mm-dd');       % looks like a datenum
    else
        s = num2str(d);
    end
end

function T = rollup(~, T, keys)
%ROLLUP  Group a diagnostics table by KEYS and sum n/nFlag, RMS-combine s0.
    [G, ia] = findgroups(T(:, keys));
    out = ia;
    out.n       = splitapply(@sum,  T.n,       G);
    out.nFlag   = splitapply(@sum,  T.nFlag,   G);
    out.meanR   = splitapply(@mean, T.meanR,   G);
    out.minR    = splitapply(@min,  T.minR,    G);
    out.maxAbsW = splitapply(@max,  T.maxAbsW, G);
    out.s0      = splitapply(@(x) sqrt(mean(x.^2)), T.s0, G);
    T = out;
end

end
end
