function demo_sensors()
%DEMO_SENSORS  The full hierarchy: satellite system / sensor / day.
%
%   Two satellite systems, each carrying several sensors, each observing over
%   several days. The common gravity field is shared; each system has its own
%   PER-SYSTEM local parameters (a scale factor, say) that all its sensors see.
%
%   THE KEY STRUCTURAL FACT
%   -----------------------
%   The v5 data struct states the whole hierarchy declaratively:
%     * a SYSTEM is an entry of data.systems (its localParams are the
%       per-system locals -- every sensor of the system sees them);
%     * a SENSOR is an observation type: one field of covariance, and one
%       (day, sensor) call of your buildDay;
%     * a DAY is an entry of data.days; a sensor that was down on a day is a
%       DATA GAP -- the arc resolver returns no arcs, nothing else changes.
%   No batching logic is yours to write: fit streams (system, sensor, day)
%   in whatever order it likes, and the solution is invariant to that order.
%
%   What the hierarchy DOES buy you is attribution: "LRI, day 14, 3 flags"
%   instead of "observation 4193" -- R.residuals() carries system, sensor and
%   day on every record.
%
%   INPUT   none
%   OUTPUT  none (prints the walk-through; ~seconds)
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.
%   Updated 2026-09-07 (CEST, UTC+02:00), R6 WP-95: the declaration is named
%   the v5 data struct (help block and the STEP 2 banner).
%   With assistance from Claude (Anthropic), model Claude Opus 5.
%   Updated 2026-08-31 (CEST, UTC+02:00), v4 SCOPED tables (Level L11,
%   WP-77): the declaration tables carry the scope column (commonParams
%   'common'|'day', localParams 'arc'|'day'|'system'); localsPerArc left the
%   USER surface. With assistance from Claude (Anthropic), model
%   Claude Fable 5.
%   Updated 2026-08-31 (CEST, UTC+02:00), v3 data contract (Level L10,
%   WP-74; v5 R5 M5): USER surface only -- the v1 combine/batch machinery is gone; the
%   multi-system fit is lsqvce.fit on the one accumulator, replaced
%   by the declarative hierarchy (systems x covariance fields x days), one
%   lsqvce.fit per question (gls for the symptom, Estimate="noise", Pool="stream", ArcFactors=false for the
%   weights), R.residuals() for the attribution. With assistance from
%   Claude (Anthropic), model Claude Fable 5.

clc;
rs = RandStream('mt19937ar', 'Seed', 2026);

%% ------------------------------------------------------------------------
% THE LAYOUT
% -------------------------------------------------------------------------
uC   = 6;                                   % common: gravity coefficients
sysN = {'GRACE-C', 'GRACE-D'};
sysId = {'C', 'D'};
sens = { {'KBR','LRI','GPS'}, {'KBR','LRI'} };
uL   = [4, 3];                              % PER-SYSTEM locals (all sensors see them)
days = 1:10;
nk   = 60;                                  % observations per (sensor, day)

% The TRUE sensor noise. We will CLAIM sigma = 1 for all of them and let the
% variance components find out otherwise -- which is the normal state of
% affairs when combining instruments whose headers disagree with reality.
sigTrue = containers.Map( ...
    {'KBR','LRI','GPS'}, ...
    {  1.0,  0.4,  2.5});     % LRI is much quieter, GPS much noisier

fprintf('LAYOUT\n');
for i = 1:numel(sysN)
    fprintf('  %-9s  sensors: %-18s  locals: %d\n', ...
            sysN{i}, strjoin(sens{i}, ', '), uL(i));
end
fprintf('  %d days, %d observations per (sensor, day)\n', numel(days), nk);

%% ------------------------------------------------------------------------
% STEP 1 -- generate the data (in real work, buildDay reads from disk)
% -------------------------------------------------------------------------
xC_true = randn(rs, uC, 1);
xL_true = cell(numel(sysN), 1);
store   = cell(numel(sysN), 1);             % store{i} maps 'sensor|day' -> batch
for i = 1:numel(sysN)
    xL_true{i} = randn(rs, uL(i), 1);
    store{i} = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for s = 1:numel(sens{i})
        for d = days
            % simulate an outage: GRACE-D's LRI is down on days 4-5
            if i == 2 && strcmp(sens{i}{s}, 'LRI') && ismember(d, [4 5])
                continue;
            end
            sg = sigTrue(sens{i}{s}) * ones(nk, 1);
            AC = randn(rs, nk, uC);
            AL = randn(rs, nk, uL(i));
            l  = AC*xC_true + AL*xL_true{i} + sg .* randn(rs, nk, 1);
            store{i}(sprintf('%s|%d', sens{i}{s}, d)) = ...
                struct('AC', AC, 'AL', AL, 'l', l);
        end
    end
end

% Plant a bad DAY: a thermal event corrupts GRACE-C's KBR on day 7.
b = store{1}('KBR|7');
b.l = b.l + 6.0;                            % a systematic offset
store{1}('KBR|7') = b;
fprintf('\n  [planted] a 6-sigma offset in GRACE-C / KBR / day 7\n');
fprintf('  [planted] GRACE-D''s LRI is down on days 4-5 (a data gap)\n');

%% ------------------------------------------------------------------------
% STEP 2 -- declare the hierarchy: the v5 data struct
% -------------------------------------------------------------------------
data = struct();
data.commonParams = [repmat({'common'}, uC, 1), num2cell((1:uC).')];  % scoped
data.days = days;
for i = 1:numel(sysN)
    s = struct();
    s.buildDay = makeBuildDay(store{i}, nk);
    for t = sens{i}                          % one covariance field per SENSOR;
        s.covariance.(t{1}) = struct('type', 'diag', 'sigma', 1);  % the CLAIM
    end
    s.localParams = [repmat({'system'}, uL(i), 1), num2cell((1:uL(i)).')];
    s.arcs = makeSchedule(store{i});         % ragged: [] on an outage day
    s.id = sysId{i};
    s.name = sysN{i};
    data.systems(i) = s;
end

%% ------------------------------------------------------------------------
% STEP 3 -- the symptom: one GLS solve under the claimed weights
% -------------------------------------------------------------------------
R = lsqvce.fit(data, lsqvce.Options(Verbose=false));
fprintf('\nSYMPTOM\n');
fprintf('  sigma0 = %.4f   (we CLAIMED sigma0 = 1)\n', R.sigma0);
fprintf('  -> far from 1, so at least one sensor is mis-weighted.\n');

%% ------------------------------------------------------------------------
% STEP 4 -- PER-SENSOR VARIANCE COMPONENTS
% -------------------------------------------------------------------------
% Both satellites' KBR are the same INSTRUMENT TYPE, so they may share one
% estimated spectral SHAPE: Pool="obsType" pools the streams of one
% observation type across systems (v5). Every stream keeps its OWN level
% regardless -- that is the variance component. (Pool="stream" gives every
% (system, sensor) pair its own shape as well -- a modelling choice the
% toolbox will not make for you.)
Rv = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="obsType", ArcFactors=false, Verbose=false));

sig = Rv.history.sigma(end, :);
fprintf('\nPER-SENSOR COMPONENTS (one level per stream)\n');
fprintf('  %-8s %12s %12s %10s\n', 'sensor', 'recovered', 'truth', 'error');
fprintf('  %s\n', repmat('-', 1, 46));
allSensors = unique([sens{:}], 'stable');
ref = [];
for gg = 1:numel(allSensors)
    k = find(contains(Rv.history.streams, allSensors{gg}), 1);
    if isempty(ref), ref = sig(k) / sigTrue(allSensors{1}); end
    got = sig(k) / ref;                      % normalise on the first sensor
    tru = sigTrue(allSensors{gg});
    fprintf('  %-8s %12.4f %12.4f %9.1f%%\n', ...
            allSensors{gg}, got, tru, 100*(got-tru)/tru);
end
fprintf(['\n  -> the components recovered the sensor noise ratios from the\n' ...
         '     data, rescaled the covariances and re-solved (sigma0 = %.3f).\n'], ...
        Rv.sigma0);
fprintf(['  RESOLUTION: with a few hundred redundancies per group you cannot\n' ...
         '  resolve a sigma ratio better than a few percent. A 5%% "correction"\n' ...
         '  is NOISE -- read the report before believing a small one.\n']);

%% ------------------------------------------------------------------------
% STEP 5 -- HUNT THE BAD DAY. This is what the hierarchy buys you.
% -------------------------------------------------------------------------
% A single day with many flags is usually an EVENT (manoeuvre, eclipse,
% thermal transient), not a set of independent blunders. Rejecting the
% observations one by one treats the symptom. R.residuals() carries system,
% sensor and day on every record, so the attribution is one loop.
fprintf('\n%s\n', repmat('=', 1, 72));
fprintf('  IS THERE A BAD DAY?\n');
fprintf('%s\n', repmat('=', 1, 72));

res = R.residuals();                        % verbosity inherited from the run
iC  = find(cellfun(@(x) isequal(x, 'C'), {res.sys.sysID}), 1);
if isempty(iC), iC = 1; end
obs = res.sys(iC).obs;
obs = obs(strcmp({obs.obsType}, 'KBR'));
fprintf('  GRACE-C / KBR, by day:\n');
fprintf('  %6s %10s\n', 'day', 'rms v');
fprintf('  %s\n', repmat('-', 1, 20));
rmsAll = zeros(1, numel(days));
for d = days
    v = vertcat(obs([obs.day] == d).v);
    rmsAll(d) = sqrt(mean(v.^2));
end
for d = days
    mark = '';
    if rmsAll(d) > 3 * median(rmsAll), mark = '   <<< EVENT'; end
    fprintf('  %6d %10.4f%s\n', d, rmsAll(d), mark);
end
fprintf(['\n  -> Day 7 stands out. That is not a set of independent blunders,\n' ...
         '     it is an EVENT. Investigate the cause; do not snoop it away.\n\n']);

end

%% ========================================================================
function f = makeBuildDay(st, nk)
%MAKEBUILDDAY  The system's buildDay over the pre-generated store.
%   Deterministic by construction: it only READS the store.
    f = @buildDay;
    function D = buildDay(day, obsType)
        b = st(sprintf('%s|%d', obsType, day));
        D.A = {b.AC};                        % v6: one cell = one arc (nk x uC)
        D.B = {b.AL};                        % nk x uL -- the per-system locals
        D.l = {b.l};
        D.t = {day*86400 + (0:nk-1).' * 5};  % seconds
    end
end

function f = makeSchedule(st)
%MAKESCHEDULE  Ragged arc resolver: one arc per (day, sensor) with data,
%   NO arcs on an outage day -- the documented data-gap mechanism.
    f = @(day, obsType) resolve(st, day, obsType);
    function ids = resolve(st, day, obsType)
        if st.isKey(sprintf('%s|%d', obsType, day))
            ids = 1;
        else
            ids = zeros(1, 0);
        end
    end
end
