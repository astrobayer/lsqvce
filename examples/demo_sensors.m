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
%   What the hierarchy DOES buy you is attribution: "GRACE-C, KBR, day 7"
%   instead of "observation 4193" -- every arc is one (system, sensor, day),
%   and so is its variance factor.
%
%   WHAT IT SHOWS, in order: the layout; the symptom (sigma0 far from 1
%   under the claimed weights); the bad day, found FIRST -- one variance
%   factor per sensor-day, read against its stream's median; then the
%   per-sensor variance components with that day left out, each stream
%   against its planted truth. The order is the lesson: a component is an
%   average over its stream, so an event inside it corrupts the level.
%
%   ACCURACY  (MEASURED 2026-10-16, R2026a, v6.25.2; the data are seeded)
%     The planted event (a 6-sigma offset on GRACE-C / KBR / day 7) comes
%     back at 15.9 x its stream's median factor; every other sensor-day lies
%     between 0.5 and 1.8. Without it, the levels are C|GPS 2.536, C|KBR
%     0.987, C|LRI 0.400, D|KBR 0.988, D|LRI 0.416 against the planted 2.5 /
%     1.0 / 0.4 / 1.0 / 0.4 -- within 4.1 %, sigma0 1.000, converged in 7
%     passes (the factor run in about 14 of its 20). WP-145 measured what the
%     earlier order cost: the levels estimated WITH the event in the data
%     came back C|KBR +46 % and D|KBR +43 % (the shared KBR shape carries
%     it), and the table, normalised on C|KBR, showed LRI and GPS 32 % low
%     while it said they were recovered. Noise.Grid 8 because an arc holds 60
%     samples: the default 60-bin law needed 168 passes on this white noise.
%
%   INPUT   none
%   OUTPUT  none (prints the walk-through; ~3 s)
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.
%   Updated 2026-10-16 (CEST, UTC+02:00), WP-145: the bad day is found
%   (per-arc factors) BEFORE the variance components, which are estimated
%   without it and shown against the truth stream by stream; Noise.Grid 8;
%   both runs converge. With assistance from Claude (Anthropic), model
%   Claude Opus 5.5.
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
% STEP 4 -- HUNT THE BAD DAY FIRST. This is what the hierarchy buys you.
% -------------------------------------------------------------------------
% A variance component is an AVERAGE over its stream: one corrupted day inside
% it is averaged in, and the level comes back wrong for every sensor read
% against it. (WP-145, measured: run on the data below with the day-7 event
% still in it, C|KBR came back 46 % high, D|KBR -- which shares the KBR shape
% -- 43 %, and a table normalised on C|KBR showed LRI and GPS 32 % low.) So the
% event is found BEFORE the levels are estimated. ArcFactors=true gives every
% arc -- here one (system, sensor, day) -- its own variance factor; read
% against the MEDIAN of its own stream (the factors are absolute against the
% claimed sigma 1, so a GPS arc sits at ~6 legitimately), a clean arc is near
% 1 and an event stands out. The run needs more passes than the toolbox's 10
% -- one factor per arc is many components -- so it gets a budget of 20 (it
% settles in about 14). Noise.Grid 8: an arc holds only nk = 60 samples, too
% few to resolve the default 60 spectral bins (bins the data barely determine
% creep toward zero and the loop needed 168 passes).
fprintf('\n%s\n', repmat('=', 1, 72));
fprintf('  IS THERE A BAD DAY?  (one variance factor per sensor-day)\n');
fprintf('%s\n', repmat('=', 1, 72));
Rh = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="obsType", ArcFactors=true, ...
        MaxPasses=20, Noise=struct('Grid', 8), Verbose=false));
factor = containers.Map(Rh.noise.arcKeys, num2cell(Rh.noise.sigma2(:)));
flagged = {};
fprintf('  %-9s %-6s %s\n', 'system', 'sensor', 'factor / stream median, days 1..10');
for i = 1:numel(sysN)
    for s = 1:numel(sens{i})
        rel = nan(1, numel(days));
        for d = days
            k = lsqvce.util.arcKey(sysId{i}, sens{i}{s}, d, 1);   % the toolbox's one arc key
            if isKey(factor, k), rel(d) = factor(k); end
        end
        rel = rel / median(rel, 'omitnan');
        fprintf('  %-9s %-6s %s\n', sysN{i}, sens{i}{s}, sprintf('%5.1f', rel));
        % a clean arc's factor scatters by about sqrt(2/60) = 0.18 around 1;
        % 5 is twenty such widths away
        for d = find(rel > 5)
            flagged(end+1, :) = {i, sens{i}{s}, d, rel(d)}; %#ok<AGROW>
        end
    end
end
for e = 1:size(flagged, 1)
    fprintf('  <<< EVENT: %s / %s / day %d, %.1f x its stream''s median\n', ...
            sysN{flagged{e,1}}, flagged{e,2}, flagged{e,3}, flagged{e,4});
end
fprintf(['  -> One sensor-day stands out. That is not a set of independent\n' ...
         '     blunders, it is an EVENT. Investigate the cause; do not snoop it\n' ...
         '     away -- and keep it OUT of the noise estimate below.\n']);

%% ------------------------------------------------------------------------
% STEP 5 -- PER-SENSOR VARIANCE COMPONENTS, without the event
% -------------------------------------------------------------------------
% The event's arc leaves the declaration the way an outage does (STEP 2): the
% arc resolver returns no arc for it, nothing else changes. Both satellites'
% KBR are the same INSTRUMENT TYPE, so they may share one estimated spectral
% SHAPE: Pool="obsType" pools the streams of one observation type across
% systems (v5). Every stream keeps its OWN level regardless -- that is the
% variance component. (Pool="stream" gives every (system, sensor) pair its own
% shape as well -- a modelling choice the toolbox will not make for you.)
for e = 1:size(flagged, 1)
    remove(store{flagged{e,1}}, sprintf('%s|%d', flagged{e,2}, flagged{e,3}));
end
Rv = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="obsType", ArcFactors=false, ...
        Noise=struct('Grid', 8), Verbose=false));

T = Rv.summary();
fprintf('\nPER-SENSOR COMPONENTS (one level per stream, the event left out)\n');
fprintf('  %-8s %12s %12s %10s\n', 'stream', 'recovered', 'truth', 'error');
fprintf('  %s\n', repmat('-', 1, 46));
for k = 1:height(T)
    tru = sigTrue(char(extractAfter(T.stream(k), '|')));
    fprintf('  %-8s %12.4f %12.4f %9.1f%%\n', T.stream(k), T.sigma(k), tru, ...
            100*(T.sigma(k) - tru)/tru);
end
fprintf(['\n  -> the components recovered every sensor''s noise from the data,\n' ...
         '     rescaled the covariances and re-solved (sigma0 = %.3f, %d passes,\n' ...
         '     converged: %d).\n'], Rv.sigma0, Rv.history.nIter, Rv.converged);
fprintf(['  RESOLUTION: with a few hundred redundancies per group you cannot\n' ...
         '  resolve a sigma ratio better than a few percent. A 5%% "correction"\n' ...
         '  is NOISE -- read the report before believing a small one.\n\n']);

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
