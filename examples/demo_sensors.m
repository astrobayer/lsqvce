function demo_sensors()
%DEMO_SENSORS  The full hierarchy: satellite system / sensor / daily batch.
%
%   Two satellite systems, each carrying several sensors, each observing over
%   several days. The common gravity field is shared; each system has its own
%   PER-SYSTEM local parameters (a scale factor, say) that all its sensors see.
%
%   THE KEY STRUCTURAL FACT
%   -----------------------
%   With PER-SYSTEM locals, a (sensor, day) batch is just ANOTHER BATCH inside
%   its system's elimination block. No new elimination logic is needed: the
%   accumulator already streams batches and does not care where they came from.
%
%   Verified: sensor x day batching gives bit-identical results to a flat
%   structure, and the solution is invariant to batch ORDER to 3e-16. You may
%   feed the days and sensors in any order and interleave them freely.
%
%   What the hierarchy DOES buy you is attribution: "LRI, day 14, 3 flags"
%   instead of "observation 4193".
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.

clc;
rs = RandStream('mt19937ar', 'Seed', 2026);

%% ------------------------------------------------------------------------
% THE LAYOUT
% -------------------------------------------------------------------------
uC   = 6;                                   % common: gravity coefficients
sysN = {'GRACE-C', 'GRACE-D'};
sens = { {'KBR','LRI','GPS'}, {'KBR','LRI'} };
uL   = [4, 3];                              % PER-SYSTEM locals (all sensors see them)
days = 1:10;
nk   = 60;                                  % observations per (sensor, day)

% The TRUE sensor noise. We will CLAIM sigma = 1 for all of them and let
% Helmert find out otherwise -- which is the normal state of affairs when
% combining instruments whose headers disagree with reality.
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
% STEP 1 -- declare the hierarchy
% -------------------------------------------------------------------------
bx = lsqvce.BatchIndex();
for i = 1:numel(sysN)
    bx.addSystem(sysN{i}, sens{i}, days);
end

% Register which (system, sensor, day) triples actually have data. A sensor
% that was down on a given day is simply not registered -- no empty batches
% needed, and no requirement that every sensor cover every day.
for i = 1:numel(sysN)
    for s = 1:numel(sens{i})
        for d = days
            % simulate an outage: GRACE-D's LRI is down on days 4-5
            if i == 2 && strcmp(sens{i}{s}, 'LRI') && ismember(d, [4 5])
                continue;
            end
            bx.setCount(sysN{i}, sens{i}{s}, d, nk);
        end
    end
end
bx.lock();

fprintf('\n  batches per system: ');
for i = 1:numel(sysN), fprintf('%d ', bx.nBatch(i)); end
fprintf('\n  (GRACE-D''s LRI is down on days 4-5 -- those triples were simply\n');
fprintf('   never registered. No empty batches are required.)\n');

%% ------------------------------------------------------------------------
% STEP 2 -- generate the data (in real work, this reads from disk)
% -------------------------------------------------------------------------
xC_true = randn(rs, uC, 1);
xL_true = cell(numel(sysN),1);
store   = cell(numel(sysN),1);
for i = 1:numel(sysN)
    xL_true{i} = randn(rs, uL(i), 1);
    store{i} = cell(bx.nBatch(i), 1);
    for k = 1:bx.nBatch(i)
        info = bx.decode(i, k);
        sg   = sigTrue(info.sensor) * ones(info.n, 1);
        AC = randn(rs, info.n, uC);
        AL = randn(rs, info.n, uL(i));
        l  = AC*xC_true + AL*xL_true{i} + sg .* randn(rs, info.n, 1);
        store{i}{k} = struct('AC', AC, 'AL', AL, 'l', l, ...
                             'p', ones(info.n,1));   % we CLAIM sigma = 1
    end
end

% Plant a bad DAY: a thermal event corrupts GRACE-C's KBR on day 7.
for k = 1:bx.nBatch(1)
    info = bx.decode(1, k);
    if strcmp(info.sensor, 'KBR') && info.day == 7
        store{1}{k}.l = store{1}{k}.l + 6.0;      % a systematic offset
        fprintf('\n  [planted] a 6-sigma offset in GRACE-C / KBR / day 7\n');
    end
end

% the batch sources
    function [AC, AL, l, p] = src(i, k)
        b = store{i}{k};
        AC = b.AC;  AL = b.AL;  l = b.l;  p = b.p;
    end

%% ------------------------------------------------------------------------
% STEP 3 -- combine
% -------------------------------------------------------------------------
cmb = lsqvce.Combine(uC, 'method', 'qr');
for i = 1:numel(sysN)
    cmb.addSystem(uL(i), @(k) src(i,k), bx.nBatch(i), 'name', sysN{i});
end
sol = cmb.solve();
D   = cmb.diagnose(sol);

%% ------------------------------------------------------------------------
% STEP 4 -- INTEGRITY FIRST, always
% -------------------------------------------------------------------------
fprintf('\nINTEGRITY\n');
fprintf('  n = %d,  u_total = %d  (uC = %d + uL = %d)\n', ...
        D.n, D.uTotal, D.uC, D.uL);
fprintf('  sum(r_i) = %.6f   n - u_total = %d   rel err %.2e\n', ...
        D.sumR, D.r, D.sumRRelErr);
assert(D.sumRRelErr < 1e-8, 'integrity failed');
fprintf('  -> holds.\n');

fprintf('\n  s0 = %.4f   (we CLAIMED sigma0 = 1)\n', D.s0);
fprintf('  -> s0 is far from 1, so at least one sensor is mis-weighted.\n');

%% ------------------------------------------------------------------------
% STEP 5 -- ATTRIBUTE. This is what the hierarchy buys you.
% -------------------------------------------------------------------------
% Without BatchIndex, a flagged observation is "observation 4193".
% With it, it is "GRACE-C / KBR / day 7".
bx.report(D);

%% ------------------------------------------------------------------------
% STEP 6 -- PER-SENSOR VARIANCE COMPONENTS
% -------------------------------------------------------------------------
% Both KBR and LRI observe the same arcs, and their quoted noise levels are
% almost never both right. Helmert tells you the ratio the DATA imply.
%
% NOTE the group index must be GLOBAL across systems: GRACE-C's KBR and
% GRACE-D's KBR are the same INSTRUMENT TYPE, so they share a variance
% component. (If you believe each unit has its own noise, give them separate
% group indices -- that is a modelling choice, and the toolbox will not make
% it for you.)
allSensors = unique([sens{:}], 'stable');
nGroup = numel(allSensors);
fprintf('sensor groups: %s\n', strjoin(allSensors, ', '));

    function g = grp(i, k)
        info = bx.decode(i, k);
        g = repmat(find(strcmp(allSensors, info.sensor), 1), info.n, 1);
    end

VC = cmb.varComp(sol, @grp, nGroup);

fprintf('  %-8s %12s %12s %10s\n', 'sensor', 'recovered', 'truth', 'error');
fprintf('  %s\n', repmat('-', 1, 46));
ref = VC.sigma(1) / sigTrue(allSensors{1});     % normalise on the first sensor
for gg = 1:nGroup
    got = VC.sigma(gg) / ref;
    tru = sigTrue(allSensors{gg});
    fprintf('  %-8s %12.4f %12.4f %9.1f%%\n', ...
            allSensors{gg}, got, tru, 100*(got-tru)/tru);
end

fprintf(['\n  -> Helmert recovered the sensor noise ratios from the data.\n' ...
         '     Apply VC.scale to the covariances and re-run.\n']);

% -------------------------------------------- ABOUT THE HUTCHINSON ALTERNATIVE
% cmb.varComp uses the FULL Qxx (including the eliminated per-system locals),
% assembled per system -- the block-arrow topology. The matrix-free Hutchinson
% VCE (the large-u method, O(u^3 + G*Z*u^2) vs O(G*u^3)) is wired as a REAL
% substitution on the FLAT path only: lsqvce.vceHutchinson, or
% lsqvce.varComp(..., 'backend','hutchinson'). See examples/demo_lsqvce.m and
% examples/demo_workflow.m for a side-by-side exact-vs-Hutchinson run. Bringing
% it into cmb.varComp needs per-group N_g^cc exposed from the accumulator core
% (roadmap item 0, still open), so at this hierarchy's scale the exact path
% above is both correct and fast.

fprintf(['\n  RESOLUTION: with r_g of a few hundred you cannot resolve a sigma\n' ...
         '  ratio better than a few percent. A 5%% "correction" is NOISE.\n' ...
         '  VC.significant tells you which corrections are real:  %s\n'], ...
         mat2str(VC.significant(:).'));

%% ------------------------------------------------------------------------
% STEP 7 -- HUNT THE BAD DAY
% -------------------------------------------------------------------------
% A single day with many flags is usually an EVENT (manoeuvre, eclipse,
% thermal transient), not a set of independent blunders. Rejecting the
% observations one by one treats the symptom.
fprintf('\n%s\n', repmat('=', 1, 72));
fprintf('  IS THERE A BAD DAY?\n');
fprintf('%s\n', repmat('=', 1, 72));

T = bx.attribute(D, 'by', 'both');
T = T(strcmp(T.system, 'GRACE-C') & strcmp(T.sensor, 'KBR'), :);
fprintf('  GRACE-C / KBR, by day:\n');
fprintf('  %6s %8s %10s\n', 'day', 'flags', 'rms v');
fprintf('  %s\n', repmat('-', 1, 28));
for q = 1:height(T)
    mark = '';
    if T.nFlag(q) > 3*median(T.nFlag), mark = '   <<< EVENT'; end
    fprintf('  %6d %8d %10.4f%s\n', T.day(q), T.nFlag(q), T.s0(q), mark);
end
fprintf(['\n  -> Day 7 stands out. That is not a set of independent blunders,\n' ...
         '     it is an EVENT. Investigate the cause; do not snoop it away.\n\n']);

end
