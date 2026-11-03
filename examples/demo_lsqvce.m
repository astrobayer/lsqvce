%% DEMO_LSQVCE  Worked examples: how to judge a least-squares adjustment.
%
%   Six scenarios, each isolating one failure mode, so you can see what the
%   report says when a specific thing is wrong -- and, just as importantly,
%   what it says when things are fine. Every run is the USER path, an in-memory
%   lsqvce.TaskDefinition and the one verb:
%
%       R = lsqvce.fit(lsqvce.TaskDefinition.fromMatrices(A, l, Weights=p));
%
%   with the diagnostics on R.diag and the verdict from R.report(). Case 5 also
%   shows the CURE: a streaming TaskDefinition with one stream per instrument
%   group, each with its own Build, and Estimate="noise", Pool="stream",
%   ArcFactors=false, so fit estimates the weights it needs.
%
%   A script: no input; run it whole or section by section. It leaves its
%   Results (R, R2, ..., R6, R5v) in the workspace and draws one figure (case 3).
%
%   ACCURACY  Not pinned: six illustrations. Their numbers come from seeded
%             draws (mt19937ar, seed 2026) and are printed, not asserted; the
%             gate runs it as a smoke test (test_examples).
%   ERRORS    none of its own; a declaration error is raised by the add-call
%             that makes it (lsqvce:TaskDefinition:*), a mistake inside a Build
%             while the day is read (lsqvce:DayContext:*).
%
%   See also LSQVCE.TASKDEFINITION, LSQVCE.FIT, LSQVCE.RESULTS, LSQVCE.REPORT,
%   DEMO_WORKFLOW, DEMO_TUTORIAL
%
%   Updated 2026-11-02 (CET, UTC+01:00), WP-161: declared as lsqvce.TaskDefinition, one Build per group. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-09-07 (CEST, UTC+02:00), R6 WP-95: the verdict is R.report()'s,
%   not R.source's, and the per-lag redundancies are Hutchinson-probed -- the
%   joint estimator has no exact route.
%   With assistance from Claude (Anthropic), model Claude Opus 5.
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.
%   Updated 2026-08-31 (CEST, UTC+02:00), v4 SCOPED tables (Level L11,
%   WP-77): the declaration tables carry the scope column (commonParams
%   'common'|'day', localParams 'arc'|'day'|'system'); localsPerArc left the
%   USER surface. With assistance from Claude (Anthropic), model
%   Claude Fable 5.
%   Updated 2026-08-31 (CEST, UTC+02:00), v3 data contract (Level L10,
%   WP-74; v5 R5 M5): USER surface only -- the flat drivers of v1 are gone; every
%   kind runs on the one accumulator through lsqvce.fit --
%   lsqvce.fit on the in-memory kind and, for the wrong-weights cure, by
%   Estimate="noise", Pool="stream", ArcFactors=false on a two-type streaming declaration. With assistance from
%   Claude (Anthropic), model Claude Fable 5.

clear; clc;
rs = RandStream('mt19937ar', 'Seed', 2026);

%% ------------------------------------------------------------------ setup
n = 2000;  u = 8;
A  = randn(rs, n, u);
xt = randn(rs, u, 1);
sg = 0.5 + rand(rs, n, 1);
p  = 1 ./ sg.^2;


%% ===================================================================== (1)
% A CLEAN adjustment. This is the baseline: everything passes, and the report
% tells you the one thing that actually matters -- that no undetectable blunder
% can move a parameter by more than a fraction of its sigma.
fprintf('\n\n#################### 1. CLEAN ADJUSTMENT ####################\n');
l = A * xt + sg .* randn(rs, n, 1);
R = lsqvce.fit(lsqvce.TaskDefinition.fromMatrices(A, l, Weights=p));
D = R.diag;


%% ===================================================================== (2)
% OUTLIERS. Three gross errors. Watch the variance factor inflate, and watch
% data snooping localise them -- with Bonferroni, so the list is short and real
% rather than long and mostly noise.
fprintf('\n\n#################### 2. OUTLIERS ####################\n');
l2 = l;
l2([137 892 1455]) = l2([137 892 1455]) + [12; -15; 9] .* sg([137 892 1455]);
R2 = lsqvce.fit(lsqvce.TaskDefinition.fromMatrices(A, l2, Weights=p));
D2 = R2.diag;

fprintf('Planted at 137, 892, 1455.  Flagged: %s\n', mat2str(D2.flagIdx(:).'));
fprintf(['\nNOTE the procedure: reject the LARGEST |w| only, re-solve, ' ...
         're-snoop.\nA single gross error smears its residual over its ' ...
         'neighbours and can\nmake innocent observations look guilty.\n']);


%% ===================================================================== (3)
% A SYSTEMATIC ERROR. This is the case that the chi-square test and data
% snooping BOTH handle badly, and it is the most common real failure.
% A smooth drift produces no single large residual, so snooping finds nothing,
% yet the parameters are corrupted.
fprintf('\n\n#################### 3. SYSTEMATIC ERROR ####################\n');
t  = linspace(0, 1, n).';
l3 = A * xt + sg .* randn(rs, n, 1) + 1.2 * sin(2*pi*t);   % unmodelled signal
R3 = lsqvce.fit(lsqvce.TaskDefinition.fromMatrices(A, l3, Weights=p));
D3 = R3.diag;

fprintf('\nOutliers flagged: %d\n', D3.nFlag);
fprintf(['\nThe variance factor is inflated but snooping finds little. That\n' ...
         'combination is the SIGNATURE of a systematic error. No scalar\n' ...
         'statistic will localise it -- you must LOOK at the residuals:\n\n']);

figure('Name', 'Systematic error: the residuals confess');
subplot(2,1,1);
plot(t, D3.w, '.', 'MarkerSize', 4); grid on;
yline( D3.kCrit, 'r--'); yline(-D3.kCrit, 'r--');
xlabel('time'); ylabel('w_i');
title('Normalised residuals vs TIME -- the structure is obvious');
subplot(2,1,2);
plot(t, D.w, '.', 'MarkerSize', 4); grid on;
yline( D.kCrit, 'r--'); yline(-D.kCrit, 'r--');
xlabel('time'); ylabel('w_i');
title('For comparison: the clean case, structureless');

fprintf(['ALWAYS plot w against every ordering variable you have: time,\n' ...
         'elevation, azimuth, temperature, satellite ID, along-track\n' ...
         'position. This is the single highest-value diagnostic and it is\n' ...
         'the one no automated test can do for you.\n']);


%% ===================================================================== (4)
% AN UNCONTROLLED OBSERVATION. A parameter that only one observation sees.
% Its residual is forced to zero: a blunder in it is INVISIBLE. Snooping will
% report a perfectly clean result while the parameter is quietly wrong.
fprintf('\n\n#################### 4. UNCONTROLLED OBSERVATION ####################\n');
A4 = A;
A4(:, 8) = 0;  A4(500, 8) = 1;        % parameter 8 <- observation 500 alone
l4 = A4 * xt + sg .* randn(rs, n, 1);
R4 = lsqvce.fit(lsqvce.TaskDefinition.fromMatrices(A4, l4, Weights=p));
D4 = R4.diag;

fprintf(['\nData snooping reports %d outliers -- a "clean" result.\n' ...
         'But observation 500 has r = %.2e and an INFINITE MDB.\n' ...
         'A gross error in it passes into parameter 8 intact and invisible.\n' ...
         'A clean snooping result does not exonerate it; it cannot see it.\n'], ...
         D4.nFlag, D4.rNum(500));


%% ===================================================================== (5)
% WRONG WEIGHTS. Two instrument groups; we mis-state the noise of the second
% by a factor of 3. The global test rejects, the residuals are clean and
% structureless -- that combination means the WEIGHTS are wrong, not the data.
fprintf('\n\n#################### 5. WRONG WEIGHTS ####################\n');
g = [ones(n/2,1); 2*ones(n/2,1)];
sgTrue = sg;  sgTrue(g==2) = sgTrue(g==2) * 3;       % group 2 is 3x noisier
l5 = A * xt + sgTrue .* randn(rs, n, 1);
% ...but we still claim the ORIGINAL weights p:
R5 = lsqvce.fit(lsqvce.TaskDefinition.fromMatrices(A, l5, Weights=p), lsqvce.Options(Verbose=false));
D5 = R5.diag;

fprintf('Global variance ratio: %.3f  (should be ~1)\n', D5.varRatio);
fprintf('Outliers flagged     : %d\n', D5.nFlag);
fprintf(['\nInflated variance factor + clean residuals + no structure\n' ...
         '   ==> the weights are wrong. Estimate them: declare the two\n' ...
         '   instrument groups as two STREAMS of one system of a task\n' ...
         '   (nl = 0, one arc) and let Estimate="noise", Pool="stream", ArcFactors=false find the weights.\n']);

vtask = lsqvce.TaskDefinition(Days=1);                % the whole set is one "day"
vtask.addCommon(compose("x%d", 1:u));                 % the same u unknowns, one name each
vsys = vtask.addSystem("groups");                     % one system, no locals (nl = 0)
% one stream per group, each with its own Build (groupDay, below): the rows of
% ONE group as ONE arc, and the CLAIMED (mis-stated) noise as a per-row sigma
vsys.addStream("g1", Build=@(day, ctx) groupDay(ctx, A, l5, t, g == 1), ...
               Noise=lsqvce.Noise.white(sg(g==1)));
vsys.addStream("g2", Build=@(day, ctx) groupDay(ctx, A, l5, t, g == 2), ...
               Noise=lsqvce.Noise.white(sg(g==2)));
R5v = lsqvce.fit(vtask, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, ...
          Verbose=false));

% a stream's key is 'system key|stream'; the system has no Id, so its key is
% 'x1' (the index 1 made a field name, LSQVCE.UTIL.STREAMKEY)
sigEnd = R5v.history.sigma(end, :);                  % converged group sigmas
ratio  = sigEnd(strcmp(R5v.history.streams, 'x1|g2')) / ...
         sigEnd(strcmp(R5v.history.streams, 'x1|g1'));
fprintf(['Recovered sigma ratio: %.3f   (truth: 3.000)\n' ...
         'The VCE loop rescaled the group covariances itself and re-solved:\n' ...
         'sigma0 = %.3f, variance ratio back at ~1.\n'], ratio, R5v.sigma0);

% The per-lag redundancies are estimated with Hutchinson probes: Options.Probes,
% the default 100, and that is what this run used -- the joint estimator has no
% exact route and no size-based switch, so u = 8 is probed like a gravity field.
% The estimate is stochastic: pin Options.Seed to replay the same probes.
%
% A FREE network needs a datum, and a datum is a DECLARED prior on the
% null space -- the exact minimum-constraint solution, whatever its weight:
%     task = lsqvce.TaskDefinition.fromMatrices(A, l, Weights=p);
%     task.addPrior("datum", L=B.', Sigma=1);
% (B spans the null space of A'PA). Without it the one accumulator refuses
% the rank-deficient system by name (lsqvce:qr:singularR); it never guesses.


%% ===================================================================== (6)
% NEAR-COLLINEAR PARAMETERS. Two parameters the data cannot separate. Their
% individual sigmas look fine; they are meaningless.
fprintf('\n\n#################### 6. COLLINEAR PARAMETERS ####################\n');
A6 = A;
A6(:, 2) = A6(:, 1) + 1e-4 * randn(rs, n, 1);     % col 2 ~ col 1
l6 = A6 * xt + sg .* randn(rs, n, 1);
R6 = lsqvce.fit(lsqvce.TaskDefinition.fromMatrices(A6, l6, Weights=p));
D6 = R6.diag;

fprintf(['\nParameters 1 and 2 are correlated at %.5f. Only their SUM is\n' ...
         'determined. Quoting a sigma for either one alone is a fiction,\n' ...
         'however small the number looks.\n'], D6.corrMax);


%% ==================================================================== recap
fprintf('\n\n################################################################\n');
fprintf('SUMMARY OF VERDICTS\n');
fprintf('################################################################\n');
names = {'clean', 'outliers', 'systematic', 'uncontrolled', ...
         'wrong weights', 'collinear'};
% the one-line verdict is report()'s (R5 M2: every kind rides the block path,
% so R.source is the sweep's solution struct, not the flat driver's verdict)
verdictOf = @(Rk) Rk.report('Verbose', false).verdict;
verds = {verdictOf(R), verdictOf(R2), verdictOf(R3), ...
         verdictOf(R4), verdictOf(R5), verdictOf(R6)};
for i = 1:numel(names)
    fprintf('  %-16s -> %s\n', names{i}, verds{i});
end
fprintf(['\nThe one to fear is #3 and #4: the adjustment can look healthy on\n' ...
         'every scalar statistic while being wrong. Redundancy numbers and\n' ...
         'residual PLOTS are what catch them.\n\n']);

%% ------------------------------------------------------------------ local
function out = groupDay(ctx, A, l, t, inGroup)
%GROUPDAY  The Build of case 5: the rows of one instrument group as one arc,
%   no local unknowns (nl = 0), the epochs in seconds.
    rows = find(inGroup);
    out = ctx.newDay();
    out.addArc(A(rows, :), l(rows), t(rows) * 86400);
end
