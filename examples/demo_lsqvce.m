%% DEMO_LSQVCE  Worked examples: how to judge a least-squares adjustment.
%
%   Six scenarios, each isolating one failure mode, so you can see what the
%   report says when a specific thing is wrong -- and, just as importantly,
%   what it says when things are fine. Every run is the USER path:
%
%       R = lsqvce.fit(struct('A', A, 'l', l, 'p', p));
%
%   with the diagnostics on R.diag and the verdict on R.source. Case 5 also
%   shows the CURE: a streaming declaration with one covariance field per
%   instrument group and Estimate="noise", Pool="stream", ArcFactors=false, so fit estimates the weights it needs.
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.
%   Updated 2026-08-31 (CEST, UTC+02:00), v4 SCOPED tables (Level L11,
%   WP-77): the declaration tables carry the scope column (commonParams
%   'common'|'day', localParams 'arc'|'day'|'system'); localsPerArc left the
%   USER surface. With assistance from Claude (Anthropic), model
%   Claude Fable 5.
%   Updated 2026-08-31 (CEST, UTC+02:00), v3 data contract (Level L10,
%   WP-74): USER surface only -- lsqvce.engine.adjust/varComp replaced by
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
R = lsqvce.fit(struct('A', A, 'l', l, 'p', p));
D = R.diag;


%% ===================================================================== (2)
% OUTLIERS. Three gross errors. Watch the variance factor inflate, and watch
% data snooping localise them -- with Bonferroni, so the list is short and real
% rather than long and mostly noise.
fprintf('\n\n#################### 2. OUTLIERS ####################\n');
l2 = l;
l2([137 892 1455]) = l2([137 892 1455]) + [12; -15; 9] .* sg([137 892 1455]);
R2 = lsqvce.fit(struct('A', A, 'l', l2, 'p', p));
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
R3 = lsqvce.fit(struct('A', A, 'l', l3, 'p', p));
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
R4 = lsqvce.fit(struct('A', A4, 'l', l4, 'p', p));
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
R5 = lsqvce.fit(struct('A', A, 'l', l5, 'p', p), lsqvce.Options(Verbose=false));
D5 = R5.diag;

fprintf('Global variance ratio: %.3f  (should be ~1)\n', D5.varRatio);
fprintf('Outliers flagged     : %d\n', D5.nFlag);
fprintf(['\nInflated variance factor + clean residuals + no structure\n' ...
         '   ==> the weights are wrong. Estimate them: declare the two\n' ...
         '   instrument groups as two OBSERVATION TYPES of one streaming\n' ...
         '   system (nl = 0, one arc) and let Estimate="noise", Pool="stream", ArcFactors=false find the weights.\n']);

vdata.commonParams = [repmat({'common'}, u, 1), num2cell((1:u).')];   % scoped, nc = u
vdata.days = 1;                                      % the whole set is one "day"
vs = struct();
vs.buildDay = @(day, ot) groupDay(A, l5, t, g, ot);  % rows of ONE group
vs.covariance.g1 = struct('type', 'diag', 'sigma', sg(g==1));   % the CLAIMED
vs.covariance.g2 = struct('type', 'diag', 'sigma', sg(g==2));   % (mis-stated) noise
vdata.systems = vs;
R5v = lsqvce.fit(vdata, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, ...
          Verbose=false));

sigEnd = R5v.history.sigma(end, :);                  % converged group sigmas
ratio  = sigEnd(strcmp(R5v.history.streams, 'x1|g2')) / ...
         sigEnd(strcmp(R5v.history.streams, 'x1|g1'));
fprintf(['Recovered sigma ratio: %.3f   (truth: 3.000)\n' ...
         'The VCE loop rescaled the group covariances itself and re-solved:\n' ...
         'sigma0 = %.3f, variance ratio back at ~1.\n'], ratio, R5v.sigma0);

% At gravity-field sizes (u in the thousands) the exact per-group traces of
% this loop dominate; the ENGINE tier has a matrix-free Hutchinson trace
% backend for that regime (doc lsqvce.engine.varComp). At u = 8 the exact path
% is both faster and exact, so that is what fit ran here.


%% ===================================================================== (6)
% NEAR-COLLINEAR PARAMETERS. Two parameters the data cannot separate. Their
% individual sigmas look fine; they are meaningless.
fprintf('\n\n#################### 6. COLLINEAR PARAMETERS ####################\n');
A6 = A;
A6(:, 2) = A6(:, 1) + 1e-4 * randn(rs, n, 1);     % col 2 ~ col 1
l6 = A6 * xt + sg .* randn(rs, n, 1);
R6 = lsqvce.fit(struct('A', A6, 'l', l6, 'p', p));
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
function D = groupDay(A, l, t, g, obsType)
% The buildDay of case 5: one instrument group as one observation type.
    rows = find(g == 1 + strcmp(obsType, 'g2'));
    D.A = A(rows, :);
    D.B = zeros(numel(rows), 0);        % no local unknowns (nl = 0, the default)
    D.l = l(rows);
    D.t = t(rows) * 86400;              % seconds
end
