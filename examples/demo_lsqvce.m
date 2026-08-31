%% DEMO_LSQVCE  Worked examples: how to judge a least-squares adjustment.
%
%   Six scenarios, each isolating one failure mode, so you can see what the
%   report says when a specific thing is wrong -- and, just as importantly,
%   what it says when things are fine.
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.

clear; clc;
rs = RandStream('mt19937ar', 'Seed', 2026);

%% ------------------------------------------------------------------ setup
% A generic batch source. In a real problem this reads from disk, or generates
% partials from an orbit integration; here we just slice a matrix.
    function src = makeSource(A, l, p, nb)
        n = numel(l);
        e = round(linspace(0, n, nb + 1));
        src = @(k) deal(A(e(k)+1:e(k+1), :), l(e(k)+1:e(k+1)), ...
                        p(e(k)+1:e(k+1)));
    end

n = 2000;  u = 8;  nb = 20;
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
[sol, D, V] = lsqvce.engine.adjust(u, makeSource(A, l, p, nb), nb);


%% ===================================================================== (2)
% OUTLIERS. Three gross errors. Watch the variance factor inflate, and watch
% data snooping localise them -- with Bonferroni, so the list is short and real
% rather than long and mostly noise.
fprintf('\n\n#################### 2. OUTLIERS ####################\n');
l2 = l;
l2([137 892 1455]) = l2([137 892 1455]) + [12; -15; 9] .* sg([137 892 1455]);
[~, D2, V2] = lsqvce.engine.adjust(u, makeSource(A, l2, p, nb), nb);

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
[~, D3, V3] = lsqvce.engine.adjust(u, makeSource(A, l3, p, nb), nb);

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
[~, D4, V4] = lsqvce.engine.adjust(u, makeSource(A4, l4, p, nb), nb);

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
[~, D5, V5] = lsqvce.engine.adjust(u, makeSource(A, l5, p, nb), nb, 'verbose', false);

fprintf('Global variance ratio: %.3f  (should be ~1)\n', D5.varRatio);
fprintf('Outliers flagged     : %d\n', D5.nFlag);
fprintf(['\nInflated variance factor + clean residuals + no structure\n' ...
         '   ==> the weights are wrong. Estimate them:\n']);

e  = round(linspace(0, n, nb + 1));
gsrc = @(k) g(e(k)+1:e(k+1));
VC = lsqvce.engine.varComp(u, makeSource(A, l5, p, nb), nb, gsrc, 2);

fprintf(['Recovered sigma ratio: %.3f   (truth: 3.000)\n' ...
         'Apply VC.scale to the covariances and re-run.\n'], ...
         VC.sigma(2) / VC.sigma(1));

% --------------------------------------------------- SAME PROBLEM, HUTCHINSON
% The EXACT path above is O(G*u^3): it forms N^-1 N_g for every group. For a
% high-degree gravity field (u in the thousands) that dominates. The matrix-free
% Hutchinson backend never forms N^-1; it estimates each trace with Z random
% probes against the Cholesky factor, costing O(u^3 + G*Z*u^2). The crossover is
% u ~ Z, so it is the method you reach for at GRAVITY-FIELD variance factors -- not here,
% where u = 8 makes the exact path both faster and exact. We run it anyway to
% show the SAME ANSWER comes back. Z is set high (4000) so the Monte-Carlo error
% (~1/sqrt(Z)) is small enough that the agreement is visibly tight; a smaller Z
% is cheaper but noisier.
VCh = lsqvce.engine.varComp(u, makeSource(A, l5, p, nb), nb, gsrc, 2, ...
                     'backend', 'hutchinson', 'probes', 4000, ...
                     'seed', 2026, 'verbose', false);

fprintf(['\n[Hutchinson backend, Z = %d probes]\n' ...
         'Recovered sigma ratio: %.3f   (exact: %.3f, truth: 3.000)\n' ...
         'Same answer, matrix-free. On small u it is SLOWER and noisier;\n' ...
         'reach for it only when u is in the thousands.\n'], ...
         VCh.probes, VCh.sigma(2) / VCh.sigma(1), VC.sigma(2) / VC.sigma(1));


%% ===================================================================== (6)
% NEAR-COLLINEAR PARAMETERS. Two parameters the data cannot separate. Their
% individual sigmas look fine; they are meaningless.
fprintf('\n\n#################### 6. COLLINEAR PARAMETERS ####################\n');
A6 = A;
A6(:, 2) = A6(:, 1) + 1e-4 * randn(rs, n, 1);     % col 2 ~ col 1
l6 = A6 * xt + sg .* randn(rs, n, 1);
[~, D6, V6] = lsqvce.engine.adjust(u, makeSource(A6, l6, p, nb), nb);

fprintf(['\nParameters 1 and 2 are correlated at %.5f. Only their SUM is\n' ...
         'determined. Quoting a sigma for either one alone is a fiction,\n' ...
         'however small the number looks.\n'], D6.corrMax);


%% ==================================================================== recap
fprintf('\n\n################################################################\n');
fprintf('SUMMARY OF VERDICTS\n');
fprintf('################################################################\n');
names = {'clean', 'outliers', 'systematic', 'uncontrolled', ...
         'wrong weights', 'collinear'};
verds = {V.verdict, V2.verdict, V3.verdict, V4.verdict, ...
         V5.verdict, V6.verdict};
for i = 1:numel(names)
    fprintf('  %-16s -> %s\n', names{i}, verds{i});
end
fprintf(['\nThe one to fear is #3 and #4: the adjustment can look healthy on\n' ...
         'every scalar statistic while being wrong. Redundancy numbers and\n' ...
         'residual PLOTS are what catch them.\n\n']);
