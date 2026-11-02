function wp122_measure(outFile)
%WP122_MEASURE  The diagnosis WP-122's first CI round needs (template/kbr_variational).
%
%   WP122_MEASURE(outFile)
%
%   A MEASUREMENT, not a test: it asserts nothing and prints. It runs on the
%   self-hosted CI runner through .github/workflows/measure.yml; the job log
%   is the record. CI round 1 of WP-122 (be0f5a9) failed
%   testKbrPartialsAgainstReintegration and testKbrClosedLoopThroughFit
%   (the latter by an error), and the gate prints only the test names:
%     PART A  the two tests through runtests (their failure diagnostics,
%             with the values and the stack, go to the log);
%     PART B  the partials case by case: |Psi * delta - d| / |d|;
%     PART C  the closed loop step by step (fixture, declaration, fit), every
%             error with its whole cause chain.
%
%   INPUT   outFile  char  optional diary file
%
%   See also kbr_variational, fixture_kbr_l1b
%
%   Created 2026-10-12 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- WP-122, CI round 1.
    if nargin > 0 && ~isempty(outFile)
        diary(outFile);  cleanDiary = onCleanup(@() diary('off')); %#ok<NASGU>
    end
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    p0 = path;  restore = onCleanup(@() path(p0)); %#ok<NASGU>
    addpath(root);  addpath(fullfile(root, 'template'));
    addpath(fullfile(root, 'examples'));  addpath(fullfile(root, 'dev', 'tests'));
    fprintf('WP-122 measurement, %s, MATLAB %s\n', char(datetime('now')), version);

    fprintf('\n===== PART A: the two failing tests through runtests =====\n');
    try
        r = runtests(fullfile(root, 'dev', 'tests', 'test_physOracle.m'), 'ProcedureName', ...
                     {'testKbrPartialsAgainstReintegration', 'testKbrClosedLoopThroughFit'});
        disp(table(r));
    catch ME
        printChain(ME);
    end

    fprintf('\n===== PART B: the partials, case by case =====\n');
    try
        partials();
    catch ME
        printChain(ME);
    end

    fprintf('\n===== PART C: the closed loop, step by step =====\n');
    try
        d = tempname;  mkdir(d);  cl = onCleanup(@() rmdir(d, 's')); %#ok<NASGU>
        t = tic;  [cfg, truth] = fixture_kbr_l1b(d);
        fprintf('fixture ok (%.1f s): %s\n', toc(t), strjoin({dir(d).name}, ' '));
        t = tic;  [task, kit, data] = kbr_variational(cfg);   % WP-161: the task first
        fprintf('declaration ok (%.1f s): nc %d, local rows %d\n', toc(t), ...
                size(data.commonParams, 1), size(data.systems.localParams, 1));
        t = tic;  R = lsqvce.fit(task, Verbose=false);
        dc = kit.unitCommon .* R.x(:);
        fprintf('fit ok (%.1f s): |dc - dc_true| / |dc_true| = %.3e\n', toc(t), ...
                norm(dc - truth.dc) / norm(truth.dc));
        disp([truth.dc, dc]);
    catch ME
        printChain(ME);
    end
end

function partials()
    gf = struct('GM', 3.986004415e14, 'R', 6378136.3, 'C', blkdiag(1, zeros(2)), 'S', zeros(3));
    [~, kit] = kbr_variational(struct('days', 60000, 'dataDir', tempdir, 'gravityField', gf, ...
                                      'trf2crf', @(t) eye(3), 'degree', [2 2]));
    GM = gf.GM;  R = gf.R;  N = 6;
    [nn, mm, cs] = kit.coeffIndex(N);
    c = zeros(numel(nn), 1);  c(nn == 0) = 1;
    c(nn == 2 & mm == 0 & cs == 0) = -4.841694573e-4;
    rs = RandStream('mt19937ar', 'Seed', 3);
    k3 = nn >= 3;  c(k3) = 1e-6 * randn(rs, nnz(k3), 1) ./ nn(k3).^2;
    est = find(nn >= 2);  np = numel(est);
    y0 = leo(GM, R, 0.4);
    run = @(cc, yy, b, s) kit.propagate(kit.makeField(cc, N, est, GM, R), @(t) eye(3), @acc, b, s, 0, yy, 5, 1080);
    [Y, Psi] = run(c, y0, zeros(3, 1), ones(3, 1));
    fprintf('Psi %d x %d, |Psi(:,1:6)| %.3e, final r %.3f km\n', size(Psi), norm(Psi(:, 1:6)), norm(Y(end, 1:3)) / 1e3);
    j = find(est == find(nn == 5 & mm == 3 & cs == 1));
    cases = {'IC x', 1, @() run(c, y0 + [1; 0; 0; 0; 0; 0], zeros(3, 1), ones(3, 1)), 1;
             'IC vy', 5, @() run(c, y0 + [0; 0; 0; 0; 1e-3; 0], zeros(3, 1), ones(3, 1)), 1e-3;
             'S53', 6 + j, @() run(c + 1e-7 * ((1:numel(c)).' == est(j)), y0, zeros(3, 1), ones(3, 1)), 1e-7;
             'bias y', 6 + np + 2, @() run(c, y0, [0; 1e-8; 0], ones(3, 1)), 1e-8;
             'scale x', 9 + np + 1, @() run(c, y0, zeros(3, 1), [1.1; 1; 1]), 0.1};
    for k = 1:size(cases, 1)
        Y2 = cases{k, 3}();
        d = (Y2(end, :) - Y(end, :)).';
        p = Psi(:, cases{k, 2}) * cases{k, 4};
        fprintf('%-8s |d| %.3e  |Psi*delta| %.3e  rel %.3e\n   d    %s\n   Psid %s\n', cases{k, 1}, ...
                norm(d), norm(p), norm(p - d) / norm(d), mat2str(d.', 4), mat2str(p.', 4));
    end
end

function y0 = leo(GM, R, u)
    a = R + 490e3;  vc = sqrt(GM / a);  i = deg2rad(89);
    y0 = [a * [cos(u); sin(u) * cos(i); sin(u) * sin(i)]; vc * [-sin(u); cos(u) * cos(i); cos(u) * sin(i)]];
end

function [Q, a] = acc(t)
    th = 2 * pi / 5650 * t + 0.2;
    Q = [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1] * [1 0 0; 0 cos(0.3) -sin(0.3); 0 sin(0.3) cos(0.3)];
    a = [-1e-7 * (1 + 0.3 * cos(th)); 2e-8 * sin(th); 5e-8 * cos(2 * th)];
end

function printChain(ME)
    fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
    for k = 1:numel(ME.cause)
        fprintf('--- cause %d ---\n', k);
        printChain(ME.cause{k});
    end
end
