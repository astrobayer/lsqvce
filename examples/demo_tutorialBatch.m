function [R, s] = demo_tutorialBatch(varargin)
%DEMO_TUTORIALBATCH  Tutorial 2: an adjustment read in batches, linear and nonlinear.
%
%   [R, s] = DEMO_TUTORIALBATCH()
%   [R, s] = DEMO_TUTORIALBATCH('verbose', true)
%
%   THE EXAMPLE OF THE MANUAL'S TUTORIAL 2 (chapter "Tutorials", part 2). Every
%   step of the chapter is a step of this file, so the chapter runs whenever
%   the test suite runs this file. Read the two side by side.
%
%   THE CASE.  Ten years of hourly sea level at a tide gauge (87 600 epochs),
%       y(t) = c + d t + a1 cos(w t) + b1 sin(w t) + a2 cos(2 w t) + b2 sin(2 w t),
%   t in years, w = 2 pi / year: an offset, a trend and the annual and
%   semi-annual tides -- six unknowns. The design has 87 600 rows: you do not
%   want it in memory, and in a real case it is spread over files. So the rows
%   are handed to the toolbox ONE MONTH AT A TIME (120 batches), by a function
%   you write. The noise is AR(1) from hour to hour (phi = 0.6, sigma 2 cm).
%   Then a nonlinear model: the post-event relaxation of a GNSS station,
%   u(t) = a (1 - exp(-t / tau)) + c, solved by Gauss-Newton from a rough start.
%
%   WHAT IT DOES (the steps of the chapter)
%     1. writes the batch function [A, l, p] = batchFcn(k): the rows of month k,
%        their observations and their weights, the same rows every call;
%     2. declares the task with lsqvce.TaskDefinition.fromBatches and fits it;
%     3. checks it against the matrix kind on the same rows (equal to 1e-12);
%     4. declares the hourly correlation: the batch function returns the full
%        m x m weight matrix of each month (the inverse of its AR(1)
%        covariance) instead of a weight vector -- the estimates hardly move,
%        their formal errors grow to honest values;
%     5. fits the nonlinear model: [J, w, p] = batchFcn(k, x) returns the
%        Jacobian at x and w = l - f(x); X0 is the start; fit runs damped
%        Gauss-Newton and records every iteration.
%
%   INPUT  (Name-Value, all optional)
%     verbose  1x1 logical  print the reports (default false; the comparison
%                           tables are always printed)
%   OUTPUT
%     R        1x1 lsqvce.Results   the fit of step 4 (correlation declared)
%     s        1x1 struct           the measured summary:
%       .truth      6x1 double   the planted [c d a1 b1 a2 b2]
%       .white      1x1 struct   step 2 (weight vector): .x .se .sigma0
%       .matrixGap  1x1 double   step 3: max|x_batched - x_matrix| / max|x|
%       .full       1x1 struct   step 4 (full weight matrix): .x .se .sigma0
%       .seRatio    6x1 double   step 4 se / step 2 se
%       .nonlinear  1x1 struct   step 5: .truth .x0 .x .se .passes .stopReason
%
%   ACCURACY
%     Measured values are printed and pinned in dev/tests/test_fit.m
%     (testDemoTutorialBatchIsTheManualsTutorial): step 4's estimates within
%     4 of its formal sigma of the truth, the batched solve
%     equal to the matrix kind to 1e-10 relative, step 4's sigma0 within 10 %
%     of 1, the trend's se ratio near sqrt((1 + phi) / (1 - phi)) = 2 (within
%     25 %), the nonlinear fit converged with its parameters within 4 formal
%     sigma of the truth.
%
%     MEASURED 2026-11-03 (R2026a, v7.0.0, MATLAB MCP bridge): batched against
%     one matrix 1.2e-15; sigma0 0.999 under both weightings, the se ratio 2.00
%     for all six parameters; the trend 0.00313 +- 0.00002 (white) and
%     +- 0.00005 (AR(1)) against the planted 0.0032 -- 3.5 against 1.4 of its
%     formal sigma; Gauss-Newton converged in 5 iterations, tau 0.3416 +-
%     0.0059 (planted 0.35).
%
%   ERRORS  none of its own; a batch function that returns rows of the wrong
%           size is refused by the engines (lsqvce:batchSystem:*).
%
%   See also LSQVCE.TASKDEFINITION, LSQVCE.FIT, LSQVCE.RESULTS, LSQVCE.OPTIONS,
%            DEMO_TUTORIALMATRIX, DEMO_TUTORIALSTREAMING, DEMO_WORKFLOW
%
%   Created 2026-11-03 (CET, UTC+01:00) with assistance from Claude (Anthropic),
%   model Claude Opus 5.5 -- WP-163: one tutorial per kind of task.

    p = inputParser;
    p.addParameter('verbose', false);
    p.parse(varargin{:});
    o = p.Results;

    % ---- the record: 10 years of hourly sea level, cut into 120 months ---------
    cfg.hours  = 87600;                                % 10 years of 8760 h
    cfg.edges  = round(linspace(0, cfg.hours, 121));   % month k = rows edges(k)+1 .. edges(k+1)
    cfg.truth  = [0.120; 0.0032; 0.085; -0.040; 0.021; 0.013];   % m, m/yr, m, m, m, m
    cfg.phi    = 0.6;                                  % hour-to-hour correlation
    cfg.sigma  = 0.02;                                 % m
    names = ["offset" "trend" "annualCos" "annualSin" "semiCos" "semiSin"];

    % ---- Steps 1 and 2: the batch function, the task, the fit ------------------
    batchFcn = @(k) monthRows(cfg, k, "white");
    [A1, l1, p1] = batchFcn(1);                        % look at one batch first
    fprintf('\nStep 1 -- batch 1: A is %d x %d, l %d x 1, p %d x 1\n', size(A1), numel(l1), numel(p1));
    task = lsqvce.TaskDefinition.fromBatches(batchFcn, NumBatches=120, NumUnknowns=6);
    Rw = lsqvce.fit(task);
    Tw = Rw.coefficients();
    Vw = Rw.report(Verbose=o.verbose);

    % ---- Step 3: the same rows as one matrix -----------------------------------
    [Aall, lall, pall] = deal(cell(120, 1));
    for k = 1:120, [Aall{k}, lall{k}, pall{k}] = batchFcn(k); end
    Rm = lsqvce.fit(lsqvce.TaskDefinition.fromMatrices(vertcat(Aall{:}), vertcat(lall{:}), ...
                                                       Weights=vertcat(pall{:})));
    matrixGap = max(abs(Rw.x - Rm.x)) / max(abs(Rm.x));
    fprintf('Step 3 -- batched against one matrix: max relative difference %.1e\n', matrixGap);

    % ---- Step 4: the hourly correlation, as a full weight matrix per month -----
    taskFull = lsqvce.TaskDefinition.fromBatches(@(k) monthRows(cfg, k, "ar1"), ...
                                                 NumBatches=120, NumUnknowns=6);
    R = lsqvce.fit(taskFull);
    T = R.coefficients();
    V = R.report(Verbose=o.verbose);
    fprintf(['\nSteps 2 and 4 -- weight vector (white) against full weight matrix (AR(1))\n' ...
             '%-10s %10s %20s %20s %8s\n'], 'param', 'truth', 'white', 'AR(1)', 'se ratio');
    for k = 1:6
        fprintf('%-10s %10.5f %10.5f +- %7.5f %10.5f +- %7.5f %8.2f\n', names(k), cfg.truth(k), ...
                Tw.estimate(k), Tw.se(k), T.estimate(k), T.se(k), T.se(k) / Tw.se(k));
    end
    fprintf('sigma0 %14s %20.3f %20.3f   verdicts %s / %s\n', '', Rw.sigma0, R.sigma0, ...
            Vw.verdict, V.verdict);

    % ---- Step 5: a nonlinear model, Gauss-Newton -------------------------------
    nl.truth = [0.042; 0.35; 0.003];                   % a (m), tau (yr), c (m)
    nl.x0    = [0.03; 0.2; 0];                         % a rough start
    nl.sigma = 0.002;                                  % m
    task5 = lsqvce.TaskDefinition.fromBatches(@(k, x) relaxRows(nl, k, x), NumBatches=2, X0=nl.x0);
    R5 = lsqvce.fit(task5, MaxPasses=20, Tolerance=1e-10);
    T5 = R5.coefficients();
    nPass = numel(R5.history);
    fprintf(['\nStep 5 -- Gauss-Newton from [%g %g %g]: %s after %d iterations\n' ...
             '%-6s %10s %20s\n'], nl.x0, R5.stopReason, nPass, 'param', 'truth', 'estimate');
    nm5 = ["a" "tau" "c"];
    for k = 1:3
        fprintf('%-6s %10.5f %10.5f +- %7.5f\n', nm5(k), nl.truth(k), T5.estimate(k), T5.se(k));
    end

    s = struct('truth', cfg.truth, ...
               'white', struct('x', Tw.estimate, 'se', Tw.se, 'sigma0', Rw.sigma0), ...
               'matrixGap', matrixGap, ...
               'full', struct('x', T.estimate, 'se', T.se, 'sigma0', R.sigma0), ...
               'seRatio', T.se ./ Tw.se, ...
               'nonlinear', struct('truth', nl.truth, 'x0', nl.x0, 'x', T5.estimate, 'se', T5.se, ...
                                   'passes', nPass, 'stopReason', string(R5.stopReason)));
end

function [A, l, P] = monthRows(cfg, k, weighting)
%MONTHROWS  The rows of month K: design A (m x 6), observations l (m x 1) and
%   weights -- a vector 1/sigma^2 ("white") or the full m x m inverse of the
%   month's AR(1) covariance ("ar1"). The noise is drawn from a stream seeded
%   on k, so every call returns the same rows (the toolbox reads every batch
%   twice: once to solve, once for the diagnostics).
    rows = (cfg.edges(k) + 1 : cfg.edges(k + 1)).';
    m = numel(rows);
    t = (rows - 1) / 8760;                             % years since the start
    w = 2 * pi;                                        % one cycle per year
    A = [ones(m, 1), t, cos(w * t), sin(w * t), cos(2 * w * t), sin(2 * w * t)];
    rs = RandStream('mt19937ar', 'Seed', 5000 + k);
    e  = filter(1, [1 -cfg.phi], cfg.sigma * sqrt(1 - cfg.phi^2) * randn(rs, m + 100, 1));
    l  = A * cfg.truth + e(101:end);                   % AR(1) noise, start-up dropped
    if weighting == "white"
        P = ones(m, 1) / cfg.sigma^2;
    else
        % the inverse of the AR(1) covariance C_ij = sigma^2 phi^|i-j| is tridiagonal
        q = cfg.phi;
        P = diag([1; (1 + q^2) * ones(m - 2, 1); 1]) - q * diag(ones(m - 1, 1), 1) ...
            - q * diag(ones(m - 1, 1), -1);
        P = P / (cfg.sigma^2 * (1 - q^2));
    end
end

function [J, w, p] = relaxRows(nl, k, x)
%RELAXROWS  Batch K (one year of daily positions) of the relaxation model at the
%   current parameters X = [a; tau; c]: the Jacobian J (m x 3), the reduced
%   observations w = l - f(x) and the weights p.
    t  = ((k - 1) * 365 + (1:365).') / 365;           % years since the event
    rs = RandStream('mt19937ar', 'Seed', 7000 + k);
    l  = nl.truth(1) * (1 - exp(-t / nl.truth(2))) + nl.truth(3) + nl.sigma * randn(rs, 365, 1);
    E  = exp(-t / x(2));
    f  = x(1) * (1 - E) + x(3);
    J  = [1 - E, -x(1) * t .* E / x(2)^2, ones(365, 1)];
    w  = l - f;
    p  = ones(365, 1) / nl.sigma^2;
end
