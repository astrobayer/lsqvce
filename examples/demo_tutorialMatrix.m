function [R, s] = demo_tutorialMatrix(varargin)
%DEMO_TUTORIALMATRIX  Tutorial 1: an adjustment whose design matrix fits in memory.
%
%   [R, s] = DEMO_TUTORIALMATRIX()
%   [R, s] = DEMO_TUTORIALMATRIX('verbose', true)
%
%   THE EXAMPLE OF THE MANUAL'S TUTORIAL 1 (chapter "Tutorials", part 1). Every
%   step of the chapter is a step of this file, so the chapter runs whenever
%   the test suite runs this file. Read the two side by side.
%
%   THE CASE.  A levelling network: six new points P1 .. P6 and a benchmark BM
%   of known height 100 m, joined by 15 levelled height differences. Each line
%   has its length L (km), and its standard deviation is 1 mm per sqrt(km):
%       h_ij = H_j - H_i + e_ij,     sigma_ij = 0.001 m * sqrt(L_ij)
%   The heights H1 .. H6 are the unknowns (u = 6, redundancy n - u = 9).
%
%   WHAT IT DOES (the steps of the chapter)
%     1. writes the design matrix A, the observations l (observed minus
%        computed) and the weights p = 1 ./ sigma.^2, and declares the task
%        with lsqvce.TaskDefinition.fromMatrices;
%     2. fits it: the estimates with their formal errors, sigma0, the verdict;
%     3. reads the reliability of every observation (redundancy number r, the
%        test statistic w, the minimal detectable bias);
%     4. plants a gross error of 20 mm in one line and lets the report find it,
%        then removes it with data snooping (R.snoop) and, as an alternative,
%        down-weights it with a robust estimator (Robust="huber");
%     5. drops the benchmark: the free network has a datum defect, fit refuses it
%        by name, and a declared datum prior (task.addPrior) gives the
%        minimum-constraint solution -- with the SAME height differences.
%
%   INPUT  (Name-Value, all optional)
%     verbose  1x1 logical  print the reports (default false; the comparison
%                           tables are always printed)
%   OUTPUT
%     R        1x1 lsqvce.Results   the fit of step 2 (benchmark fixed, clean data)
%     s        1x1 struct           the measured summary:
%       .truth      6x1 double   the planted heights H1 .. H6 (m)
%       .x          6x1 double   step 2's estimates (m)
%       .se         6x1 double   their formal errors (m)
%       .sigma0     1x1 double   step 2's sigma0
%       .verdict    1x1 string   step 2's verdict
%       .rMin       1x1 double   the smallest redundancy number of step 3
%       .rSum       1x1 double   sum of the redundancy numbers (= n - u = 9)
%       .blunder    1x1 struct   step 4: .row (the planted line), .flagged (the
%                                rows the report flags), .verdict, .rejected (the
%                                rows snooping rejects), .xSnoop, .factor (the
%                                robust weight factor of the planted row), .xRobust
%       .free       1x1 struct   step 5: .refusal (the identifier without a
%                                datum), .x (7x1 heights incl. BM), .dhGap (the
%                                largest difference of the height differences
%                                against step 2, m), .datumSum (sum of the
%                                free heights minus the prior's value, m)
%
%   ACCURACY
%     Measured values are printed and pinned in dev/tests/test_examples.m
%     (testTutorialMatrixIsTheManualsTutorial): the estimates within 3 formal
%     sigma of the planted heights, sigma0 within the chi-square band of 9
%     degrees of freedom, sum(r) = 9 to 1e-10, the planted blunder flagged and
%     rejected (and only it), the robust weight factor of its row below 0.5, the
%     free network's height differences equal to the fixed network's to 1e-10 m.
%
%   ERRORS  none of its own. Without a datum, fit refuses the free network with
%           the engines' identifier; the tutorial catches and prints it.
%
%   See also LSQVCE.TASKDEFINITION, LSQVCE.FIT, LSQVCE.RESULTS, LSQVCE.OPTIONS,
%            DEMO_TUTORIALBATCH, DEMO_TUTORIALSTREAMING, DEMO_LSQVCE
%
%   Created 2026-11-03 (CET, UTC+01:00) with assistance from Claude (Anthropic),
%   model Claude Opus 5.5 -- WP-163: one tutorial per kind of task.

    p = inputParser;
    p.addParameter('verbose', false);
    p.parse(varargin{:});
    o = p.Results;

    % ---- the network: 15 lines between BM (point 0) and P1 .. P6 -------------
    from = [0 1 2 3 1 4 5 6 4 3 0 2 1 6 5].';
    to   = [1 2 3 0 4 5 6 2 6 5 4 5 3 0 1].';
    rs   = RandStream('mt19937ar', 'Seed', 163);
    L    = 1 + 3 * rand(rs, 15, 1);                    % line lengths, 1 .. 4 km
    sigma = 0.001 * sqrt(L);                           % 1 mm / sqrt(km)
    Hbm   = 100;                                       % the benchmark (m)
    truth = [101.234; 102.871; 100.456; 99.312; 98.760; 101.905];   % H1 .. H6 (m)
    Hall  = [Hbm; truth];                              % point 0 is BM
    h     = Hall(to + 1) - Hall(from + 1) + sigma .* randn(rs, 15, 1);   % the levelled values

    % ---- Step 1: the design matrix, the observations, the weights -------------
    n = numel(h);  u = 6;
    A = zeros(n, u);
    for i = 1:n
        if to(i)   > 0, A(i, to(i))   =  1; end        % + H_to
        if from(i) > 0, A(i, from(i)) = -1; end        % - H_from
    end
    l = h + Hbm * (from == 0) - Hbm * (to == 0);       % observed minus computed: BM moved to l
    w = 1 ./ sigma.^2;                                 % WEIGHTS, not standard deviations
    task = lsqvce.TaskDefinition.fromMatrices(A, l, Weights=w);

    % ---- Step 2: fit -----------------------------------------------------------
    R = lsqvce.fit(task);
    T = R.coefficients();
    V = R.report(Verbose=o.verbose);
    fprintf('\nStep 2 -- the heights (m)\n%-6s %10s %10s %10s\n', 'point', 'truth', 'estimate', 'se [mm]');
    for k = 1:u
        fprintf('P%-5d %10.4f %10.4f %10.2f\n', k, truth(k), T.estimate(k), 1e3 * T.se(k));
    end
    fprintf('sigma0 = %.3f, verdict %s\n', R.sigma0, V.verdict);

    % ---- Step 3: how well is every observation controlled? ---------------------
    Q = R.reliability();
    [rMin, iMin] = min(Q.r);
    fprintf(['\nStep 3 -- sum(r) = %.4f (= n - u = %d); the weakest line is %d (BM/P%d - P%d), ' ...
             'r = %.2f, MDB = %.1f mm\n'], sum(Q.r), n - u, iMin, from(iMin), to(iMin), rMin, ...
            1e3 * Q.mdb(iMin));

    % ---- Step 4: a gross error -------------------------------------------------
    bad = 7;                                           % line P5 - P6
    lb = l;  lb(bad) = lb(bad) + 0.020;                % a 20 mm blunder
    Rb = lsqvce.fit(lsqvce.TaskDefinition.fromMatrices(A, lb, Weights=w));
    Vb = Rb.report(Verbose=o.verbose);
    flagged = Rb.diag.flagIdx(:).';
    [~, REJ, solSnoop] = Rb.snoop();
    rejected = REJ.idx(:).';
    Rr = lsqvce.fit(lsqvce.TaskDefinition.fromMatrices(A, lb, Weights=w), Robust="huber");
    factor = Rr.history(end).factor(bad);
    fprintf(['\nStep 4 -- a 20 mm blunder in line %d: verdict %s, flagged %s, snooping rejects %s;\n' ...
             '          the robust weight factor of line %d is %.3f\n'], bad, Vb.verdict, ...
            mat2str(flagged), mat2str(rejected), bad, factor);
    fprintf('%-6s %10s %10s %10s %10s\n', 'point', 'clean', 'blunder', 'snooped', 'huber');
    for k = 1:u
        fprintf('P%-5d %10.4f %10.4f %10.4f %10.4f\n', k, R.x(k), Rb.x(k), solSnoop.x(k), Rr.x(k));
    end

    % ---- Step 5: a free network and its datum ----------------------------------
    Af = zeros(n, u + 1);                              % BM is an unknown too: H0 .. H6
    for i = 1:n
        Af(i, to(i) + 1)   =  1;
        Af(i, from(i) + 1) = -1;
    end
    refusal = "";
    try
        lsqvce.fit(lsqvce.TaskDefinition.fromMatrices(Af, h, Weights=w));
    catch ME
        refusal = string(ME.identifier);
    end
    free = lsqvce.TaskDefinition.fromMatrices(Af, h, Weights=w);
    free.addPrior("datum", L=ones(1, u + 1) / (u + 1), Value=100.9, Sigma=1);   % the mean height
    Rf = lsqvce.fit(free);
    dhFixed = [R.x; Hbm] - Hbm;                         % heights relative to BM, step 2
    dhFree  = [Rf.x(2:end); Rf.x(1)] - Rf.x(1);         % the same differences, free network
    dhGap   = max(abs(dhFree - dhFixed));
    datumSum = mean(Rf.x) - 100.9;
    fprintf(['\nStep 5 -- without a datum fit refuses: %s\n' ...
             '          with the datum prior: mean height %.4f m (declared 100.9), the height\n' ...
             '          differences agree with step 2 to %.1e m\n'], refusal, mean(Rf.x), dhGap);

    s = struct('truth', truth, 'x', T.estimate, 'se', T.se, 'sigma0', R.sigma0, ...
               'verdict', string(V.verdict), 'rMin', rMin, 'rSum', sum(Q.r), ...
               'blunder', struct('row', bad, 'flagged', flagged, 'verdict', string(Vb.verdict), ...
                                 'rejected', rejected, 'xSnoop', solSnoop.x(:), ...
                                 'factor', factor, 'xRobust', Rr.x(:)), ...
               'free', struct('refusal', refusal, 'x', Rf.x(:), 'dhGap', dhGap, ...
                              'datumSum', datumSum));
end
