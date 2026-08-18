function results = run_all_lsqvce(root, skip)
%RUN_ALL_LSQVCE  Run every lsqvce test suite file-by-file and print a summary.
%   results = run_all_lsqvce()          uses the folder this file lives in
%   results = run_all_lsqvce(root)      root is the toolbox root (contains +lsqvce)
%   results = run_all_lsqvce(root, skip)  additionally skips the named suites
%
%   WHY FILE-BY-FILE, NOT runtests('tests','IncludeSubfolders',true)
%     Each accum test file defines a local helper makeCfg(); running the whole
%     tree in one runtests call is fine, but running file-by-file keeps every
%     invocation short enough to finish under a 4-minute MCP timeout and makes a
%     hang attributable to a single file. It also isolates the parallel suites,
%     which spin up a pool.
%
%   INPUT
%     root  1x1 char/string  toolbox root (contains +lsqvce). Default: the
%                            parent of the folder this file lives in.
%     skip  1xS cell/string  suite names (file names WITHOUT extension, e.g.
%                            'test_parallel') to SKIP. Skipped suites are
%                            printed in the summary AND recorded in the log
%                            file, so a skip is always visible, never silent.
%                            Default: {} (skip nothing). The POLICY of what
%                            to skip and when lives in CI_GATE (one place per
%                            decision); this argument is only the mechanism.
%
%   OUTPUT
%     results  1xF struct  one row per suite file with fields
%              .name    char     the test file (no extension)
%              .passed  1x1      number of passing test points
%              .total   1x1      total test points
%              .failed  1x1      number of failing test points
%              .incomplete 1x1   number that did not complete
%              .secs    1x1      wall-clock seconds for that suite. Compare
%                                these between runs: a test that quietly starts
%                                running a full VCE loop costs seconds and shows
%                                up NOWHERE in a pass count.
%
%   `run_lsqvce('all')` delegates here, so it is the same run under a different
%   name -- one gate, one implementation. Use RUN_LSQVCE for a single suite or a
%   name pattern while iterating.
%
%   ERRORS
%     run_all_lsqvce:noToolbox   root does not contain a +lsqvce package.
%     run_all_lsqvce:badSkip     skip is not a cell array of char / string
%                                array of suite names.
%
%   See also RUN_LSQVCE, CI_GATE
%
%   Created 2026-07-16 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.
%   Updated 2026-08-18 (CEST, UTC+02:00): added the `skip` argument (WP-30;
%   the CI gate skips pool suites on runners without PCT — visibly, in the
%   summary and the log). With assistance from Claude (Anthropic), model
%   Claude Fable 5.
    if nargin < 1 || isempty(root)
        here = fileparts(mfilename('fullpath'));   % .../tests
        root = fileparts(here);                    % toolbox root
    end
    if nargin < 2 || isempty(skip)
        skip = {};
    end
    if isstring(skip), skip = cellstr(skip); end
    if ~iscellstr(skip)                                        %#ok<ISCLSTR>
        error('run_all_lsqvce:badSkip', ...
            'skip must be a cell array of char or a string array of suite names.');
    end
    if exist(fullfile(root, '+lsqvce'), 'dir') ~= 7
        error('run_all_lsqvce:noToolbox', ...
            '%s does not contain a +lsqvce package.', root);
    end
    addpath(root);

    suites = { ...
        fullfile(root,'tests','test_lsqvce.m'), ...
        fullfile(root,'tests','test_lsqvce_ext.m'), ...
        fullfile(root,'tests','test_lsqvce_opts.m'), ...
        fullfile(root,'tests','test_lsqvce_combine.m'), ...
        fullfile(root,'tests','test_util.m'), ...
        fullfile(root,'tests','test_traceEstimate.m'), ...
        fullfile(root,'tests','test_ellmer.m'), ...
        fullfile(root,'tests','test_lsqAccum.m'), ...
        fullfile(root,'tests','test_vceMath.m'), ...
        fullfile(root,'tests','test_parallel.m'), ...
        fullfile(root,'tests','test_vceFunction.m'), ...
        fullfile(root,'tests','test_signalCore.m'), ...
        fullfile(root,'tests','test_streamUtils.m'), ...
        fullfile(root,'tests','test_precision.m'), ...
        fullfile(root,'tests','test_reliability.m'), ...
        fullfile(root,'tests','test_adjustBlock.m'), ...
        fullfile(root,'tests','test_adjustVCE.m'), ...
        fullfile(root,'tests','test_contribution.m'), ...
        fullfile(root,'tests','test_ci_gate.m') };
    % Note: test_speed.m is a benchmark, not a pass/fail suite; run it manually.

    results = struct('name', {}, 'passed', {}, 'total', {}, ...
                     'failed', {}, 'incomplete', {}, 'secs', {});
    grandP = 0; grandT = 0; grandF = 0;
    % A file that does not LOAD is the worst outcome there is: it hides an
    % unknown number of tests. It used to be printed and then dropped, so
    % the headline still read '0 failed' while 41 tests never ran (a
    % duplicated function name in test_util, 2026-07-19). Counted now.
    loadErrors = {};
    capture = {};
    skipped = {};
    for i = 1:numel(suites)
        f = suites{i};
        [~, nm] = fileparts(f);
        if ismember(nm, skip)
            % A skip is a DECISION, so it is printed like a result line and
            % written to the log below - a suite that silently vanishes is
            % the same failure class as one that fails to load (§0a-2).
            fprintf('  %-20s : SKIPPED (requested via the skip argument)\n', nm);
            skipped{end+1} = nm;                               %#ok<AGROW>
            continue
        end
        if exist(f, 'file') ~= 2
            fprintf('  %-20s : MISSING (%s)\n', nm, f);
            continue
        end
        try
            % Capture the console output so warnings and failure diagnostics
            % land in the log file rather than scrolling past. The summary
            % printed below is unaffected.
            tSuite = tic;
            capt = evalc('r = runtests(f);');
            secs = toc(tSuite);
            capture{end+1} = sprintf('\n===== %s =====\n%s', nm, capt); %#ok<AGROW>
            p = sum([r.Passed]); t = numel(r);
            fl = sum([r.Failed]); ic = sum([r.Incomplete]);
        catch ME
            fprintf('  %-20s : LOAD ERROR %s   <-- NOTHING IN THIS FILE RAN\n', ...
                    nm, ME.message);
            loadErrors{end+1} = sprintf('%s: %s', nm, ME.message); %#ok<AGROW>
            results(end+1) = struct('name',nm,'passed',0,'total',NaN, ...
                                    'failed',NaN,'incomplete',NaN,'secs',NaN); %#ok<AGROW>
            capture{end+1} = sprintf(['\n===== %s =====\nLOAD ERROR: %s\n' ...
                'NOTHING IN THIS FILE RAN. The pass counts below EXCLUDE it, ' ...
                'so a green total here does not mean the suite is green.\n'], ...
                nm, ME.message); %#ok<AGROW>
            continue
        end
        results(end+1) = struct('name',nm,'passed',p,'total',t, ...
                                'failed',fl,'incomplete',ic,'secs',secs); %#ok<AGROW>
        grandP = grandP + p; grandT = grandT + t; grandF = grandF + fl;
        flag = '';
        if fl > 0 || ic > 0, flag = '   <-- CHECK'; end
        fprintf('  %-20s : %3d/%-3d pass, %d fail, %d incomplete  %6.1fs%s\n', ...
                nm, p, t, fl, ic, secs, flag);
    end

    % WALL CLOCK. Recorded because a slowdown is otherwise unattributable: a
    % test that quietly starts running a full VCE loop costs seconds and shows up
    % nowhere in a pass count. Compare the per-suite seconds between runs.
    totalSecs = sum([results.secs], 'omitnan');
    fprintf('\n==== GRAND TOTAL : %d / %d passed, %d failed  (%.1fs) ====\n', ...
            grandP, grandT, grandF, totalSecs);
    if ~isempty(skipped)
        fprintf('==== %d SUITE(S) SKIPPED ON REQUEST: %s ====\n', ...
                numel(skipped), strjoin(skipped, ', '));
    end
    if ~isempty(loadErrors)
        fprintf(['==== %d FILE(S) FAILED TO LOAD - THE TOTAL ABOVE IS ' ...
                 'INCOMPLETE ====\n'], numel(loadErrors));
        for i = 1:numel(loadErrors)
            fprintf('       %s\n', loadErrors{i});
        end
    end

    % Full log, same convention as RUN_LSQVCE, so a run can be inspected or
    % shared verbatim afterwards instead of scrolling the command window.
    logFile = fullfile(fileparts(mfilename('fullpath')), 'last_run_all.txt');
    fid = fopen(logFile, 'w');
    if fid > 0
        c = onCleanup(@() fclose(fid));
        fprintf(fid, '=== run_all_lsqvce log ===\n');
        fprintf(fid, 'when   : %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS')); %#ok<TNOW1,DATST>
        fprintf(fid, 'total  : %.1f s wall clock\n', totalSecs);
        fprintf(fid, 'result : %d / %d passed, %d failed\n', ...
                grandP, grandT, grandF);
        if ~isempty(skipped)
            % On the HEADLINE, same reason as the load errors: the first
            % line anyone reads must say the run was not the full set.
            fprintf(fid, 'result : %d SUITE(S) SKIPPED ON REQUEST: %s\n', ...
                    numel(skipped), strjoin(skipped, ', '));
        end
        if isempty(loadErrors)
            fprintf(fid, '\n');
        else
            % Said on the HEADLINE, because that is the line anyone reads first.
            fprintf(fid, ['result : %d FILE(S) FAILED TO LOAD - the count ' ...
                          'above EXCLUDES them\n'], numel(loadErrors));
            for i = 1:numel(loadErrors)
                fprintf(fid, '         %s\n', loadErrors{i});
            end
            fprintf(fid, '\n');
        end
        for i = 1:numel(results)
            fprintf(fid, '  %-20s : %3d/%-3d pass, %d fail, %d incomplete  %6.1fs\n', ...
                    results(i).name, results(i).passed, results(i).total, ...
                    results(i).failed, results(i).incomplete, results(i).secs);
        end
        for i = 1:numel(capture)
            fprintf(fid, '%s', capture{i});
        end
        fprintf('full log: %s\n', logFile);
    end
end
