function ci_gate()
%CI_GATE  The ONE CI entry point: full suite + coverage + log, or error.
%
%   ci_gate()
%
%   Called by the CI (`matlab-actions/run-command`: `cd('dev/tests'); ci_gate`)
%   and by the local postflight (MATLAB MCP bridge / `matlab -batch`). It is
%   the gate the D8 loop protects `main` with - NEVER a bare `runtests`:
%   a bare runtests reports failures but says nothing about coverage, and a
%   suite that fails to LOAD would pass unseen (§0a-2, the 2026-07-19
%   test_util lesson).
%
%   WHAT IT DOES
%     1. Runs RUN_ALL_LSQVCE on the toolbox root (the folder above tests/).
%     2. Under CI (environment variable LSQVCE_CI=1):
%        - skips the pool suite test_parallel (hosted runners have no
%          Parallel Computing Toolbox licence until `products:` requests it
%          in ci.yml) - the skip is PRINTED by run_all and recorded in the
%          log, never silent;
%        - caps any parallel pool at 2 workers should one be opened anyway
%          (2 = the hosted-runner core count; more workers only queue);
%        - creates figures invisible (headless Linux runner; MATLAB -batch
%          is display-less, this guards tests that open figures).
%     3. Globs tests/test_*.m and hands run + skip + disk + the documented
%        benchmark exclusion (test_speed) to GATEVERDICT.
%     4. Prints the verdict; on red it ERRORS so the CI job (and any
%        -batch invocation) exits non-zero.
%
%   INPUT
%     (none)      policy lives here, mechanism in RUN_ALL_LSQVCE: the set
%                 of CI skips is decided in THIS file only (rule 3).
%
%   OUTPUT
%     (none)      side effects: console summary, dev/tests/last_run_all.txt
%                 (written by RUN_ALL_LSQVCE, uploaded as the CI artifact),
%                 and an error on red.
%
%   ERRORS
%     ci_gate:red   the run is not green: >0 failed, >0 incomplete, a suite
%                   failed to load, a test_*.m file on disk was neither run
%                   nor skipped/excluded, a phantom suite ran, or the log
%                   was not written. The message lists every reason.
%
%   EXAMPLE
%     % local postflight, exactly what the CI runs:
%     %   matlab -batch "cd('dev/tests'); ci_gate"
%     ci_gate();
%
%   See also RUN_ALL_LSQVCE, GATEVERDICT, LSQVCE_SETUP
%
%   Created 2026-08-18 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.

    here = fileparts(mfilename('fullpath'));       % .../dev/tests (L8, WP-69)
    root = fileparts(fileparts(here));             % toolbox root
    % L8 (WP-69): temporary path, restored when the gate returns -- the CI
    % invocation is `cd('dev/tests'); ci_gate`, no lsqvce_setup involved.
    pathBefore  = path();
    restorePath = onCleanup(@() path(pathBefore)); %#ok<NASGU>

    % --- CI policy (the ONE place it is decided) --------------------------
    isCI = strcmp(getenv('LSQVCE_CI'), '1');
    skip = {};
    if isCI
        skip = {'test_parallel'};   % no PCT licence on the hosted runner;
                                    % enable products: in ci.yml to lift.
        try                          % cap a pool if one is opened anyway
            c = parcluster('Processes');
            c.NumWorkers = min(c.NumWorkers, 2);
        catch
            % no PCT at all - nothing to cap.
        end
        set(groot, 'DefaultFigureVisible', 'off');  % headless runner
    end

    % --- run --------------------------------------------------------------
    results = run_all_lsqvce(root, skip);

    % --- gather the coverage evidence ------------------------------------
    d = dir(fullfile(here, 'test_*.m'));
    onDisk = cellfun(@(f) f(1:end-2), {d.name}, 'UniformOutput', false);
    excluded = {'test_speed'};      % benchmark, not pass/fail (see the
                                    % note in run_all_lsqvce's suite list).
    logFile = fullfile(here, 'last_run_all.txt');

    % --- judge ------------------------------------------------------------
    v = gateVerdict(results, skip, onDisk, excluded, logFile);
    fprintf('\nci_gate: %s\n', v.summary);
    if ~v.ok
        error('ci_gate:red', 'ci_gate RED:\n  %s', ...
              strjoin(v.reasons, sprintf('\n  ')));
    end
end
