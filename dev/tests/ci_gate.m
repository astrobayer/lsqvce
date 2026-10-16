function ci_gate()
%CI_GATE  The ONE CI entry point: full suite + coverage + log, or error.
%
%   ci_gate()
%
%   Called by the CI (`matlab-actions/run-command`: `cd('dev/tests'); ci_gate`
%   on the self-hosted Windows runner, tier B) and interactively over the
%   MATLAB MCP bridge / `matlab -batch` (tier A). It is
%   the gate the D8 loop protects `main` with - NEVER a bare `runtests`:
%   a bare runtests reports failures but says nothing about coverage, and a
%   suite that fails to LOAD would pass unseen (§0a-2, the 2026-07-19
%   test_util lesson).
%
%   WHAT IT DOES
%     1. Runs RUN_ALL_LSQVCE on the toolbox root (the folder above tests/).
%     2. Under CI (environment variable LSQVCE_CI=1):
%        - runs EVERY suite, the pool suite test_parallel included (R6
%          WP-103: the self-hosted runner is the user's PC and carries the
%          Parallel Computing Toolbox; the hosted-runner skip and its
%          2-worker pool cap are gone, so the CI run IS the full gate on
%          the target platform -- one gate, no postflight debt);
%        - creates figures invisible (MATLAB -batch is display-less; this
%          guards tests that open figures).
%        A skip, should one ever be reinstated here, is PRINTED by run_all
%        and recorded in the log, never silent.
%     3. Globs dev/tests/test_*.m and hands run + skip + disk + the
%        documented benchmark exclusion (test_speed) + the time the run
%        started (WP-145: the log must be THIS run's) to GATEVERDICT.
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
%                   nor skipped/excluded, a phantom suite ran, the log
%                   was not written by this run, a new warning identifier
%                   was inventoried or a probe suite's capture is blind.
%                   The message lists every reason.
%
%   EXAMPLE
%     % local postflight, exactly what the CI runs:
%     %   matlab -batch "cd('dev/tests'); ci_gate"
%     ci_gate();
%
%   See also RUN_ALL_LSQVCE, GATEVERDICT, LSQVCE_SETUP
%
%   Updated 2026-10-16 (CEST, UTC+02:00), WP-145: the run's start time goes
%   to GATEVERDICT, so a stale log cannot pass as this run's. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-09-09 (CEST, UTC+02:00), R6 WP-103: the CI skip of
%   test_parallel and the 2-worker pool cap removed -- the self-hosted
%   runner has the PCT (user, 2026-09-09). With assistance from Claude
%   (Anthropic), model Claude Fable 5.1.
%   Created 2026-08-18 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.

    here = fileparts(mfilename('fullpath'));       % .../dev/tests (L8, WP-69)
    root = fileparts(fileparts(here));             % toolbox root
    % L8 (WP-69): temporary path, restored when the gate returns -- the CI
    % invocation is `cd('dev/tests'); ci_gate`, no lsqvce_setup involved.
    pathBefore  = path();
    restorePath = onCleanup(@() path(pathBefore)); %#ok<NASGU>

    % --- CI policy (the ONE place it is decided) --------------------------
    % R6 WP-103: nothing is skipped and no pool is capped -- the self-hosted
    % runner is the target platform with the PCT, so the CI run and the
    % interactive run judge the SAME suites (the hosted-runner skip of
    % test_parallel lived here 2026-08-18 .. 2026-09-09).
    isCI = strcmp(getenv('LSQVCE_CI'), '1');
    skip = {};
    if isCI
        set(groot, 'DefaultFigureVisible', 'off');  % -batch is display-less
    end

    % --- run --------------------------------------------------------------
    runStart = datenum(datetime('now'));           % WP-145 (A9): rule 5's reference
    results = run_all_lsqvce(root, skip);

    % --- gather the coverage evidence ------------------------------------
    d = dir(fullfile(here, 'test_*.m'));
    onDisk = cellfun(@(f) f(1:end-2), {d.name}, 'UniformOutput', false);
    excluded = {'test_speed'};      % benchmark, not pass/fail (see the
                                    % note in run_all_lsqvce's suite list).
    logFile = fullfile(here, 'last_run_all.txt');

    % --- judge ------------------------------------------------------------
    v = gateVerdict(results, skip, onDisk, excluded, logFile, runStart);
    fprintf('\nci_gate: %s\n', v.summary);
    if ~v.ok
        error('ci_gate:red', 'ci_gate RED:\n  %s', ...
              strjoin(v.reasons, sprintf('\n  ')));
    end
end
