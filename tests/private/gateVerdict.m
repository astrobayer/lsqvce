function verdict = gateVerdict(results, skipped, onDisk, excluded, logFile)
%GATEVERDICT  Decide green/red for one full-suite run. Pure logic, no I/O.
%
%   verdict = gateVerdict(results, skipped, onDisk, excluded, logFile)
%
%   The ONE place the CI verdict is decided (rule 3). CI_GATE gathers the
%   inputs (runs the suite, globs the disk) and this function judges them;
%   keeping the judgment PURE makes it unit-testable without recursively
%   re-running the suite (test_ci_gate.m feeds it synthetic runs).
%
%   A run is GREEN if and only if ALL of:
%     1. no suite failed to LOAD (a NaN total in results - §0a-2: a load
%        error hides an unknown number of tests and must trip the gate),
%     2. 0 failed test points,
%     3. 0 incomplete test points,
%     4. suites RUN + suites SKIPPED + documented EXCLUSIONS account for
%        every test_*.m file ON DISK, and nothing was run/skipped that is
%        not on disk (coverage in BOTH directions: a forgotten new suite
%        file and a stale runner entry both trip),
%     5. the log file exists and is non-empty (the artifact the CI uploads
%        and the postflight quotes must never be missing).
%
%   INPUT
%     results   1xF struct   as returned by RUN_ALL_LSQVCE: fields .name
%                            (char), .passed, .total, .failed, .incomplete,
%                            .secs (1x1 double each; .total is NaN for a
%                            suite that failed to load).
%     skipped   1xS cell     suite names (char, no extension) that were
%                            skipped ON REQUEST (the skip argument of
%                            RUN_ALL_LSQVCE). May be {}.
%     onDisk    1xD cell     suite names (char, no extension) of every
%                            test_*.m file found on disk in tests/.
%     excluded  1xE cell     suite names DOCUMENTED as never gate-relevant
%                            (currently the benchmark test_speed). May be {}.
%     logFile   1x1 char     absolute path of the run log
%                            (tests/last_run_all.txt).
%
%   OUTPUT
%     verdict   1x1 struct   .ok       1x1 logical  true = green
%                            .reasons  1xR cell     one char row per red
%                                                   reason; {} when green
%                            .summary  1x1 char     one headline line, e.g.
%                                      'GREEN: 512/512 passed, 19/19 suites
%                                       accounted for (1 skipped, 1 excluded)'
%
%   ERRORS
%     gateVerdict:badInput   any input has the wrong type or shape.
%
%   EXAMPLE
%     r = struct('name','test_util','passed',41,'total',41, ...
%                'failed',0,'incomplete',0,'secs',1.2);
%     v = gateVerdict(r, {}, {'test_util','test_speed'}, {'test_speed'}, f);
%     assert(v.ok)
%
%   See also CI_GATE, RUN_ALL_LSQVCE
%
%   Created 2026-08-18 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.

    if ~isstruct(results) || ~all(isfield(results, ...
            {'name','passed','total','failed','incomplete'}))
        error('gateVerdict:badInput', ...
            'results must be the struct array returned by run_all_lsqvce.');
    end
    skipped  = local_cellstr(skipped,  'skipped');
    onDisk   = local_cellstr(onDisk,   'onDisk');
    excluded = local_cellstr(excluded, 'excluded');
    if ~(ischar(logFile) || (isstring(logFile) && isscalar(logFile)))
        error('gateVerdict:badInput', 'logFile must be a char row or string.');
    end
    logFile = char(logFile);

    reasons = {};

    % 1. load errors: run_all records them with total = NaN.
    loadFail = arrayfun(@(r) any(isnan(r.total)), results);
    if any(loadFail)
        reasons{end+1} = sprintf( ...
            'LOAD ERROR in: %s - nothing in those files ran.', ...
            strjoin({results(loadFail).name}, ', '));
    end

    % 2./3. failures and incompletes (NaN-safe: load errors counted above).
    nFail = sum([results.failed],     'omitnan');
    nInc  = sum([results.incomplete], 'omitnan');
    if nFail > 0
        reasons{end+1} = sprintf('%d test point(s) FAILED.', nFail);
    end
    if nInc > 0
        reasons{end+1} = sprintf('%d test point(s) INCOMPLETE.', nInc);
    end

    % 4. coverage, both directions. A count of failures is not a statement
    %    about coverage (§0a-2), so the set identity is checked exactly.
    ranNames  = {results.name};
    accounted = [ranNames, skipped, excluded];
    missing   = setdiff(onDisk, accounted);        % on disk, never touched
    phantom   = setdiff(ranNames, onDisk);         % ran, but not on disk
    if ~isempty(missing)
        reasons{end+1} = sprintf( ...
            'suite file(s) on disk but NOT run/skipped/excluded: %s.', ...
            strjoin(missing, ', '));
    end
    if ~isempty(phantom)
        reasons{end+1} = sprintf( ...
            'suite(s) ran that have no test_*.m file on disk: %s.', ...
            strjoin(phantom, ', '));
    end

    % 5. the log artifact.
    d = dir(logFile);
    if isempty(d) || d.bytes == 0
        reasons{end+1} = sprintf('log file missing or empty: %s.', logFile);
    end

    verdict.ok      = isempty(reasons);
    verdict.reasons = reasons;
    nP = sum([results.passed], 'omitnan');
    nT = sum([results.total],  'omitnan');
    if verdict.ok
        verdict.summary = sprintf( ...
            'GREEN: %d/%d passed, %d/%d suites accounted for (%d skipped, %d excluded)', ...
            nP, nT, numel(unique(accounted)), numel(onDisk), ...
            numel(skipped), numel(excluded));
    else
        verdict.summary = sprintf('RED: %s', strjoin(reasons, ' | '));
    end
end

function c = local_cellstr(c, name)
%LOCAL_CELLSTR  Normalise to a cell array of char rows or error.
    if isstring(c), c = cellstr(c); end
    if isempty(c), c = {}; end
    if ~iscellstr(c)                                           %#ok<ISCLSTR>
        error('gateVerdict:badInput', ...
            '%s must be a cell array of char or a string array.', name);
    end
    c = reshape(c, 1, []);
end
