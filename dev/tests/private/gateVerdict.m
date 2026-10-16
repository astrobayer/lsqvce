function verdict = gateVerdict(results, skipped, onDisk, excluded, logFile, runStart)
%GATEVERDICT  Decide green/red for one full-suite run. Pure logic, no I/O.
%
%   verdict = gateVerdict(results, skipped, onDisk, excluded, logFile)
%   verdict = gateVerdict(results, skipped, onDisk, excluded, logFile, runStart)
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
%        and the postflight quotes must never be missing), and -- when
%        runStart is given -- it was written by THIS run: its modification
%        time is not before runStart (WP-145, A9: run_all_lsqvce skips the
%        log when fopen fails, so a read-only log, or one held open by an
%        editor on Windows, used to pass as this run's; reproduced),
%     6. every suite EXECUTED exactly the points it DECLARES (WP-47,
%        results.declared vs .total; a suite that runs fewer points than
%        it declares is hiding tests),
%     7. the warning inventory holds NO identifier other than
%        lsqvce:internal:testProbe, and every PROBE SUITE that ran must
%        have reported it: test_ci_gate (the suite capture) and, since
%        WP-145, test_examples (the examples' inner capture, which until
%        then discarded every warning an example raised). An empty
%        inventory there means the capture is blind, not that the code is
%        warning-free. Rules 6/7 apply only when the fields exist --
%        pre-WP-47 result structs stay judgeable.
%
%   INPUT
%     results   1xF struct   as returned by RUN_ALL_LSQVCE: fields .name
%                            (char), .passed, .total, .failed, .incomplete,
%                            .secs (1x1 double each; .total is NaN for a
%                            suite that failed to load); since WP-47 also
%                            .declared (1x1 double, NaN on load failure)
%                            and .warnIds (1xW cell of char identifiers) --
%                            both OPTIONAL, rules 6/7 engage when present.
%     skipped   1xS cell     suite names (char, no extension) that were
%                            skipped ON REQUEST (the skip argument of
%                            RUN_ALL_LSQVCE). May be {}.
%     onDisk    1xD cell     suite names (char, no extension) of every
%                            test_*.m file found on disk in tests/.
%     excluded  1xE cell     suite names DOCUMENTED as never gate-relevant
%                            (currently the benchmark test_speed). May be {}.
%     logFile   1x1 char     absolute path of the run log
%                            (dev/tests/last_run_all.txt).
%     runStart  1x1 double   OPTIONAL: the datenum at which the run began
%                            (CI_GATE takes it before it calls the runner).
%                            Omitted or [], rule 5 checks presence only.
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
%   See also CI_GATE, RUN_ALL_LSQVCE, ECHOWARNINGS
%
%   Updated 2026-10-16 (CEST, UTC+02:00), WP-145: rule 5 binds the log to
%   this run (runStart), rule 7 requires the probe of test_examples too.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-08-20 (CEST, UTC+02:00), WP-47: rules 6 (declared vs
%   executed, external handover 2.6.4) and 7 (single-allowed-identifier
%   warning gate + capture-health probe, 2.5). With assistance from
%   Claude (Anthropic), model Claude Fable 5.
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
    if nargin < 6, runStart = []; end
    if ~(isempty(runStart) || (isnumeric(runStart) && isscalar(runStart) ...
            && isfinite(runStart)))
        error('gateVerdict:badInput', 'runStart must be a finite datenum scalar or [].');
    end

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
    elseif ~isempty(runStart) && d.datenum < runStart - 2/86400
        % 2 s of slack for a file system that stores coarse times (FAT);
        % a full run takes minutes, so a stale log is never inside it
        reasons{end+1} = sprintf(['log file was NOT written by this run: %s ' ...
            '(last written %s, the run started %s) -- read-only or held open?'], ...
            logFile, datestr(d.datenum, 'yyyy-mm-dd HH:MM:SS'), ...
            datestr(runStart, 'yyyy-mm-dd HH:MM:SS')); %#ok<DATST>
    end

    % 6. declared vs executed, per suite (WP-47, external handover 2.6.4):
    %    a suite that EXECUTES fewer points than it DECLARES is hiding
    %    tests, and the pass count alone cannot see it. Optional field for
    %    backward compatibility with pre-WP-47 result structs.
    if isfield(results, 'declared')
        for k = 1:numel(results)
            if ~isnan(results(k).total) && ~isnan(results(k).declared) ...
                    && results(k).declared ~= results(k).total
                reasons{end+1} = sprintf( ...
                    '%s DECLARES %d test point(s) but RAN %d.', ...
                    results(k).name, results(k).declared, ...
                    results(k).total); %#ok<AGROW>
            end
        end
    end

    % 7. the warning inventory (WP-47, external handover 2.5): EXACTLY ONE
    %    identifier may appear in a clean run, lsqvce:internal:testProbe --
    %    no allow-list to maintain, which is the point. And when
    %    test_ci_gate ran, the probe MUST have been seen: an empty
    %    inventory there means the capture plumbing is broken, which is
    %    indistinguishable from warning-free code otherwise.
    if isfield(results, 'warnIds')
        allW = unique([results.warnIds]);
        foreign = setdiff(allW, {'lsqvce:internal:testProbe'});
        if ~isempty(foreign)
            reasons{end+1} = sprintf( ...
                'new warning identifier(s) in the inventory: %s.', ...
                strjoin(foreign, ', '));
        end
        for probeSuite = {'test_ci_gate', 'test_examples'}
            ci = strcmp({results.name}, probeSuite{1});
            if any(ci) && ~isnan(results(ci).total) ...
                    && ~ismember('lsqvce:internal:testProbe', results(ci).warnIds)
                reasons{end+1} = sprintf(['warning capture is BLIND: %s ran ' ...
                    'but its probe identifier never reached the inventory.'], ...
                    probeSuite{1}); %#ok<AGROW>
            end
        end
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
