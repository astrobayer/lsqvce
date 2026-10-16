function tests = test_ci_gate
%TEST_CI_GATE  Contract tests for the CI verdict (tests/private/gateVerdict).
%
%   The verdict is PURE logic, so it is tested here with synthetic run
%   records - the suite never re-runs itself recursively. CI_GATE itself is
%   only the gatherer (run + glob + log path); its judgment is entirely this
%   contract. §6.1 fail-on-old-code: on the pre-WP-30 tree neither
%   gateVerdict nor ci_gate exists, so every test point here fails there.
%
%   Run:
%     results = runtests('tests/test_ci_gate');  disp(table(results))
%
%   See also CI_GATE, GATEVERDICT, RUN_ALL_LSQVCE, ECHOWARNINGS
%
%   Created 2026-08-18 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.
%   Updated 2026-10-16 (CEST, UTC+02:00), WP-145: rule 5 bound to the run's
%   start, the probe of test_examples, and the lock on uncaught captures.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.
    tests = functiontests(localfunctions);
end

%% ---------------------------------------------------------------- helpers
function r = greenRun()
%GREENRUN  Two healthy suites, the minimal green fixture.
    r(1) = struct('name','test_a','passed',10,'total',10, ...
                  'failed',0,'incomplete',0,'secs',1.0);
    r(2) = struct('name','test_b','passed', 5,'total', 5, ...
                  'failed',0,'incomplete',0,'secs',0.5);
end

function f = liveLog(tc)
%LIVELOG  A temporary, non-empty log file, deleted on teardown.
    f = [tempname '.txt'];
    fid = fopen(f, 'w');  fprintf(fid, 'log\n');  fclose(fid);
    tc.addTeardown(@() delete(f));
end

%% ------------------------------------------------------------------ green
function testGreenRunIsGreen(tc)
% categories: robustness | tags: diagnostics
    v = gateVerdict(greenRun(), {}, {'test_a','test_b'}, {}, liveLog(tc));
    tc.verifyTrue(v.ok);
    tc.verifyEmpty(v.reasons);
    tc.verifySubstring(v.summary, 'GREEN');
end

function testSkippedAndExcludedAccountedFor(tc)
% categories: contract | tags: diagnostics
% The CI case: pool suite skipped, benchmark excluded - still green, and the
% summary says so (a skip must be visible in the verdict, never silent).
    onDisk = {'test_a','test_b','test_parallel','test_speed'};
    v = gateVerdict(greenRun(), {'test_parallel'}, onDisk, ...
                    {'test_speed'}, liveLog(tc));
    tc.verifyTrue(v.ok);
    tc.verifySubstring(v.summary, '1 skipped');
    tc.verifySubstring(v.summary, '1 excluded');
end

%% -------------------------------------------------------------------- red
function testFailedTestPointIsRed(tc)
% categories: contract | tags: diagnostics
    r = greenRun();  r(2).failed = 1;  r(2).passed = 4;
    v = gateVerdict(r, {}, {'test_a','test_b'}, {}, liveLog(tc));
    tc.verifyFalse(v.ok);
    tc.verifySubstring(strjoin(v.reasons, ' '), 'FAILED');
end

function testIncompleteTestPointIsRed(tc)
% categories: contract | tags: diagnostics
    r = greenRun();  r(1).incomplete = 2;
    v = gateVerdict(r, {}, {'test_a','test_b'}, {}, liveLog(tc));
    tc.verifyFalse(v.ok);
    tc.verifySubstring(strjoin(v.reasons, ' '), 'INCOMPLETE');
end

function testLoadErrorIsRed(tc)
% categories: robustness | tags: diagnostics
% §0a-2: run_all records a suite that failed to LOAD with total = NaN. That
% must trip the gate even though the failed-count is 0 (the 2026-07-19
% test_util lesson: 41 hidden tests behind a green headline).
    r = greenRun();
    r(3) = struct('name','test_c','passed',0,'total',NaN, ...
                  'failed',NaN,'incomplete',NaN,'secs',NaN);
    v = gateVerdict(r, {}, {'test_a','test_b','test_c'}, {}, liveLog(tc));
    tc.verifyFalse(v.ok);
    tc.verifySubstring(strjoin(v.reasons, ' '), 'LOAD ERROR');
    tc.verifySubstring(strjoin(v.reasons, ' '), 'test_c');
end

function testForgottenSuiteFileIsRed(tc)
% categories: contract | tags: diagnostics
% A new test_*.m on disk that the runner list does not know: the exact
% regression the coverage identity exists for.
    v = gateVerdict(greenRun(), {}, {'test_a','test_b','test_new'}, {}, ...
                    liveLog(tc));
    tc.verifyFalse(v.ok);
    tc.verifySubstring(strjoin(v.reasons, ' '), 'test_new');
end

function testPhantomSuiteIsRed(tc)
% categories: contract | tags: diagnostics
% The other direction: the runner ran something that has no file on disk
% (a stale list entry after a rename would look like this).
    v = gateVerdict(greenRun(), {}, {'test_a'}, {}, liveLog(tc));
    tc.verifyFalse(v.ok);
    tc.verifySubstring(strjoin(v.reasons, ' '), 'test_b');
end

function testMissingLogIsRed(tc)
% categories: contract | tags: diagnostics
    ghost = fullfile(tempdir, 'no_such_lsqvce_log_file.txt');
    v = gateVerdict(greenRun(), {}, {'test_a','test_b'}, {}, ghost);
    tc.verifyFalse(v.ok);
    tc.verifySubstring(strjoin(v.reasons, ' '), 'log file');
end

function testEmptyLogIsRed(tc)
% categories: robustness | tags: diagnostics
    f = [tempname '.txt'];  fid = fopen(f, 'w');  fclose(fid);
    tc.addTeardown(@() delete(f));
    v = gateVerdict(greenRun(), {}, {'test_a','test_b'}, {}, f);
    tc.verifyFalse(v.ok);
end

function testAllReasonsAreListed(tc)
% categories: contract | tags: diagnostics
% Red causes must not shadow each other: a failure AND a coverage gap AND a
% missing log yield three reasons, so the CI log names everything at once.
    r = greenRun();  r(1).failed = 3;  r(1).passed = 7;
    ghost = fullfile(tempdir, 'no_such_lsqvce_log_file.txt');
    v = gateVerdict(r, {}, {'test_a','test_b','test_new'}, {}, ghost);
    tc.verifyFalse(v.ok);
    tc.verifyNumElements(v.reasons, 3);
end

%% ------------------------------------------------------------- validation
function testBadResultsErrors(tc)
% categories: contract | tags: diagnostics
    tc.verifyError(@() gateVerdict(42, {}, {}, {}, 'x'), ...
                   'gateVerdict:badInput');
end

function testBadSkipTypeErrors(tc)
% categories: contract | tags: diagnostics
    tc.verifyError(@() gateVerdict(greenRun(), 5, {'test_a','test_b'}, ...
                   {}, 'x'), 'gateVerdict:badInput');
end

function testStringInputsAreAccepted(tc)
% categories: robustness | tags: diagnostics
% String arrays must behave exactly like cellstr (no second way, but the
% MATLAB string type is normalised, not rejected).
    v = gateVerdict(greenRun(), string.empty, ["test_a","test_b"], ...
                    string.empty, string(liveLog(tc)));
    tc.verifyTrue(v.ok);
end

%% ------------------------------------------- runner mechanism (skip path)
function testRunAllRejectsBadSkip(tc)
% categories: contract
% The mechanism half in run_all_lsqvce: a malformed skip must error before
% anything runs, with the documented identifier.
    tc.verifyError(@() run_all_lsqvce([], 42), 'run_all_lsqvce:badSkip');
end

%% -------------------------------------------------- WP-47: declared count
function testDeclaredMismatchIsRed(tc)
% categories: contract, robustness | tags: diagnostics
% WP-47 rule 6 (external handover 2.6.4): a suite that EXECUTES fewer points
% than it DECLARES is hiding tests and must trip the gate -- and the reason
% must name both numbers, because "which points vanished" is the first
% question the reader asks.
    r = greenRun47();
    r(2).declared = 7;                       % declares 7, ran 5
    v = gateVerdict(r, {}, {'test_a','test_b'}, {}, liveLog(tc));
    tc.verifyFalse(v.ok);
    tc.verifyTrue(any(contains(v.reasons, 'DECLARES 7')));
    tc.verifyTrue(any(contains(v.reasons, 'RAN 5')));
end

function testLegacyResultsWithoutNewFieldsStayJudgeable(tc)
% categories: contract, robustness
% Pre-WP-47 result structs carry neither .declared nor .warnIds; rules 6/7
% must simply not engage -- backward compatibility is part of the contract.
    v = gateVerdict(greenRun(), {}, {'test_a','test_b'}, {}, liveLog(tc));
    tc.verifyTrue(v.ok);
end

function testDeclaredPointsCounter(tc)
% categories: contract, robustness
% The MATLAB-side static counter (tests/private/declaredPoints): counts test
% signatures, skips the entry point and %{ %} blocks -- the block-comment
% case IS the mismatch a naive grep cannot see. Also the documented
% rejection.
    d = tempname;  mkdir(d);  tc.addTeardown(@() rmdir(d, 's'));
    f = fullfile(d, 'test_fixture.m');
    fid = fopen(f, 'w');
    fprintf(fid, [ ...
        'function tests = test_fixture\n' ...
        'tests = functiontests(localfunctions);\n' ...
        'end\n' ...
        'function testAlpha(tc)\n', 'end\n' ...
        '%%{\n' ...
        'function testCommentedOut(tc)\n', 'end\n' ...
        '%%}\n' ...
        'function testBeta(tc)\n', 'end\n']);
    fclose(fid);
    tc.verifyEqual(declaredPoints(f), 2, ...
        'the block-commented test must NOT be counted');
    caught = MException.empty;
    try
        declaredPoints(fullfile(d, 'nope.m'));
    catch e
        caught = e;
    end
    tc.verifyNotEmpty(caught);
    tc.verifyEqual(caught.identifier, 'lsqvce:declaredPoints:noFile');
end

%% ------------------------------------------------ WP-47: warning inventory
function testForeignWarningIdentifierIsRed(tc)
% categories: contract, robustness | tags: diagnostics
% WP-47 rule 7 (external handover 2.5): exactly ONE identifier may appear in
% a clean run's inventory. Any other is new and fails the gate -- no
% allow-list to maintain, which is the point.
    r = greenRun47();
    r(1).warnIds = {'lsqvce:foo:surprise'};
    v = gateVerdict(r, {}, {'test_a','test_b'}, {}, liveLog(tc));
    tc.verifyFalse(v.ok);
    tc.verifyTrue(any(contains(v.reasons, 'lsqvce:foo:surprise')));
end

function testUnidentifiedWarningIsRed(tc)
% categories: contract, robustness
% A warning line without an identifier is inventoried as '<unidentified>'
% so it can never hide behind the parse.
    r = greenRun47();
    r(2).warnIds = {'<unidentified>'};
    v = gateVerdict(r, {}, {'test_a','test_b'}, {}, liveLog(tc));
    tc.verifyFalse(v.ok);
end

function testProbeAllowedAndRequiredWhenCiGateRan(tc)
% categories: contract, robustness | tags: diagnostics
% The single allowed identifier passes; and when test_ci_gate itself ran,
% its probe MUST have reached the inventory -- an empty inventory there
% means the capture plumbing is BLIND, which the pass count cannot see.
    r = greenRun47();
    r(3) = struct('name','test_ci_gate','passed',3,'total',3, ...
                  'failed',0,'incomplete',0,'secs',0.2, ...
                  'declared',3,'warnIds',{{'lsqvce:internal:testProbe'}});
    v = gateVerdict(r, {}, {'test_a','test_b','test_ci_gate'}, {}, liveLog(tc));
    tc.verifyTrue(v.ok, 'the probe identifier alone must stay green');

    r(3).warnIds = {};                       % capture went blind
    v = gateVerdict(r, {}, {'test_a','test_b','test_ci_gate'}, {}, liveLog(tc));
    tc.verifyFalse(v.ok);
    tc.verifyTrue(any(contains(v.reasons, 'BLIND')));
end

function testWarningProbeFires(tc)
% categories: robustness | tags: diagnostics
% THE LIVE PROBE (external handover 2.5): this warning is the one identifier
% the gate allows, fired un-suppressed so the run_all capture has something
% to see. If gateVerdict ever reports the capture BLIND, this test is where
% the signal was supposed to enter. It must NOT be silenced.
    warning('lsqvce:internal:testProbe', ...
        'gate probe: this warning proves the inventory capture works');
    tc.verifyTrue(true);
end

%% --------------------------------------- WP-145: rule 5 and rule 7 (A9, A1)
function testStaleLogIsRed(tc)
% categories: contract, robustness | tags: diagnostics
% A9: run_all_lsqvce skips the log when fopen fails, so a read-only log (or
% one an editor holds open on Windows) passed rule 5 as this run's --
% reproduced on main e483142 with the log made read-only. With the run's
% start given, a log older than the run is red; one written after it is not.
    f = liveLog(tc);
    v = gateVerdict(greenRun(), {}, {'test_a','test_b'}, {}, f, datenum(datetime('now')) + 1);
    tc.verifyFalse(v.ok);
    tc.verifyTrue(any(contains(v.reasons, 'NOT written by this run')));
    v = gateVerdict(greenRun(), {}, {'test_a','test_b'}, {}, f, datenum(datetime('now')) - 1);
    tc.verifyTrue(v.ok, 'a log written after the run started is this run''s');
end

function testBadRunStartErrors(tc)
% categories: contract
    tc.verifyError(@() gateVerdict(greenRun(), {}, {'test_a','test_b'}, {}, liveLog(tc), 'now'), ...
                   'gateVerdict:badInput');
end

function testExamplesProbeIsRequired(tc)
% categories: contract, robustness | tags: diagnostics
% A1: test_examples carries a probe since WP-145, because its inner capture
% discarded every warning the examples raised. When it ran without
% reporting the probe, its capture is deaf again -- red, named.
    r = greenRun47();
    r(3) = struct('name','test_examples','passed',3,'total',3, ...
                  'failed',0,'incomplete',0,'secs',20, ...
                  'declared',3,'warnIds',{{'lsqvce:internal:testProbe'}});
    onDisk = {'test_a','test_b','test_examples'};
    tc.verifyTrue(gateVerdict(r, {}, onDisk, {}, liveLog(tc)).ok);
    r(3).warnIds = {};
    v = gateVerdict(r, {}, onDisk, {}, liveLog(tc));
    tc.verifyFalse(v.ok);
    tc.verifyTrue(any(contains(v.reasons, 'BLIND: test_examples')));
end

function testEveryCaptureInTheSuitesEchoesItsWarnings(tc)
% categories: contract | tags: diagnostics
% The LOCK of WP-145 (A1): a suite that captures code with evalc and drops
% the text drops its warnings too, unless the capture goes to echoWarnings
% on the same or the next line. A new capture without it is named here.
    here = fileparts(mfilename('fullpath'));
    d = dir(fullfile(here, 'test_*.m'));
    call = ['evalc' '('];                       % split: this file is scanned too
    bad = {};
    for k = 1:numel(d)
        L = splitlines(fileread(fullfile(here, d(k).name)));
        for i = 1:numel(L)
            if startsWith(strtrim(L{i}), '%') || ~contains(L{i}, call), continue, end
            nxt = '';
            if i < numel(L), nxt = L{i+1}; end
            if ~contains([L{i} nxt], 'echoWarnings(')
                bad{end+1} = sprintf('%s:%d', d(k).name, i); %#ok<AGROW>
            end
        end
    end
    tc.verifyEmpty(bad, sprintf('captures whose warnings are dropped: %s', strjoin(bad, ', ')));
end

function r = greenRun47()
%GREENRUN47  greenRun with the WP-47 fields present and healthy.
    r = greenRun();
    [r.declared] = deal(10, 5);
    [r.warnIds]  = deal({}, {});
end
