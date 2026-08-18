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
%   See also CI_GATE, GATEVERDICT, RUN_ALL_LSQVCE
%
%   Created 2026-08-18 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.
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
    v = gateVerdict(greenRun(), {}, {'test_a','test_b'}, {}, liveLog(tc));
    tc.verifyTrue(v.ok);
    tc.verifyEmpty(v.reasons);
    tc.verifySubstring(v.summary, 'GREEN');
end

function testSkippedAndExcludedAccountedFor(tc)
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
    r = greenRun();  r(2).failed = 1;  r(2).passed = 4;
    v = gateVerdict(r, {}, {'test_a','test_b'}, {}, liveLog(tc));
    tc.verifyFalse(v.ok);
    tc.verifySubstring(strjoin(v.reasons, ' '), 'FAILED');
end

function testIncompleteTestPointIsRed(tc)
    r = greenRun();  r(1).incomplete = 2;
    v = gateVerdict(r, {}, {'test_a','test_b'}, {}, liveLog(tc));
    tc.verifyFalse(v.ok);
    tc.verifySubstring(strjoin(v.reasons, ' '), 'INCOMPLETE');
end

function testLoadErrorIsRed(tc)
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
% A new test_*.m on disk that the runner list does not know: the exact
% regression the coverage identity exists for.
    v = gateVerdict(greenRun(), {}, {'test_a','test_b','test_new'}, {}, ...
                    liveLog(tc));
    tc.verifyFalse(v.ok);
    tc.verifySubstring(strjoin(v.reasons, ' '), 'test_new');
end

function testPhantomSuiteIsRed(tc)
% The other direction: the runner ran something that has no file on disk
% (a stale list entry after a rename would look like this).
    v = gateVerdict(greenRun(), {}, {'test_a'}, {}, liveLog(tc));
    tc.verifyFalse(v.ok);
    tc.verifySubstring(strjoin(v.reasons, ' '), 'test_b');
end

function testMissingLogIsRed(tc)
    ghost = fullfile(tempdir, 'no_such_lsqvce_log_file.txt');
    v = gateVerdict(greenRun(), {}, {'test_a','test_b'}, {}, ghost);
    tc.verifyFalse(v.ok);
    tc.verifySubstring(strjoin(v.reasons, ' '), 'log file');
end

function testEmptyLogIsRed(tc)
    f = [tempname '.txt'];  fid = fopen(f, 'w');  fclose(fid);
    tc.addTeardown(@() delete(f));
    v = gateVerdict(greenRun(), {}, {'test_a','test_b'}, {}, f);
    tc.verifyFalse(v.ok);
end

function testAllReasonsAreListed(tc)
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
    tc.verifyError(@() gateVerdict(42, {}, {}, {}, 'x'), ...
                   'gateVerdict:badInput');
end

function testBadSkipTypeErrors(tc)
    tc.verifyError(@() gateVerdict(greenRun(), 5, {'test_a','test_b'}, ...
                   {}, 'x'), 'gateVerdict:badInput');
end

function testStringInputsAreAccepted(tc)
% String arrays must behave exactly like cellstr (no second way, but the
% MATLAB string type is normalised, not rejected).
    v = gateVerdict(greenRun(), string.empty, ["test_a","test_b"], ...
                    string.empty, string(liveLog(tc)));
    tc.verifyTrue(v.ok);
end

%% ------------------------------------------- runner mechanism (skip path)
function testRunAllRejectsBadSkip(tc)
% The mechanism half in run_all_lsqvce: a malformed skip must error before
% anything runs, with the documented identifier.
    tc.verifyError(@() run_all_lsqvce([], 42), 'run_all_lsqvce:badSkip');
end
