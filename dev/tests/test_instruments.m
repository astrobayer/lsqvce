function tests = test_instruments
%TEST_INSTRUMENTS  The suite's own instruments: the timing helper timePair.
%   tests = test_instruments
%
%   R5 Level M0 (WP-78): "no budget may be written before its instrument
%   exists" -- and an instrument nobody has tested is a claim. This suite
%   pins the CONTRACT (fields, rejection ids), the PRECISION on the identity
%   pair, the METAMORPHIC properties (swapping sides inverts the ratio; a
%   doubled workload doubles it) and the ROBUSTNESS to the Inner setting of
%   dev/tests/private/timePair.m. Workloads are ~ms so the points cost
%   seconds on a hosted runner and never touch the toolbox numerics.
%
%   Categories per point are declared on the line after each signature
%   (audit check 16). Timing assertions are BANDS, never equalities: a
%   band that a bad minute can leave is a flake, not a test, so every bound
%   here is wide (0.5x .. 2x) and states why.
%
%   Since WP-145 also the two instruments of the capture and the iteration
%   loop: ECHOWARNINGS (an inner evalc's warnings reach the gate) and
%   RUN_LSQVCE (its shortcuts find dev/tests, LogFile '' means no log).
%
%   See also TIMEPAIR, ECHOWARNINGS, RUN_LSQVCE, RUN_ALL_LSQVCE, TEST_SPEED
%
%   Created 2026-09-01 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M0, WP-78.
%   Updated 2026-10-16 (CEST, UTC+02:00), WP-145: echoWarnings and
%   run_lsqvce. With assistance from Claude (Anthropic), model Claude Opus 5.5.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    % a ~2 ms deterministic workload: a matrix product on fixed data
    rs = RandStream('mt19937ar', 'Seed', 7);
    tc.TestData.M = randn(rs, 220, 220);
end

function f = work(tc, times)
    M = tc.TestData.M;
    f = @() workBody(M, times);
end

function workBody(M, times)
    for k = 1:times
        X = M * M;  %#ok<NASGU>
    end
end

%% ============================================================ contract
function testRecordFieldsAndSettings(tc)
% categories: contract
    rec = timePair(work(tc, 1), work(tc, 1), 'Repeats', 5, 'Inner', 2, 'Label', 'idem');
    want = {'ratio','band','tA','tB','ratios','repeats','inner','label','machine'};
    tc.verifyEqual(sort(fieldnames(rec)).', sort(want), 'the record has exactly the documented fields');
    tc.verifyEqual(rec.repeats, 5);
    tc.verifyEqual(rec.inner, 2);
    tc.verifyEqual(rec.label, 'idem');
    tc.verifySize(rec.ratios, [5 1]);
    tc.verifyEqual(rec.band, [min(rec.ratios) max(rec.ratios)], 'band is the per-repeat min..max');
    tc.verifyEqual(rec.ratio, median(rec.ratios), 'the statistic is the median of per-repeat ratios');
    tc.verifyTrue(ischar(rec.machine) && contains(rec.machine, version('-release')));
end

function testRejectsNonHandles(tc)
% categories: contract
    tc.verifyError(@() timePair(1, work(tc, 1)), 'lsqvce:timePair:badHandle');
    tc.verifyError(@() timePair(work(tc, 1), 'x'), 'lsqvce:timePair:badHandle');
end

function testRejectsBadSettings(tc)
% categories: contract
    % typed arguments: Repeats < 3 and Inner < 1 are rejected by the
    % signature (MATLAB validator ids), before any call is timed
    tc.verifyError(@() timePair(work(tc, 1), work(tc, 1), 'Repeats', 2), ...
                   'MATLAB:validators:mustBeGreaterThanOrEqual');
    tc.verifyError(@() timePair(work(tc, 1), work(tc, 1), 'Inner', 0), ...
                   'MATLAB:validators:mustBePositive');
end

%% ============================================================ precision
function testIdentityPairIsNearOne(tc)
% categories: precision
    rec = timePair(work(tc, 1), work(tc, 1), 'Repeats', 15);
    % tol: identical callables; the ratio's median must sit near 1 -- the
    % 0.5..2 band is 5x the measured 10 % spread at 15 repeats and covers a
    % hosted runner's jitter without ever accepting a real 2x asymmetry
    tc.verifyGreaterThanOrEqual(rec.ratio, 0.5, sprintf('identity ratio %.2f band [%.2f %.2f]', rec.ratio, rec.band));
    tc.verifyLessThanOrEqual(rec.ratio, 2.0, sprintf('identity ratio %.2f band [%.2f %.2f]', rec.ratio, rec.band));
    tc.verifyGreaterThan(rec.tA, 0);
    tc.verifyGreaterThan(rec.tB, 0);
end

%% ============================================================ metamorphic
function testSwappingSidesInvertsTheRatio(tc)
% categories: metamorphic
    % eps-level is NOT the claim (timing); the claim is r(A,B) * r(B,A) ~ 1
    a = work(tc, 1);  b = work(tc, 3);
    r1 = timePair(a, b, 'Repeats', 9);
    r2 = timePair(b, a, 'Repeats', 9);
    p  = r1.ratio * r2.ratio;
    % tol: product of two independent medians of the same 1:3 workload pair;
    % each carries ~15 % spread at 9 repeats -> 0.5..2 is >3x that
    tc.verifyGreaterThanOrEqual(p, 0.5, sprintf('r1 %.2f r2 %.2f', r1.ratio, r2.ratio));
    tc.verifyLessThanOrEqual(p, 2.0, sprintf('r1 %.2f r2 %.2f', r1.ratio, r2.ratio));
end

function testTripledWorkloadTriplesTheRatio(tc)
% categories: metamorphic
    r = timePair(work(tc, 3), work(tc, 1), 'Repeats', 9);
    % tol: 3x the work on the numerator; call overhead of a 2 ms body is
    % < 5 %, so the ratio must land in 1.5..6 (a 2x band either side)
    tc.verifyGreaterThanOrEqual(r.ratio, 1.5, sprintf('ratio %.2f band [%.2f %.2f]', r.ratio, r.band));
    tc.verifyLessThanOrEqual(r.ratio, 6.0, sprintf('ratio %.2f band [%.2f %.2f]', r.ratio, r.band));
end

%% ============================================================ robustness
function testInnerBatchingDoesNotMoveTheRatio(tc)
% categories: robustness
    r1 = timePair(work(tc, 2), work(tc, 1), 'Repeats', 9, 'Inner', 1);
    r3 = timePair(work(tc, 2), work(tc, 1), 'Repeats', 9, 'Inner', 3);
    % tol: Inner changes only how many calls sit inside one timing; the
    % ratio is a property of the workloads -- agreement to a factor 2
    tc.verifyGreaterThanOrEqual(r3.ratio / r1.ratio, 0.5, sprintf('inner1 %.2f inner3 %.2f', r1.ratio, r3.ratio));
    tc.verifyLessThanOrEqual(r3.ratio / r1.ratio, 2.0, sprintf('inner1 %.2f inner3 %.2f', r1.ratio, r3.ratio));
end

%% ------------------------------------------- WP-145: echoWarnings (A1)
function testEchoWarningsReemitsWhatAnInnerCaptureSwallowed(tc)
% categories: contract, robustness | tags: diagnostics
    % The defect, then the repair, in one capture-inside-a-capture -- the
    % shape of test_examples inside run_all_lsqvce. Without the echo the
    % identifier never reaches the outer text; with it, both warnings do,
    % framed as MATLAB frames them, and the line that merely STARTS with
    % "Warning:" (demo_messaging replays its log sink like that) does not.
    probe = 'lsqvce:instruments:echoProbe';
    vb = warning('query', 'verbose');  warning('on', 'verbose');
    tc.addTeardown(@() warning(vb.state, 'verbose'));
    raise = sprintf(['warning(''%s'', ''with an id''); warning(''without an id''); ' ...
                     'fprintf(''Warning: printed text, not a warning\\n'');'], probe);
    deaf  = evalc('inner = evalc(raise);');
    echoWarnings(deaf);   % the lock's rule; deaf holds nothing to echo
    tc.verifyFalse(contains(deaf, probe), 'an inner capture alone swallows the warning');
    heard = evalc('inner = evalc(raise); blocks = echoWarnings(inner);');
    tc.verifyTrue(contains(heard, ['warning off ' probe]), ...
        'the identifier hint reaches the outer capture, where the inventory parses it');
    tc.verifyNumElements(blocks, 2, 'both warnings, the one without an id included');
    tc.verifyFalse(contains(heard, 'printed text'), 'printed text is not a warning');
end

function testEchoWarningsRejectsANonText(tc)
% categories: contract
    tc.verifyError(@() echoWarnings(42), 'lsqvce:echoWarnings:badCapture');
    tc.verifyEmpty(echoWarnings(''), 'an empty capture holds no warning');
end

%% --------------------------------------- WP-145: run_lsqvce (A5, A8)
function testRunLsqvceShortcutRunsItsSuite(tc)
% categories: contract | tags: diagnostics
    % A5: every shortcut resolved to the root tests/ folder that L8 removed,
    % so run_lsqvce('streams') ran NOTHING and printed "RUN ERRORED"
    % (reproduced on main e483142). One cheap point through a shortcut, its
    % log in a temporary file, so the session's own logs stay untouched.
    f = [tempname '.txt'];
    tc.addTeardown(@() deleteIfThere(f));
    out = evalc('r = run_lsqvce(''streams'', ''Name'', ''*MaxLagCorrShortSeries*'', ''LogFile'', f);');
    echoWarnings(out);
    tc.verifyNumElements(r, 1, 'the shortcut finds its suite in dev/tests');
    tc.verifyTrue(all([r.Passed]), out);
    tc.verifyTrue(isfile(f), 'the chosen log is written');
end

function testRunLsqvceEmptyLogFileWritesNoLog(tc)
% categories: contract, robustness
    % A8: '' was replaced by the default path, so it could not switch the
    % log off and the no-log branch (the only place Quiet acts) never ran.
    % Not given and given empty are two requests.
    here = fileparts(mfilename('fullpath'));
    dflt = fullfile(here, 'last_run_streams.txt');
    before = dir(dflt);
    out = evalc('r = run_lsqvce(''streams'', ''Name'', ''*MaxLagCorrShortSeries*'', ''LogFile'', '''');');
    echoWarnings(out);
    tc.verifyNumElements(r, 1);
    after = dir(dflt);
    tc.verifyEqual(numel(after), numel(before), 'no default log appears');
    if ~isempty(after)
        tc.verifyEqual(after.datenum, before.datenum, 'the default log is not rewritten');
    end
    tc.verifyFalse(contains(out, 'full log:'), 'the summary names no log');
end

function testRunLsqvceSaysWhenNoTestMatched(tc)
% categories: contract
    % A filter that matches nothing is not an error, and is not reported as one.
    out = evalc('r = run_lsqvce(''streams'', ''Name'', ''*noSuchTestName*'', ''LogFile'', '''');');
    echoWarnings(out);
    tc.verifyEmpty(r);
    tc.verifyTrue(contains(out, 'no test matched'), out);
    tc.verifyFalse(contains(out, 'ERRORED'), out);
end

function deleteIfThere(f)
%DELETEIFTHERE  Teardown helper: a log the test may or may not have written.
    if isfile(f), delete(f); end
end
