function tests = test_util
%TEST_UTIL  Direct tests for the internal helpers lsqvce.util.getfielddef and
%   lsqvce.util.mergeStruct. Both had ~10 internal callers but no assertion on
%   their own behaviour (roadmap item 2). Mirrored first in
%   validation/v8verify.py.
%
%   Run:
%     results = runtests('tests/test_util');  disp(table(results))
%
%   See also LSQVCE.UTIL.GETFIELDDEF, LSQVCE.UTIL.MERGESTRUCT
    tests = functiontests(localfunctions);
end

%% ================================================= getfielddef
function testGetFieldPresent(tc)
    s = struct('arOrder', 30, 'type', 'psd');
    tc.verifyEqual(lsqvce.util.getfielddef(s, 'arOrder', 99), 30);
end

function testGetFieldMissingReturnsDefault(tc)
    s = struct('arOrder', 30);
    tc.verifyEqual(lsqvce.util.getfielddef(s, 'nope', 99), 99);
end

function testGetFieldEmptyIsAbsent(tc)
% The distinction the function exists to make: an empty field means "use the
% default", NOT "use []". This is why 'BatchOverlap',[] selects the safe
% default overlap instead of zero (zero overlap silently destroys accuracy).
    s = struct('BatchOverlap', []);
    tc.verifyEqual(lsqvce.util.getfielddef(s, 'BatchOverlap', 7), 7);
end

function testGetFieldNonStruct(tc)
    tc.verifyEqual(lsqvce.util.getfielddef([], 'x', 5), 5);
    tc.verifyEqual(lsqvce.util.getfielddef('not a struct', 'x', 5), 5);
end

function testGetFieldNonScalarStructUsesDefault(tc)
    s(1).x = 1; s(2).x = 2;                 % 1x2 struct array -> not scalar
    tc.verifyEqual(lsqvce.util.getfielddef(s, 'x', 42), 42);
end

%% ================================================= mergeStruct
function testMergeOverrides(tc)
    a = struct('type', 'psd', 'arOrder', 30);
    b = struct('arOrder', 60);
    c = lsqvce.util.mergeStruct(a, b);
    tc.verifyEqual(c.arOrder, 60);          % overridden
    tc.verifyEqual(c.type, 'psd');          % kept
end

function testMergeAddsNewField(tc)
    a = struct('type', 'psd');
    c = lsqvce.util.mergeStruct(a, struct('sigma', 2.5));
    tc.verifyEqual(c.type, 'psd');
    tc.verifyEqual(c.sigma, 2.5);
end

function testMergeEmptyOverlayIsNoOp(tc)
    a = struct('type', 'psd', 'arOrder', 30);
    tc.verifyEqual(lsqvce.util.mergeStruct(a, struct()), a);
    tc.verifyEqual(lsqvce.util.mergeStruct(a, []), a);
end

function testMergeNonStructBase(tc)
    b = struct('x', 1);
    tc.verifyEqual(lsqvce.util.mergeStruct([], b), b);
end

function testMergeCovOverridePattern(tc)
% The real use: overlay a runtime covariance override onto a static spec.
    covStatic = struct('type', 'diag', 'sigma', 1.0);
    covRun    = struct('sigma', 2.5);
    cm = lsqvce.util.mergeStruct(covStatic, covRun);
    tc.verifyEqual(cm.type, 'diag');
    tc.verifyEqual(cm.sigma, 2.5);
end

%% ================================================= lsqvce.util.log
function testLogDefaultsAreOneAndTrue(tc)
    lsqvce.util.log.reset();
    [fid, screen] = lsqvce.util.log.get();
    tc.verifyEqual(fid, 1);
    tc.verifyTrue(screen);
end

function testLogMsgWritesToFile(tc)
    % msg must write to logfid; with a file fid and screen=false, the text lands
    % in the file and NOT on screen.
    f = [tempname, '.log'];
    fid = fopen(f, 'w');
    c1 = onCleanup(@() cleanupLog(fid, f));
    lsqvce.util.log.set('fid', fid, 'screen', false);
    lsqvce.util.log.msg('hello %d\n', 42);
    lsqvce.util.log.reset();
    fclose(fid);
    txt = fileread(f);
    tc.verifyTrue(contains(txt, 'hello 42'));
end

function testLogWarnRaisesAndMirrors(tc)
    % warn must raise a real MATLAB warning (catchable) AND mirror to the file.
    f = [tempname, '.log'];
    fid = fopen(f, 'w');
    c1 = onCleanup(@() cleanupLog(fid, f));
    lsqvce.util.log.set('fid', fid, 'screen', false);
    ws = warning('off', 'lsqvce:test:logwarn');
    c2 = onCleanup(@() warning(ws));
    lastwarn('');
    lsqvce.util.log.warn('lsqvce:test:logwarn', 'bad value %g', 3.5);
    [~, id] = lastwarn();
    lsqvce.util.log.reset();
    fclose(fid);
    tc.verifyEqual(id, 'lsqvce:test:logwarn');   % real warning was raised
    txt = fileread(f);
    tc.verifyTrue(contains(txt, 'bad value 3.5')); % and mirrored to the log file
end

function testLogSetPartial(tc)
    lsqvce.util.log.reset();
    lsqvce.util.log.set('screen', false);        % set only one field
    [fid, screen] = lsqvce.util.log.get();
    tc.verifyEqual(fid, 1);                       % fid untouched
    tc.verifyFalse(screen);
    lsqvce.util.log.reset();
end

function testLogStartParallelSafeWithoutPool(tc)
    % startParallel must be a safe no-op when the PCT / a pool is unavailable:
    % it returns [] and does NOT disturb the (fid, screen) state. With a pool it
    % returns a DataQueue; either way stopParallel must clean up without error.
    lsqvce.util.log.reset();
    q = lsqvce.util.log.startParallel();
    hasPar = ~isempty(ver('parallel')) && ...
             exist('parallel.pool.DataQueue', 'class') == 8;
    if hasPar
        % A queue may or may not be created depending on whether a pool opens;
        % if one was created it must be the right type.
        if ~isempty(q)
            tc.verifyClass(q, 'parallel.pool.DataQueue');
        end
    else
        tc.verifyEmpty(q, 'no PCT -> startParallel must return []');
    end
    [fid, screen] = lsqvce.util.log.get();        % state untouched
    tc.verifyEqual(fid, 1);
    tc.verifyTrue(screen);
    lsqvce.util.log.stopParallel();               % idempotent, no error
    lsqvce.util.log.stopParallel();
end

function testLogSendToQueueFalseOnClient(tc)
    % On the CLIENT (no worker task), sendToQueue must return false even after a
    % sink is installed, so client-side msg/warn still emit locally as before.
    lsqvce.util.log.reset();
    lsqvce.util.log.startParallel();
    c = onCleanup(@() lsqvce.util.log.stopParallel());
    tc.verifyFalse(lsqvce.util.log.sendToQueue('msg', 'x'));
    % and a client msg still reaches the file
    f = [tempname, '.log'];  fid = fopen(f, 'w');
    c2 = onCleanup(@() cleanupLog(fid, f));
    lsqvce.util.log.set('fid', fid, 'screen', false);
    lsqvce.util.log.msg('client %d\n', 7);
    lsqvce.util.log.reset();
    tc.verifyTrue(contains(fileread(f), 'client 7'));
end

function testLogDrainWritesToFid(tc)
    % The client-side drain callback must write a {kind,txt} worker message to
    % the captured fid, prefixing warnings like warn does.
    f = [tempname, '.log'];  fid = fopen(f, 'w');
    c = onCleanup(@() cleanupLog(fid, f));
    lsqvce.util.log.drain(fid, false, {'msg',  'from worker 3\n'});
    lsqvce.util.log.drain(fid, false, {'warn', 'r_g low'});
    fclose(fid);
    txt = fileread(f);
    tc.verifyTrue(contains(txt, 'from worker 3'));
    tc.verifyTrue(contains(txt, 'Warning: r_g low'));
    lsqvce.util.log.reset();
end

function cleanupLog(fid, f)
    try, if ~isempty(fopen(fid)), fclose(fid); end, catch, end
    lsqvce.util.log.reset();
    if exist(f, 'file'), delete(f); end
end

%% ================================================= showExport
function testShowExportEmptyIsNoOp(tc)
    % A hook that returned [] must produce no file and no error.
    out = lsqvce.util.showExport([], 'nope', 0, tempdir);
    tc.verifyEqual(out, '');
end

function testShowExportInvalidHandleIsNoOp(tc)
    % A non-figure value must be a silent no-op (never break the solve).
    out = lsqvce.util.showExport(42, 'nope', 0, tempdir);
    tc.verifyEqual(out, '');
end

function testShowExportWritesPng(tc)
    % With a real figure and an export dir, a PNG must appear with the sanitised
    % name (the 'satKey|obsType' pipe becomes '_'). Skipped if no display/figure.
    d = fullfile(tempdir, sprintf('lsqvce_figtest_%d', randi(1e6)));
    f = [];
    try
        f = figure('Visible','off');  plot(1:10);
    catch
        tc.assumeFail('no figure support in this environment');
    end
    out = lsqvce.util.showExport(f, 'postSolve_x1_posC_pass1', 0, d);
    tc.verifyTrue(~isempty(out) && exist(out,'file')==2, 'PNG must be written');
    tc.verifyTrue(contains(out, 'postSolve_x1_posC_pass1.png'));
    % cleanup
    if exist(out,'file'), delete(out); end
    if exist(d,'dir'), rmdir(d, 's'); end
end

function testShowExportSanitisesPipe(tc)
    % A raw stream key with '|' must be sanitised in the file name.
    d = fullfile(tempdir, sprintf('lsqvce_figtest_%d', randi(1e6)));
    f = [];
    try
        f = figure('Visible','off');  plot(1:10);
    catch
        tc.assumeFail('no figure support in this environment');
    end
    out = lsqvce.util.showExport(f, 'postModel_x1|posC_pass3', 0, d);
    tc.verifyFalse(contains(out, '|'), 'the pipe must be sanitised out');
    if exist(out,'file'), delete(out); end
    if exist(d,'dir'), rmdir(d, 's'); end
end

%% ================================================= verbosity levels
% WHY THESE EXIST: until 2026-07-18 not one test asserted a single line of the
% toolbox's own output -- every test ran with Verbose off. That is how a literal
% '\n' bug survived until the user hit it by hand (HANDOVER 1, item 5). These
% tests assert the LEVEL LOGIC and the RENDERED TEXT, so messaging regressions
% fail here instead of in front of a user.

function testVerbosityLogicalMapsToOldBehaviour(tc)
    % The boolean must keep its exact former meaning, or every existing call
    % site silently changes what it prints.
    tc.verifyEqual(lsqvce.util.verbosity.resolve(false), 0);
    tc.verifyEqual(lsqvce.util.verbosity.resolve(true),  1);
end

function testVerbosityNumericLevelsPassThrough(tc)
    for lvl = 0:2
        tc.verifyEqual(lsqvce.util.verbosity.resolve(lvl), lvl);
    end
end

function testVerbosityRejectsRetiredLevelThree(tc)
    % Level 3 was merged into 2 on 2026-07-19. Erroring is deliberate: someone
    % asking for 3 wants MORE output, and silently clamping to 2 would leave
    % them wondering why nothing changed.
    tc.verifyError(@() lsqvce.util.verbosity.resolve(3), 'lsqvce:verbosity:badLevel');
end

function testVerbosityEmptyDefaultsToProgress(tc)
    tc.verifyEqual(lsqvce.util.verbosity.resolve([]), 1);
end

function testVerbosityRejectsBadLevel(tc)
    tc.verifyError(@() lsqvce.util.verbosity.resolve(4), 'lsqvce:verbosity:badLevel');
    tc.verifyError(@() lsqvce.util.verbosity.resolve(-1), 'lsqvce:verbosity:badLevel');
    tc.verifyError(@() lsqvce.util.verbosity.resolve(1.5), 'lsqvce:verbosity:badLevel');
    tc.verifyError(@() lsqvce.util.verbosity.resolve([1 2]), 'lsqvce:verbosity:badLevel');
    tc.verifyError(@() lsqvce.util.verbosity.resolve(9),     'lsqvce:verbosity:badLevel');
end

function testVerbosityPredicatesAreCumulative(tc)
    % Each level INCLUDES those below it. A regression to equality testing
    % (v == 2) would silence detail at a higher level, so pin the whole grid.
    want = [false false      % v = 0
            true  false      % v = 1
            true  true ];    % v = 2
    for lvl = 0:2
        got = [lsqvce.util.verbosity.progress(lvl), ...
               lsqvce.util.verbosity.detail(lvl)];
        tc.verifyEqual(got, want(lvl+1, :), ...
            sprintf('predicate grid wrong at level %d', lvl));
    end
end

function testVerbosityBannersAllShareOneWidth(tc)
    % The point of centralising: the 2026-07-18 log showed '[resid ]' 499 times
    % and '[resid   ]' 63 times -- the same tag at two widths in ONE run.
    tags = {'Variance components','Adjustment','Stochastic model','Convergence', ...
            'Contribution','Joint estimation','residuals','redundancy', ...
            'normal equations','conditioning','data snooping','group pooling'};
    w = cellfun(@(t) numel(lsqvce.util.verbosity.banner(t)), tags);
    tc.verifyEqual(numel(unique(w)), 1, ...
        'every stage label must render to one width so log columns line up');
    tc.verifyEqual(unique(w), lsqvce.util.verbosity.labelWidth(), ...
        'labelWidth must report the width callers should indent continuations by');
end

function testVerbosityBannerIsPlainWordsNotBrackets(tc)
    % Labels are human-readable words, NOT bracketed abbreviations: the old
    % '[SCALE   ]' invented a term ('scale' is not what anyone calls
    % variance-component estimation) and then repeated the message text.
    b = lsqvce.util.verbosity.banner('Variance components');
    tc.verifyFalse(contains(b, '['), 'labels must not be bracketed');
    tc.verifyTrue(startsWith(b, 'Variance components'));
end

function testVerbosityBannerPreservesCase(tc)
    % BANNER must not normalise case: the capitalisation convention is enforced
    % at the call sites (all labels are Capitalised since 2026-07-18), and a
    % function that silently re-cased its input would hide a call site that
    % broke the convention.
    tc.verifyTrue(startsWith(lsqvce.util.verbosity.banner('Adjustment'), 'Adjustment'));
    tc.verifyTrue(startsWith(lsqvce.util.verbosity.banner('residuals'), 'residuals'));
end

function testAllStageLabelsAreCapitalised(tc)
    % The convention itself, checked where it actually lives -- the call sites.
    % A lowercase label would break the column's visual consistency and is the
    % kind of thing that creeps back in one call site at a time.
    root = fileparts(fileparts(mfilename('fullpath')));
    files = dir(fullfile(root, '+lsqvce', '**', '*.m'));
    bad = {};
    for i = 1:numel(files)
        txt = fileread(fullfile(files(i).folder, files(i).name));
        if endsWith(files(i).name, 'verbosity.m'), continue, end
        tok = regexp(txt, 'verbosity\.banner\(''([a-z][^'']*)''\)', 'tokens');
        for j = 1:numel(tok)
            bad{end+1} = sprintf('%s: %s', files(i).name, tok{j}{1}); %#ok<AGROW>
        end
    end
    tc.verifyEmpty(bad, sprintf( ...
        'these stage labels are not capitalised:\n  %s', strjoin(bad, '\n  ')));
end

function testVerbosityBannerKeepsOverlongTagWhole(tc)
    % A label wider than the column loses alignment rather than losing letters --
    % the visible break is the cue to widen WIDTH or shorten the label.
    s = lsqvce.util.verbosity.banner('an extremely long stage label');
    tc.verifyTrue(contains(s, 'an extremely long stage label'));
end

function testVerbosityRejectsBadTag(tc)
    tc.verifyError(@() lsqvce.util.verbosity.banner(''),  'lsqvce:verbosity:badTag');
    tc.verifyError(@() lsqvce.util.verbosity.banner(42),  'lsqvce:verbosity:badTag');
end

function testVerbosityDescribe(tc)
    tc.verifyEqual(lsqvce.util.verbosity.describe(0), 'silent');
    tc.verifyEqual(lsqvce.util.verbosity.describe(false), 'silent');
    tc.verifyEqual(lsqvce.util.verbosity.describe(true), 'progress');
    tc.verifyEqual(lsqvce.util.verbosity.describe(2), 'detail');
end

function testVerbosityBannerRendersThroughLog(tc)
    % END-TO-END on the real output path: a banner written via log.msg must
    % arrive in the log file with its brackets and padding INTACT and with a
    % real newline -- not the two characters backslash-n. This is precisely the
    % bug class that went undetected because no test read the output.
    f = [tempname, '.log'];  fid = fopen(f, 'w');
    c = onCleanup(@() cleanupLog(fid, f));
    lsqvce.util.log.set('fid', fid, 'screen', false);
    lsqvce.util.log.msg('%sdone: %d group(s)\n', ...
        lsqvce.util.verbosity.banner('Variance components'), 3);
    lsqvce.util.log.reset();
    txt = fileread(f);
    tc.verifyTrue(contains(txt, 'Variance components  done: 3 group(s)'), ...
        'the label must render padded, with its arguments substituted');
    tc.verifyFalse(contains(txt, '\n'), ...
        'a literal backslash-n in the output means the format string was escaped twice');
    tc.verifyTrue(endsWith(txt, newline), 'the line must end with a real newline');
end

function testLogSingleArgExpandsEscapes(tc)
    % REGRESSION GUARD (2026-07-18). Callers write log.msg('TEXT\n') expecting
    % FPRINTF semantics. render() used to return a lone argument untouched and
    % emit() writes with fprintf(fid,'%s',txt), which expands nothing -- so the
    % two characters backslash-n reached the log verbatim. ~20 call sites are
    % written that way (report, snoop, robust, BatchIndex, Results, ftest,
    % systematic, adjustVCE) and NOT ONE was covered, so it stayed latent.
    f = [tempname, '.log'];  fid = fopen(f, 'w');
    c = onCleanup(@() cleanupLog(fid, f));
    lsqvce.util.log.set('fid', fid, 'screen', false);
    lsqvce.util.log.msg('  PER-SENSOR SUMMARY\n');
    lsqvce.util.log.reset();
    txt = fileread(f);
    tc.verifyFalse(contains(txt, '\n'), ...
        'a lone format string must have its escapes expanded, like fprintf');
    tc.verifyTrue(endsWith(txt, newline));
    tc.verifyTrue(contains(txt, 'PER-SENSOR SUMMARY'));
end

function testLogSingleArgKeepsEscapedPercent(tc)
    % The flip side of expanding escapes: a caller written against fprintf spells
    % a literal percent '%%', and it must still arrive as one '%'.
    f = [tempname, '.log'];  fid = fopen(f, 'w');
    c = onCleanup(@() cleanupLog(fid, f));
    lsqvce.util.log.set('fid', fid, 'screen', false);
    lsqvce.util.log.msg('share 61.2%% of the information\n');
    lsqvce.util.log.reset();
    txt = fileread(f);
    tc.verifyTrue(contains(txt, 'share 61.2% of the information'));
end

%% ================================================= stage-detail layout
function testDurationFormatsAsClock(tc)
    tc.verifyEqual(lsqvce.util.verbosity.duration(0.519), '00:00:00.519');
    tc.verifyEqual(lsqvce.util.verbosity.duration(1.204), '00:00:01.204');
    tc.verifyEqual(lsqvce.util.verbosity.duration(0),     '00:00:00.000');
    tc.verifyEqual(lsqvce.util.verbosity.duration(5025.006), '01:23:45.006');
end

function testDurationRoundingCarries(tc)
    % 59.9996 s must not render as ':60.000'. The millisecond rounding has to
    % carry into seconds/minutes/hours or the clock shows an impossible time.
    tc.verifyEqual(lsqvce.util.verbosity.duration(59.9996), '00:01:00.000');
    tc.verifyEqual(lsqvce.util.verbosity.duration(3599.9999), '01:00:00.000');
end

function testDurationHandlesNonFinite(tc)
    tc.verifyEqual(lsqvce.util.verbosity.duration(NaN), '--:--:--.---');
    tc.verifyEqual(lsqvce.util.verbosity.duration(Inf), '--:--:--.---');
    % A negative elapsed time (clock skew) clamps rather than printing '-1'.
    tc.verifyEqual(lsqvce.util.verbosity.duration(-3), '00:00:00.000');
end

function testStageDetailRendersTheAgreedLayout(tc)
    % The layout is a user-facing SPECIFICATION, not an implementation detail.
    tc.verifyEqual(lsqvce.util.verbosity.heading('Compute the residuals'), ...
        sprintf('\nCompute the residuals\n'), ...
        'a heading opens with a blank line and carries no colon');
    tc.verifyEqual(lsqvce.util.verbosity.dayLine(54466), sprintf(' Day 54466\n'));

    t = lsqvce.util.verbosity.typeLine(1, 2, 'posD', 'demo-CD-pair', '', 0.368);
    tc.verifyTrue(contains(t, 'Type {1,2}: "posD" of "demo-CD-pair"'));
    tc.verifyTrue(contains(t, '00:00:00.368'));
end

function testTypeLineCarriesTheSystemIndex(tc)
    % The first index is the SYSTEM. It was hardcoded to 1, so every system's
    % types printed as {1,k} and system 2's lines claimed to be system 1's.
    a = lsqvce.util.verbosity.typeLine(1, 1, 'posC', 'pair', '', 0.5);
    b = lsqvce.util.verbosity.typeLine(2, 1, 'posE', 'E',    '', 0.5);
    tc.verifyTrue(contains(a, '{1,1}'));
    tc.verifyTrue(contains(b, '{2,1}'), 'the system index must follow the system');
end

function testTypeLineAlignsItsColumns(tc)
    % Names differ in length; the count and time columns must still line up.
    a = lsqvce.util.verbosity.typeLine(1, 1, 'posC', 'demo-CD-pair', '1400 obs', 0.5);
    b = lsqvce.util.verbosity.typeLine(2, 1, 'rr',   'E',            '900 obs',  1.5);
    tc.verifyEqual(strfind(a, '00:00:'), strfind(b, '00:00:'), ...
        'the time column must sit at one offset regardless of name length');
end

function testStreamLineIsSharedByEveryPerStreamBlock(tc)
    % Stochastic model, Contribution and Convergence all use this one shape, so
    % the blocks of a pass read as one document instead of three formats.
    a = lsqvce.util.verbosity.streamLine('posC', 'demo-CD-pair', 'AR(12)   white');
    b = lsqvce.util.verbosity.streamLine('rr',   'demo-CD-pair', 'sigma = 1.2044');
    tc.verifyTrue(contains(a, '"posC" of "demo-CD-pair"'));
    tc.verifyEqual(numel(a) - numel('AR(12)   white'), ...
                   numel(b) - numel('sigma = 1.2044'), ...
                   'the trailing text must start at the same column in both');
end

function testFieldBlockIsSetApartByBlankLines(tc)
    % A system summary is followed by the NEXT system's work, so it needs a gap
    % on both sides or the two run together.
    s = lsqvce.util.verbosity.fieldBlock('System 1 "demo-CD-pair"', ...
            {'observations','50400'; 'local unknowns','8'});
    tc.verifyTrue(startsWith(s, newline));
    tc.verifyTrue(endsWith(s, newline));
    tc.verifyTrue(contains(s, 'observations     50400'));
end

function testEveryPerStreamBlockSharesOneColumn(tc)
    % Stochastic model, Contribution and Convergence must line up. Convergence
    % used to indent its rows by the 21-character stage-label column while the
    % other two started at 2, so the three tables of one pass stepped in and
    % out. Routing all of them through labelledRow is what fixes that, and this
    % pins it: same start column, same value column.
    a = lsqvce.util.verbosity.streamLine('posC', 'demo-CD-pair', 'AR(12)');
    b = lsqvce.util.verbosity.labelledRow('posC', '41.2%');
    c = lsqvce.util.verbosity.streamLine('rr', 'pair', 'sigma = 1.2044');

    starts = cellfun(@(s) find(s ~= ' ', 1), {a, b, c});
    tc.verifyEqual(starts, repmat(starts(1), 1, 3), ...
        'every per-stream row must begin at the same column');

    tc.verifyEqual(strfind(a, 'AR(12)'), strfind(c, 'sigma'), ...
        'the value column must be at one offset in every block');
end

function testLabelledRowKeepsOverlongLabelWhole(tc)
    % A label wider than the column pushes the value right rather than being
    % truncated -- losing the name is worse than losing the alignment, and the
    % visible break is the cue to widen STREAMCOL.
    s = lsqvce.util.verbosity.labelledRow( ...
            '"aVeryLongObservationType" of "aVeryLongSystemName"', 'x');
    tc.verifyTrue(contains(s, 'aVeryLongObservationType'));
    tc.verifyTrue(contains(s, 'aVeryLongSystemName'));
    tc.verifyTrue(endsWith(strtrim(s), 'x'));
end
