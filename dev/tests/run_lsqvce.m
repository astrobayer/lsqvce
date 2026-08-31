function results = run_lsqvce(what, varargin)
%RUN_LSQVCE  Run a SUBSET of the lsqvce test suite, for fast iteration.
%
%   run_lsqvce                     % list the available shortcuts and exit
%   results = run_lsqvce(what)
%   results = run_lsqvce(what, 'Name', pattern)
%
%   RUN_ALL_LSQVCE runs everything (498 tests, ~130 s) and is the gate before
%   shipping; `run_lsqvce('all')` is the same run under this function's name.
%   While iterating on one area that is cumbersome, so this runner takes a
%   shortcut name (or a file path, or a name pattern) and runs only that,
%   printing a one-line summary instead of a wall of output.
%
%   INPUT
%     what  char   one of the shortcuts below, OR a path to a test .m file, OR
%                  a '*pattern*' matched against test NAMES across the suite:
%
%       'joint'      the joint/vceMode end-to-end tests (tests/test_adjustVCE)
%       'ellmer'     the Ellmer primitives    (tests/test_ellmer.m, ~47 tests)
%       'contrib'    the contribution tests   (tests/test_contribution.m)
%       'adjustvce'  the whole adjustVCE file (tests/test_adjustVCE.m)
%       'messaging'  the verbosity wiring tests (levels 0-3, captured text)
%       'accum'      the accumulation layer   (tests/test_lsqAccum.m)
%       'block'      adjustBlock/Results/makeLoader (tests/test_adjustBlock.m)
%       'parallel'   mapreduce + PSD iteration (tests/test_parallel.m)
%       'vce'        lsqvce.engine.vce contract      (tests/test_vceFunction.m)
%       'vcemath'    VCE mathematics by hand  (tests/test_vceMath.m)
%       'signal'     outliers/varWhiten/ar2arma (tests/test_signalCore.m)
%       'streams'    per-stream helpers       (tests/test_streamUtils.m)
%       'all'        the whole suite -- delegates to RUN_ALL_LSQVCE (see OUTPUT)
%
%   OPTIONS (Name/Value)
%     'Name'    char     an extra name filter, e.g. '*degenerate*'. Combined with
%                        the shortcut, so run_lsqvce('ellmer','Name','*psdVCE*')
%                        runs only the psdVCE tests inside test_ellmer.
%     'Quiet'   logical  suppress warnings during the run (default true). The
%                        suite deliberately provokes warnings in failure-mode
%                        tests, which otherwise bury the result.
%     'LogFile' char     ALSO write the full console output -- warnings, failure
%                        diagnostics, stack traces -- to this file. Default is
%                        'tests/last_run_<target>.txt', so runs of DIFFERENT
%                        targets do not overwrite each other and several can be
%                        inspected after the fact. Pass '' to disable, or a path
%                        to choose your own.
%                        With a LogFile, 'Quiet' is IGNORED for the file: the
%                        warnings are captured there even though the console
%                        stays clean.
%
%   OUTPUT
%     results  1xN matlab.unittest.TestResult   the raw results (also summarised
%                                               to the command window)
%
%              EXCEPT for 'all', which delegates to RUN_ALL_LSQVCE and therefore
%              returns THAT function's 1xS per-suite summary struct array
%              (.name .passed .total .failed .incomplete .secs). The two runners
%              answer different questions -- one test vs one suite -- so they
%              summarise at different granularity. Check the target before
%              indexing the result programmatically.
%
%   EXAMPLES
%     run_lsqvce('joint')                        % the four joint tests
%     run_lsqvce('ellmer', 'Name', '*psdVCE*')   % just the PSD VCE tests
%     run_lsqvce('*degenerate*')                 % that test wherever it lives
%
%   Each run writes tests/last_run_<target>.txt, so successive runs of different
%   targets can all be inspected afterwards.
%
%   See also RUN_ALL_LSQVCE (the full gate), RUNTESTS
%
%   Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.

    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % .../dev/tests -> toolbox root (L8)

    if nargin < 1 || isempty(what)
        printShortcuts();
        results = [];
        return
    end

    p = inputParser;  p.FunctionName = 'run_lsqvce';
    p.addParameter('Name' , '', @(s) ischar(s) || isstring(s));
    p.addParameter('Quiet', true, @(v) isscalar(v));
    p.addParameter('LogFile', '', @(s) ischar(s) || isstring(s));
    p.parse(varargin{:});
    o = p.Results;
    logFile = char(o.LogFile);

    % ---- resolve the shortcut to a target + optional name filter ----------
    nameFilter = char(o.Name);
    switch lower(char(what))
        case 'joint'
            target = fullfile(root,'tests','test_adjustVCE.m');
            if isempty(nameFilter), nameFilter = '*joint*'; end
        case 'ellmer'
            target = fullfile(root,'tests','test_ellmer.m');
        case 'contrib'
            target = fullfile(root,'tests','test_contribution.m');
        case 'adjustvce'
            target = fullfile(root,'tests','test_adjustVCE.m');
        case 'messaging'
            % The verbosity wiring: levels 0-3 asserted on captured log text.
            target = fullfile(root,'tests','test_adjustVCE.m');
            if isempty(nameFilter), nameFilter = '*Messaging*'; end
        case 'accum'
            target = fullfile(root,'tests','test_lsqAccum.m');
        case 'block'
            target = fullfile(root,'tests','test_adjustBlock.m');
        case 'parallel'
            target = fullfile(root,'tests','test_parallel.m');
        case 'vce'
            target = fullfile(root,'tests','test_vceFunction.m');
        case 'vcemath'
            target = fullfile(root,'tests','test_vceMath.m');
        case 'signal'
            target = fullfile(root,'tests','test_signalCore.m');
        case 'streams'
            target = fullfile(root,'tests','test_streamUtils.m');
        case 'all'
            % DELEGATES to RUN_ALL_LSQVCE -- it does not re-implement it.
            %
            % This used to be `runtests(fullfile(root,'tests'))`, i.e. runtests
            % on a FOLDER, which is the trap HISTORY_v1.md 1.7b removed: a file that
            % fails to LOAD is silently skipped and the remaining files still
            % report green, so the suite can look clean while an entire file
            % never ran. RUN_ALL_LSQVCE iterates the registered list, counts a
            % load failure AS a failure, keeps the warning inventory (a signal,
            % not noise -- 4d) and records per-suite seconds.
            %
            % One gate, one implementation. `run_lsqvce('all')` and
            % `run_all_lsqvce` are the same run; this is a spelling, not a
            % second engine.
            if ~isempty(nameFilter)
                error('lsqvce:runLsqvce:allTakesNoFilter', ...
                    ['''all'' runs the whole gate and cannot be name-filtered. ' ...
                     'Use a pattern on its own, e.g. run_lsqvce(''%s'').'], ...
                     nameFilter);
            end
            results = run_all_lsqvce();
            return
        otherwise
            w = char(what);
            if exist(w, 'file') == 2
                target = w;                       % an explicit file path
            elseif contains(w, '*')
                target = fullfile(root,'tests');  % a bare pattern: search all
                if isempty(nameFilter), nameFilter = w; end
            else
                printShortcuts();
                error('lsqvce:run_lsqvce:unknownTarget', ...
                      'Unknown target ''%s''. See the shortcuts above.', w);
            end
    end

    args = {target};
    if ~isempty(nameFilter), args = [args, {'Name', nameFilter}]; end
    if exist(target, 'dir') == 7
        % A folder target: recurse. tests/ is FLAT since 2026-07-18, but
        % IncludingSubfolders costs nothing and keeps this correct if a
        % subfolder is ever reintroduced.
        args = [args, {'IncludeSubfolders', true}];
    end

    % Default log name is PER TARGET, so run_lsqvce('ellmer') and
    % run_lsqvce('joint') do not overwrite each other and both can be read
    % afterwards. Non-alphanumerics in a pattern target are squashed so the
    % name stays a valid filename.
    if isempty(logFile)
        tag = regexprep(lower(char(what)), '[^a-z0-9]+', '_');
        tag = regexprep(tag, '^_+|_+$', '');
        if isempty(tag), tag = 'run'; end
        logFile = fullfile(root, 'tests', sprintf('last_run_%s.txt', tag));
    end

    % ---- run, capturing everything when a log file is wanted --------------
    %   With a LogFile we deliberately do NOT mute warnings: the whole point of
    %   the file is to preserve the failure diagnostics and the warnings that
    %   explain them. The console still gets only the summary.
    if ~isempty(logFile)
        try
            captured = evalc('results = runtests(args{:});');
        catch ME
            captured = sprintf('RUN ERRORED: %s\n%s\n', ME.identifier, ME.message);
            results = matlab.unittest.TestResult.empty;
        end
        writeLog(logFile, what, nameFilter, captured, results);
    else
        ws = [];
        if logical(o.Quiet)
            ws = warning('off','all');
        end
        try
            results = runtests(args{:});
        catch ME
            if ~isempty(ws), warning(ws); end   % never leave warnings muted
            rethrow(ME);
        end
        if ~isempty(ws), warning(ws); end
    end

    % ---- one-line summary -------------------------------------------------
    if isempty(results)
        fprintf('\n==== %s : run ERRORED, see the log ====\n', char(what));
        if ~isempty(logFile), fprintf('full log: %s\n', logFile); end
        return
    end
    nP = sum([results.Passed]);
    nF = sum([results.Failed]);
    nI = sum([results.Incomplete]);
    fprintf('\n==== %s : %d/%d passed', char(what), nP, numel(results));
    if nF > 0, fprintf(', %d FAILED', nF); end
    if nI > 0, fprintf(', %d incomplete', nI); end
    fprintf(' ====\n');
    if ~isempty(logFile)
        fprintf('full log: %s\n', logFile);
    end

    if nF > 0
        fprintf('\nFailed:\n');
        bad = results([results.Failed]);
        for k = 1:numel(bad)
            fprintf('  %s\n', bad(k).Name);
        end
        fprintf(['\nRe-run one of them with its diagnostics:\n' ...
                 '  run_lsqvce(''*%s*'')\n'], lastNamePart(bad(1).Name));
    end
end

%==========================================================================
function writeLog(logFile, what, nameFilter, captured, results)
%WRITELOG  Persist the full run -- summary, console output, and the per-failure
%   diagnostics -- so it can be read or shared verbatim afterwards.
    fid = fopen(logFile, 'w');
    if fid < 0
        warning('lsqvce:run_lsqvce:noLog', ...
                'Could not open %s for writing; skipping the log.', logFile);
        return
    end
    c = onCleanup(@() fclose(fid));

    fprintf(fid, '=== run_lsqvce log ===\n');
    fprintf(fid, 'when   : %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS')); %#ok<TNOW1,DATST>
    fprintf(fid, 'target : %s\n', char(what));
    fprintf(fid, 'filter : %s\n', nameFilter);
    if ~isempty(results)
        fprintf(fid, 'result : %d passed, %d failed, %d incomplete of %d\n\n', ...
            sum([results.Passed]), sum([results.Failed]), ...
            sum([results.Incomplete]), numel(results));
    else
        fprintf(fid, 'result : (no results -- the run errored)\n\n');
    end

    % The per-failure diagnostics carry the assertion text and any warning the
    % code emitted to explain itself; they are the most useful part of the log.
    if ~isempty(results) && any([results.Failed])
        bad = results([results.Failed]);
        fprintf(fid, '=== FAILURES (%d) ===\n', numel(bad));
        for k = 1:numel(bad)
            fprintf(fid, '\n--- %s ---\n', bad(k).Name);
            try
                d = bad(k).Details.DiagnosticRecord;
                for j = 1:numel(d)
                    if ~isempty(d(j).Report)
                        fprintf(fid, '%s\n', d(j).Report);
                    end
                end
            catch
                fprintf(fid, '(diagnostics unavailable)\n');
            end
        end
        fprintf(fid, '\n');
    end

    fprintf(fid, '=== CONSOLE OUTPUT (warnings included) ===\n%s\n', captured);
end

%==========================================================================
function printShortcuts()
    fprintf(['run_lsqvce -- run part of the suite\n\n' ...
             '  run_lsqvce(''joint'')      joint / vceMode end-to-end tests\n' ...
             '  run_lsqvce(''ellmer'')     Ellmer primitives (~47)\n' ...
             '  run_lsqvce(''contrib'')    contribution tests\n' ...
             '  run_lsqvce(''adjustvce'')  the whole adjustVCE file\n' ...
             '  run_lsqvce(''messaging'')  verbosity wiring, levels 0-3\n' ...
             '  run_lsqvce(''accum'')      accumulation layer\n' ...
             '  run_lsqvce(''block'')      adjustBlock / Results / makeLoader\n' ...
             '  run_lsqvce(''parallel'')   mapreduce + PSD iteration\n' ...
             '  run_lsqvce(''vce'')        lsqvce.engine.vce contract\n' ...
             '  run_lsqvce(''vcemath'')    VCE mathematics by hand\n' ...
             '  run_lsqvce(''signal'')     outliers / varWhiten / ar2arma\n' ...
             '  run_lsqvce(''streams'')    per-stream helpers\n' ...
             '  run_lsqvce(''all'')        the whole suite (= run_all_lsqvce)\n\n' ...
             '  run_lsqvce(''ellmer'', ''Name'', ''*psdVCE*'')   narrow further\n' ...
             '  run_lsqvce(''*degenerate*'')                 by name, anywhere\n\n' ...
             'Use run_all_lsqvce for the full check before shipping.\n']);
end

%==========================================================================
function s = lastNamePart(fullName)
%LASTNAMEPART  'test_adjustVCE/test_joint_x' -> 'test_joint_x', for the hint.
    parts = strsplit(fullName, '/');
    s = parts{end};
end
