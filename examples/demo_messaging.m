function demo_messaging(varargin)
%DEMO_MESSAGING  See what each verbosity level actually prints, side by side.
%
%   DEMO_MESSAGING()
%   DEMO_MESSAGING('levels', [0 2])
%   DEMO_MESSAGING('toFile', '/tmp/levels.txt')
%
%   Runs the SAME small GRACE-FO problem once per verbosity level and shows the
%   output of each, so the levels can be compared on real output rather than on
%   a description of it.
%
%   WHY THIS EXISTS
%     The verbosity levels and the stage labels were agreed from samples pasted
%     into a chat. Nothing in the toolbox let you SEE them: every test runs with
%     verbose off (deliberately -- a test suite that printed a full run would
%     bury its own failures), and DEMO_SEPARABLE runs at one level only. So the
%     one thing a user most wants to review before adopting the toolbox -- what
%     it will say while a long job runs -- was the one thing not demonstrable.
%
%   THE LEVELS (see LSQVCE.UTIL.VERBOSITY). Each INCLUDES the ones below it.
%     0  SILENT      nothing. Warnings still appear -- see the note below.
%     1  PROGRESS    outer passes, each stage's start/done with timing, and the
%                    convergence report (PER STREAM, not just max|sigma-1|).
%     2  DETAIL      everything else: per-day accumulation and residuals, per-
%                    stream AR order and whiteness, per-group shares and sigmas,
%                    system conditioning, inner VCE iterations, fallbacks.
%
%   WARNINGS ARE NOT VERBOSITY-CONTROLLED, by design. A run that failed to
%   converge says so even at level 0 -- otherwise a silent batch job would hide
%   its own failure. You will see warnings under level 0 below; that is correct.
%
%   READING THE OUTPUT
%     Stage labels sit in a fixed column so the log reads as a table. CASE
%     CARRIES MEANING: 'Variance components' (capitalised) is an outer-loop
%     STAGE, 'residuals' (lowercase) is detail beneath one. Continuation lines
%     and tables are indented to the label column, so they read as belonging to
%     the stage above them rather than starting a new one.
%
%   INPUT (Name/Value)
%     'levels'  1xK double   which levels to show          (default [0 1 2])
%     'toFile'  1xN char     also write everything here    (default '', off)
%
%   OUTPUT
%     none. Everything is printed.
%
%   EXAMPLE
%     demo_messaging('levels', [1 2]);      % compare the default with full detail
%
%   See also LSQVCE.UTIL.VERBOSITY, DEMO_SEPARABLE, DEMO_SEPARABLE_DATA
%
%   Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.

    p = inputParser;  p.FunctionName = 'demo_messaging';
    p.addParameter('levels', [0 1 2], @(v) isnumeric(v) && all(v >= 0 & v <= 2));
    p.addParameter('toFile', '', @(s) ischar(s) || isstring(s));
    p.parse(varargin{:});
    o = p.Results;

    outFid = -1;
    if ~isempty(char(o.toFile))
        outFid = fopen(char(o.toFile), 'w');
        if outFid < 0
            error('demo_messaging:badFile', ...
                  'Could not open %s for writing.', char(o.toFile));
        end
        cleanOut = onCleanup(@() fclose(outFid));   %#ok<NASGU>
    end

    for lvl = o.levels(:).'
        txt = runAtLevel(lvl);
        header = sprintf(['\n' repmat('=', 1, 78) '\n' ...
                          '  verbose = %d   (%s)\n' ...
                          repmat('=', 1, 78) '\n'], ...
                          lvl, lsqvce.util.verbosity.describe(lvl));
        fprintf('%s', header);
        if isempty(strtrim(txt))
            fprintf('  (no output)\n');
        else
            fprintf('%s', txt);
        end
        if outFid > 0
            fprintf(outFid, '%s', header);
            fprintf(outFid, '%s', txt);
        end
    end

    fprintf(['\n' repmat('-', 1, 78) '\n' ...
             'Each level is a strict superset of the one below it, so a line you\n' ...
             'see at level 1 is still there at level 3. Pick the LOWEST level that\n' ...
             'answers the question you actually have.\n']);
end

%==========================================================================
function txt = runAtLevel(lvl)
%RUNATLEVEL  Run the demo problem at one verbosity level, capturing all output.
%   The toolbox routes every status line AND every warning through one sink
%   (LSQVCE.UTIL.LOG), so pointing that sink at a temp file captures exactly what
%   the user would have seen on screen -- which is the point of the demo.
%
%   INPUT   lvl  1x1 double   0..3
%   OUTPUT  txt  1xN char     everything the run logged
    f   = [tempname, '.log'];
    fid = fopen(f, 'w');
    c   = onCleanup(@() cleanupLevelLog(fid, f));   %#ok<NASGU>

    lsqvce.util.log.set('fid', fid, 'screen', false);
    try
        [cfg, L] = demo_grace_dataLoader();
        % MaxOuter 1 / MaxIter 2 keep the demo quick: the point is WHAT is
        % printed at each level, not how long the loop takes to converge. A
        % short run does mean it stops before convergence, which is why the
        % level-0 block below still shows warnings.
        P = lsqvce.adjustVCE(L, cfg.days, 'verbose', lvl, 'maxOuter', 1);
        P = P.groupComponentsBy('obsType').varianceComponents('MaxIter', 2);
        P.run();
    catch ME
        lsqvce.util.log.reset();
        rethrow(ME);
    end
    lsqvce.util.log.reset();
    txt = fileread(f);
end

%==========================================================================
function cleanupLevelLog(fid, f)
%CLEANUPLEVELLOG  Restore the global log state and remove the temp file.
    lsqvce.util.log.reset();
    if fid > 0
        try, fclose(fid); catch, end   %#ok<NOCOM>
    end
    if exist(f, 'file'), delete(f); end
end
