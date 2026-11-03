function demo_messaging(varargin)
%DEMO_MESSAGING  See what each verbosity level actually prints, side by side.
%
%   DEMO_MESSAGING()
%   DEMO_MESSAGING('levels', [0 1])
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
%     bury its own failures), and DEMO_GRACE runs at one level only. So the
%     one thing a user most wants to review before adopting the toolbox -- what
%     it will say while a long job runs -- was the one thing not demonstrable.
%
%   THE LEVELS (see LSQVCE.UTIL.VERBOSITY). Each INCLUDES the ones below it.
%     0  SILENT      nothing. Warnings still appear -- see the note below.
%     1  PROGRESS    MEASURED on this demo's own run (2026-09-07): the stage
%                    label of each swept system ('Joint estimation  swept
%                    system 1'), and the per-pass convergence line
%                    max|alpha2-1| with the arc count and gammaBar. No timing
%                    accompanies these labels on this path. The PER-STREAM view
%                    is recorded, not printed -- read it after the run with
%                    R.summary() (one row per stream) or
%                    R.convergence('Plot', true).
%     2  DETAIL      NOTHING BEYOND LEVEL 1 THROUGH LSQVCE.FIT, measured: the
%                    two runs below print identical text. LSQVCE.OPTIONS/
%                    engineArgs hands every driver 'Verbose',
%                    logical(Verbose), so 1 and 2 both arrive as PROGRESS; and
%                    the accumulator that does keep the numeric level
%                    (LSQVCE.ENGINE.BLOCKACCUMULATOR, which is where a per-day
%                    accumulation block or a cond(Nll) line would come from) is
%                    built silent on the joint path. To see that detail today
%                    you drive the ENGINE tier yourself and hand it a numeric
%                    level. This is a real limitation of the level, honestly
%                    printed by this demo rather than described as working.
%
%   WARNINGS ARE NOT VERBOSITY-CONTROLLED, by design. A run that failed to
%   converge says so even at level 0 -- otherwise a silent batch job would hide
%   its own failure. You will see warnings under level 0 below; that is correct.
%
%   READING THE OUTPUT
%     Stage labels sit in a fixed column so the log reads as a table. EVERY
%     LABEL IS CAPITALISED (user decision, 2026-07-18) -- 'Joint estimation' is
%     the one this demo's path emits; 'Residuals' and 'Variance components' come
%     from the sweep and the component step of other runs. Case carries no
%     meaning; the stage/detail
%     hierarchy is carried by INDENTATION alone. Continuation lines, per-day
%     blocks and per-stream tables are indented beneath the label they belong
%     to, so they read as belonging to the stage above them rather than starting
%     a new one. The only lowercase text in that column is quoted DATA -- an
%     observation type, a system name -- never a label.
%
%   INPUT (Name/Value)
%     'levels'  1xK double   which levels to show          (default [0 1 2])
%     'toFile'  1xN char     also write everything here    (default '', off)
%
%   OUTPUT
%     none. Everything is printed.
%
%   EXAMPLE
%     demo_messaging('levels', [0 1]);      % silence against the default log
%
%   See also LSQVCE.UTIL.VERBOSITY, LSQVCE.UTIL.LOG, LSQVCE.FIT,
%            LSQVCE.OPTIONS, DEMO_GRACE, DEMO_GRACE_DATALOADER
%
%   Updated 2026-11-02 (CET, UTC+01:00), WP-161: the run fits the loader's
%   lsqvce.TaskDefinition. With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-09-07 (CEST, UTC+02:00), R6 WP-95: level 2 is documented as
%   not reachable through LSQVCE.FIT, the level-1 line names what is actually
%   printed, the retired case convention is replaced by the capitalised-label
%   rule, and the dead demo_separable cross-references are re-pointed.
%   With assistance from Claude (Anthropic), model Claude Opus 5.
%   Updated 2026-08-31 (CEST, UTC+02:00), v4 SCOPED tables (Level L11,
%   WP-77): the declaration tables carry the scope column (commonParams
%   'common'|'day', localParams 'arc'|'day'|'system'); localsPerArc left the
%   USER surface. With assistance from Claude (Anthropic), model
%   Claude Fable 5.
%   Updated 2026-08-31 (CEST, UTC+02:00), v3 data contract (Level L10,
%   WP-74): the run is lsqvce.fit over the demo's v3 data struct;
%   Options.Verbose carries the level. With assistance from Claude
%   (Anthropic), model Claude Fable 5.
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
             'see at level 1 is still there at level 2. Pick the LOWEST level that\n' ...
             'answers the question you actually have. Level 2 currently prints the\n' ...
             'same as level 1 through LSQVCE.FIT -- see the help for why.\n']);
end

%==========================================================================
function txt = runAtLevel(lvl)
%RUNATLEVEL  Run the demo problem at one verbosity level, capturing all output.
%   The toolbox routes every status line AND every warning through one sink
%   (LSQVCE.UTIL.LOG), so pointing that sink at a temp file captures exactly what
%   the user would have seen on screen -- which is the point of the demo.
%
%   INPUT   lvl  1x1 double   0..2
%   OUTPUT  txt  1xN char     everything the run logged
    f   = [tempname, '.log'];
    fid = fopen(f, 'w');
    c   = onCleanup(@() cleanupLevelLog(fid, f));   %#ok<NASGU>

    lsqvce.util.log.set('fid', fid, 'screen', false);
    try
        [~, task] = demo_grace_dataLoader();      % its lsqvce.TaskDefinition
        % MaxPasses 1 keeps the demo quick: the point is WHAT is
        % printed at each level, not how long the loop takes to converge. A
        % short run does mean it stops before convergence, which is why the
        % level-0 block below still shows warnings. Options.Verbose carries the
        % level 0..2 (v3, L10 WP-74), but engineArgs hands the drivers a
        % LOGICAL, so 1 and 2 reach them alike -- see the help block.
        lsqvce.fit(task, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, Verbose=lvl, ...
                MaxPasses=1));
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
