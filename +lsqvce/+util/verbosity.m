classdef verbosity
%LSQVCE.UTIL.VERBOSITY  Message levels and banner formatting, in one place.
%
%   The toolbox used a BOOLEAN 'verbose', which forced one decision on two
%   different audiences: someone watching a long run wants a progress line per
%   outer pass, while someone debugging a stream that will not whiten wants the
%   per-stream detail behind it. One flag cannot serve both, so the output was
%   either too sparse to diagnose or too noisy to watch. This class replaces the
%   flag with LEVELS (three today; four when it was written, see below), and
%   owns the banner vocabulary so every line in the
%   toolbox is laid out the same way.
%
%   THE LEVELS
%     0  SILENT    nothing at all.
%     1  PROGRESS  the shape of the run: outer passes, each stage's start/done
%                  with timings, and the convergence report. The default, and
%                  what you want when watching a long job.
%     2  DETAIL    everything else: per-day accumulation and residuals, the
%                  estimated spectral shape and each stream's level and
%                  whiteness, per-group shares and sigmas, the conditioning of
%                  each system, inner VCE iterations, fallbacks.
%
%   THERE USED TO BE A LEVEL 3 (STAGE / DIAGNOSTIC were separate). It was
%   dropped on 2026-07-19: once the inner accumulator stopped repeating its
%   per-day block on every Helmert iteration, the only thing level 3 still added
%   over level 2 was the per-group sigmas -- too thin a difference to justify a
%   tier, and a user then has to guess which of two "verbose" levels they want.
%   One detail level, or none.
%
%   Each level INCLUDES the ones below it, so a check is always "is the current
%   level at least L", never an equality test. Use the predicates:
%
%     if lsqvce.util.verbosity.atLeast(v, 2), ... end     % explicit
%     if lsqvce.util.verbosity.detail(v),     ... end     % same thing, named
%
%   BACKWARD COMPATIBILITY (this is why nothing broke)
%     Every public 'verbose' option still accepts a LOGICAL, and RESOLVE maps it:
%     false -> 0 (SILENT) and true -> 1 (PROGRESS). Those are exactly the two
%     behaviours the boolean used to select, so existing calls, tests and demos
%     keep their meaning. Levels 2 and 3 are new capability, reachable only by
%     passing a number.
%
%   THE STAGE LABELS (why they are centralised, and why they are words)
%     Labels were written by hand at ~30 call sites and had drifted. The MATLAB
%     log of 2026-07-18 showed '[resid ]' 499 times and '[resid   ]' 63 times --
%     the SAME tag at two widths in ONE run. Misaligned labels make a log hard to
%     scan, which matters because reading the log is how defects are caught here
%     (HISTORY_v1.md 4d). BANNER puts every label in one column.
%
%     They are PLAIN WORDS, not abbreviations. The abbreviated form was dropped
%     because the codes were invented rather than geodetic -- 'scale' and 'shape'
%     are not what anyone calls variance-component estimation and stochastic
%     modelling -- and because a tag that abbreviates the message merely repeats
%     it. Compare:
%
%       [SCALE   ] variance-component estimation between groups   (says it twice)
%       Variance components   estimating between groups           (says it once)
%
%     EVERY LABEL IS CAPITALISED (user decision, 2026-07-18). An earlier version
%     used case to mark stage-vs-detail ('Variance components' against
%     'residuals'), but the hierarchy is already carried by INDENTATION -- detail
%     sits under the stage it belongs to -- so case was a second, weaker signal
%     for the same thing, and a log with mixed capitalisation just reads as
%     inconsistent. BANNER still preserves whatever case it is given; the
%     convention is enforced by the call sites, not by this function.
%
%   EXAMPLE
%     v = lsqvce.util.verbosity.resolve(2);          % or true / false / 0..2
%     if lsqvce.util.verbosity.progress(v)
%         lsqvce.util.log.msg('%s estimating between groups ...\n', ...
%             lsqvce.util.verbosity.banner('Variance components'));
%     end
%     if lsqvce.util.verbosity.detail(v)
%         lsqvce.util.log.msg('%s %s: AR(%d)\n', ...
%             lsqvce.util.verbosity.banner('Stochastic model'), key, p);
%     end
%
%   ERRORS
%     lsqvce:verbosity:badLevel  not a scalar, or outside 0..2 (level 3 was
%                                merged into 2; a silent clamp would leave a
%                                caller wondering why nothing changed)
%     lsqvce:verbosity:badTag    an unknown message tag
%
%   See also LSQVCE.UTIL.LOG, LSQVCE.ENGINE.ADJUSTVCEJOINT
%
%   Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.

    properties (Constant)
        SILENT   = 0   % nothing
        PROGRESS = 1   % outer passes, stage start/done, the convergence report
        DETAIL   = 2   % + everything the run can say about what it is doing
    end

    properties (Constant, Access = private)
        DESCCOL   = 46   % width of the 'Type {i,j}: "x" of "y"' column
        EXTRACOL  = 14   % width of the optional middle column (counts)
        STREAMCOL = 38   % width of the '"type" of "system"' column
        WIDTH = 21       % label column width; widen freely if a label needs it
                         % (longest today: 'variance components')
    end

    methods (Static)

        function lvl = resolve(v)
        %RESOLVE  Turn a user-supplied verbosity into a level 0..2.
        %   Accepts the historical LOGICAL as well as a numeric level, so a
        %   caller written against the old boolean keeps its exact behaviour.
        %
        %   INPUT
        %     v    1x1 logical | 1x1 numeric | []   false/true, or 0..2.
        %                                           [] means "unset" -> PROGRESS.
        %   OUTPUT
        %     lvl  1x1 double   0 (SILENT) | 1 (PROGRESS) | 2 (DETAIL)
        %
        %   MAPPING
        %     false -> 0, true -> 1. Those are the two behaviours the boolean
        %     used to select; level 2 needs an explicit number.
        %
        %   ERRORS
        %     lsqvce:verbosity:badLevel   not a scalar, or outside 0..2. Level 3
        %                                 is rejected with its own message: it
        %                                 was merged into 2 on 2026-07-19, and a
        %                                 silent clamp would leave a caller who
        %                                 asked for 3 wondering why nothing
        %                                 changed.
            if nargin < 1 || isempty(v)
                lvl = lsqvce.util.verbosity.PROGRESS;
                return
            end
            if ~isscalar(v)
                error('lsqvce:verbosity:badLevel', ...
                    ['verbose must be a scalar: a logical (false/true) or a ' ...
                     'level 0..2 (0 silent, 1 progress, 2 detail).']);
            end
            if islogical(v)
                lvl = double(v);            % false->0 SILENT, true->1 PROGRESS
                return
            end
            if isnumeric(v) && isfinite(v) && v == 3
                % Level 3 existed until 2026-07-19 and was merged into 2. Say so
                % rather than silently clamping: a caller who asks for 3 has a
                % reason, and a clamp would leave them wondering why the output
                % did not change.
                error('lsqvce:verbosity:badLevel', ...
                    ['verbose level 3 no longer exists: it was merged into ' ...
                     'level 2 on 2026-07-19, which now shows every detail. ' ...
                     'Use 2.']);
            end
            if ~isnumeric(v) || ~isfinite(v) || v ~= fix(v) || v < 0 || v > 2
                error('lsqvce:verbosity:badLevel', ...
                    ['verbose must be a logical or an integer level 0..2 ' ...
                     '(0 silent, 1 progress, 2 detail); got %s.'], ...
                    strtrim(evalc('disp(v)')));
            end
            lvl = double(v);
        end

        function tf = atLeast(v, lvl)
        %ATLEAST  True when the current verbosity reaches a level.
        %   The levels are cumulative, so this is the ONLY comparison the
        %   toolbox should make -- never an equality test like `v == 1`, which
        %   would silence a line at every HIGHER level and is how detail gets
        %   lost when a level is added.
        %
        %   INPUT
        %     v    1x1 logical|numeric  current verbosity (resolved internally)
        %     lvl  1x1 double           the level being asked about, 0..2
        %   OUTPUT
        %     tf   1x1 logical
            tf = lsqvce.util.verbosity.resolve(v) >= lvl;
        end

        function tf = progress(v)
        %PROGRESS  True at level >= 1: outer passes, stage start/done, verdict.
        %   INPUT  v 1x1 logical|numeric   OUTPUT tf 1x1 logical
            tf = lsqvce.util.verbosity.atLeast(v, lsqvce.util.verbosity.PROGRESS);
        end

        function tf = detail(v)
        %DETAIL  True at level >= 2: everything beyond the progress skeleton.
        %   The ONLY detail predicate. There were two (stage/diagnostic) while
        %   there were two detail levels; keeping both after the merge would be
        %   two names for one test, which is the ambiguity this toolbox has been
        %   removing everywhere else.
        %   INPUT  v 1x1 logical|numeric   OUTPUT tf 1x1 logical
            tf = lsqvce.util.verbosity.atLeast(v, lsqvce.util.verbosity.DETAIL);
        end

        function s = banner(tag)
        %BANNER  A stage label, padded to one column so the log reads as a table.
        %   Written in PLAIN WORDS, not an abbreviated code. The earlier form was
        %   a bracketed 8-character tag ('[SCALE   ]', '[shape   ]'), which had
        %   two problems: the abbreviations were invented rather than geodetic
        %   ('scale' and 'shape' are not what anyone calls variance-component
        %   estimation and stochastic modelling), and the tag then repeated what
        %   the message already said. Words in a fixed column carry the same
        %   structure without either fault, and the column can simply be widened
        %   when a longer label is genuinely needed.
        %
        %   CASE STILL CARRIES MEANING and is preserved, not normalised:
        %     Capitalised  an outer-loop stage    'Variance components'
        %     lowercase    detail beneath a stage 'residuals', 'redundancy'
        %
        %   INPUT
        %     tag  1xN char | string   e.g. 'Variance components', 'residuals'
        %   OUTPUT
        %     s    1xM char            the label padded to the common width
        %
        %   A label longer than the width is kept WHOLE (the line loses alignment
        %   rather than the reader losing the word) -- that visible break is the
        %   cue to widen WIDTH or shorten the label.
        %
        %   ERRORS
        %     lsqvce:verbosity:badTag   empty or not text
            if nargin < 1 || isempty(tag) || ~(ischar(tag) || isstring(tag))
                error('lsqvce:verbosity:badTag', ...
                    'banner needs a non-empty char or string tag.');
            end
            t = char(tag);
            w = lsqvce.util.verbosity.WIDTH;
            if numel(t) < w
                t(end+1:w) = ' ';
            end
            s = t;
        end

        function w = labelWidth()
        %LABELWIDTH  The width of the stage-label column, for aligning under it.
        %   Callers that print a continuation or a table BELONGING to a stage
        %   indent by this much, so the rows line up beneath the label rather
        %   than beside it. Exposed so the column width lives in ONE place: a
        %   hardcoded indent elsewhere would silently drift the day WIDTH changes.
        %
        %   OUTPUT  w  1x1 double  characters
            w = lsqvce.util.verbosity.WIDTH;
        end

        function s = describe(v)
        %DESCRIBE  The level as a short human-readable word, for reports.
        %   INPUT  v 1x1 logical|numeric   OUTPUT s 1xN char
            switch lsqvce.util.verbosity.resolve(v)
                case 0, s = 'silent';
                case 1, s = 'progress';
                otherwise, s = 'detail';
            end
        end

        % ---- multi-line stage detail (level 2/3 layout) -------------------
        % The banner style above is a FLAT log: one bracketed tag per line, good
        % for scanning a long run. The helpers below render the INDENTED,
        % NUMBERED form the user actually reads when following a single day:
        %
        %   Accumulate the normal matrix:
        %    Day 54466
        %     1. Preparing partial derivatives for observations of type "satellite":
        %      -> Group {1,1}: type "kinematicPosition" for "GraceA"   - 00:00:00.519
        %      -> Group {1,2}: type "kinematicPosition" for "GraceB"   - 00:00:00.368
        %     2. Accumulate the normal matrix:
        %      -> Group {1,1}: type "kinematicPosition" for "GraceA"   - 00:00:01.204
        %
        % Indentation encodes the hierarchy (activity > day > step > group), and
        % every duration is HH:MM:SS.mmm so columns line up and a slow group is
        % visible at a glance. These return TEXT, so a caller still decides the
        % level at which to emit it.

        function s = heading(text)
        %HEADING  Top-level activity line, preceded by a blank line.
        %   'Accumulate the normal matrix', 'Compute the residuals',
        %   'Stochastic model', 'Contribution', 'Convergence'.
        %
        %   THE BLANK LINE IS THE POINT. One pass emits four or five of these
        %   blocks back to back; without a gap they merge into one wall and the
        %   reader has to parse the text to find where one ends. No trailing
        %   colon: the block's own rows say what it is, and a colon on some
        %   headings but not others is exactly the inconsistency that makes a
        %   log look unmaintained.
        %
        %   INPUT   text 1xN char|string   the activity
        %   OUTPUT  s    1xM char          blank line + text + newline
            s = sprintf('\n%s\n', char(text));
        end

        function s = dayLine(day)
        %DAYLINE  The day under an activity heading, e.g. ' Day 54466'.
        %   INPUT   day 1x1 numeric | char | string   day identifier (e.g. MJD)
        %   OUTPUT  s   1xN char
            s = sprintf(' Day %s\n', char(string(day)));
        end


        function s = typeLine(sysIdx, typeIdx, obsType, sysName, extra, seconds)
        %TYPELINE  One observation type's line, under a day header.
        %   '  Type {1,2}: "posD" of "demo-CD-pair"    1400 obs    00:00:00.368'
        %
        %   ONE LINE PER TYPE, in fixed columns so the eye can run down the
        %   times and the counts. Replaces an earlier two-line form (a numbered
        %   "Accumulate observations of type X" step above a "-> Group ..." line)
        %   which restated the type twice and doubled the height of every block.
        %
        %   The index is {SYSTEM, TYPE}. Both matter: with several systems the
        %   type number alone is ambiguous, and a hardcoded system index made
        %   system 2's lines claim to belong to system 1.
        %
        %   INPUT
        %     sysIdx   1x1 double     the system's ordinal
        %     typeIdx  1x1 double     the type's ordinal within that system
        %     obsType  1xN char|string  observation type, quoted in the output
        %     sysName  1xN char|string  the system's readable name, quoted
        %     extra    1xN char       optional middle column (e.g. '1400 obs');
        %                             '' leaves the column blank
        %     seconds  1x1 double | []  elapsed time; [] omits the time column
        %   OUTPUT
        %     s        1xM char       one line, newline-terminated
        %
        %   Created 2026-07-19 (CEST, UTC+02:00) with assistance from Claude
        %   (Anthropic), model Claude Opus 4.8.
            body = sprintf('  Type {%d,%d}: "%s" of "%s"', ...
                           sysIdx, typeIdx, char(obsType), char(sysName));
            if nargin < 5 || isempty(extra), extra = ''; end
            % sprintf's '*' width, not pad(): pad returns a STRING scalar and
            % is a needless dependency for something sprintf already does.
            if nargin < 6 || isempty(seconds)
                s = sprintf('%-*s%s\n', lsqvce.util.verbosity.DESCCOL, body, ...
                            char(extra));
            else
                s = sprintf('%-*s%-*s%s\n', ...
                    lsqvce.util.verbosity.DESCCOL,  body, ...
                    lsqvce.util.verbosity.EXTRACOL, char(extra), ...
                    lsqvce.util.verbosity.duration(seconds));
            end
        end

        function s = streamLine(obsType, sysName, rest)
        %STREAMLINE  One stream's row in a per-stream block.
        %   '  "posC" of "demo-CD-pair"     sigma = 1.0342   white = yes'
        %
        %   The shared shape for every per-stream report -- Stochastic model,
        %   Contribution, Convergence -- so the four blocks of one pass read as
        %   one document rather than four formats. Only the trailing text
        %   differs between them.
        %
        %   INPUT
        %     obsType 1xN char|string  observation type, quoted
        %     sysName 1xN char|string  system name, quoted
        %     rest    1xN char         whatever that block reports
        %   OUTPUT
        %     s       1xM char         one line, newline-terminated
        %
        %   Created 2026-07-19 (CEST, UTC+02:00) with assistance from Claude
        %   (Anthropic), model Claude Opus 4.8.
            lbl = sprintf('"%s" of "%s"', char(obsType), char(sysName));
            s   = lsqvce.util.verbosity.labelledRow(lbl, rest);
        end

        function s = labelledRow(label, rest)
        %LABELLEDROW  One row of a per-stream/per-group block: label then value.
        %   '  "posC" of "demo-CD-pair"          sigma = 1.0342   white = yes'
        %
        %   The single indent-and-column rule behind EVERY such block --
        %   Stochastic model, Contribution and Convergence. They used to differ:
        %   Convergence indented its rows by the 21-character stage-label column
        %   while the other two started at 2, so the three tables of one pass
        %   stepped in and out. Routing them all through here is what makes them
        %   line up, and makes it impossible for one to drift again.
        %
        %   STREAMLINE builds the '"type" of "system"' label and calls this;
        %   callers whose label is not a stream (a VCE group can span systems)
        %   call this directly with whatever label they have.
        %
        %   INPUT
        %     label  1xN char|string   the row's identity
        %     rest   1xN char|string   what this block reports about it
        %   OUTPUT
        %     s      1xM char          one row, newline-terminated
        %
        %   Created 2026-07-19 (CEST, UTC+02:00) with assistance from Claude
        %   (Anthropic), model Claude Opus 4.8.
            s = sprintf('  %-*s%s\n', lsqvce.util.verbosity.STREAMCOL, ...
                        char(label), char(rest));
        end

        function s = fieldBlock(title, fields)
        %FIELDBLOCK  A titled, one-field-per-line block with a leading blank line.
        %   'System 1 "demo-CD-pair"' / '  observations   50400' / ...
        %
        %   INPUT
        %     title   1xN char     the header line
        %     fields  Nx2 cell     {name, value} pairs, both char, in order
        %   OUTPUT
        %     s       1xM char     the block, newline-terminated
        %
        %   Created 2026-07-19 (CEST, UTC+02:00) with assistance from Claude
        %   (Anthropic), model Claude Opus 4.8.
            w = 0;
            for i = 1:size(fields,1)
                w = max(w, numel(fields{i,1}));
            end
            s = sprintf('\n%s\n', char(title));
            for i = 1:size(fields,1)
                s = [s, sprintf('  %-*s   %s\n', w, fields{i,1}, fields{i,2})]; %#ok<AGROW>
            end
            % TRAILING blank line too: whatever follows a system's summary is
            % the NEXT system's work, and without the gap the two run together.
            s = [s, newline];
        end


        function s = duration(seconds)
        %DURATION  Seconds as HH:MM:SS.mmm, so timings align in a column.
        %   Chosen over a bare '%.2f s' because a run spans milliseconds (one
        %   small group) to hours (a full accumulation), and a fixed-width
        %   clock format keeps both readable in the same column.
        %
        %   INPUT   seconds 1x1 double   elapsed time, >= 0 (negative is clamped
        %                                to 0; non-finite renders as --:--:--.---)
        %   OUTPUT  s       1xN char     e.g. '00:00:00.519', '01:23:45.006'
            if ~isscalar(seconds) || ~isnumeric(seconds) || ~isfinite(seconds)
                s = '--:--:--.---';
                return
            end
            t  = max(0, double(seconds));
            hh = floor(t/3600);
            mm = floor((t - 3600*hh)/60);
            ss = floor(t - 3600*hh - 60*mm);
            ms = round((t - floor(t))*1000);
            if ms == 1000        % rounding carried into the next second
                ms = 0;  ss = ss + 1;
                if ss == 60, ss = 0; mm = mm + 1; end
                if mm == 60, mm = 0; hh = hh + 1; end
            end
            s = sprintf('%02d:%02d:%02d.%03d', hh, mm, ss, ms);
        end
    end
end
