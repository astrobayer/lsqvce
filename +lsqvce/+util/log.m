classdef log
%LSQVCE.UTIL.LOG  Toolbox-wide logging: one file identifier, one screen flag.
%
%   The whole lsqvce toolbox routes its status text and its warnings through
%   this one place, so a caller can capture everything to a log file and decide
%   independently whether it also appears on screen.
%
%   R5 M1 (2026-09-02): warnOnce(id, ...) fires an identifier once per run;
%   inventory() returns every identifier warn/warnOnce raised since
%   inventory('reset') -- lsqvce.fit resets at start and stores the list in
%   R.meta.warnings, so a run's warnings are a RESULT, not scrollback.
%
%   TWO SETTINGS (persistent for the MATLAB session)
%     logfid   a FILE IDENTIFIER that receives ALL output ("standard output").
%              Defaults to 1 (the MATLAB command window). Set it to the fid
%              returned by fopen(...,'a') to capture to a file, e.g.
%                  fid = fopen('run.log','a');
%                  lsqvce.util.log.set('fid', fid, 'screen', false);
%              Everything the toolbox prints, and every warning it raises, is
%              then written to run.log.
%     logflag  a LOGICAL. When true (default), output ALSO appears on screen
%              (fid 1). When false, output goes ONLY to logfid -- nothing on
%              screen. If logfid is itself 1, screen output happens regardless
%              (writing to 1 IS the screen), so set a real file fid when you
%              want a quiet run.
%
%   USAGE INSIDE THE TOOLBOX (replaces bare fprintf / warning)
%     lsqvce.util.log.msg('processed %d arcs\n', n);         % status line
%     lsqvce.util.log.warn('lsqvce:foo:bar', 'r_g = %g low', r);  % warning
%
%   The msg/warn calls take printf-style format + args, exactly like fprintf and
%   warning, so converting a call site is a one-token change. warn ALSO raises a
%   real MATLAB warning (so lastwarn, warning state, and 'error'-on-warning all
%   keep working) and MIRRORS the text into logfid; it appears on screen only if
%   logflag is true (when false, warn suppresses the on-screen display but still
%   registers the warning and writes it to the log file).
%
%   CONFIGURE
%     lsqvce.util.log.set('fid', logfid, 'screen', logflag)   set either/both
%     [fid, screen] = lsqvce.util.log.get()                   read current state
%     lsqvce.util.log.reset()                                 back to (1, true)
%
%   THREAD SAFETY / PARALLEL
%     The (fid, screen) state is persistent PER PROCESS. On parallel workers
%     (parfor/spmd) each worker has its own copy at its default (1, true) unless
%     you set it inside the worker; without help a worker's output prints on the
%     worker's own stream and is LOST to the client's log file.
%
%     WORKER-AWARE LOGGING (parfor-safe sink). On the CLIENT, before a parfor
%     pass, call
%         q = lsqvce.util.log.startParallel();
%     to install a parallel.pool.DataQueue whose callback drains every worker
%     message into the client's current (fid, screen). A worker does not see
%     the client's persistent state, so the queue must be HANDED to it: the
%     parfor body starts with lsqvce.util.log.attachWorker(q) (R7-41 -- without
%     it the sink was inert). Worker calls to msg/warn then SEND their text to
%     the queue; the client's callback writes it to the one log file, in
%     arrival order, a warnOnce id once per run across all workers. The task
%     returns inventory() with its result and the client merges it
%     (mergeInventory), so worker warnings reach R.meta.warnings in task
%     order. After the pass, call
%         lsqvce.util.log.stopParallel();
%     If the Parallel Computing Toolbox is absent or no pool exists, startParallel
%     is a safe no-op (q = []) and attachWorker still collects the inventory.
%     LSQVCE.ENGINE.MAPREDUCE does all of this for the toolbox's own passes.
%
%   EXAMPLE
%     fid = fopen('vce.log', 'w');
%     lsqvce.util.log.set('fid', fid, 'screen', true);   % file AND screen
%     R = lsqvce.fit(task, Estimate="noise");            % all output captured
%     lsqvce.util.log.reset();
%     fclose(fid);
%
%   EXAMPLE (parfor-safe)
%     fid = fopen('run.log','w');  lsqvce.util.log.set('fid', fid);
%     q = lsqvce.util.log.startParallel();               % client installs sink
%     ids = cell(nTask, 1);
%     parfor k = 1:nTask
%         lsqvce.util.log.attachWorker(q);                % worker: hand it the sink
%         lsqvce.util.log.msg('worker did task %d\n', k);  % -> client's run.log
%         ids{k} = lsqvce.util.log.inventory();          % this task's warnings
%     end
%     for k = 1:nTask, lsqvce.util.log.mergeInventory(ids{k}); end
%     lsqvce.util.log.stopParallel();  lsqvce.util.log.reset();  fclose(fid);
%
%   ACCURACY
%     Not applicable: it routes text and warnings (file identifier, screen
%     flag, the worker queue); no number is formed. WP-139.
%
%   ERRORS
%     lsqvce:util:log:badKey        an unknown log option key
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:log:noParallelQueue                 [warning] Worker log messages cannot be forwarded to the client (X); they are p...
%
%   See also FPRINTF, WARNING, PARALLEL.POOL.DATAQUEUE, LSQVCE.ENGINE.ADJUSTVCEJOINT,
%            LSQVCE.ENGINE.MAPREDUCE
%
%   Updated 2026-11-03 (CET, UTC+01:00), WP-162: the examples fit a
%   lsqvce.TaskDefinition (lsqvce.fit takes no struct since 7.0.0). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-139: the help states its
%   accuracy, "Not applicable" with its reason. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-118 part 2 (R7-41): the
%   parallel sink reaches the workers (attachWorker), worker warnings reach
%   the run inventory (mergeInventory), and a warnOnce id is shown once per
%   run, not once per worker. With assistance from Claude (Anthropic), model
%   Claude Opus 5.5.
%   Created 2026-07-16 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.

    methods (Static)

        function set(varargin)
        %SET  Configure the log file identifier and/or the screen flag.
        %   lsqvce.util.log.set('fid', logfid, 'screen', logflag)
        %   Either name/value pair may be given alone.
            p = inputParser;  p.FunctionName = 'lsqvce.util.log.set';
            p.addParameter('fid',    [], @(v) isnumeric(v) && isscalar(v));
            p.addParameter('screen', [], @(v) islogical(v) && isscalar(v));
            p.parse(varargin{:});
            o = p.Results;
            if ~isempty(o.fid),    lsqvce.util.log.state('fid',    o.fid);    end
            if ~isempty(o.screen), lsqvce.util.log.state('screen', o.screen); end
        end

        function [fid, screen] = get()
        %GET  Return the current file identifier and screen flag.
            fid    = lsqvce.util.log.state('fid');
            screen = lsqvce.util.log.state('screen');
        end

        function reset()
        %RESET  Restore defaults: fid = 1 (command window), screen = true.
        %   Does NOT remove a parallel sink; call stopParallel for that.
            lsqvce.util.log.state('fid',    1);
            lsqvce.util.log.state('screen', true);
        end

        function q = startParallel()
        %STARTPARALLEL  Install a DataQueue sink so parfor workers' log output is
        %   drained into the CLIENT's current (fid, screen). Call on the CLIENT
        %   before a parfor pass. Returns the queue handle (or [] if unavailable).
        %   Safe no-op without the Parallel Computing Toolbox or an open pool.
            q = [];
            % Need the PCT and the DataQueue class to exist.
            if isempty(ver('parallel')) || ...
               exist('parallel.pool.DataQueue', 'class') ~= 8
                return
            end
            % Capture the client's CURRENT sink; the callback writes there.
            [cfid, cscreen] = lsqvce.util.log.get();
            try
                q = parallel.pool.DataQueue;
                afterEach(q, @(m) lsqvce.util.log.drain(cfid, cscreen, m));
            catch ME
                % No queue: worker messages are emitted locally instead, and
                % the client is told once why (R5 M1).
                warning('lsqvce:log:noParallelQueue', ...
                    ['Worker log messages cannot be forwarded to the client ' ...
                     '(%s); they are printed on the workers.'], ME.identifier);
                q = [];  return
            end
            lsqvce.util.log.state('queue', q);
            % R7-41: the display dedup of this pass starts from what the client
            % already fired (a SNAPSHOT: the drain runs asynchronously, possibly
            % after mergeInventory, so it must not read the live inventory --
            % it would then suppress a worker warning's only display)
            lsqvce.util.log.state('seen', []);
            fired = lsqvce.util.log.state('fired');
            for k = 1:numel(fired), lsqvce.util.log.state('seen', fired{k}); end
        end

        function onWorker = attachWorker(q)
        %ATTACHWORKER  Worker side of the parallel sink (R7-41). Call at the
        %   START of every parfor task / spmd block with the queue that
        %   startParallel returned (may be []): the worker forwards its log
        %   text to the client and starts a fresh warning inventory for THIS
        %   task, which the task returns (inventory()) and the client merges
        %   (mergeInventory) -- so worker warnings reach R.meta.warnings. Before,
        %   the queue lived only in the client's persistent state, workers never
        %   saw it, and their warnings were missing from the inventory.
        %   No-op on the client (a parfor without a pool runs there), so the
        %   run's own inventory is never cleared by it.
        %   OUTPUT  onWorker  1x1 logical  true when running on a worker
            onWorker = exist('getCurrentTask', 'file') == 2 && ~isempty(getCurrentTask());
            if ~onWorker, return, end
            lsqvce.util.log.state('queue', q);
            lsqvce.util.log.state('fired', []);
        end

        function mergeInventory(ids)
        %MERGEINVENTORY  Client side (R7-41): add the warning identifiers a
        %   worker task fired to this run's inventory, first-fired order kept.
        %   INPUT  ids  1xK cell  as returned by inventory() on the worker
            for k = 1:numel(ids)
                lsqvce.util.log.state('fired', ids{k});
            end
        end

        function stopParallel()
        %STOPPARALLEL  Remove the DataQueue sink installed by startParallel.
        %   Call on the CLIENT after the parfor pass. Idempotent.
            lsqvce.util.log.state('queue', []);
        end

        function msg(varargin)
        %MSG  Write a status line to logfid, and to the screen if logflag.
        %   lsqvce.util.log.msg(fmt, ...)   -- printf-style, like fprintf.
        %   If a parallel sink is installed (startParallel) and this call runs on
        %   a worker, the text is SENT to the client's queue instead of printed.
            txt = lsqvce.util.log.render(varargin{:});
            if lsqvce.util.log.sendToQueue('msg', txt), return, end
            [fid, screen] = lsqvce.util.log.get();
            lsqvce.util.log.emit(fid, screen, txt);
        end

        function fired = warnOnce(id, varargin)
        %WARNONCE  Raise a warning the FIRST time its identifier fires in this
        %   run; later calls with the same id are silent (the count belongs in
        %   the skip record, not in the warning stream -- R5 M1 doctrine
        %   "warnings once per run"). Returns true when the warning was raised.
        %   lsqvce.util.log.warnOnce(id, fmt, ...)
            fired = ~any(strcmp(lsqvce.util.log.state('fired'), char(id)));
            if fired
                lsqvce.util.log.raise(id, true, varargin{:});
            end
        end

        function ids = inventory(reset)
        %INVENTORY  The warning identifiers raised through warn/warnOnce since
        %   the last reset, first-fired order, each once (R5 M1).
        %   ids = lsqvce.util.log.inventory()        read
        %   lsqvce.util.log.inventory('reset')       clear (fit does this at start)
            if nargin > 0 && strcmpi(char(reset), 'reset')
                lsqvce.util.log.state('fired', []);
            end
            ids = lsqvce.util.log.state('fired');
        end

        function warn(id, varargin)
        %WARN  Register a MATLAB warning, mirror its text to logfid, and show it
        %   on screen only if logflag is true.
        %   lsqvce.util.log.warn(id, fmt, ...)   -- id is a warning identifier.
            lsqvce.util.log.raise(id, false, varargin{:});
        end

    end

    methods (Static, Hidden)

        function raise(id, once, varargin)
        %RAISE  warn's body; once = true when called by warnOnce (R7-41: the
        %   client shows a once-per-run warning ONCE across all workers).
        %
        %   The warning is always RAISED in the warning system (so lastwarn, the
        %   warning state, and error-on-warning keep working) and always written
        %   to logfid. Whether it also appears on SCREEN follows logflag, exactly
        %   like msg: logflag=false => file only, nothing on screen.
            txt = lsqvce.util.log.render(varargin{:});
            lsqvce.util.log.state('fired', id);       % R5 M1: the run inventory

            % On a worker with a sink installed, forward to the client (tagged as
            % a warning, with its id) and still raise a local warning so
            % worker-side state and error-on-warning keep working -- with its
            % DISPLAY off: parfor relays a displayed worker warning to the
            % client itself, which printed a once-per-run warning once per
            % worker (R7-41). The client's drain shows it, deduplicated.
            kind = 'warn';  if once, kind = 'warnOnce'; end
            if lsqvce.util.log.sendToQueue(kind, txt, id)
                st = warning('off', id);
                cw = onCleanup(@() warning(st));
                warning(id, '%s', txt);
                return
            end
            [fid, screen] = lsqvce.util.log.get();

            % Raise the real warning so lastwarn / state / error-on-warning work.
            % MATLAB prints it to the command window when warnings are displayed;
            % to honour logflag we temporarily turn DISPLAY off (verbose/backtrace
            % state), which does NOT affect lastwarn or an 'error' state.
            showOnScreen = screen || (fid == 1);   % fid==1 IS the screen
            if showOnScreen
                warning(id, '%s', txt);            % normal: goes to screen/stderr
            else
                st = warning('off', id);           % suppress THIS id's display
                cu = onCleanup(@() warning(st));    % restore on exit
                warning(id, '%s', txt);            % still sets lastwarn & state
            end

            % Mirror to the log file (unless the file IS the screen, fid==1,
            % where the warning() above already wrote it).
            %
            % SET APART DELIBERATELY. A warning is the one thing in the log that
            % must not be skimmed past, and in a long run it was a single line
            % among hundreds of status lines of the same shape. It now gets a
            % blank line either side, a marker, and WRAPPED continuation text
            % indented under that marker -- so it breaks the column rhythm of the
            % surrounding output and the eye lands on it.
            %
            % No colour: the sink may be a file, and ANSI escapes written to a
            % file are noise. On screen MATLAB's own warning() above already
            % renders in the command window's warning colour.
            if fid ~= 1
                try
                    fprintf(fid, '\n%s\n\n', ...
                        lsqvce.util.log.wrapWarning(txt));
                catch ME
                    % The log FILE is unwritable: the warning is on screen
                    % already; say the file lost it (R5 M1, no empty catch).
                    fprintf(1, '[lsqvce.util.log: warning not written to the log file: %s]\n', ...
                            ME.identifier);
                end
            end
        end

        function sent = sendToQueue(kind, txt, id)
        %SENDTOQUEUE  If a DataQueue sink is installed, send {kind, txt} to it and
        %   return true; otherwise return false. On the CLIENT the queue exists
        %   too, but afterEach drains on the client's event loop -- to avoid a
        %   message looping back, only SEND when running on a worker (spmdIndex or
        %   getCurrentTask is non-empty). Robust to the PCT being absent.
            sent = false;
            q = lsqvce.util.log.state('queue');
            if isempty(q), return, end
            % getCurrentTask exists only with PCT: a CONDITION, not a catch (R5 M1)
            onWorker = exist('getCurrentTask', 'file') == 2 && ~isempty(getCurrentTask());
            if ~onWorker, return, end
            if nargin < 3, id = ''; end
            try
                send(q, {kind, txt, char(id)});
                sent = true;
            catch ME
                sent = false;                  % the caller emits locally; the
                fprintf(1, '[lsqvce.util.log: queue send failed, %s]\n', ...
                        ME.identifier);        % reason is printed first (R5 M1)
            end
        end

        function drain(fid, screen, m)
        %DRAIN  Client-side afterEach callback: write a worker message {kind,txt}
        %   to the captured (fid, screen). Warnings are prefixed like warn does.
            id = '';
            if iscell(m) && numel(m) >= 2
                kind = m{1};  txt = m{2};
                if numel(m) >= 3, id = m{3}; end
            else
                kind = 'msg'; txt = char(m);
            end
            if strcmp(kind, 'warnOnce')
                % R7-41: once per RUN, not once per worker -- skip it when the
                % client had fired it before this pass or another worker's copy
                % was shown ('seen', seeded by startParallel)
                if any(strcmp(lsqvce.util.log.state('seen'), id))
                    return
                end
                lsqvce.util.log.state('seen', id);
                kind = 'warn';
            end
            if strcmp(kind, 'warn')
                lsqvce.util.log.emit(fid, screen, sprintf('Warning: %s\n', txt));
            else
                lsqvce.util.log.emit(fid, screen, txt);
            end
        end

        function out = wrapWarning(txt)
        %WRAPWARNING  'Warning: ...' with continuation lines indented under it.
        %   Toolbox warnings are long by design -- they say what happened AND
        %   what to do about it -- so as one unbroken line they run off the
        %   window and are hard to read back. Wrapped at a fixed width with the
        %   continuations indented, the whole message stays visible and the
        %   block reads as one unit.
        %
        %   INPUT   txt  1xN char   the warning message
        %   OUTPUT  out  1xM char   'Warning: <first line>' + indented rest
            width  = 74;
            indent = '         ';                 % under the word 'Warning:'
            words  = strsplit(strtrim(txt));
            line   = 'Warning:';
            out    = '';
            for k = 1:numel(words)
                cand = [line ' ' words{k}];
                if numel(cand) > width && ~isempty(line)
                    out  = [out line newline];    %#ok<AGROW>
                    line = [indent words{k}];
                else
                    line = cand;
                end
            end
            out = [out line];
        end

        function txt = render(varargin)
        %RENDER  sprintf the args, tolerating a pre-formatted single string.
        %   THE SINGLE-ARGUMENT CASE IS THE SUBTLE ONE. Callers write
        %       log.msg('  PER-SENSOR SUMMARY\n')
        %   expecting fprintf semantics, where the trailing \n is an ESCAPE that
        %   becomes a newline. An earlier version returned that lone argument
        %   untouched, and EMIT writes with fprintf(fid,'%s',txt), which does not
        %   expand escapes either -- so the two characters backslash-n reached
        %   the log verbatim. ~20 call sites are written that way; none was
        %   covered by a test, so the defect stayed latent (it is the same bug
        %   class the user hit by hand in the joint path).
        %
        %   sprintf(s) with ONE argument expands escapes and, having no
        %   conversion specifiers to consume, leaves the rest of the text alone
        %   -- exactly fprintf's behaviour. A literal percent must still be
        %   written '%%', which is already true for every caller because they
        %   were written against fprintf.
            if isempty(varargin)
                txt = '';
            elseif numel(varargin) == 1
                txt = sprintf(varargin{1});                 % expand escapes
            else
                txt = sprintf(varargin{:});
            end
        end

        function emit(fid, screen, txt)
        %EMIT  Write txt to fid, and to screen (fid 1) if requested and distinct.
            if isempty(txt), return, end
            wroteScreen = false;
            if ~isempty(fid) && fid >= 1
                try
                    fprintf(fid, '%s', txt);
                    wroteScreen = (fid == 1);
                catch ME
                    % bad fid -> the screen, WITH the reason, so nothing is lost
                    % and nothing is silent (R5 M1)
                    fprintf(1, '[lsqvce.util.log: file write failed, %s]\n%s', ...
                            ME.identifier, txt);
                    wroteScreen = true;
                end
            end
            if screen && ~wroteScreen
                fprintf(1, '%s', txt);
            end
        end

        function out = state(key, val)
        %STATE  Persistent (fid, screen, queue) store. state(key) reads;
        %   state(key,val) sets. 'queue' holds the parallel DataQueue or [].
            persistent FID SCREEN QUEUE FIRED SEEN
            if isempty(FID),    FID    = 1;    end
            if isempty(SCREEN), SCREEN = true; end
            if isempty(FIRED),  FIRED  = {};   end
            if isempty(SEEN),   SEEN   = {};   end
            switch key
                case 'fired'
                    % R5 M1: the run's warning INVENTORY -- every identifier
                    % lsqvce.util.log.warn raised since the last reset, in
                    % first-fired order, each once. Read into R.meta.warnings.
                    if nargin > 1
                        if isempty(val), FIRED = {};                % reset
                        elseif ~any(strcmp(FIRED, val)), FIRED{end+1} = char(val);
                        end
                    end
                    out = FIRED;
                case 'seen'
                    % R7-41: once-per-run warning ids a parallel pass already
                    % displayed (the drain's dedup; reset by startParallel)
                    if nargin > 1
                        if isempty(val), SEEN = {};
                        elseif ~any(strcmp(SEEN, val)), SEEN{end+1} = char(val);
                        end
                    end
                    out = SEEN;
                case 'fid'
                    if nargin > 1, FID = val; end
                    out = FID;
                case 'screen'
                    if nargin > 1, SCREEN = val; end
                    out = SCREEN;
                case 'queue'
                    if nargin > 1, QUEUE = val; end
                    out = QUEUE;
                otherwise
                    error('lsqvce:util:log:badKey', 'Unknown log state ''%s''.', key);
            end
        end

    end
end
