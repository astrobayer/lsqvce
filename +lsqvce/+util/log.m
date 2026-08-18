classdef log
%LSQVCE.UTIL.LOG  Toolbox-wide logging: one file identifier, one screen flag.
%
%   The whole lsqvce toolbox routes its status text and its warnings through
%   this one place, so a caller can capture everything to a log file and decide
%   independently whether it also appears on screen.
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
%     message into the client's current (fid, screen). Inside the parfor body,
%     worker calls to msg/warn detect the queue and SEND their text to it instead
%     of printing locally; the client's callback writes it to the one log file,
%     in arrival order. After the pass, call
%         lsqvce.util.log.stopParallel();
%     If the Parallel Computing Toolbox is absent or no pool exists, startParallel
%     is a safe no-op and logging stays per-process.
%
%   EXAMPLE
%     fid = fopen('vce.log', 'w');
%     lsqvce.util.log.set('fid', fid, 'screen', true);   % file AND screen
%     R = lsqvce.adjustVCE(L, days).run();         % all output captured
%     lsqvce.util.log.reset();
%     fclose(fid);
%
%   EXAMPLE (parfor-safe)
%     fid = fopen('run.log','w');  lsqvce.util.log.set('fid', fid);
%     lsqvce.util.log.startParallel();                   % client installs sink
%     parfor k = 1:nTask
%         lsqvce.util.log.msg('worker did task %d\n', k);  % -> client's run.log
%     end
%     lsqvce.util.log.stopParallel();  lsqvce.util.log.reset();  fclose(fid);
%
%   ERRORS
%     lsqvce:util:log:badKey        an unknown log option key
%
%   See also FPRINTF, WARNING, PARALLEL.POOL.DATAQUEUE, LSQVCE.ADJUSTVCE
%
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
            catch
                q = [];  return                       % pool problem -> no-op
            end
            lsqvce.util.log.state('queue', q);
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

        function warn(id, varargin)
        %WARN  Register a MATLAB warning, mirror its text to logfid, and show it
        %   on screen only if logflag is true.
        %   lsqvce.util.log.warn(id, fmt, ...)   -- id is a warning identifier.
        %
        %   The warning is always RAISED in the warning system (so lastwarn, the
        %   warning state, and error-on-warning keep working) and always written
        %   to logfid. Whether it also appears on SCREEN follows logflag, exactly
        %   like msg: logflag=false => file only, nothing on screen.
            txt = lsqvce.util.log.render(varargin{:});

            % On a worker with a sink installed, forward to the client (tagged as
            % a warning) and still raise a local warning so worker-side state and
            % error-on-warning keep working; skip the local file/screen mirror.
            if lsqvce.util.log.sendToQueue('warn', txt)
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
                catch
                end
            end
        end

    end

    methods (Static, Hidden)

        function sent = sendToQueue(kind, txt)
        %SENDTOQUEUE  If a DataQueue sink is installed, send {kind, txt} to it and
        %   return true; otherwise return false. On the CLIENT the queue exists
        %   too, but afterEach drains on the client's event loop -- to avoid a
        %   message looping back, only SEND when running on a worker (spmdIndex or
        %   getCurrentTask is non-empty). Robust to the PCT being absent.
            sent = false;
            q = lsqvce.util.log.state('queue');
            if isempty(q), return, end
            onWorker = false;
            try
                t = getCurrentTask();          %#ok<NASGU>  % [] on the client
                onWorker = ~isempty(t);
            catch
                onWorker = false;              % no PCT on this process
            end
            if ~onWorker, return, end
            try
                send(q, {kind, txt});
                sent = true;
            catch
                sent = false;                  % fall back to local emit
            end
        end

        function drain(fid, screen, m)
        %DRAIN  Client-side afterEach callback: write a worker message {kind,txt}
        %   to the captured (fid, screen). Warnings are prefixed like warn does.
            if iscell(m) && numel(m) == 2
                kind = m{1};  txt = m{2};
            else
                kind = 'msg'; txt = char(m);
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
                catch
                    % bad fid -> fall back to the screen so nothing is lost
                    fprintf(1, '%s', txt);  wroteScreen = true;
                end
            end
            if screen && ~wroteScreen
                fprintf(1, '%s', txt);
            end
        end

        function out = state(key, val)
        %STATE  Persistent (fid, screen, queue) store. state(key) reads;
        %   state(key,val) sets. 'queue' holds the parallel DataQueue or [].
            persistent FID SCREEN QUEUE
            if isempty(FID),    FID    = 1;    end
            if isempty(SCREEN), SCREEN = true; end
            switch key
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
