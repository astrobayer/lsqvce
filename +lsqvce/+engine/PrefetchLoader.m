classdef PrefetchLoader < handle
%LSQVCE.ENGINE.PREFETCHLOADER  Overlap loader I/O with accumulation (ROADMAP 3.3).
%
%   pf  = lsqvce.engine.PrefetchLoader(dataFcn, dayList);
%   acc.addSystem(sysID, nl, pf.handle(), dayList);
%
%   While the accumulator consumes day d, a background worker
%   (parfeval + backgroundPool) already loads day d+1 -- the whole day, all
%   obsTypes and arcs, in one future. Worth it exactly when the loader is
%   I/O-BOUND (real L1B reads, RINEX parsing); for a cheap in-memory loader
%   the copy overhead loses. Profile first.
%
%   MEMORY
%     Exactly ONE extra day is ever resident (the prefetched one). The
%     depth-one policy is deliberate: deeper pipelines multiply peak memory
%     for no additional overlap once loading is slower than accumulating.
%
%   CORRECTNESS CONTRACT
%     * The wrapped loader must be DETERMINISTIC (the toolbox already
%       requires this for the residuals pass) and worker-safe -- it runs on a
%       background worker. persistent caches inside it are worker-local and
%       fine.
%     * Served data is byte-identical to a direct call: the future's result
%       IS a dataFcn result, cached and handed out. Nothing is recomputed or
%       approximated. Verified serially in dev/tests/test_v4.m by comparing full
%       accumulations through the wrapper against the bare loader.
%     * Out-of-order access (the redundancy pass under parfor visits arcs in
%       task order) simply MISSES the cache and falls through to a direct
%       call -- correct, just not accelerated. This wrapper targets the
%       ACCUMULATION pass, which is strictly day-ordered.
%
%   GRACEFUL DEGRADATION
%     No Parallel Computing Toolbox, no backgroundPool, or 'Enable' false ->
%     pure passthrough. The wrapper never becomes a hard dependency.
%
%   NOT tested on a live pool in the build environment (no MATLAB+PCT
%   there): the passthrough and cache-serving logic is what the test suite
%   pins down; confirm the actual overlap with a profile run on real
%   hardware.
%
%   OPTIONS (name-value)
%     'Enable'   'auto' (default: on iff backgroundPool is available) |
%                true | false
%
%   ERRORS
%     lsqvce:prefetch:badLoader     the wrapped loader is not a struct of closures
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:prefetch:badEnable                  'Enable' must be 'auto', true or false
%     lsqvce:prefetch:disabled                   [warning] Background prefetch could not be started (X
%     lsqvce:prefetch:futureFailed               [warning] The prefetched day X failed on the background pool (X
%     lsqvce:prefetch:noPool                     [warning] Prefetch requested but no backgroundPool is available; running as a p...
%
%   See also LSQVCE.ENGINE.BLOCKACCUMULATOR, LSQVCE.ENGINE.SYNTHLOADER, LSQVCE.ENGINE.MAPREDUCE
%
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (c) c2 (F-07):
%   typed `arguments` signatures, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   Updated 2026-08-18 (CEST, UTC+02:00): pool availability now probed by
%   the ONE probe lsqvce.util.hasBackgroundPool (probe-by-USE; the private
%   exist()-based poolOK was false-negative under `matlab -batch` and made
%   prefetch silently degrade to passthrough in batch runs - WP-30
%   postflight finding). With assistance from Claude (Anthropic), model
%   Claude Fable 5.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    properties (SetAccess = private)
        dataFcn                 % the wrapped loader
        dayList                 % day sequence, defines "next day"
        enabled = false         % background prefetch active?
        hits    = 0             % served from prefetch
        misses  = 0             % fell through to a direct call
    end

    properties (Access = private)
        fut      = []           % pending parfeval future ([] if none)
        futDay   = []           % day the future is loading
        cacheDay = []           % day currently cached
        cache    = []           % containers.Map arcKey -> arc struct
    end

    methods
        function o = PrefetchLoader(dataFcn, dayList, opt)
        %PREFETCHLOADER  Construct a loader that prefetches the next day in the
        %   background (parfeval) while the current one is consumed. See classdef
        %   header for options. INPUT dataFcn 4-arg loader; dayList 1xD days.
        arguments
            dataFcn
            dayList
            opt.Enable = 'auto'
        end
            narginchk(2, inf);
            if ~isa(dataFcn, 'function_handle')
                error('lsqvce:prefetch:badLoader', ...
                      'dataFcn must be a function handle.');
            end
            % WP-39 c2: typed `arguments` signature. 'Enable' keeps its
            % v1 contract "auto" | true | false, so it stays UNTYPED in
            % the block (a mustBeMember would kill the logical form) and
            % is validated here exactly as the body always consumed it.
            en = opt.Enable;
            if isstring(en), en = char(en); end
            if ~(islogical(en) && isscalar(en)) && ~strcmp(en, 'auto')
                error('lsqvce:prefetch:badEnable', ...
                      '''Enable'' must be ''auto'', true or false.');
            end

            o.dataFcn = dataFcn;
            o.dayList = dayList(:).';

            if islogical(en)
                o.enabled = en && lsqvce.util.hasBackgroundPool();
                if en && ~o.enabled
                    lsqvce.util.log.warn('lsqvce:prefetch:noPool', ...
                        ['Prefetch requested but no backgroundPool is ' ...
                         'available; running as a passthrough.']);
                end
            else
                o.enabled = lsqvce.util.hasBackgroundPool();
            end
        end

        function f = handle(o)
        %HANDLE  The @(d,ot,a,uo) function handle to give the accumulator.
            f = @(d, ot, a, uo) o.fetch(d, ot, a, uo);
        end

        function out = fetch(o, d, ot, a, uo)
            % meta calls: pass through, and use the moment to trigger the
            % prefetch of the NEXT day (we now know a new day is starting)
            if (ischar(d) || isstring(d)) && strcmp(d, 'meta')
                out = o.dataFcn(d, ot, a, uo);
                if o.enabled
                    o.kick(ot, uo);            % ot carries the day here
                end
                return
            end

            if ~o.enabled
                out = o.dataFcn(d, ot, a, uo);
                return
            end

            % day switch: harvest the future if it holds this day
            if ~isequal(o.cacheDay, d)
                o.harvest(d);
            end

            key = sprintf('%s|%s', string(ot), string(a));
            if isequal(o.cacheDay, d) && ~isempty(o.cache) && o.cache.isKey(key)
                out = o.cache(key);
                o.hits = o.hits + 1;
            else
                out = o.dataFcn(d, ot, a, uo);   % miss: direct, still correct
                o.misses = o.misses + 1;
            end
        end

        function delete(o)
        %DELETE  Destructor: cancel any in-flight prefetch future.
            if ~isempty(o.fut)
                try cancel(o.fut); catch, end
            end
        end
    end

    methods (Access = private)
        function kick(o, dNow, uo)
        %KICK  Submit the future for the day AFTER dNow, unless already done.
            k = find(o.dayList == dNow, 1);
            if isempty(k) || k == numel(o.dayList), return, end
            dNext = o.dayList(k+1);
            if isequal(o.futDay, dNext) || isequal(o.cacheDay, dNext)
                return
            end
            if ~isempty(o.fut)
                try cancel(o.fut); catch, end
            end
            try
                o.fut = parfeval(backgroundPool, ...
                                 @lsqvce.engine.PrefetchLoader.loadWholeDay, 1, ...
                                 o.dataFcn, dNext, uo);
                o.futDay = dNext;
            catch ME
                o.fut = []; o.futDay = [];
                o.enabled = false;             % degrade once, and SAY WHY (R5 M1)
                lsqvce.util.log.warn('lsqvce:prefetch:disabled', ...
                    ['Background prefetch could not be started (%s: %s); ' ...
                     'the loader continues serially.'], ME.identifier, ME.message);
            end
        end

        function harvest(o, d)
        %HARVEST  Collect a completed prefetch for day d into the cache (drops any
        %   stale future/cache first); silent no-op if the future is missing/mismatched.
            o.cacheDay = [];  o.cache = [];
            if isempty(o.fut) || ~isequal(o.futDay, d), return, end
            try
                o.cache    = fetchOutputs(o.fut);
                o.cacheDay = d;
            catch ME
                o.cache = []; o.cacheDay = [];  % the day is re-read serially
                lsqvce.util.log.warn('lsqvce:prefetch:futureFailed', ...
                    ['The prefetched day %s failed on the background pool ' ...
                     '(%s: %s); it is read again in the foreground.'], ...
                    string(d), ME.identifier, ME.message);
            end
            o.fut = []; o.futDay = [];
        end
    end

    methods (Static, Hidden)
        function C = loadWholeDay(dataFcn, d, uo)
        %LOADWHOLEDAY  Runs ON THE BACKGROUND WORKER: one whole day.
            C = containers.Map('KeyType','char','ValueType','any');
            meta = dataFcn('meta', d, [], uo);
            for ot = 1:numel(meta.obsTypes)
                uo.obsType = meta.obsTypes{ot};
                arcs = meta.arcs{ot};
                for a = arcs(:).'
                    uo.day = d;  uo.arc = a;
                    C(sprintf('%s|%s', string(uo.obsType), string(a))) = ...
                        dataFcn(d, uo.obsType, a, uo);
                end
            end
        end
    end
end
