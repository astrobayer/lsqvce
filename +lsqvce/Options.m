classdef Options
%LSQVCE.OPTIONS  The one settings object for LSQVCE.FIT (v5: 14 knobs + 5 groups).
%
%   opts = lsqvce.Options()                       every default
%   opts = lsqvce.Options(Name, Value, ...)       override some
%   opts.Name = Value                             change one later (validated)
%   R = lsqvce.fit(data, Name, Value, ...)        the same names, straight on fit
%
%   WHAT TO ESTIMATE, HOW TO RUN. Everything ABOUT A STREAM (its law, whether
%   its law / level / arc factors are estimated) is declared on the stream
%   (data.systems(k).covariance.<type>: type, sigma, estimate, arcFactors);
%   this object holds only RUN policy. R5 M5 (WP-87, decision R5-D2) cut it
%   from 62 flat knobs to the 14 below plus five typed groups; the migration
%   from every v4 name is in RELEASE_NOTES.md and HANDOVER §12.10.
%
%   TOP LEVEL
%     Estimate   "none" | "noise"   "none": one GLS solve under the declared
%                                   laws. "noise": the joint (Ellmer) variance
%                                   component estimator on every stream whose
%                                   spec allows it (estimate ~= "none"), with
%                                   the policy below. Replaces Model: the v4
%                                   spellings are Estimate="noise" plus
%                                     "ellmer" -> Pool as set, ArcFactors true
%                                     "vce"    -> Pool "stream", ArcFactors false
%                                     "pooled" -> Pool "all",    ArcFactors false
%     Pool       "all" | "obsType" | "stream"   which streams SHARE one estimated
%                                   spectral shape (every stream keeps its own
%                                   level regardless); default "all"
%     ArcFactors logical            per-arc variance factors on (Ellmer 6.5.4),
%                                   default true; false = one level per stream
%     Robust     "none" | "huber" | "igg3" | "bisquare"   IRLS on the same
%                                   loop, every kind (details in RobustOpts)
%     MaxPasses  double  10         the ONE pass budget (v4: MaxOuter / JointMaxIter)
%     MinPasses  double  1          passes run even if the tolerance is met
%                                   (GRAVFIRE's minimumIterations)
%     Tolerance  double  1e-2       stop when max|alpha2 - 1| over frequencies
%                                   and arcs is below it (v4: JointTol)
%     Backend    "qr" | "normal"    square-root vs normal equations
%     Probes     double  100        Hutchinson probes for the per-lag
%                                   redundancies (0 = the exact route)
%     Seed       [] | double        RNG seed for the probes (recorded)
%     CacheDays  double  1          built days kept resident per system; Inf =
%                                   the whole calendar (one build per run)
%     Verbose    logical | 0..2     screen messages (0 silent, 1 progress,
%                                   2 detail; level 3 was merged into 2 in
%                                   2026-07 and lsqvce.util.verbosity refuses it)
%     LogFile    string  ""         mirror the messages to this file ("" = none)
%     OnPass     [] | @(ctx)        YOUR inspection hook, called at the end of
%                                   EVERY pass of every iterative run (the
%                                   variance-component loop, IRLS, Gauss-Newton)
%                                   with ONE ctx struct: .driver .pass .x .Cxx
%                                   .sigma0, the loop state (noise: .C .Cgroups
%                                   .sigma2 .arcKeys .streams .sigma .white
%                                   .whitePval .maxDev .skips; IRLS: .factor
%                                   .keys .scale .dx; Gauss-Newton: .vPv .dxRel
%                                   .lambda), .acc (the LIVE accumulator: call
%                                   lsqvce.engine.residuals(ctx.acc, ...) or any
%                                   engine function inside the hook), .data (your
%                                   declaration as passed to fit) and
%                                   .processingOptions (per system, as declared).
%                                   Plot, drawnow, save wherever you like; return
%                                   nothing to continue, false to stop after this
%                                   pass (stopReason "userStop"). A throwing hook
%                                   is NOT swallowed (lsqvce:callUser:threw with
%                                   your error as cause). Replaces v4's
%                                   plotPostModel / plotPostSolve seams.
%
%   GROUPS (structs; set fields one at a time or give a whole struct)
%     Noise        Grid (60 | "auto")  frequencies of the estimated shape
%                  Whitener "toeplitz"|"ar"|"arma"  how a PSD law whitens (v4 PsdMode)
%                  MeanType "robust"|"plain"        arc-factor normalisation
%                  Smooth 0                         moving average on the PSD
%                  FreezeRedundancy 1               components below it hold
%                  SkipThreshold 0.2                skip RATE that refuses the pass
%     RobustOpts   Tuning [], MaxIter 30, Tol 1e-8
%     Diagnostics  Alpha 1e-3, Bonferroni true, Beta 0.2, RMin 0.01, RWeak 0.1,
%                  Test "auto"|"pope"|"baarda"|"mad", RobustCheck true,
%                  Store "all"|"flagged"|"none", NWorst 1000
%     Parallel     Mode "serial"|"parfor"|"spmd", DistributeDays false, Threshold 8
%     GaussNewton  TolF 1e-10, Lambda0 1e-6     (x0 problems; MaxPasses / Tolerance
%                  are the budget and the step tolerance)
%
%   REGULARISATION is not an option: it is DECLARED on the data, as
%   pseudo-observations in the user's own columns -- data.priors(k) on the
%   common parameters (L in the columns of A) and data.systems(k).priorDay on
%   the local ones (L in the columns of B), each with x0, sigma, from, group
%   (R5-D8). The v4 Reg / RegLambda / RegMatrix / DatumB / Sigma0 / Sparse
%   knobs are gone (HANDOVER §12.10).
%
%   METHODS
%     engineArgs(target)   the Name-Value list a driver takes (validated)
%     toStruct()           every public property as a struct (provenance)
%     disp()               grouped display, engine defaults marked *
%
%   ERRORS
%     MATLAB:validators:*          a typed property receives a value outside
%                                  its type / set / range (at assignment)
%     lsqvce:Options:badGrid       Noise.Grid is neither "auto" nor an integer >= 2
%     lsqvce:Options:badVerbose    Verbose is neither logical nor an integer 0..2
%     lsqvce:Options:badOnPass     OnPass is neither [] nor a function handle
%     lsqvce:Options:badGroupField a group struct carries an unknown field
%     lsqvce:Options:badTarget     engineArgs target is unknown
%
%   EXAMPLE
%     opts = lsqvce.Options(Estimate="noise", Pool="obsType", Robust="huber");
%     opts.Noise.Grid = 2160;  opts.Parallel.Mode = "parfor";
%     R = lsqvce.fit(data, opts);
%
%   See also LSQVCE.FIT, LSQVCE.RESULTS, the Manual, chapter 3, RELEASE_NOTES.md
%
%   Created 2026-09-03 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M5, WP-87 (v5 interface).

    properties
        Estimate   (1,1) string {mustBeMember(Estimate, ["none","noise"])} = "none"
        Pool       (1,1) string {mustBeMember(Pool, ["all","obsType","stream"])} = "all"
        ArcFactors (1,1) logical = true
        Robust     (1,1) string {mustBeMember(Robust, ["none","huber","igg3","bisquare"])} = "none"
        MaxPasses  (1,1) double {mustBeInteger, mustBePositive} = 10
        MinPasses  (1,1) double {mustBeInteger, mustBePositive} = 1
        Tolerance  (1,1) double {mustBePositive} = 1e-2
        Backend    (1,1) string {mustBeMember(Backend, ["qr","normal"])} = "qr"
        Probes     (1,1) double {mustBeInteger, mustBeNonnegative} = 100
        Seed       double = []
        CacheDays  (1,1) double {mustBePositive} = 1
        Verbose    = false
        LogFile    (1,1) string = ""
        OnPass     = []

        Noise       (1,1) struct = lsqvce.Options.noiseDefaults()
        RobustOpts  (1,1) struct = struct('Tuning', [], 'MaxIter', 30, 'Tol', 1e-8)
        Diagnostics (1,1) struct = lsqvce.Options.diagnosticsDefaults()
        Parallel    (1,1) struct = struct('Mode', "serial", 'DistributeDays', false, 'Threshold', 8)
        GaussNewton (1,1) struct = struct('TolF', 1e-10, 'Lambda0', 1e-6)

    end

    methods
        function obj = Options(varargin)
            for k = 1:2:numel(varargin)
                name = char(string(varargin{k}));
                obj.(name) = varargin{k+1};
            end
            obj.validate();
        end

        % group setters MERGE a partial struct onto the defaults, so
        % Options(Noise=struct('Grid', 2160)) and opts.Noise.Grid = 2160 both work
        function obj = set.Noise(obj, s),       obj.Noise       = lsqvce.Options.mergeGroup(lsqvce.Options.noiseDefaults(), s, 'Noise'); end
        function obj = set.RobustOpts(obj, s),  obj.RobustOpts  = lsqvce.Options.mergeGroup(struct('Tuning', [], 'MaxIter', 30, 'Tol', 1e-8), s, 'RobustOpts'); end
        function obj = set.Diagnostics(obj, s), obj.Diagnostics = lsqvce.Options.mergeGroup(lsqvce.Options.diagnosticsDefaults(), s, 'Diagnostics'); end
        function obj = set.Parallel(obj, s),    obj.Parallel    = lsqvce.Options.mergeGroup(struct('Mode', "serial", 'DistributeDays', false, 'Threshold', 8), s, 'Parallel'); end
        function obj = set.GaussNewton(obj, s), obj.GaussNewton = lsqvce.Options.mergeGroup(struct('TolF', 1e-10, 'Lambda0', 1e-6), s, 'GaussNewton'); end

        function validate(obj)
            g = obj.Noise.Grid;
            okAuto = (ischar(g) || (isstring(g) && isscalar(g))) && strcmp(char(g), 'auto');
            okNum  = isnumeric(g) && isscalar(g) && isfinite(g) && g >= 2 && g == round(g);
            if ~(okAuto || okNum)
                error('lsqvce:Options:badGrid', ...
                    'Noise.Grid must be "auto" or an integer >= 2 (got %s).', mat2str(g));
            end
            v = obj.Verbose;
            % R6 WP-95: the range is the one lsqvce.util.verbosity has enforced
            % since 2026-07-19 (level 3 merged into 2). Options accepted 3 and
            % the refusal then came from a driver, by another id, mid-run.
            if isnumeric(v) && ~(isscalar(v) && isfinite(v) && v >= 0 && v <= 2 && v == round(v))
                error('lsqvce:Options:badVerbose', ...
                    'Verbose must be logical or an integer level 0..2 (got %g).', v);
            end
            if ~isempty(obj.OnPass) && ~isa(obj.OnPass, 'function_handle')
                error('lsqvce:Options:badOnPass', 'OnPass must be [] or a function handle @(ctx).');
            end
            groups = {'Noise', 'RobustOpts', 'Diagnostics', 'Parallel', 'GaussNewton'};
            d = lsqvce.Options.groupDefaults();
            for k = 1:numel(groups)
                extra = setdiff(fieldnames(obj.(groups{k})), fieldnames(d.(groups{k})));
                if ~isempty(extra)
                    error('lsqvce:Options:badGroupField', ...
                        '%s has no field ''%s'' (fields: %s).', groups{k}, extra{1}, ...
                        strjoin(fieldnames(d.(groups{k})).', ', '));
                end
            end
        end

        function args = engineArgs(obj, target)
            obj.validate();
            target = char(string(target));
            P = obj.Parallel;  N = obj.Noise;
            par = {'Parallel', string(P.Mode), 'DistributeDays', logical(P.DistributeDays), ...
                   'ParallelThreshold', P.Threshold};
            switch target
                case 'adjustBlock'
                    args.nv = [{'Backend', obj.Backend}, par, ...
                               {'Verbose', logical(obj.Verbose), 'Diagnostics', obj.diagnoseNV()}];
                case 'adjustNonlinear'
                    args.nv = [{'MaxIter', obj.MaxPasses, 'TolX', obj.Tolerance, ...
                               'TolF', obj.GaussNewton.TolF, 'Lambda0', obj.GaussNewton.Lambda0, ...
                               'Backend', obj.Backend}, par, {'Verbose', logical(obj.Verbose), 'Diagnostics', obj.diagnoseNV(), ...
                               'OnPass', obj.OnPass}];
                case 'adjustRobust'
                    args.nv = [{'Backend', obj.Backend, 'Influence', obj.Robust, ...
                               'Tuning', obj.RobustOpts.Tuning, 'MaxIter', obj.RobustOpts.MaxIter, ...
                               'Tol', obj.RobustOpts.Tol}, par, {'Verbose', logical(obj.Verbose), 'Diagnostics', obj.diagnoseNV(), ...
                               'OnPass', obj.OnPass}];
                case 'adjustVCEjoint'
                    args.nv = {'NumLags', N.Grid, 'MaxIter', obj.MaxPasses, ...
                               'Tol', obj.Tolerance, 'Probes', obj.Probes, 'Seed', obj.Seed, ...
                               'SkipThreshold', N.SkipThreshold, 'Pool', obj.Pool, ...
                               'MeanType', N.MeanType, 'Verbose', logical(obj.Verbose), ...
                               'OnPass', obj.OnPass};
                case 'diagnose'
                    args.nv = obj.diagnoseNV();
                otherwise
                    error('lsqvce:Options:badTarget', ...
                        ['engineArgs target must be one of adjustBlock, adjustVCEjoint, ' ...
                         'adjustRobust, adjustNonlinear, diagnose (got ''%s'').'], target);
            end
        end

        function s = toStruct(obj)
            mc = metaclass(obj);
            s  = struct();
            for k = 1:numel(mc.PropertyList)
                pr = mc.PropertyList(k);
                if strcmp(pr.GetAccess, 'public') && ~pr.Constant && ~pr.Hidden
                    s.(pr.Name) = obj.(pr.Name);
                end
            end
        end

        function disp(obj)
            groups = lsqvce.Options.groupTable();
            d = lsqvce.Options();                       % defaults for the marker
            fprintf('  lsqvce.Options  (* = default)\n');
            for g = 1:size(groups, 1)
                fprintf('\n  %s\n', groups{g, 1});
                names = groups{g, 2};
                for k = 1:numel(names)
                    n = names{k};
                    v = obj.(n);
                    if isequal(v, d.(n)), mark = '*'; else, mark = ' '; end
                    fprintf('   %s %-12s %s\n', mark, n, lsqvce.Options.fmt(v));
                end
            end
            fprintf('\n');
        end
    end

    methods (Access = private)
        function nv = diagnoseNV(obj)
            D = obj.Diagnostics;
            nv = {'Alpha', D.Alpha, 'Bonferroni', logical(D.Bonferroni), 'Beta', D.Beta, ...
                  'RMin', D.RMin, 'RWeak', D.RWeak, 'Test', string(D.Test), ...
                  'RobustCheck', logical(D.RobustCheck), 'Store', string(D.Store), ...
                  'NWorst', D.NWorst};
        end
    end

    methods (Static)
        function s = noiseDefaults()
            s = struct('Grid', 60, 'Whitener', "toeplitz", 'MeanType', "robust", ...
                       'Smooth', 0, 'FreezeRedundancy', 1, 'SkipThreshold', 0.2);
        end
        function s = diagnosticsDefaults()
            s = struct('Alpha', 1e-3, 'Bonferroni', true, 'Beta', 0.2, 'RMin', 0.01, ...
                       'RWeak', 0.1, 'Test', "auto", 'RobustCheck', true, ...
                       'Store', "all", 'NWorst', 1000);
        end
        function g = mergeGroup(d, s, name)
        %MERGEGROUP  Partial struct s onto defaults d; unknown fields refused.
            if ~isstruct(s) || ~isscalar(s)
                error('lsqvce:Options:badGroupField', '%s must be a scalar struct.', name);
            end
            g = d;
            f = fieldnames(s);
            for k = 1:numel(f)
                if ~isfield(d, f{k})
                    error('lsqvce:Options:badGroupField', ...
                        '%s has no field ''%s'' (fields: %s).', name, f{k}, strjoin(fieldnames(d).', ', '));
                end
                g.(f{k}) = s.(f{k});
            end
        end
        function d = groupDefaults()
            d = struct('Noise', lsqvce.Options.noiseDefaults(), ...
                       'RobustOpts', struct('Tuning', [], 'MaxIter', 30, 'Tol', 1e-8), ...
                       'Diagnostics', lsqvce.Options.diagnosticsDefaults(), ...
                       'Parallel', struct('Mode', "serial", 'DistributeDays', false, 'Threshold', 8), ...
                       'GaussNewton', struct('TolF', 1e-10, 'Lambda0', 1e-6));
        end
    end

    methods (Static, Access = private)
        function t = groupTable()
            t = {'What to estimate',   {'Estimate', 'Pool', 'ArcFactors', 'Robust'}; ...
                 'The loop',           {'MaxPasses', 'MinPasses', 'Tolerance', 'Probes', 'Seed'}; ...
                 'Run control',        {'Backend', 'CacheDays', 'Verbose', 'LogFile', 'OnPass'}; ...
                 'Groups',             {'Noise', 'RobustOpts', 'Diagnostics', 'Parallel', 'GaussNewton'}};
        end
        function s = fmt(v)
            if isa(v, 'function_handle')
                s = func2str(v);
                if s(1) ~= '@', s = ['@' s]; end
            elseif isstring(v) || ischar(v)
                s = ['"' char(v) '"'];
            elseif isstruct(v)
                f = fieldnames(v);
                parts = cell(1, numel(f));
                for k = 1:numel(f)
                    parts{k} = sprintf('%s=%s', f{k}, lsqvce.Options.fmt(v.(f{k})));
                end
                s = ['struct(' strjoin(parts, ', ') ')'];
            elseif islogical(v)
                if v, s = 'true'; else, s = 'false'; end
            elseif isempty(v)
                s = '[]';
            elseif isnumeric(v) && isscalar(v)
                s = num2str(v);
            else
                s = mat2str(v);
            end
        end
    end
end
