classdef Options
%LSQVCE.OPTIONS  The one settings object for LSQVCE.FIT (typed, with engine defaults).
%
%   opts = lsqvce.Options()                       every default
%   opts = lsqvce.Options(Name, Value, ...)       override some
%   opts.Name = Value                             change one later (validated)
%
%   WHY ONE OBJECT (R4 findings G-02/G-03). The estimation options used to
%   live in four conventions: per-driver Name-Value names, the fluent
%   configurator chain of LSQVCE.ENGINE.ADJUSTVCE, the loader UserOpt, and the joint
%   translation. Options is the single place every knob is NAMED, TYPED, has
%   its DEFAULT stated, and is checked -- at construction through an
%   arguments block (unknown names and bad values are rejected by the
%   signature), and on later assignment through property validation. The
%   defaults are EXACTLY the engines' own defaults (pinned by the reference
%   test in dev/tests/test_options.m), so a run through fit with a default
%   Options is the same run as a direct driver call: no behaviour change.
%
%   The MATLAB model is optimoptions/statset: one object, canonical
%   UpperCamelCase names, string enums with case-sensitive lowercase VALUES.
%   It is a CLASS named Options (not a function options) because the user's
%   Windows working copy cannot hold options.m and Options.m side by side
%   (case-insensitive file system; audit check 18b guards this tree-wide).
%
%   PROPERTIES  (group -- name  type  default  meaning)
%
%   What to do
%     Model      string   "gls"   "gls" one GLS solve (adjustBlock) |
%                                 "vce" alternate stochastic model <-> variance
%                                 components (adjustVCE 'separate') |
%                                 "pooled" pooled spectrum, no per-arc factors
%                                 (adjustVCE 'joint') | "ellmer" the ten-step
%                                 joint co-estimation with keyed arc factors
%                                 (adjustVCEjoint)
%     Robust     string   "none"  "none" | "irls" M-estimation (batch path,
%                                 LSQVCE.ENGINE.ROBUST)
%
%   Solver
%     Backend    string   "qr"    "qr" | "normal" (the engines' default is
%                                 "qr"; changing the DEFAULT is a decision,
%                                 not an option -- HANDOVER §11.2 item 2)
%     Reg        string   "none"  "none"|"tikhonov"|"tsvd"|"lm"|"auto"
%     RegLambda  double   0       regularisation weight, >= 0
%     RegMatrix  double   []      p x nc regularisation matrix
%     Sparse     logical  false   batch path: sparse normal matrix
%     DatumB     double   []      batch paths: datum constraint matrix
%     Sigma0     double   1       batch paths: a-priori sigma0 (> 0)
%
%   Run control
%     Parallel   string   "serial"  "serial" | "parfor" | "spmd"
%     DistributeDays logical false  parallelise over the DAYS of a system
%                                 (one raw partial per day, merged in day
%                                 order) -- the axis that scales for a
%                                 2-satellite mission; needs Parallel ~=
%                                 "serial" and >= ParallelThreshold days
%     ParallelThreshold double 8  below this many tasks a pool is not used
%     Verbose    logical|double  true   logical, or a level 0..3. PINNED
%                                 DEVIATION: adjust and adjustVCE default to
%                                 true, adjustBlock to false; Options follows
%                                 the majority, so fit's gls path talks by
%                                 default -- set Verbose=false for silence
%     LogFid     double   []      log file id (see LSQVCE.UTIL.LOG); [] = leave
%     LogFlag    logical  []      log to screen; [] = leave
%     ShowSeconds double  0       plot hooks: seconds to show each figure
%     ExportDir  string   ""      plot hooks: PNG folder ("" = none)
%
%   Variance components (between groups, Helmert)
%     GroupBy    char|string|fh 'satObs'  'satObs'|'obsType'|'system'|'arc'
%                                 or @(sysID,obsType,arc,day) -> key
%     VceMaxIter double   10      inner VCE iterations
%     VceTol     double   5e-3    convergence on max|sigma-1|
%     MinRedundancy double 50     below this a group cannot be estimated
%     ComponentsPerEpoch struct struct('pos',3)  obsType -> components/epoch
%
%   Stochastic model (temporal shape, within a group)
%     MaxOrder   double   60      highest AR order considered
%     AcvTol     double   1e-3    autocovariance fit tolerance
%     MaxLag     double   []      lags for the ACV estimate ([] = auto)
%     WhiteAlpha double   0.01    Ljung-Box significance (whiteness gate)
%
%   Outer loop
%     MaxOuter   double   6       outer passes
%     OuterTol   double   5e-3    convergence on max|sigma-1| and AR drift
%
%   Pooled / Ellmer joint spectrum
%     NumLags    double|"auto" 60 frequencies in the noise model. "auto"
%                                 (pooled only) sizes it from the data. For
%                                 "ellmer" the model reaches the loaders as an
%                                 'acf' that is zero beyond its last lag, so
%                                 NumLags below the longest arc length gives a
%                                 truncated (in general indefinite) covariance:
%                                 every arc is skipped as nonPD and the driver
%                                 reports noArcs. Use NumLags >= longest arc.
%     Probes     double   100     Hutchinson probes
%     SkipThreshold double 0.2    R5 M1: fraction of a stream's arcs that may be
%                                 SKIPPED in one pass (non-PD covariance, a
%                                 whitener that threw, ...) before the run
%                                 REFUSES instead of updating the stream on the
%                                 rest (lsqvce:adjustVCE:skipThreshold)
%     Seed       double   []      RNG seed for the probes ([] = unseeded)
%     Norm       string   "ortho" DCT normalisation "ortho" | "ellmer"
%     MeanType   string   "robust" arc-factor mean "robust" | "plain"
%     PsdMode    string   "toeplitz" "toeplitz" | "ar" | "arma"   (pooled)
%     Smooth     double   0       PSD moving-average half-width (pooled)
%     C0         double   []      starting covariance function (pooled)
%     JointMaxIter double 10      inner joint iterations (ellmer)
%     JointTol   double   1e-2    joint convergence (ellmer)
%     StepSchedule string "none"  "none" | "decay" (ellmer)
%     AvgTail    double   0       passes averaged at the end (ellmer)
%
%   Robust IRLS (batch path)
%     Influence  string   "huber" "huber" | "igg3" | "bisquare"
%     Tuning     double   []      tuning constant ([] = the influence default)
%     RobustMaxIter double 30     IRLS iterations
%     RobustTol  double   1e-8    IRLS convergence
%
%   Nonlinear (batch path, LSQVCE.ENGINE.SOLVENONLINEAR)
%     NonlinearMaxIter double 30  Gauss-Newton iterations
%     TolX       double   1e-3    step tolerance
%     TolF       double   1e-10   cost tolerance
%     Lambda0    double   1e-6    initial damping
%
%   Diagnostics (forwarded to LSQVCE.ENGINE.DIAGNOSE; replaces KeepUnmatched)
%     Alpha      double   1e-3    test significance
%     Bonferroni logical  true    family-wise correction
%     Beta       double   0.20    power for the MDB
%     RMin       double   1e-6    redundancy floor
%     RWeak      double   0.10    weak-redundancy threshold
%     Test       string   "auto"  "auto" | "pope" | "baarda" | "mad"
%     RobustCheck logical true    run the robust cross-check
%     Store      string   "all"   "all" | "flagged" | "none"
%     NWorst     double   1000    worst observations kept
%     Verify     logical  true    second-pass integrity verification
%     Progress   logical  false   diagnostics progress lines
%
%   METHODS
%     args = opts.engineArgs(target)   the mapping onto ONE existing driver:
%         target  string  "adjustBlock" | "adjustVCE" | "adjustVCEjoint" |
%                         "adjust" | "robust" | "solveNonlinear" | "diagnose"
%         args    1x1 struct  .nv  1x2k cell  Name-Value pairs for that call
%                 for "adjustVCE" additionally the configurator groups, in the
%                 engine's own field order and value classes:
%                 .shape .scale .outer .joint (structs), .axes (struct
%                 obsType -> count), .mode (char 'separate' | 'joint')
%     opts.validate()                  cross-property checks (called by the
%                                      constructor and by engineArgs)
%     s = opts.toStruct()              every public property as a struct field
%                                      (what fit records in Results.meta.options)
%
%   ACCURACY
%     No numerics of its own. The claim is EXACTNESS of the mapping: for a
%     default Options, engineArgs("adjustVCE") returns structs isequal to a
%     fresh adjustVCE's shapeOpt/scaleOpt/outerOpt/jointOpt/blockOpt/plotOpt
%     (dev/tests/test_options.m, reference), and every other target's .nv
%     carries the literal defaults read from the driver signatures on
%     2026-08-31 (pinned in the same suite; a driver default change breaks
%     that test on purpose).
%
%   EXAMPLE
%     opts = lsqvce.Options(Model="vce", GroupBy="satObs", Parallel="parfor");
%     opts.MaxOuter = 4;
%     a = opts.engineArgs("adjustVCE");
%     P = lsqvce.engine.adjustVCE(L, days, a.nv{:}, 'Shape', a.shape, ...
%             'Scale', a.scale, 'Outer', a.outer, 'Joint', a.joint, 'Axes', a.axes);
%     R = P.run();                              % (fit does exactly this, WP-56/64b)
%
%   ERRORS
%     Your value is wrong (signature / property validation, MATLAB ids):
%       MATLAB:validators:mustBeMember     enum VALUE not in the list (values
%                                          are case-sensitive lowercase)
%       MATLAB:validators:mustBePositive   and the other typed constraints
%       MATLAB:InputParser:...             / MATLAB:TooManyInputs  unknown
%                                          option name
%     Your combination is wrong (body, lsqvce ids):
%       lsqvce:Options:badGroupBy   GroupBy is not a char/string or handle
%       lsqvce:Options:badNumLags   NumLags is neither "auto" nor an integer >= 2
%       lsqvce:Options:badVerbose   Verbose is neither logical nor a level 0..3
%       lsqvce:Options:badTarget    engineArgs target unknown
%
%   See also LSQVCE.ENGINE.ADJUSTBLOCK, LSQVCE.ENGINE.ADJUSTVCE, LSQVCE.ENGINE.ADJUSTVCEJOINT,
%            LSQVCE.ENGINE.ADJUST, LSQVCE.ENGINE.ROBUST, LSQVCE.ENGINE.SOLVENONLINEAR,
%            LSQVCE.ENGINE.DIAGNOSE, LSQVCE.RESULTS
%
%   Updated 2026-08-31 (CEST, UTC+02:00), R4 WP-62 (Level L4): DistributeDays
%   and ParallelThreshold -- the parallel policy of the streaming engines,
%   reachable through fit. With assistance from Claude (Anthropic), model
%   Claude Fable 5.
%   Updated 2026-08-31 (CEST, UTC+02:00), R4 WP-56b (Level L4): toStruct()
%   for the provenance record of fit's in-memory / batched kinds. With
%   assistance from Claude (Anthropic), model Claude Fable 5.
%   Created 2026-08-31 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5 -- R4 Level L1, WP-55.

    properties
        % ---- what to do
        Model  (1,1) string {mustBeMember(Model,  ["gls","vce","pooled","ellmer"])} = "gls"
        Robust (1,1) string {mustBeMember(Robust, ["none","irls"])} = "none"
        % ---- solver
        Backend   (1,1) string {mustBeMember(Backend, ["qr","normal"])} = "qr"
        Reg       (1,1) string {mustBeMember(Reg, ["none","tikhonov","tsvd","lm","auto"])} = "none"
        RegLambda (1,1) double {mustBeNonnegative} = 0
        RegMatrix double = []
        Sparse    (1,1) logical = false
        DatumB    double = []
        Sigma0    (1,1) double {mustBePositive} = 1
        % ---- run control
        Parallel    (1,1) string {mustBeMember(Parallel, ["serial","parfor","spmd"])} = "serial"
        DistributeDays    (1,1) logical = false
        ParallelThreshold (1,1) double {mustBeInteger, mustBeNonnegative} = 8
        Verbose     (1,1) {mustBeA(Verbose, ["logical","double"])} = true
        LogFid      double  {mustBeScalarOrEmpty} = []
        LogFlag     logical = logical([])
        ShowSeconds (1,1) double = 0
        ExportDir   (1,1) string = ""
        % ---- variance components
        GroupBy = 'satObs'
        VceMaxIter    (1,1) double {mustBeInteger, mustBePositive} = 10
        VceTol        (1,1) double {mustBePositive} = 5e-3
        MinRedundancy (1,1) double {mustBeNonnegative} = 50
        ComponentsPerEpoch (1,1) struct = struct('pos', 3)
        % ---- stochastic model
        MaxOrder   (1,1) double {mustBeInteger, mustBeNonnegative} = 60
        AcvTol     (1,1) double {mustBePositive} = 1e-3
        MaxLag     double {mustBeScalarOrEmpty} = []
        WhiteAlpha (1,1) double {mustBePositive, mustBeLessThan(WhiteAlpha, 1)} = 0.01
        % ---- outer loop
        MaxOuter (1,1) double {mustBeInteger, mustBePositive} = 6
        OuterTol (1,1) double {mustBePositive} = 5e-3
        % ---- pooled / ellmer joint spectrum
        NumLags = 60
        Probes   (1,1) double {mustBeInteger, mustBePositive} = 100
        SkipThreshold (1,1) double {mustBeInRange(SkipThreshold, 0, 1)} = 0.2
        Seed     double {mustBeScalarOrEmpty} = []
        Norm     (1,1) string {mustBeMember(Norm, ["ortho","ellmer"])} = "ortho"
        MeanType (1,1) string {mustBeMember(MeanType, ["robust","plain"])} = "robust"
        PsdMode  (1,1) string {mustBeMember(PsdMode, ["toeplitz","ar","arma"])} = "toeplitz"
        Smooth   (1,1) double {mustBeInteger, mustBeNonnegative} = 0
        C0       double = []
        JointMaxIter (1,1) double {mustBeInteger, mustBePositive} = 10
        JointTol     (1,1) double {mustBePositive} = 1e-2
        StepSchedule (1,1) string {mustBeMember(StepSchedule, ["none","decay"])} = "none"
        AvgTail      (1,1) double {mustBeInteger, mustBeNonnegative} = 0
        % ---- robust IRLS
        Influence     (1,1) string {mustBeMember(Influence, ["huber","igg3","bisquare"])} = "huber"
        Tuning        double = []
        RobustMaxIter (1,1) double {mustBeInteger, mustBePositive} = 30
        RobustTol     (1,1) double {mustBePositive} = 1e-8
        % ---- nonlinear
        NonlinearMaxIter (1,1) double {mustBeInteger, mustBePositive} = 30
        TolX    (1,1) double {mustBePositive} = 1e-3
        TolF    (1,1) double {mustBePositive} = 1e-10
        Lambda0 (1,1) double {mustBeNonnegative} = 1e-6
        % ---- diagnostics
        Alpha       (1,1) double {mustBePositive, mustBeLessThan(Alpha, 1)} = 1e-3
        Bonferroni  (1,1) logical = true
        Beta        (1,1) double {mustBePositive, mustBeLessThan(Beta, 1)} = 0.20
        RMin        (1,1) double {mustBeNonnegative} = 1e-6
        RWeak       (1,1) double {mustBeNonnegative, mustBeLessThanOrEqual(RWeak, 1)} = 0.10
        Test        (1,1) string {mustBeMember(Test, ["auto","pope","baarda","mad"])} = "auto"
        RobustCheck (1,1) logical = true
        Store       (1,1) string {mustBeMember(Store, ["all","flagged","none"])} = "all"
        NWorst      (1,1) double {mustBeNonnegative} = 1000
        Verify      (1,1) logical = true
        Progress    (1,1) logical = false
    end

    methods
        function obj = Options(opt)
        %OPTIONS  Construct from Name-Value pairs; every name is a property.
        arguments
            opt.?lsqvce.Options
        end
            f = fieldnames(opt);
            for k = 1:numel(f)
                obj.(f{k}) = opt.(f{k});
            end
            obj.validate();
        end

        function validate(obj)
        %VALIDATE  Cross-property checks the property validators cannot express.
            g = obj.GroupBy;
            if ~(ischar(g) || (isstring(g) && isscalar(g)) || isa(g, 'function_handle'))
                error('lsqvce:Options:badGroupBy', ...
                    ['GroupBy must be a char/string key (''satObs'', ''obsType'', ' ...
                     '''system'', ''arc'') or a function handle @(sysID,obsType,arc,day).']);
            end
            n = obj.NumLags;
            okAuto = (ischar(n) || (isstring(n) && isscalar(n))) && strcmp(char(n), 'auto');
            okNum  = isnumeric(n) && isscalar(n) && isfinite(n) && n >= 2 && n == round(n);
            if ~(okAuto || okNum)
                error('lsqvce:Options:badNumLags', ...
                    'NumLags must be "auto" or an integer >= 2 (got %s).', mat2str(n));
            end
            v = obj.Verbose;
            if isnumeric(v) && ~(isfinite(v) && v >= 0 && v <= 3 && v == round(v))
                error('lsqvce:Options:badVerbose', ...
                    'Verbose must be logical or an integer level 0..3 (got %g).', v);
            end
        end

        function args = engineArgs(obj, target)
        %ENGINEARGS  Map this object onto ONE existing driver's inputs.
        %   args = opts.engineArgs(target), target as listed in the class help.
        %   The .nv cell is meant for a comma-list call: driver(..., args.nv{:}).
        %   Value classes match what the drivers store (char for the adjustVCE
        %   configurator groups), so a default Options reproduces the drivers'
        %   own defaults byte for byte (reference test).
            obj.validate();
            target = char(string(target));
            switch target
                case 'adjustBlock'
                    % adjustBlock takes a LOGICAL Verbose (no levels); a level
                    % >= 1 means "talk", so logical() is the faithful map.
                    args.nv = {'Backend', obj.Backend, 'Reg', obj.Reg, ...
                               'RegLambda', obj.RegLambda, 'RegMatrix', obj.RegMatrix, ...
                               'Parallel', obj.Parallel, 'DistributeDays', obj.DistributeDays, ...
                               'ParallelThreshold', obj.ParallelThreshold, ...
                               'Verbose', logical(obj.Verbose), ...
                               'ShowSeconds', obj.ShowSeconds, 'ExportDir', obj.ExportDir};
                case 'adjustNonlinear'
                    % R5 M2 2c': Gauss-Newton / Levenberg on the block path
                    args.nv = {'MaxIter', obj.NonlinearMaxIter, 'TolX', obj.TolX, ...
                               'TolF', obj.TolF, 'Lambda0', obj.Lambda0, ...
                               'Backend', obj.Backend, 'Parallel', obj.Parallel, ...
                               'DistributeDays', obj.DistributeDays, ...
                               'ParallelThreshold', obj.ParallelThreshold, ...
                               'Verbose', logical(obj.Verbose)};
                case 'adjustRobust'
                    % R5 M2 2b: IRLS on the block path (every kind)
                    args.nv = {'Backend', obj.Backend, 'Influence', obj.Influence, ...
                               'Tuning', obj.Tuning, 'MaxIter', obj.RobustMaxIter, ...
                               'Tol', obj.RobustTol, ...
                               'Parallel', obj.Parallel, 'DistributeDays', obj.DistributeDays, ...
                               'ParallelThreshold', obj.ParallelThreshold, ...
                               'Verbose', logical(obj.Verbose)};
                case 'adjustVCE'
                    if strcmp(obj.Model, "pooled")
                        mode = 'joint';
                    else
                        mode = 'separate';
                    end
                    args.nv = {'GroupBy', obj.GroupBy, 'VceMode', string(mode), ...
                               'MaxOuter', obj.MaxOuter, 'Backend', obj.Backend, ...
                               'Reg', obj.Reg, 'Parallel', obj.Parallel, ...
                               'DistributeDays', obj.DistributeDays, ...
                               'ParallelThreshold', obj.ParallelThreshold, ...
                               'Verbose', obj.Verbose, 'LogFid', obj.LogFid, ...
                               'LogFlag', obj.LogFlag, 'ShowSeconds', obj.ShowSeconds, ...
                               'ExportDir', obj.ExportDir};
                    args.mode  = mode;
                    % Field ORDER and value CLASSES follow adjustVCE's constructor
                    % (adjustVCE.m, "obj.shapeOpt = struct(...)" and neighbours).
                    args.shape = struct('MaxOrder', obj.MaxOrder, 'AcvTol', obj.AcvTol, ...
                                        'MaxLag', obj.MaxLag, 'WhiteAlpha', obj.WhiteAlpha);
                    args.scale = struct('MaxIter', obj.VceMaxIter, 'Tol', obj.VceTol, ...
                                        'MinRedundancy', obj.MinRedundancy);
                    args.outer = struct('MaxOuter', obj.MaxOuter, 'Tol', obj.OuterTol);
                    nl = obj.NumLags;
                    if ~isnumeric(nl), nl = char(nl); end
                    args.joint = struct('NumLags', nl, 'Probes', obj.Probes, 'SkipThreshold', obj.SkipThreshold, ...
                                        'Seed', obj.Seed, 'Norm', char(obj.Norm), ...
                                        'MeanType', char(obj.MeanType), 'C0', obj.C0, ...
                                        'PsdMode', char(obj.PsdMode), 'Smooth', obj.Smooth);
                    args.axes  = obj.ComponentsPerEpoch;
                case 'adjustVCEjoint'
                    nl = obj.NumLags;
                    if ~isnumeric(nl)
                        error('lsqvce:Options:badNumLags', ...
                            'adjustVCEjoint needs a numeric NumLags ("auto" is pooled-only).');
                    end
                    args.nv = {'NumLags', nl, 'MaxIter', obj.JointMaxIter, ...
                               'Tol', obj.JointTol, 'Probes', obj.Probes, 'Seed', obj.Seed, 'SkipThreshold', obj.SkipThreshold, ...
                               'Norm', obj.Norm, 'MeanType', obj.MeanType, ...
                               'StepSchedule', obj.StepSchedule, 'AvgTail', obj.AvgTail, ...
                               'Verbose', logical(obj.Verbose)};
                case 'adjust'
                    args.nv = [{'Method', char(obj.Backend), 'Sparse', obj.Sparse, ...
                                'DatumB', obj.DatumB, 'Sigma0', obj.Sigma0, ...
                                'Verbose', obj.Verbose}, obj.diagnoseNV()];
                case 'robust'
                    args.nv = [{'Influence', char(obj.Influence), 'Tuning', obj.Tuning, ...
                                'MaxIter', obj.RobustMaxIter, 'Tol', obj.RobustTol, ...
                                'DatumB', obj.DatumB, 'Method', char(obj.Backend), ...
                                'Sigma0', obj.Sigma0, 'Verbose', obj.Verbose}, obj.diagnoseNV()];
                case 'solveNonlinear'
                    args.nv = [{'MaxIter', obj.NonlinearMaxIter, 'TolX', obj.TolX, ...
                                'TolF', obj.TolF, 'Lambda0', obj.Lambda0, ...
                                'DatumB', obj.DatumB, 'Method', char(obj.Backend), ...
                                'Sigma0', obj.Sigma0, 'Verbose', obj.Verbose}, obj.diagnoseNV()];
                case 'diagnose'
                    args.nv = obj.diagnoseNV();
                otherwise
                    error('lsqvce:Options:badTarget', ...
                        ['engineArgs target must be one of adjustBlock, adjustVCE, ' ...
                         'adjustVCEjoint, adjust, robust, solveNonlinear, diagnose (got ''%s'').'], ...
                        target);
            end
        end

        function s = toStruct(obj)
        %TOSTRUCT  Every public property as a plain struct field (for Results.meta).
        %   s = opts.toStruct()   the RESOLVED option set a run actually used,
        %   in the shape runMeta expects. Read through the metaclass so a new
        %   property is picked up without a second list.
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
        %DISP  Grouped display: one line per property, engine defaults marked.
            groups = lsqvce.Options.groupTable();
            d = lsqvce.Options();                       % defaults for the marker
            fprintf('  lsqvce.Options  (* = engine default)\n');
            for g = 1:size(groups, 1)
                fprintf('\n  %s\n', groups{g, 1});
                names = groups{g, 2};
                for k = 1:numel(names)
                    n = names{k};
                    v = obj.(n);
                    if isequal(v, d.(n)), mark = '*'; else, mark = ' '; end
                    fprintf('   %s %-18s %s\n', mark, n, lsqvce.Options.fmt(v));
                end
            end
            fprintf('\n');
        end
    end

    methods (Access = private)
        function nv = diagnoseNV(obj)
        %DIAGNOSENV  The LSQVCE.ENGINE.DIAGNOSE pass-through pairs (replaces KeepUnmatched).
            nv = {'Alpha', obj.Alpha, 'Bonferroni', obj.Bonferroni, 'Beta', obj.Beta, ...
                  'RMin', obj.RMin, 'RWeak', obj.RWeak, 'Test', obj.Test, ...
                  'RobustCheck', obj.RobustCheck, 'Store', obj.Store, ...
                  'NWorst', obj.NWorst, 'Verify', obj.Verify, 'Progress', obj.Progress};
        end
    end

    methods (Static, Access = private)
        function t = groupTable()
        %GROUPTABLE  Display groups (name, property list), in help-text order.
            t = {'What to do',            {'Model', 'Robust'}; ...
                 'Solver',                {'Backend', 'Reg', 'RegLambda', 'RegMatrix', 'Sparse', 'DatumB', 'Sigma0'}; ...
                 'Run control',           {'Parallel', 'DistributeDays', 'ParallelThreshold', 'Verbose', 'LogFid', 'LogFlag', 'ShowSeconds', 'ExportDir'}; ...
                 'Variance components',   {'GroupBy', 'VceMaxIter', 'VceTol', 'MinRedundancy', 'ComponentsPerEpoch'}; ...
                 'Stochastic model',      {'MaxOrder', 'AcvTol', 'MaxLag', 'WhiteAlpha'}; ...
                 'Outer loop',            {'MaxOuter', 'OuterTol'}; ...
                 'Joint spectrum',        {'NumLags', 'Probes', 'SkipThreshold', 'Seed', 'Norm', 'MeanType', 'PsdMode', 'Smooth', 'C0', 'JointMaxIter', 'JointTol', 'StepSchedule', 'AvgTail'}; ...
                 'Robust IRLS',           {'Influence', 'Tuning', 'RobustMaxIter', 'RobustTol'}; ...
                 'Nonlinear',             {'NonlinearMaxIter', 'TolX', 'TolF', 'Lambda0'}; ...
                 'Diagnostics',           {'Alpha', 'Bonferroni', 'Beta', 'RMin', 'RWeak', 'Test', 'RobustCheck', 'Store', 'NWorst', 'Verify', 'Progress'}};
        end

        function s = fmt(v)
        %FMT  One-line rendering of a property value for disp.
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
            elseif isempty(v)
                s = '[]';
            elseif islogical(v) && isscalar(v)
                if v, s = 'true'; else, s = 'false'; end
            elseif isnumeric(v) && isscalar(v)
                s = num2str(v, '%.6g');
            else
                s = sprintf('[%s %s]', mat2str(size(v)), class(v));
            end
        end
    end
end
