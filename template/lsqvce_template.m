function L = lsqvce_template(procOpt)
%LSQVCE_TEMPLATE  Starting point for your own least-squares adjustment case.
%
%   L = LSQVCE_TEMPLATE(procOpt)
%
%   Copy this file, rename it, and fill in the two stubs at the bottom. It
%   returns the cell array of Loaders that the toolbox consumes:
%
%       L = myCase(procOpt);
%       R = lsqvce.engine.adjustVCE(L, procOpt.days, ...
%               'verbose', procOpt.loggingLevel, ...
%               'logFid',  procOpt.loggingOverviewFileID) ...
%             .groupComponentsBy('obsType') ...
%             .run();
%
%   THIS FILE HAS TWO REGISTERS -- READ THIS BEFORE RUNNING IT
%     * buildLoaders (below) is REAL, RUNNABLE code. It derives one Loader per
%       entry of procOpt.observationGroups. You edit it only if your topology
%       differs from one-system-per-group (see MERGING SYSTEMS below).
%     * myBuildDay and myEpochCov are STUBS. They ERROR until you write them.
%       That is deliberate: only you know how to read your data, and a plausible
%       fake would be worse than an honest stop.
%   So: running this as shipped gets you Loaders that error on first use. That
%   is the expected state until the two stubs are filled.
%
%   WHAT YOU MUST WRITE                                       (the stubs)
%     myBuildDay   one day of design matrix + observations    REQUIRED
%     myEpochCov   the fixed per-epoch covariance D           only for
%                                                             multi-axis streams
%
%   WHAT YOU MUST DECLARE                                     (in your procOpt)
%     .observationGroups   one entry per system: its streams and identity
%     .days                the day list
%     .estimation          your common/local parametrisation
%
%   WHAT THE TOOLBOX DERIVES -- DO NOT STATE IT TWICE
%     nc  from commonParams rows   |   nl  from localParams rows
%   The common-parameter count is stated ONCE, where the parametrisation is
%   defined. There is no nc argument on adjustVCE or adjustBlock.
%
%   INPUT
%     procOpt  1x1  struct  your processing configuration (see PROCOPT below).
%                           Piped UNCHANGED to myBuildDay as its 5th argument
%                           and surfaced to plot hooks as ctx.procOpt.
%
%   OUTPUT
%     L        1xS  cell    one Loader per system, ready for LSQVCE.ENGINE.ADJUSTVCE
%                           or LSQVCE.ENGINE.ADJUSTBLOCK.
%
%   PROCOPT -- IT IS IMMUTABLE INPUT, NOT SHARED STATE
%     procOpt is captured when makeLoader is called and piped through unchanged
%     for the whole run. So it must NOT carry live cursors: a field like
%     'iteridx' would freeze at its value on the first pass while the loop
%     iterates around it, and your builder would silently read the wrong pass
%     forever. The engine tracks these itself -- use ctx.pass and ctx.sysID.
%
%     Group your settings into nested structs; only the fields your builder
%     needs have to be there:
%
%       procOpt.days            1xD   day list             (fed to adjustVCE)
%       procOpt.loggingLevel    1x1   0 silent | 1 progress | 2 detail
%                                     NOTE 3 is REJECTED (merged into 2), so a
%                                     pipeline that can emit 3 must clamp first.
%       procOpt.loggingOverviewFileID 1x1  fid             (fed as 'logFid')
%       procOpt.showAndSaveGraphics   1x1  gates hook attachment (an IF in your
%                                     driver, not a toolbox option)
%       procOpt.observationGroups 1xS cell of structs      (drives the loop below)
%       procOpt.estimation      1x1   struct  your parametrisation tables
%       procOpt.data            1x1   struct  paths, maximumDegree, ...
%       procOpt.epochs          1x1   struct  mEpoch, dt, arcsPerDay, ...
%       procOpt.forces .constants .ephemeris .earthRotation .metadata
%                                     your own sub-structs, opaque to the toolbox
%
%     THE FIVE SHARED VALUES. days, loggingLevel, loggingOverviewFileID,
%     showAndSaveGraphics and observationGroups are needed by your pipeline
%     BEFORE any toolbox object exists, so they live here -- and the toolbox is
%     FED from them (see the call at the top of this help). Never type the same
%     value into both: that is exactly the drift this interface was built to
%     remove.
%
%   ERRORS
%     lsqvce_template:notImplemented   a stub was called before you wrote it
%     lsqvce_template:badProcOpt       procOpt is missing a required field
%
%   EXAMPLE
%     procOpt = myPipelineConfig();
%     L = lsqvce_template(procOpt);
%     R = lsqvce.engine.adjustVCE(L, procOpt.days).groupComponentsBy('obsType').run();
%     R.report();
%
%   See also LSQVCE.ENGINE.MAKELOADER, LSQVCE.ENGINE.ADJUSTVCE, LSQVCE.ENGINE.ADJUSTBLOCK,
%            DEMO_SEPARABLE_DATA
%
%   Created 2026-07-19 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.

    narginchk(1, 1);
    requireFields(procOpt, {'observationGroups', 'days'});
    L = buildLoaders(procOpt);
end

%% ========================================================== RUNNABLE: loaders
function L = buildLoaders(procOpt)
%BUILDLOADERS  One Loader per observation group. REAL CODE -- edit to taste.
%
%   Each entry of procOpt.observationGroups describes ONE system: which
%   observation types it carries and how it is identified. This loop turns each
%   into a Loader via LSQVCE.ENGINE.MAKELOADER.
%
%   MERGING SYSTEMS -- the one topology decision you cannot avoid
%     Two satellites must share ONE system whenever an observation's design
%     couples both of their LOCAL parameters (an inter-satellite range or
%     range-rate does exactly this). The block-arrow model has one local block
%     per system, so a coupling observation cannot span two of them. In that
%     case give both satellites' streams the SAME group entry and stack their
%     locals: [y_A ; y_B]. Streams that are merely correlated do NOT need
%     merging -- VCE handles that through the variance components.
%
%   INPUT
%     procOpt  1x1 struct  the processing configuration (see the file header)
%   OUTPUT
%     L        1xS cell    one Loader per group
    groups = procOpt.observationGroups;
    L = cell(1, numel(groups));

    for k = 1:numel(groups)
        g = groups{k};

        % --- what this group observes -------------------------------------
        % obsTypes are the STREAM names. They are the unit of variance-component
        % estimation when you call groupComponentsBy('obsType'), so name them
        % after what they physically are, not after the file they came from.
        obsTypes = cellstr(g.observationTypes);

        % --- a-priori covariance, one entry per obsType --------------------
        % Deliberately allowed to be WRONG: the loop's job is to discover the
        % truth. A plain white start is honest and lets you see it work.
        cov = aprioriCovariance(obsTypes, g);

        % --- the day-builder, bound over your configuration ----------------
        % NOTE the 5-argument form. The 5th argument is procOpt, and your
        % builder must read its settings FROM IT -- not from a variable captured
        % in this workspace. Two copies of one value is the bug this interface
        % exists to prevent.
        buildDay = @(day, obsType, common, local, po) ...
                       myBuildDay(day, obsType, common, local, po);

        % --- the per-epoch covariance hook (multi-axis streams only) -------
        % Supply it only when a stream has more than one axis per epoch (e.g.
        % 3-D positions). Leave it empty for scalar streams.
        if anyMultiAxis(g)
            epochCov = @(day, obsType) myEpochCov(day, obsType, procOpt);
        else
            epochCov = [];
        end

        % --- build the Loader ---------------------------------------------
        % systemId MUST be numeric or a SINGLE char: the accumulator indexes
        % systems by [store.sysID], which flattens a multi-char string into
        % separate elements. Carry a readable label in 'name' instead.
        L{k} = lsqvce.engine.makeLoader(buildDay, ...
                 'observationTypes',  obsTypes, ...
                 'covariance',        cov, ...
                 'systemId',          g.systemId, ...
                 'name',              groupName(g, k), ...
                 'arcs',              g.arcsPerDay, ...
                 'days',              procOpt.days, ...
                 'commonParams',      procOpt.estimation.commonParams, ...
                 'localParams',       g.localParams, ...
                 'processingOptions', procOpt, ...
                 'epochCovariance',   epochCov);
    end
end

%% ============================================================= STUBS: your data
function D = myBuildDay(day, obsType, common, local, procOpt)  %#ok<STOUT,INUSD>
%MYBUILDDAY  Read ONE day of ONE observation type. YOU WRITE THIS.
%
%   It must be DETERMINISTIC. The accumulator reads every day a SECOND time for
%   the residual pass and must get byte-identical data back; a random draw or a
%   changing file would corrupt the residuals silently, not loudly.
%
%   INPUT
%     day      1x1   double   the day identifier from procOpt.days
%     obsType  1xk   char     which stream to build
%     common   Rc x14 cell    the common parametrisation (nc = Rc rows)
%     local    Rl x14 cell    this system's local parametrisation (nl = Rl rows)
%     procOpt  1x1   struct   your processing configuration, piped unchanged
%
%   OUTPUT
%     D  1x1 struct with fields
%       .A       m x nc   common design matrix for the whole day
%       .B       m x nl   local design matrix (may be sparse)
%       .l       m x 1    observations
%       .t       m x 1    epoch timestamps [s]
%       .arcIdx  1 x nArc cell of row-index vectors, one per arc
%                         (optional; omitted -> the whole day is ONE arc)
%       .nAxis   1x1      axes per epoch, e.g. 3 for 3-D positions
%                         (optional; omit for scalar streams)
%
%   ROW ORDER FOR MULTI-AXIS STREAMS: rows are EPOCH-MAJOR, i.e.
%   [x1 y1 z1 | x2 y2 z2 | ...]. The temporal correlation model strides by
%   nAxis, so a stream stored axis-major comes out with its time law stretched
%   by the axis count -- and it will not announce itself.
%
%   GETTING nc AND nl: read them from the tables you were handed,
%   nc = size(common,1) and nl = size(local,1). Do not hardcode them; that is
%   the duplication the parametrisation pipe removes.
    error('lsqvce_template:notImplemented', ...
        ['myBuildDay is a stub -- write it for your data.\n' ...
         'Asked for day %s, obsType ''%s''. Return .A .B .l .t (and .arcIdx, ' ...
         '.nAxis if multi-axis); see `help lsqvce_template`.'], ...
         num2str(day), obsType);
end

function dd = myEpochCov(day, obsType, procOpt)  %#ok<STOUT,INUSD>
%MYEPOCHCOV  The FIXED per-epoch covariance D. YOU WRITE THIS (multi-axis only).
%
%   This is the spatial half of the separable model
%       C = sigma^2 * D (R_t kron I_d) D
%   where D is known per epoch (e.g. from a PPP solution) and R_t is the
%   temporal law the toolbox estimates. Supplying D lets the temporal law be
%   recovered by exact de-scaling instead of from whitened residuals.
%
%   INPUT
%     day      1x1  double  the day identifier
%     obsType  1xk  char    which stream
%     procOpt  1x1  struct  your processing configuration
%
%   OUTPUT
%     dd  1x1 struct with fields
%       .Cep          nEpochDay x d x d  per-epoch covariance
%       .arc(i).epochs  1xE  day-epoch indices of arc i
%       .arc(i).rows    1xR  row indices of arc i (R = E*d)
%
%   Scalar streams do not need this -- return [] from the hook instead.
    error('lsqvce_template:notImplemented', ...
        ['myEpochCov is a stub -- write it, or pass no ''epochCovariance'' ' ...
         'if every stream is scalar.\nAsked for day %s, obsType ''%s''.'], ...
         num2str(day), obsType);
end

%% ====================================================== OPTIONAL: plot hooks
% A hook is a function handle YOU write, @(ctx) -> figureHandle. The engine
% calls it, takes the figure back, and owns the show/export/close entirely
% (showSeconds / exportDir), so your hook only BUILDS a figure. A hook that
% throws is logged and skipped -- plots never break the solve.
%
% Attach via makeLoader above: a bare handle is broadcast to every obsType of
% that loader, or pass a containers.Map keyed by obsType for per-stream control.
%       'plotPostSolve', @myPostSolvePlot, ...
%       'plotPostModel', @myPostModelPlot
%
% WHEN THEY FIRE
%   postSolve  after EACH solve, before the residual pass  (ctx has the estimate)
%   postModel  ONCE, after the stochastic model is fitted  (ctx also has the
%              spectrum: .psd .psdFreq .acv .arOrder -- these do not exist
%              earlier, because the model has not been fitted yet)
%
% WHAT ctx CARRIES (schema 1; fields are only ever ADDED, never renamed)
%   .schema .stage .stream .satKey .obsType .sysID .nAxes .pass
%   .R        the lsqvce.Results  (.x, .sol, .acc, .residuals(...))
%   .x .sd    common estimates and their standard deviations ([] if unavailable)
%   .common .procOpt   YOUR parametrisation and configuration, coming back
%   .model .sigmaG     the fitted noise model and this group's VCE sigma
%   .exportDir .showSeconds
%
% function fh = myPostModelPlot(ctx)
%     fh = figure('Visible', 'off');
%     semilogy(ctx.psdFreq, ctx.psd);
%     title(sprintf('%s  (AR order %d)', ctx.stream, ctx.arOrder), ...
%           'Interpreter', 'none');
%     xlabel('frequency');  ylabel('PSD');  grid on;
% end

%% ========================================================== small local helpers
function cov = aprioriCovariance(obsTypes, g)
%APRIORICOVARIANCE  A-priori covariance spec per obsType, for LSQVCE.ENGINE.COVMODEL.
%   White and unit-variance by default: the VCE loop discovers the truth, so a
%   deliberately plain start makes the discovery visible. Replace with your best
%   prior if you have one.
%
%   INPUT
%     obsTypes  1xT cellstr  the stream names
%     g         1x1 struct   the observation-group entry (may carry .nAxis)
%   OUTPUT
%     cov       1x1 struct   .(obsType) = a covariance spec
    cov = struct();
    for k = 1:numel(obsTypes)
        ot = obsTypes{k};
        if axesOf(g, ot) > 1
            cov.(ot) = struct('type','separable', 'ar',1, 'sigma',1, ...
                              'd', axesOf(g, ot));
        else
            cov.(ot) = struct('type','diag', 'sigma',1);
        end
    end
end

function n = axesOf(g, obsType)
%AXESOF  Axes per epoch for one stream (3 for 3-D positions, 1 for scalars).
%   Reads g.nAxis.(obsType) when present, else assumes scalar. Adapt to your own
%   convention if your group struct describes axes differently.
%
%   INPUT   g 1x1 struct; obsType 1xk char
%   OUTPUT  n 1x1 double  axes per epoch (>= 1)
    n = 1;
    if isfield(g, 'nAxis') && isstruct(g.nAxis) && isfield(g.nAxis, obsType)
        n = g.nAxis.(obsType);
    end
end

function tf = anyMultiAxis(g)
%ANYMULTIAXIS  True when any stream of this group has more than one axis.
%   INPUT   g 1x1 struct   the observation-group entry
%   OUTPUT  tf 1x1 logical
    tf = false;
    ots = cellstr(g.observationTypes);
    for k = 1:numel(ots)
        if axesOf(g, ots{k}) > 1, tf = true; return; end
    end
end

function s = groupName(g, k)
%GROUPNAME  A human-readable label for a group (display only, never indexed).
%   INPUT   g 1x1 struct; k 1x1 double  the group index
%   OUTPUT  s 1xn char
    if isfield(g, 'name') && ~isempty(g.name)
        s = char(g.name);
    else
        s = sprintf('system-%d', k);
    end
end

function requireFields(s, names)
%REQUIREFIELDS  Error clearly when procOpt is missing something required.
%   Fails at the TOP with the field name, rather than deep inside makeLoader
%   with a message about an argument the user never typed.
%
%   INPUT   s 1x1 struct; names 1xN cellstr
%   OUTPUT  none (throws lsqvce_template:badProcOpt)
    missing = names(~isfield(s, names));
    if ~isempty(missing)
        error('lsqvce_template:badProcOpt', ...
            'procOpt is missing required field(s): %s.', strjoin(missing, ', '));
    end
end
