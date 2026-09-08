function L = lsqvce_template(procOpt)
%LSQVCE_TEMPLATE  Starting point for your own least-squares adjustment case.
%
%   L = LSQVCE_TEMPLATE(procOpt)
%
%   Copy this file, rename it, and fill in the one stub at the bottom. It
%   returns the cell array of Loaders that the toolbox consumes:
%
%       L = myCase(procOpt);
%       lsqvce.util.log.set('fid', procOpt.loggingOverviewFileID, ...
%               'screen', procOpt.loggingLevel > 0);    % one sink for the text
%       R = lsqvce.engine.adjustBlock(L, procOpt.days, ...
%               'Verbose', procOpt.loggingLevel > 0);   % one GLS solve on the Loaders
%
%   ENGINE-TIER TEMPLATE (since 2.0). The user path needs no Loader at all:
%   describe the same case as a plain `data` struct (the Manual, chapter 3) and
%   call lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false)) -- examples/demo_fit.m
%   shows it end to end. Keep this template when you want the Loaders in
%   hand (procOpt-driven cases, direct engine work, 1.x code).
%
%   THIS FILE HAS TWO REGISTERS -- READ THIS BEFORE RUNNING IT
%     * buildLoaders (below) is REAL, RUNNABLE code. It derives one Loader per
%       entry of procOpt.observationGroups. You edit it only if your topology
%       differs from one-system-per-group (see MERGING SYSTEMS below).
%     * myBuildDay is a STUB. It ERRORS until you write it. That is
%       deliberate: only you know how to read your data, and a plausible fake
%       would be worse than an honest stop.
%   So: running this as shipped gets you Loaders that error on first use. That
%   is the expected state until the stub is filled.
%
%   WHAT YOU MUST WRITE                                       (the stub)
%     myBuildDay   one day of design matrix + observations    REQUIRED
%                  + their error information (D.Cep: a sigma per epoch, or a
%                  d x d covariance per epoch for multi-axis streams -- R6
%                  WP-90; the toolbox composes it with the temporal law)
%
%   WHAT YOU MUST DECLARE                                     (in your procOpt)
%     .observationGroups   one entry per system: its streams and identity
%     .days                the day list
%     .estimation          your common/local parametrisation
%
%   WHAT THE TOOLBOX DERIVES -- DO NOT STATE IT TWICE
%     nc  from commonParams rows   |   nl  from localParams rows
%   The common-parameter count is stated ONCE, where the parametrisation is
%   defined. There is no nc argument on adjustBlock.
%
%   INPUT
%     procOpt  1x1  struct  your processing configuration (see PROCOPT below).
%                           Piped UNCHANGED to myBuildDay as its 5th argument
%                           and kept on every Loader as L{k}.procOpt. Only
%                           LSQVCE.FIT surfaces it to a hook as
%                           ctx.processingOptions; at engine tier hand it over
%                           yourself -- 'OnPassContext',
%                           struct('processingOptions', {{procOpt}}) -- and
%                           ADJUSTBLOCK has no OnPass at all.
%
%   OUTPUT
%     L        1xS  cell    one Loader per system, ready for LSQVCE.ENGINE.ADJUSTBLOCK
%                           or LSQVCE.ENGINE.ADJUSTROBUST.
%
%   PROCOPT -- IT IS IMMUTABLE INPUT, NOT SHARED STATE
%     procOpt is captured when makeLoader is called and piped through unchanged
%     for the whole run. So it must NOT carry live cursors: a field like
%     'iteridx' would freeze at its value on the first pass while the loop
%     iterates around it, and your builder would silently read the wrong pass
%     forever. The engine tracks these itself -- an OnPass hook is handed the
%     live pass number as ctx.pass and the driver's name as ctx.driver, and the
%     system identity rides on the Loader as L{k}.sysID.
%
%     Group your settings into nested structs; only the fields your builder
%     needs have to be there:
%
%       procOpt.days            1xD   day list             (fed to adjustBlock)
%       procOpt.loggingLevel    1x1   0 silent | 1 progress | 2 detail
%                                     NOTE 3 is REJECTED (merged into 2), so a
%                                     pipeline that can emit 3 must clamp first.
%       procOpt.loggingOverviewFileID 1x1  fid
%                                     NOT an engine argument: route it to the
%                                     toolbox-wide logger before the solve,
%                                     lsqvce.util.log.set('fid', ...) (user
%                                     tier: lsqvce.util.log's file sink).
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
%     R = lsqvce.engine.adjustBlock(L, procOpt.days);
%     R.report();
%
%   See also LSQVCE.ENGINE.MAKELOADER, LSQVCE.ENGINE.ADJUSTBLOCK, LSQVCE.FIT,
%            DEMO_GRACE_DATALOADER
%
%   Updated 2026-09-07 (CEST, UTC+02:00), R6 WP-95: one stub, not two; the
%   overview fid goes to lsqvce.util.log.set (there is no 'logFid'); procOpt
%   reaches a hook only through 'OnPassContext' (ctx.processingOptions is a
%   LSQVCE.FIT service) and there is no ctx.sysID; the parametrisation tables
%   are Rc x C / Rl x C, the column count your own convention.
%   With assistance from Claude (Anthropic), model Claude Opus 5.
%   Updated 2026-09-06 (CEST, UTC+02:00), R6 WP-90: the error information
%   comes back from myBuildDay with the day (D.Cep) -- the myEpochCov stub,
%   the epochCovariance hook and the retired plot-hook / adjustVCE text are
%   gone. With assistance from Claude (Anthropic), model Claude Fable 5.1.
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
        % estimation, so name them after what they physically are, not after the
        % file they came from.
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
                 'processingOptions', procOpt);
    end
end

%% ============================================================== STUB: your data
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
%     common   Rc x C cell    the common parametrisation, C columns of YOUR
%                             convention (nc = Rc rows)
%     local    Rl x C cell    this system's local parametrisation, C columns of
%                             YOUR convention (nl = Rl rows)
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
%       .Cep     nEpoch x 1 | nEpoch x d | nEpoch x d x d   THE ERROR INFORMATION
%                         of the day (R6 WP-90): the per-epoch standard deviation
%                         of a scalar stream, the per-axis standard deviations, or
%                         the full per-epoch covariance (e.g. from a PPP solution)
%                         of a multi-axis stream. The toolbox composes it with the
%                         temporal law: C = D (R_t kron I_d) D. The axis count is
%                         derived from it. (optional; omit when every observation
%                         is equally accurate)
%       .nAxis   1x1      axes per epoch of a multi-axis stream WITHOUT Cep
%                         (optional; derived from Cep when present)
%
%   ROW ORDER FOR MULTI-AXIS STREAMS: rows are EPOCH-MAJOR, i.e.
%   [x1 y1 z1 | x2 y2 z2 | ...]. The temporal correlation model strides by
%   the axis count, so a stream stored axis-major comes out with its time law
%   stretched by the axis count -- and it will not announce itself.
%
%   GETTING nc AND nl: read them from the tables you were handed,
%   nc = size(common,1) and nl = size(local,1). Do not hardcode them; that is
%   the duplication the parametrisation pipe removes.
    error('lsqvce_template:notImplemented', ...
        ['myBuildDay is a stub -- write it for your data.\n' ...
         'Asked for day %s, obsType ''%s''. Return .A .B .l .t (and .arcIdx, ' ...
         '.Cep for the error information); see `help lsqvce_template`.'], ...
         num2str(day), obsType);
end

%% ========================================================== small local helpers
function cov = aprioriCovariance(obsTypes, ~)
%APRIORICOVARIANCE  A-priori covariance spec per obsType, for LSQVCE.ENGINE.COVMODEL.
%   The TEMPORAL law only -- white and unit-level by default: the VCE loop
%   discovers the truth, so a deliberately plain start makes the discovery
%   visible. Replace with your best prior if you have one. The per-epoch
%   accuracy of the observations is not declared here: myBuildDay returns it
%   with the day (D.Cep), for scalar and multi-axis streams alike.
%
%   INPUT
%     obsTypes  1xT cellstr  the stream names
%     g         1x1 struct   the observation-group entry (unused here)
%   OUTPUT
%     cov       1x1 struct   .(obsType) = a covariance spec
    cov = struct();
    for k = 1:numel(obsTypes)
        cov.(obsTypes{k}) = struct('type','diag', 'sigma',1);
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
