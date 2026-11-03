function task = gravfire_data(commonData, myOptions, systemsData, opt)
%GRAVFIRE_DATA  Declare a GRAVFIRE preprocessing bundle as an lsqvce.TaskDefinition.
%   task = gravfire_data(commonData, myOptions, systemsData)
%   task = gravfire_data(..., 'DesignFcn', @myDesign)      % your A/B builder
%   task = gravfire_data(..., 'DryRun', true)              % random A/B: pipeline check
%
%   TEMPLATE (R5 M0 close-out, HANDOVER §12.8): the GRAVFIRE bundle
%   (commonData / myOptions / systemsData, as saved in testdata.mat) already
%   contains everything a task declaration asks for EXCEPT the design
%   matrices. This function is the declaration half; the design half is the
%   'DesignFcn' you supply (the variational-equation partials of GRAVFIRE).
%   The file is meant to be READ as the step-by-step guide
%   dev/GRAVFIRE_INTERFACE.md: the STEP numbers in the comments are its steps.
%
%   THE DECLARATION IS A LSQVCE.TASKDEFINITION, filled by flat add-calls. Each
%   call checks its arguments at once and names the one at fault, before any
%   file is read; R = lsqvce.fit(task, opts) then reads the days.
%     lsqvce.TaskDefinition(Days=...)    the calendar, myOptions.days
%     task.addCommon(names, Info=T)      the Stokes rows, columns of A in this order;
%                                        every row needs a UNIQUE name, so the name is
%                                        type_degree_order ('c_2_0', 's_2_1'), and the
%                                        14 GRAVFIRE columns travel as Info (one table
%                                        row per name, variables col1 .. col14;
%                                        ctx.Common.Info in Build)
%     sys = task.addSystem("GRACE", ...) ONE system: GraceA + GraceB + KBR -- the
%                                        range-rate couples the locals of both
%                                        satellites, so they cannot be two systems
%     sys.addLocal(names, Scope="day")   ONE day's local templates; the toolbox
%                                        instantiates them on every calendar day
%     sys.addStream(type, Build=@b, Noise=lsqvce.Noise.white(sigma))
%                                        posA, posB, kbr: one Build per stream, the
%                                        STARTING temporal law (the estimator updates it)
%     task.addPrior / sys.addPrior       column 11 (Tikhonov) of the two tables
%   There is NO arc schedule: the pieces Build adds with out.addArc ARE the
%   arcs of the day (the locals are 'day'-scoped, so no plan needs the piece
%   count before a read).
%
%   INPUT
%     commonData   1x1 struct   .parameter  R x 14 cell, the common table
%     myOptions    1x1 struct   the GRAVFIRE run options (.days, .observationGroups, ...)
%     systemsData  1xS cell     one bundle per system; bundle{1}.observations is a
%                               (days x groups) cell of per-day structs, .parameter the
%                               EXPANDED local table (one row per day and satellite),
%                               .numberOfPieces (days x groups), .variance (1 x groups)
%   NAME-VALUE
%     'DesignFcn'  function_handle   [A, B] = f(o, day, ctx)
%                  o = the (day, group) bundle struct, ctx = the LSQVCE.DAYCONTEXT of
%                  the stream's Build (ctx.Stream the observation type, ctx.Settings
%                  the options bound in STEP 6, ctx.Common / ctx.Local the templates
%                  with their GRAVFIRE columns as Info). Must return A (m x
%                  ctx.NumCommon, columns in addCommon order) and B (m x ctx.NumLocal,
%                  every column from ctx.cols(names) -- the locals are 'day'-scoped,
%                  so a name gives its column on this day; never an offset by hand),
%                  rows EPOCH-MAJOR for multi-axis streams, m = EVERY epoch of the day
%                  times its axes -- the flagged epochs included (their rows are
%                  zeroed and weighted out, whatever you return). Default: a function
%                  that ERRORS with gravfire:designMissing (so a forgotten design step
%                  fails loudly, wrapped by lsqvce:makeLoader:buildDay).
%     'DryRun'     logical, default false   replace the design step by DETERMINISTIC
%                  random A/B (seeded per day and observation type) so that the
%                  declaration, the cut into pieces, the per-epoch covariance, the
%                  column counts and the whole fit pipeline can be exercised BEFORE the
%                  design step exists. The numbers of a dry run mean nothing.
%   OUTPUT
%     task  1x1 lsqvce.TaskDefinition   the declaration: R = lsqvce.fit(task, opts);
%           disp(task) lists it, task.layout() names every column of A and B
%           (task.declaration() is the engines' compiled struct, for
%           engine-tier work and for inspecting the compiled tables)
%
%   ACCURACY
%     Not applicable: bookkeeping -- the bundle's tables, flags and per-epoch
%     covariances are placed into the declaration and the days unchanged; the
%     one conversion is the local spans and knots from [MJD sec] / MJD.frac to
%     seconds of day (STEP 4). A dry run's numbers are random and mean nothing.
%
%   ERRORS
%     gravfire:designMissing        DryRun false and no DesignFcn given (raised at read time)
%     gravfire:unknownQuantity      an observation group's observationQuantity is not
%                                   'kinematicPosition' or 'range-rate'
%     gravfire:designShape          DesignFcn returned A/B with the wrong row/column count
%     gravfire:badConstraint        column 11 is malformed, or the rows of one label
%                                   disagree on sigma or on the activation pass, or
%                                   (local rows) on the a-priori value, column 10
%     lsqvce:TaskDefinition:duplicateName   two common rows with one type, degree and
%                                   order, or two local templates of one name (raised by
%                                   the add-call); the other declaration checks are
%                                   lsqvce:TaskDefinition:* too (LSQVCE.TASKDEFINITION)
%     lsqvce:DayContext:*           a mistake of Build while a day is read (a
%                                   DesignFcn asking ctx.cols for an unknown name)
%
%   EXAMPLE
%     S    = load('testdata.mat');
%     task = gravfire_data(S.commonData, S.myOptions, S.systemsData, 'DryRun', true);
%     R    = lsqvce.fit(task, lsqvce.Options(Verbose=1));   % pipeline check
%     task = gravfire_data(S.commonData, S.myOptions, S.systemsData, 'DesignFcn', @gravfireDesign);
%     R    = lsqvce.fit(task, lsqvce.Options(Estimate="noise", Noise=struct('Grid', 2160), MaxPasses=7));
%
%   See also LSQVCE.TASKDEFINITION, LSQVCE.SYSTEMDEFINITION, LSQVCE.DAYCONTEXT,
%            LSQVCE.DAYDATA, LSQVCE.NOISE, LSQVCE.FIT, LSQVCE.OPTIONS,
%            dev/GRAVFIRE_INTERFACE.md
%
%   Updated 2026-11-03 (CET, UTC+01:00), WP-162: returns the task only (the
%   compiled struct is task.declaration()); a position stream declares level 1,
%   because iSpp is absolute and travels as Sigma (user decision 2026-11-03: the
%   bundle's sqrt(variance) on top counted the level twice in an Estimate "none"
%   solve). With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-11-02 (CET, UTC+01:00), WP-161: declared as a lsqvce.TaskDefinition
%   (addCommon / addSystem / addLocal / addStream / addPrior, the GRAVFIRE columns as
%   Info); one Build per stream adds the pieces as arcs; DesignFcn is f(o, day, ctx);
%   returns the task, and data = task.declaration() as a second output. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-09-09 (CEST, UTC+02:00), R6 WP-102 (v6): THE CELLS ARE THE
%   ARCS -- buildGravfireDay returns one cell per piece (A, B, l, t, Cep),
%   the arc schedule (STEP 5, arcsFcn) is gone: the locals are 'day'-scoped,
%   so no plan needs the piece count before a read, and the toolbox learns
%   the pieces of a day from its cells. With assistance from Claude
%   (Anthropic), model Claude Fable 5.1.
%
%   VERIFIED 2026-09-02 on the user's testdata.mat (R2026a): the declaration
%   passes every declaration check of fit; DryRun + Estimate="none" accumulates 2 days x 16
%   pieces (4 posA + 4 posB + 8 kbr), nc = 117, nl = 144 derived from the 72
%   templates, Cep accepted, in 2.0 s. Estimate="noise", Pool="stream", ArcFactors=false runs
%   (62 s, 2 passes). R5 M3 made the 3-axis pieces fast on the joint estimator (38.6 s for both
%   days at 2160 lags, HANDOVER §12.5 round 3c).
%
%   Updated 2026-09-07 (CEST, UTC+02:00), R6 WP-95: the v5 data contract is
%   named as such, STEP 6 emits ONE prior struct per label with the member
%   rows stacked into L (no 'group' beside L), and STEP 7 documents the
%   declared data.systems.priors instead of the retired priorDay handle.
%   With assistance from Claude (Anthropic), model Claude Opus 5.
%   Updated 2026-09-06 (CEST, UTC+02:00), R6 WP-90: the positions declare
%   the temporal law only and return D.Cep with the day (no 'separable'
%   type, no axis count in the spec). With assistance from Claude
%   (Anthropic), model Claude Fable 5.1.
%   Created 2026-09-02 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M0 close-out (§12.8).
%
%   Updated 2026-10-17 (CEST, UTC+02:00), WP-148: the flagged epochs stay on
%   the grid with weight 0 (D.excluded) instead of being removed -- a removed
%   epoch left a gap a 'psd' law took for one step. DesignFcn returns the rows
%   of every epoch. With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-17 (CEST, UTC+02:00), WP-144 (A4): the members of a
%   labelled LOCAL group must share fromPass and the a-priori value too (a
%   later member's were dropped). With assistance from Claude (Anthropic),
%   model Claude Opus 5.5.
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-124: references to the data
%   contract point at the API reference. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.

arguments
    commonData  (1,1) struct
    myOptions   (1,1) struct
    systemsData (1,:) cell
    opt.DesignFcn = @designMissing
    opt.DryRun (1,1) logical = false
end
    sd   = systemsData{1};
    days = double(myOptions.days(:).');                       % 1 x D, integers (MJD)

    % ---- STEP 2: the calendar -- the task is made with it (STEP 1 needs the task)
    task = lsqvce.TaskDefinition(Days=days);

    % ---- STEP 1: the common parameters -- one unique name per row, in table
    %      order (the columns of A); your 14 columns ride along as Info
    P = commonData.parameter;                                 % R x 14
    task.addCommon(commonNames(P), Info=gravfireInfo(P));

    % ---- STEP 3: observation types are the STREAM names
    groups = myOptions.observationGroups;                     % 1 x G cell of structs
    G      = numel(groups);
    obsTypes = cell(1, G);
    for g = 1:G
        obsTypes{g} = obsTypeOf(groups{g});
    end

    % ---- STEP 4: the local templates -- ONE day's rows, all 'day'-scoped,
    %      times as seconds-of-day (the day itself is NOT in a template)
    T = templateFromExpanded(sd.parameter, days(1));          % n x 14

    % ---- STEP 5: the arc schedule -- NONE (v6). The pieces of a (day, group)
    %      are the arcs Build adds, one addArc per piece, and the toolbox takes
    %      them as the arcs of the day. numberOfPieces is not restated.
    %      (A schedule, addSystem's Arcs, is required only by Scope="arc" locals.)

    % ---- STEP 6: what Build needs, bound once -- the system's Settings,
    %      read in Build as ctx.Settings
    po = struct();
    po.days      = days;
    po.obsTypes  = obsTypes;                                  % 1 x G cell
    po.bundles   = sd.observations;                           % (days x groups) cell
    po.designFcn = opt.DesignFcn;
    po.dryRun    = opt.DryRun;
    po.myOptions = myOptions;

    % ---- STEP 7: the system, its local templates and its streams. Id is
    %      numeric or ONE character; the label is the name. The stream's Noise
    %      is the TEMPORAL law only (white to start; the estimator updates it):
    %      the per-epoch error information of the positions -- the 3x3
    %      covariance from iSpp -- comes back from Build with each arc (Sigma=,
    %      STEP 9) and the toolbox composes C = D (R_t kron I_3) D; the axis
    %      count is derived from it, so no Axes is declared. iSpp is ABSOLUTE
    %      (user decision 2026-11-03), so a position stream's law has level 1:
    %      Sigma is relative to it, and sqrt(variance) on top would count the
    %      level twice. A stream without iSpp (kbr) takes its level from the
    %      bundle's variance.
    sys = task.addSystem("GRACE", Settings=po, Id=1);
    if ~isempty(T)
        sys.addLocal(string(T(:, 1)), Scope="day", Info=gravfireInfo(T));
    end
    for g = 1:G
        level = sqrt(sd.variance(g));                         % kbr: the bundle's variance
        if strcmp(char(groups{g}.observationQuantity), 'kinematicPosition')
            level = 1;                                        % positions: Sigma = iSpp, absolute
        end
        sys.addStream(obsTypes{g}, Build=@buildGravfireDay, Noise=lsqvce.Noise.white(level));
    end

    % ---- STEP 6 (R5-D8 round 1): column 11 on COMMON rows -> task.addPrior
    % Your constraining factor {'Tikhonov', k, sigma[, label]} on a COMMON row
    % j constrains that parameter to x0 = column 10 (the a-priori, [] -> 0)
    % with standard deviation sigma, from pass k on. ONE prior group per
    % label: without a label the row is its own group, L = e_j' and sigma
    % fixed (Estimate="none" -- one row has redundancy <= 1); a 4th entry is a
    % group LABEL, and its rows become ONE group whose L stacks their unit
    % rows, with one shared sigma estimated by VCE (Estimate="scale",
    % Koch-Kusche). A group is declared by L OR by a name pattern (On=), never
    % both; these are declared by L, in the columns of A (addCommon order).
    cp = priorsFromColumn11(P);
    for k = 1:numel(cp)
        task.addPrior(cp(k).name, L=cp(k).L, Value=cp(k).x0, Sigma=cp(k).sigma, ...
                      Estimate=cp(k).estimate, From=cp(k).from);
    end
    % ---- STEP 7 (R5-D8 round 6): column 11 on LOCAL rows -> sys.addPrior
    % Regularisation of local parameters is DECLARED, never computed: one
    % prior group per label on the system, given EITHER by On= -- a regexp on
    % the local NAMES, from which the toolbox builds L with one row per matched
    % template per INSTANCE, so a 'day' template contributes one row per day
    % and the group spans the calendar -- OR by L=, written in the columns of B
    % as task.layout() announces them. Never both. Rows sharing a label (4th
    % entry of column 11, e.g. 'accscale') are ONE group with one sigma and
    % Estimate="scale"; without a label each constrained row is its own group
    % with a fixed sigma. How many groups your bundle yields therefore depends
    % on how many templates carry a constraint and how they are labelled. The
    % helper below emits the On= form.
    lp = gravfireLocalPriors(T);
    for k = 1:numel(lp)
        sys.addPrior(lp(k).name, On=lp(k).group, Value=lp(k).x0, Sigma=lp(k).sigma, ...
                     Estimate=lp(k).estimate, From=lp(k).from);
    end
end

% ======================================================================
function out = buildGravfireDay(day, ctx)
%BUILDGRAVFIREDAY  One (day, observation type): raw A, B, l, t, Sigma -- ONE ARC PER PIECE.
%   The Build of every stream (addStream): ctx.Stream names the observation
%   type, ctx.Settings holds what STEP 6 bound. Returns the LSQVCE.DAYDATA of
%   ctx.newDay(), one addArc per piece in piece order.
    po      = ctx.Settings;
    obsType = char(ctx.Stream);                               % 'posA' | 'posB' | 'kbr'
    g  = groupIndex(po.obsTypes, obsType);
    d  = find(po.days == day, 1);
    o  = po.bundles{d, g};

    nAx = o.numberOfAxes;
    % The flagged epochs STAY on the grid with weight 0 (WP-148): removing one
    % leaves a gap inside its piece that a 'psd' law, counting lags in
    % samples, takes for ONE step (lsqvce:covModel:unevenGrid). Their rows are
    % zeroed below and passed as Excluded; the toolbox deletes them from the
    % law's covariance, so their values -- NaN included -- never matter.
    flag  = logical(o.nanIndexOfPiece{1}(:));                 % nEp x 1
    lobs  = o.lobs{1};                                        % nEp x nAx
    lobs(flag, :) = 0;
    tEp   = o.time{1}(:, 2);                                  % seconds of day
    piece = o.indexOfPiece{1}(:);                             % nEp x 1
    nEp   = size(lobs, 1);
    m     = nEp * nAx;
    ex    = reshape(repmat(flag.', nAx, 1), [], 1);           % m x 1, epoch-major

    % STEP 8: epoch-major rows (x1 y1 z1 x2 y2 z2 ...), the whole day first
    l = reshape(lobs.', [], 1);
    t = repelem(tEp(:), nAx);

    % STEP 9: the ERROR INFORMATION -- the per-epoch 3x3 covariance from
    % iSpp as nEpoch x 3 x 3 (the axis count is derived from it, WP-90). A
    % scalar stream with a per-epoch sigma would carry it as nEpoch x 1.
    Cep = [];
    if nAx > 1
        C   = o.iSpp{1};                                      % nEp x 9 (3x3 per epoch)
        C(flag, :) = repmat(mean(C(~flag, :), 1), nnz(flag), 1);   % deleted: any valid 3x3
        Cep = permute(reshape(C.', nAx, nAx, []), [3 1 2]);   % nEp x 3 x 3
    end

    % STEP 10: the design matrices -- yours, or the dry-run stand-in (whole
    % day, m x ctx.NumCommon / m x ctx.NumLocal; the cut into pieces follows)
    nc = ctx.NumCommon;
    nl = ctx.NumLocal;
    if po.dryRun
        sd = mod(7919 * day + sum(double(obsType)), 2^31 - 1);
        rs = RandStream('mt19937ar', 'Seed', sd);             % deterministic per call
        A = randn(rs, m, nc);
        B = randn(rs, m, nl);
    else
        [A, B] = po.designFcn(o, day, ctx);
        if ~isequal(size(A), [m nc]) || ~isequal(size(B), [m nl])
            error('gravfire:designShape', ...
                  ['DesignFcn returned A %s / B %s; expected A [%d %d] and B [%d %d] ' ...
                   '(day %d, %s).'], mat2str(size(A)), mat2str(size(B)), ...
                   m, nc, m, nl, day, obsType);
        end
    end
    A(ex, :) = 0;  B(ex, :) = 0;                              % weight 0: zero rows (API reference)

    % STEP 11: THE PIECES ARE THE ARCS (v6) -- one addArc per piece, in piece
    % order; the toolbox takes them as this day's arcs (ids 1:nPiece). B is
    % the DesignFcn's m x nl matrix (B=), Sigma the piece's per-epoch
    % covariance, Excluded its flagged rows.
    ids = unique(piece(:)).';
    out = ctx.newDay();
    for k = 1:numel(ids)
        ep   = find(piece == ids(k));                          % this piece's epochs
        rows = (ep(:) - 1) * nAx + (1:nAx);                   % nEp_k x nAx, epoch per row
        rows = reshape(rows.', [], 1);                         % their epoch-major rows
        sigma = [];
        if nAx > 1, sigma = Cep(ep, :, :); end
        out.addArc(A(rows, :), l(rows), t(rows), B=B(rows, :), Sigma=sigma, Excluded=ex(rows));
    end
end

% ======================================================================
function [A, B] = designMissing(varargin)
    error('gravfire:designMissing', ...
          ['No DesignFcn given: the design matrices are not part of the ' ...
           'preprocessing bundle. Pass ''DesignFcn'', @yourDesign (or ''DryRun'', true ' ...
           'for a pipeline check).']);
end

function name = obsTypeOf(group)
%OBSTYPEOF  'GraceA' + kinematicPosition -> 'posA'; range-rate -> 'kbr'.
    q = char(group.observationQuantity);
    switch q
        case 'kinematicPosition'
            name = ['pos' upper(group.name(end))];
        case 'range-rate'
            name = 'kbr';
        otherwise
            error('gravfire:unknownQuantity', ...
                  'observation group %s has quantity %s (expected kinematicPosition or range-rate).', ...
                  group.name, q);
    end
end

function g = groupIndex(obsTypes, obsType)
    g = find(strcmp(obsTypes, char(obsType)), 1);
end

function names = commonNames(P)
%COMMONNAMES  The unique name of every common row, type_degree_order ('c_2_0').
%   INPUT   P      R x 14 cell   the common table (columns 1-3: 'c'|'s', degree, order)
%   OUTPUT  names  1 x R string
    names = strings(1, size(P, 1));
    for j = 1:size(P, 1)
        names(j) = commonName(P, j);
    end
end

function nm = commonName(P, j)
%COMMONNAME  The name of common row j -- the ONE rule, used by addCommon and the priors.
    nm = sprintf('%s_%d_%d', char(string(P{j, 1})), double(P{j, 2}), double(P{j, 3}));
end

function T = gravfireInfo(P)
%GRAVFIREINFO  The GRAVFIRE columns of a parameter table as an Info table.
%   One table row per parameter; variable colK holds GRAVFIRE column K as a
%   cell column, its entries as the bundle has them (ctx.Common.Info /
%   ctx.Local.Info in Build give one struct per parameter with fields col1 ..).
%   INPUT   P  n x C cell
%   OUTPUT  T  n x C table
    vars = cell(1, size(P, 2));
    for c = 1:size(P, 2)
        vars{c} = P(:, c);
    end
    T = table(vars{:}, 'VariableNames', compose('col%d', 1:size(P, 2)));
end

function P = templateFromExpanded(P, day0)
%TEMPLATEFROMEXPANDED  One day's rows of the EXPANDED local table, as templates.
%   Selects the rows whose flag column (6) names day0, converts the absolute
%   [MJD sec] spans (columns 2, 3) to seconds-of-day, converts the knot table
%   (column 4: [MJD-fraction value]) to seconds-of-day, and drops the day from
%   the flags. The 14 GRAVFIRE columns stay in position (column 1 the name);
%   the scope is addLocal's Scope="day". The spans are NOT declared as
%   validity (addLocal's Valid): the templates of day0 serve every calendar
%   day. When GRAVFIRE emits templates itself, this helper disappears.
%   INPUT   P     N x 14 cell   the expanded local table
%           day0  1x1 double    the day whose rows are the templates (MJD)
%   OUTPUT  P     n x 14 cell   one day's templates
    flags = P(:, 6);
    sel   = cellfun(@(f) ~isempty(f) && f(1) == day0, flags);
    P     = P(sel, :);
    n     = size(P, 1);
    for k = 1:n
        for c = [2 3]                                         % [MJD sec] -> sec
            v = P{k, c};
            if numel(v) == 2, P{k, c} = v(2); end
        end
        v = P{k, 4};                                          % knots [MJD.frac value]
        if ~isempty(v) && size(v, 2) == 2
            P{k, 4} = [(v(:, 1) - day0) * 86400, v(:, 2)];
        end
        f = P{k, 6};
        if ~isempty(f), P{k, 6} = f(2:end); end               % drop the day
    end
end

%==========================================================================
function priors = priorsFromColumn11(P)
%PRIORSFROMCOLUMN11  Common-row constraints -> the prior groups of task.addPrior (R5-D8 round 1).
%   P is the 14-column common table: column 10 = a-priori value ([] -> 0),
%   column 11 = 'none' | [] | {'Tikhonov', fromPass, sigma[, label]}.
%   ONE struct per GROUP, regrouped for R5-D8 round 6: rows sharing a label
%   become one group whose L STACKS their unit rows e_j', with one shared
%   sigma and one activation pass, estimated by VCE (estimate "scale",
%   Koch-Kusche -- several rows are what make a sigma estimable); a row
%   without a label is its own group, L = e_j' and sigma fixed (estimate
%   "none": one row has redundancy <= 1). Every group is declared by L in
%   the columns of A (addCommon order). Members of one label must agree on
%   sigma and fromPass.
%   OUTPUT  priors  1xG struct (.L .x0 .sigma .name .estimate .from) | []
    priors = struct('L', {}, 'x0', {}, 'sigma', {}, 'name', {}, 'estimate', {}, 'from', {});
    nc = size(P, 1);
    groups = containers.Map('KeyType', 'char', 'ValueType', 'any');  order = {};
    for j = 1:nc
        c = P{j, 11};
        if isNoneConstraint(c), continue, end
        if ~iscell(c) || numel(c) < 3 || ~ischar(c{1})
            error('gravfire:badConstraint', ...
                  'row %d: column 11 must be ''none'' or {''Tikhonov'', fromPass, sigma[, label]}.', j);
        end
        if ~strcmpi(c{1}, 'Tikhonov')
            error('gravfire:badConstraint', ...
                  'row %d: unknown constraint kind ''%s'' (only Tikhonov is mapped).', j, c{1});
        end
        x0 = 0;
        if ~isempty(P{j, 10}), x0 = double(P{j, 10}); end
        lab = '';
        if numel(c) >= 4 && ~isempty(c{4}), lab = char(string(c{4})); end
        nm  = commonName(P, j);
        key = lab;  if isempty(key), key = nm; end
        if ~isKey(groups, key)
            groups(key) = struct('rows', [], 'x0', [], 'sigma', double(c{3}), ...
                                 'from', double(c{2}), 'labelled', ~isempty(lab));
            order{end+1} = key; %#ok<AGROW>
        end
        g = groups(key);
        if g.sigma ~= double(c{3})
            error('gravfire:badConstraint', ...
                  'group ''%s'': members must share one sigma (row %d has %g, the group %g).', ...
                  key, j, double(c{3}), g.sigma);
        end
        if g.from ~= double(c{2})
            error('gravfire:badConstraint', ...
                  'group ''%s'': members must share one fromPass (row %d has %g, group %g).', ...
                  key, j, double(c{2}), g.from);
        end
        g.rows(end+1) = j;  g.x0(end+1, 1) = x0;
        groups(key) = g;
    end
    for k = 1:numel(order)
        g = groups(order{k});
        L = zeros(numel(g.rows), nc);
        for r = 1:numel(g.rows), L(r, g.rows(r)) = 1; end
        est = 'none';  if g.labelled, est = 'scale'; end
        priors(end+1) = struct('L', L, 'x0', g.x0, 'sigma', g.sigma, 'name', order{k}, ...
                               'estimate', est, 'from', g.from); %#ok<AGROW>
    end
    if isempty(priors), priors = []; end
end

%==========================================================================
function P = gravfireLocalPriors(T)
%GRAVFIRELOCALPRIORS  Column 11 of the LOCAL templates -> the prior groups of sys.addPrior (R5-D8 round 6).
%   One struct per GROUP in the PATTERN form (On=): templates sharing a label
%   (4th entry of column 11) become one group whose pattern is the regexp of
%   their names ('^(GraceA_xaccscale|GraceB_xaccscale)$'); templates without
%   a label are one group each. The toolbox builds L from the names (column 1)
%   for every instance -- per day for a 'day' template. x0 = column 10, sigma
%   and fromPass from the entry; estimate "scale" for labelled groups
%   (Koch-Kusche), "none" otherwise.
%   INPUT   T  n x 14 cell   one day's local templates (templateFromExpanded)
%   OUTPUT  P  1xG struct (.group .x0 .sigma .from .name .estimate) | []
    P = struct('group', {}, 'x0', {}, 'sigma', {}, 'from', {}, 'name', {}, 'estimate', {});
    groups = containers.Map('KeyType', 'char', 'ValueType', 'any');  order = {};
    for j = 1:size(T, 1)
        c = T{j, 11};
        if isNoneConstraint(c), continue, end
        if ~iscell(c) || numel(c) < 3 || ~ischar(c{1}) || ~strcmpi(c{1}, 'Tikhonov')
            error('gravfire:badConstraint', ...
                  'local row %d: column 11 must be ''none'' or {''Tikhonov'', fromPass, sigma[, label]}.', j);
        end
        x0 = 0;  if ~isempty(T{j, 10}), x0 = double(T{j, 10}); end
        lab = '';  if numel(c) >= 4 && ~isempty(c{4}), lab = char(string(c{4})); end
        nm = char(string(T{j, 1}));
        key = lab;  if isempty(key), key = nm; end
        if ~isKey(groups, key)
            groups(key) = struct('names', {{}}, 'x0', x0, 'sigma', double(c{3}), 'from', double(c{2}), ...
                                 'estimate', "none", 'labelled', ~isempty(lab));
            order{end+1} = key; %#ok<AGROW>
        end
        g = groups(key);
        if g.sigma ~= double(c{3})
            error('gravfire:badConstraint', ...
                  'group ''%s'': members must share one sigma (row %d has %g, the group %g).', ...
                  key, j, double(c{3}), g.sigma);
        end
        % WP-144 (A4): a later member's fromPass and a-priori value were
        % DROPPED -- the group kept its first member's (reproduced). fromPass
        % is refused as on the common rows; x0 too, because the pattern form
        % builds L per member and per day in the toolbox's order, so one value
        % per member cannot be written here without restating that order.
        if g.from ~= double(c{2})
            error('gravfire:badConstraint', ...
                  'group ''%s'': members must share one fromPass (row %d has %g, the group %g).', ...
                  key, j, double(c{2}), g.from);
        end
        if g.x0 ~= x0
            error('gravfire:badConstraint', ['group ''%s'': members of a local group must share one a-priori ' ...
                  'value, column 10 (row %d has %g, the group %g); give them separate labels, or ' ...
                  'declare the group by L in the columns task.layout() announces.'], key, j, x0, g.x0);
        end
        g.names{end+1} = regexptranslate('escape', nm);
        groups(key) = g;
    end
    for k = 1:numel(order)
        g = groups(order{k});
        est = "none";  if g.labelled, est = "scale"; end
        P(end+1) = struct('group', ['^(' strjoin(unique(g.names), '|') ')$'], 'x0', g.x0, ...
                          'sigma', g.sigma, 'from', g.from, 'name', order{k}, 'estimate', est); %#ok<AGROW>
    end
    if isempty(P), P = []; end
end
%==========================================================================
function tf = isNoneConstraint(c)
%ISNONECONSTRAINT  'none' in every spelling the tables use: [], 'none', {'none'}.
    if iscell(c) && isscalar(c), c = c{1}; end
    tf = isempty(c) || (ischar(c) && strcmpi(c, 'none')) || (isstring(c) && strcmpi(c, "none"));
end
