function [cfg, data, task] = demo_grace_dataLoader()
%DEMO_GRACE_DATALOADER  Synthetic GRACE-FO data for the VCE demo, declared as a lsqvce.TaskDefinition.
%
%   [cfg, data, task] = demo_grace_dataLoader()
%   cfg = demo_grace_dataLoader()
%
%   Returns a ready-to-run configuration CFG and the declaration TASK, so the
%   VCE loop on multi-axis and scalar streams can be exercised end to end
%   without real data:
%
%       [cfg, ~, task] = demo_grace_dataLoader();
%       R = lsqvce.fit(task, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false));
%
%   DATA is the same declaration compiled to the engine-tier struct,
%   task.declaration(). It stays the second output because the ENGINE-tier
%   fixtures and the tests that edit a declaration take the struct
%   (fixture_grace, test_vcePolicy, test_fit, test_speedBudget); lsqvce.fit
%   compiles a task to exactly that struct, so both give the same fit.
%
%   TOPOLOGY (real GRACE-FO, "Option B", plus an independent third satellite).
%     System 1  -- the C/D PAIR.  TWO satellites C, D each with a 3-D position
%       stream (posC, posD, covariance C = D (R_t kron I_3) D with a FIXED
%       per-epoch 3x3 D that each Build passes with every arc as Sigma= --
%       the toolbox composes it with the temporal law), and ONE inter-satellite
%       range-rate stream (rr, scalar) whose design couples BOTH satellites'
%       local parameters. Because rr touches two satellites' locals and the
%       block-arrow model has one local block per system, C and D are MERGED
%       into one pair-system with a stacked local block [y_C | y_D].
%     System 2  -- satellite E, INDEPENDENT.  A third satellite with ONLY a 3-D
%       position stream (posE) and NO inter-satellite link, so it is a clean,
%       separate block-arrow system with its own local parameters y_E. It
%       observes the SAME common gravity field x_c as the pair (that is what
%       makes combining the systems meaningful), but its locals are its own.
%
%   Both systems share the common (gravity) parameters; VCE pools by stream,
%   giving FOUR variance components (posC, posD, rr, posE). Each stream starts
%   from a deliberately WRONG a-priori noise, lsqvce.Noise.white(1).
%
%   THE DECLARATION (graceTask, the TEMPLATE for your own case)
%       task = lsqvce.TaskDefinition(Days=cfg.days);
%       task.addCommon(names, Info=...)                      6 gravity coefficients
%       pair = task.addSystem("demo-CD-pair", Id=1, Settings=cfg.procOpt);
%       pair.addLocal(names, Scope="system", Info=...)       the block [y_C | y_D]
%       pair.addStream("posC", Build=@..., Noise=lsqvce.Noise.white(1), Axes=3)
%       ... posD, rr; then system "demo-E" (Id=2) with its stream posE
%   ONE Build per stream, out = build(day, ctx): it draws the stream's rows,
%   puts each local partial in the column ctx.cols(names) announces (never an
%   offset of its own) and adds the day's arcs with out.addArc. No Arcs
%   schedule is declared: the arcs of a day are those Build adds (two here),
%   and the locals are Scope="system", so no plan needs the arc count before
%   a day is read. The per-epoch 3x3 covariance of a position stream is not
%   part of its law: it travels with each arc (addArc Sigma=, mEpoch x 3 x 3)
%   -- there is no second place to state it.
%
%   THE PARAMETER CATALOGUE.  cfg.common, cfg.localPair and cfg.localE hold one
%   ROW per unknown in 15 columns -- this demo's convention, yours may differ:
%
%     col  1  scope           'common' | 'system'; the add-call states it
%                             as Scope=                     (char)
%     col  2  kind            'gravity' | 'arcBias' | ...   (char)
%     col  3  name            the parameter's NAME, unique in the task (char)
%     col  4  system          system id this unknown belongs to (numeric; 0 =
%                             common/shared)
%     col  5  satellite       'C'|'D'|'E'|'' ('' = shared)  (char)
%     col  6  degree  n       SH degree   (gravity) or []   (numeric)
%     col  7  order   m       SH order    (gravity) or []   (numeric)
%     col  8  cs              'C'|'S' coeff (gravity) or ''  (char)
%     col  9  arc             arc index (local) or []        (numeric)
%     col 10  component       'x'|'y'|'z'|'' axis or ''      (char)
%     col 11  unit            e.g. '-'|'m'|'m/s'             (char)
%     col 12  a-priori value  x0                             (numeric)
%     col 13  a-priori sigma  or Inf if free                 (numeric)
%     col 14  estimate flag   true/false                     (logical)
%     col 15  group tag       free grouping label            (char)
%
%   Column 3 becomes the name, columns 2 and 4..15 the parameter's Info table
%   (catalogueInfo), which every Build reads as ctx.Common.Info and
%   ctx.Local.Info: this demo's Builds pick a satellite's locals by
%   Info.satellite. The toolbox reads the names and scopes; the other
%   MEANINGS above are this demo's own.
%
%   INPUT
%     (none)
%   OUTPUT
%     cfg   1x1 struct  the configuration and planted truth: .nc .nlSat
%                       .nlPair .nlE .days .arcsPerDay .mEpoch .dt .seed,
%                       .x_c (nc x 1) .y_C .y_D .y_E (nlSat x 1), .trueAR
%                       .trueSigma (per stream) .Dep (3x3), the catalogue
%                       .common (nc x 15 cell) .localPair (2 nlSat x 15)
%                       .localE (nlSat x 15), .procOpt (the systems' Settings),
%                       .obsTypes .obsTypesE .allObsTypes (stream names)
%     data  1x1 struct  task.declaration(): the engine-tier declaration
%                       struct (commonParams, days, systems(1x2))
%     task  1x1 lsqvce.TaskDefinition   6 common parameters; system
%                       "demo-CD-pair" (8 locals, streams posC posD rr) and
%                       system "demo-E" (4 locals, stream posE)
%
%   ACCURACY
%     Not applicable: a deterministic generator (one seeded RandStream per day
%     and stream), no floating-point result of its own. The WP-161 migration
%     keeps every seed, the order of every draw and the DENSE local design of
%     the v6 builder, so the rows are those of the v6 declaration by
%     construction (not measured in the migrating session, which had no MATLAB).
%
%   ERRORS
%     None of its own. A declaration mistake raises lsqvce:TaskDefinition:* at
%     the add-call that makes it; a mistake inside a Build raises
%     lsqvce:DayContext:* while the day is read. (The v6 builder's
%     demo_grace_dataLoader:obsType left with its switch: every stream has its
%     own Build.)
%
%   See also LSQVCE.TASKDEFINITION, LSQVCE.SYSTEMDEFINITION, LSQVCE.NOISE,
%            LSQVCE.DAYCONTEXT, LSQVCE.DAYDATA, LSQVCE.FIT, LSQVCE.OPTIONS, DEMO_GRACE,
%            DEMO_GRAVITY_DATALOADER
%
%   Updated 2026-11-02 (CET, UTC+01:00), WP-161: declared as a lsqvce.TaskDefinition, one Build per stream. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-09-09 (CEST, UTC+02:00), R6 WP-102 (v6): THE CELLS ARE THE
%   ARCS -- demoBuildDay returns A, B, l, t and Cep as one cell per arc and
%   the declaration carries no `arcs` (the cells are the schedule; the
%   locals are 'system'-scoped). With assistance from Claude (Anthropic),
%   model Claude Fable 5.1.
%   Updated 2026-09-07 (CEST, UTC+02:00), R6 WP-95: the contract is named v5,
%   the column catalogue is renumbered onto the 15 columns the tables really
%   carry (the scope token is column 1, checked by the toolbox), and the dead
%   DEMO_SEPARABLE See also is replaced by live units.
%   With assistance from Claude (Anthropic), model Claude Opus 5.
%   Updated 2026-09-06 (CEST, UTC+02:00), R6 WP-90: the error information
%   travels with the day -- the positions declare the temporal law only
%   ('diag', deliberately wrong) and return D.Cep; the epochCovariance hook
%   and the 'separable' type are gone. With assistance from Claude
%   (Anthropic), model Claude Fable 5.1.
%   Updated 2026-08-31 (CEST, UTC+02:00), v4 SCOPED tables (Level L11,
%   WP-77): the declaration tables carry the scope column (commonParams
%   'common'|'day', localParams 'arc'|'day'|'system'); localsPerArc left the
%   USER surface. With assistance from Claude (Anthropic), model
%   Claude Fable 5.
%   Updated 2026-08-31 (CEST, UTC+02:00), v3 data contract (Level L10,
%   WP-74): USER surface only -- returns [cfg, data], the declaration
%   (obsTypes from the covariance fieldnames, localParams tables, id);
%   the ENGINE-tier loader construction moved to dev/tests/fixture_grace.m.
%   With assistance from Claude (Anthropic), model Claude Fable 5.
%   Created 2026-07-16 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.
%
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-124: references to the data
%   contract point at the API reference. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-29 (CET, UTC+01:00), WP-159: the declaration is the data
%   struct of the API reference (it was called "v5"; since v6 a day is one
%   cell per arc). With assistance from Claude (Anthropic), model Claude
%   Opus 5.5.
    cfg  = demoCfg();
    task = graceTask(cfg);
    data = task.declaration();       % the engine-tier struct of the same task
end

%==========================================================================
function task = graceTask(cfg)
%GRACETASK  The declaration: the six gravity coefficients, the merged C/D
%   pair (posC, posD, rr) and independent E (posE). The streams are added
%   in processing order, each with its own Build and a deliberately WRONG
%   a-priori noise (white, sigma 1): the loop discovers the colour.
%   INPUT   cfg   1x1 struct  see demoCfg
%   OUTPUT  task  1x1 lsqvce.TaskDefinition
    task = lsqvce.TaskDefinition(Days=cfg.days);
    task.addCommon(string(cfg.common(:, 3)), Scope="common", Info=catalogueInfo(cfg.common));
    nAx   = size(cfg.Dep, 1);                 % 3: the axes of a position epoch
    wrong = lsqvce.Noise.white(1);            % the a-priori law of every stream

    % Id: a number or a single character -- it must survive [store.sysID]
    % concatenation as ONE element, never a multi-char string like 'CD'.
    pair = task.addSystem("demo-CD-pair", Id=1, Settings=cfg.procOpt);
    pair.addLocal(string(cfg.localPair(:, 3)), Scope="system", Info=catalogueInfo(cfg.localPair));
    pair.addStream("posC", Build=@(day, ctx) buildPosC(day, ctx, cfg), Noise=wrong, Axes=nAx);
    pair.addStream("posD", Build=@(day, ctx) buildPosD(day, ctx, cfg), Noise=wrong, Axes=nAx);
    pair.addStream("rr",   Build=@(day, ctx) buildRr(day, ctx, cfg),   Noise=wrong);

    sysE = task.addSystem("demo-E", Id=2, Settings=cfg.procOpt);
    sysE.addLocal(string(cfg.localE(:, 3)), Scope="system", Info=catalogueInfo(cfg.localE));
    sysE.addStream("posE", Build=@(day, ctx) buildPosE(day, ctx, cfg), Noise=wrong, Axes=nAx);
end

function info = catalogueInfo(T)
%CATALOGUEINFO  The catalogue's descriptive columns (2 and 4..15) as the Info
%   table of an add-call, one row per parameter (column 3 is the name).
%   INPUT   T     n x 15 cell   a parameter catalogue (commonParamTable,
%                               localParamTable)
%   OUTPUT  info  n x 13 table  readable in Build as ctx.Common.Info /
%                               ctx.Local.Info
    info = cell2table(T(:, [2, 4:15]), 'VariableNames', ...
        {'kind', 'system', 'satellite', 'degree', 'order', 'cs', 'arc', 'component', ...
         'unit', 'apriori', 'aprioriSigma', 'estimate', 'group'});
end

function cfg = demoCfg()
%DEMOCFG  Config: C/D pair (posC,posD,rr) + independent satellite E (posE).
    cfg.nc         = 6;               % common (gravity) parameters
    cfg.nlSat      = 4;               % local (arc) parameters per satellite
    cfg.nlPair     = 2*cfg.nlSat;     % merged C/D pair-system local block
    cfg.nlE        = cfg.nlSat;       % satellite E local block
    cfg.obsTypes   = {'posC','posD','rr'};   % system 1 (pair) streams
    cfg.obsTypesE  = {'posE'};               % system 2 (E) streams
    cfg.allObsTypes = [cfg.obsTypes, cfg.obsTypesE];  % all streams, both systems
    cfg.days       = 1:4;
    cfg.arcsPerDay = 2;
    cfg.mEpoch     = 300;             % epochs per arc
    cfg.dt         = 5;
    cfg.seed       = 7;
    cfg.x_c        = (1:cfg.nc).';
    cfg.y_C        = 0.10 * (1:cfg.nlSat).';
    cfg.y_D        = 0.15 * (1:cfg.nlSat).';
    cfg.y_E        = 0.12 * (1:cfg.nlSat).';
    % TRUE noise. Positions: separable C = D (R_t kron I_3) D with a FIXED
    % per-epoch 3x3 D and a shared temporal AR; rr: scalar white. Streams differ
    % in level so VCE has real work; a-priori start is 1.
    cfg.trueAR.posC = [1 -0.5 0.1];
    cfg.trueAR.posD = [1 -0.4 0.05];
    cfg.trueAR.rr   = 1;
    cfg.trueAR.posE = [1 -0.3 0.08];
    cfg.trueSigma   = struct('posC',0.8, 'posD',1.3, 'rr',1.6, 'posE',1.1);
    % a representative per-epoch 3x3 covariance factor (same each epoch here for
    % simplicity; a real PPP solution varies it per epoch)
    Ld = [1 0 0; 0.3 1.0 0; 0.2 0.4 2.5];   % lower-tri "Cholesky" of D
    cfg.Dep = Ld * Ld.';                     % 3x3 SPD per-epoch covariance

    % ---- the parameter catalogue (15 columns, one row per unknown) ---------
    cfg.common = commonParamTable(cfg);      % nc rows -> nc derived
    cfg.localPair = localParamTable(cfg, {'C','D'}, [cfg.y_C, cfg.y_D]); % 2*nlSat
    cfg.localE    = localParamTable(cfg, {'E'},     cfg.y_E);            % nlSat
    % the systems' Settings: what every Build reads as ctx.Settings
    cfg.procOpt   = struct('mEpoch',cfg.mEpoch, 'dt',cfg.dt, ...
                           'arcsPerDay',cfg.arcsPerDay, 'seed',cfg.seed);
end

%% ===================================================== the Builds, one per stream
function out = buildPosC(day, ctx, cfg)
%BUILDPOSC  Satellite C's 3-D positions on DAY: they see C's locals only.
%   INPUT day 1x1 double; ctx 1x1 lsqvce.DayContext; cfg 1x1 struct.
%   OUTPUT out 1x1 lsqvce.DayData.
    out = graceDay(day, ctx, cfg, size(cfg.Dep, 1), {'C', 1});
end

function out = buildPosD(day, ctx, cfg)
%BUILDPOSD  Satellite D's 3-D positions on DAY: they see D's locals only.
%   INPUT day 1x1 double; ctx 1x1 lsqvce.DayContext; cfg 1x1 struct.
%   OUTPUT out 1x1 lsqvce.DayData.
    out = graceDay(day, ctx, cfg, size(cfg.Dep, 1), {'D', 1});
end

function out = buildRr(day, ctx, cfg)
%BUILDRR  The C-D range-rate on DAY (scalar): it couples BOTH satellites'
%   locals, weakly (partials scaled by 0.1).
%   INPUT day 1x1 double; ctx 1x1 lsqvce.DayContext; cfg 1x1 struct.
%   OUTPUT out 1x1 lsqvce.DayData.
    out = graceDay(day, ctx, cfg, 1, {'C', 0.1; 'D', 0.1});
end

function out = buildPosE(day, ctx, cfg)
%BUILDPOSE  Satellite E's 3-D positions on DAY: its own system, its own locals.
%   INPUT day 1x1 double; ctx 1x1 lsqvce.DayContext; cfg 1x1 struct.
%   OUTPUT out 1x1 lsqvce.DayData.
    out = graceDay(day, ctx, cfg, size(cfg.Dep, 1), {'E', 1});
end

function out = graceDay(day, ctx, cfg, nAx, sees)
%GRACEDAY  One synthetic stream-day, shared by the four Builds.
%   The design is random, the truth is cfg's, and the RandStream is seeded on
%   (day, stream), so every call returns the same rows -- Build must be
%   deterministic. The locals' partials are placed in the columns
%   ctx.cols(names) announces, in ONE dense B per day (as the v6 builder
%   held it; Local= would build it sparse).
%   INPUT   day   1x1 double            the day (an entry of cfg.days)
%           ctx   1x1 lsqvce.DayContext  Settings = cfg.procOpt (mEpoch, dt,
%                                       arcsPerDay, seed), NumCommon, NumLocal
%           cfg   1x1 struct            the planted truth (demoCfg)
%           nAx   1x1 double            axes per epoch (3 positions, 1 range-rate)
%           sees  k x 2 cell            {satellite, scale} rows: whose locals the
%                                       stream sees, in the order their partials
%                                       are drawn
%   OUTPUT  out   1x1 lsqvce.DayData     arcsPerDay arcs of mEpoch epochs, rows
%                                       epoch-major; a position arc carries its
%                                       per-epoch 3x3 covariance as Sigma=
    po   = ctx.Settings;
    ot   = char(ctx.Stream);            % posC | posD | rr | posE (seed, truth)
    nArc = po.arcsPerDay;
    mE   = po.mEpoch;                   % epochs per arc
    m    = nArc * mE * nAx;             % rows per day (epoch-major)

    s = RandStream('mt19937ar', 'Seed', ...
            mod(po.seed*1e5 + day*997 + sum(double(ot))*31, 2^31-1));

    A = A_design(s, m, ctx.NumCommon);  % the commons' columns, in addCommon order
    B = zeros(m, ctx.NumLocal);         % this system's local block
    y = zeros(ctx.NumLocal, 1);         % the planted locals this stream sees
    for k = 1:size(sees, 1)
        cols = ctx.cols(satelliteLocals(ctx, sees{k, 1}));
        B(:, cols) = sees{k, 2} * B_design(s, m, numel(cols), nArc, mE*nAx);
        y(cols)    = cfg.(['y_' sees{k, 1}]);
    end

    l0 = A * cfg.x_c + B * y;
    if nAx > 1
        sigT = cfg.trueSigma.(ot);
        Ld   = chol(cfg.Dep, 'lower');       % 3x3
        e    = zeros(m, 1);
        for a = 1:nArc
            base = (a-1) * mE * nAx;
            W = zeros(nAx, mE);
            for ax = 1:nAx
                w = randn(s, mE, 1);
                W(ax,:) = filter(1, cfg.trueAR.(ot), w).';  % time colour
            end
            W = sigT * (Ld * W);             % cross-axis mix + scale
            e(base + (1:mE*nAx)) = W(:);      % epoch-major [x1 y1 z1 ...]
        end
    else
        e = colouredNoise(s, m, cfg.trueAR.(ot), cfg.trueSigma.(ot), nAx, nArc);
    end
    l = l0 + e;
    t = reshape(repmat((0:mE*nArc-1)*po.dt, nAx, 1), [], 1);  % epoch-major

    out = ctx.newDay();
    for a = 1:nArc
        rows = (a-1)*mE*nAx + (1:mE*nAx);
        if nAx > 1
            % the ERROR INFORMATION of the arc: the per-epoch 3x3 covariance
            % (mE x 3 x 3), relative to the stream's law
            out.addArc(A(rows, :), l(rows), t(rows), B=B(rows, :), ...
                       Sigma=repmat(reshape(cfg.Dep, 1, nAx, nAx), mE, 1, 1));
        else
            out.addArc(A(rows, :), l(rows), t(rows), B=B(rows, :));
        end
    end
end

function names = satelliteLocals(ctx, sat)
%SATELLITELOCALS  This system's local parameters of satellite SAT, read from
%   the catalogue's satellite column (ctx.Local.Info), in declaration order.
%   INPUT   ctx    1x1 lsqvce.DayContext;  sat  1x1 char  'C' | 'D' | 'E'
%   OUTPUT  names  1 x n string
    L = ctx.Local;
    isSat = cellfun(@(info) strcmp(info.satellite, sat), L.Info);
    names = reshape(L.Name(isSat), 1, []);
end

%% ===================================================== the parameter catalogue
function C = commonParamTable(cfg)
%COMMONPARAMTABLE  nc rows x 15 cols: the scope token 'common' first, then
%   14 cols describing the shared gravity unknowns (column 3 the unique
%   name addCommon takes, the rest its Info). Synthetic degree/order
%   just enumerate the nc coefficients; a real table would carry the true SH
%   (n,m,C/S). Shared -> system 0, satellite ''.
    nc = cfg.nc;
    C  = cell(nc, 15);
    for i = 1:nc
        n = floor(sqrt(i-1));  m = (i-1) - n*n;      % toy (n,m) enumeration
        C(i,:) = { 'common', 'gravity', sprintf('Cnm_%d_%d', n, m), 0, '', ...
                   n, m, 'C', [], '', '-', cfg.x_c(i), Inf, true, 'gravity' };
    end
end

function T = localParamTable(cfg, sats, yTrue)
%LOCALPARAMTABLE  (nlSat*numel(sats)) rows x 15 cols: the scope token
%   'system' first (this demo hand-stacks the pair block [y_C | y_D], so its
%   rows stay hand-numbered -- exactly what 'system' scope means), then 14
%   cols (column 3 the unique name addLocal takes, column 5 the satellite the
%   Builds select by). yTrue is nlSat x numel(sats) (or nlSat x 1) of
%   a-priori truth.
    nlSat = cfg.nlSat;  ns = numel(sats);
    T = cell(nlSat*ns, 15);
    r = 0;
    for is = 1:ns
        sat = sats{is};
        for k = 1:nlSat
            r = r + 1;
            T(r,:) = { 'system', 'arcBias', sprintf('bias_%s_%d', sat, k), NaN, sat, ...
                       [], [], '', k, '', 'm', yTrue(k, is), Inf, true, ...
                       sprintf('loc_%s', sat) };
        end
    end
end

%% ===================================================== local helpers
function A = A_design(rs, m, nc)
    A = randn(rs, m, nc);
    A(:,1) = 1;                    % keep it well-conditioned
end

function B = B_design(rs, m, nl, nArc, mArc)
    B = zeros(m, nl);
    for a = 1:nArc
        r = (a-1)*mArc + (1:mArc);
        B(r,:) = 0.3 * randn(rs, mArc, nl);
    end
end

function e = colouredNoise(rs, m, ar, sig, nAx, nArc)
%COLOUREDNOISE  Synthetic noise matching the separable model: the temporal AR
%   law runs DOWN EPOCHS, independently per axis and per arc (rows are
%   epoch-major: [ax1 ax2 ax3] for epoch 1, then epoch 2, ...).
    if nargin < 5 || isempty(nAx),  nAx  = 1; end
    if nargin < 6 || isempty(nArc), nArc = 1; end
    nE = m / (nAx * nArc);                       % epochs per arc
    e  = zeros(m, 1);
    for a = 1:nArc
        base = (a-1) * nE * nAx;
        for ax = 1:nAx
            idx = base + (ax:nAx:(nE*nAx));       % this axis, down epochs
            w   = sig * randn(rs, nE, 1);
            if isscalar(ar) && ar == 1
                e(idx) = w;                       % white
            else
                e(idx) = filter(1, ar, w);        % AR(p) colour down epochs
            end
        end
    end
end

