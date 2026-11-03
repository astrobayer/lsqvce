function [cfg, task] = demo_gravity_dataLoader()
%DEMO_GRAVITY_DATALOADER  GRACE-like gravity-recovery data for DEMO_GRAVITY / DEMO_FIT, as a lsqvce.TaskDefinition.
%
%   [cfg, task] = demo_gravity_dataLoader()
%   cfg = demo_gravity_dataLoader()
%
%   Returns a configuration CFG and the declaration TASK (what LSQVCE.FIT
%   takes), so a PHYSICALLY STRUCTURED gravity recovery can be run end to end
%   without any external data:
%
%       [cfg, task] = demo_gravity_dataLoader();
%       R = lsqvce.fit(task, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false));
%
%   The ENGINE tier (fixture_gravity, the tests) takes the compiled struct,
%   task.declaration(); lsqvce.fit(task) fits exactly it.
%
%   INPUT
%     (none)
%   OUTPUT
%     cfg   1x1 struct  the configuration (days, nc, sys, cov, arcsOfDay,
%                       truth) -- CFG IN DETAIL below
%     task  1x1 lsqvce.TaskDefinition   77 common parameters named "C2_0",
%                       "S2_1", ... with Info degree / order / term; three
%                       systems "GRACE-C", "GRACE-D", "KBR" (Id 1, 2, 3) on a
%                       ragged arc schedule, each with its own per-arc locals
%                       and ONE stream (potC, potD, kbr)
%
%   WHAT MAKES THIS DIFFERENT FROM DEMO_GRACE_DATALOADER
%     That loader builds its design from randn -- fine for exercising the
%     algebra, but a Gaussian design is BENIGN: near-orthogonal columns,
%     well conditioned. Here every design column is a real spherical-harmonic
%     potential partial evaluated along a sampled orbit
%     (DEMO_GRAVITY_DESIGN), so the normal matrix carries what production
%     carries: order-block correlation, polar-gap conditioning, (1/r)^n
%     degree attenuation. Nothing in the design is random.
%
%   TOPOLOGY (GRACE-like: two satellites, three streams)
%     System 1 'GRACE-C'  stream 'potC'  potential along GRACE-C's orbit
%     System 2 'GRACE-D'  stream 'potD'  potential along GRACE-D's orbit, which
%                                        trails C along-track by a fixed lag
%     System 3 'KBR'      stream 'kbr'   the inter-satellite POTENTIAL
%                                        DIFFERENCE, A_kbr = A_D - A_C. This is
%                                        the energy-balance reading of the KBR
%                                        observable and it is what makes the
%                                        demo GRACE-like: one stream whose
%                                        design couples both orbits and whose
%                                        noise is far quieter than the two
%                                        absolute streams.
%
%   THE THREE THINGS THIS FIXTURE EXERCISES THAT THE OTHERS DO NOT
%     1. A DIFFERENT NUMBER OF LOCAL UNKNOWNS PER SATELLITE. Locals are
%        per-arc empirical parameters and each system declares its own:
%        C carries 2 per arc ("bias", "drift"), D carries 3 ("bias", "drift",
%        "quadratic"), the KBR stream carries 1 (a range "bias"). nl is
%        therefore 32 / 48 / 16 -- the block-arrow solver must handle three
%        different local block sizes in one adjustment.
%     2. A DIFFERENT NUMBER OF ARCS PER DAY. The schedule is ragged (3, 4, 2,
%        4, 3 arcs), and the missing arcs are missing BY ID, not by count:
%        day 3 has arcs 1 and 3, day 5 has arcs 2, 3 and 4. This is what a
%        real arc schedule looks like once manoeuvres and data gaps have
%        taken their share.
%     3. A DIFFERENT ARC LENGTH FOR EVERY ARC (140 to 220 epochs). Nothing in
%        the toolbox may assume a fixed arc length.
%        Points 1 to 3 are DECLARED, not hand-wired (WP-52): each system's
%        Arcs (addSystem) is a per-day resolver of the arc ids -- required
%        here, because Scope="arc" locals (addLocal) make nl a function of
%        the schedule -- and the toolbox does the local bookkeeping: it
%        numbers every (day, arc), DERIVES nl (ctx.NumLocal), and announces
%        the columns of an arc's locals by name, ctx.cols(names, Arc=id).
%        Each Build adds one arc per scheduled id (ctx.ArcIds) and never
%        computes an offset of its own.
%
%   PLANTED TRUTH (deterministic -- the loop must recover it)
%     cfg.xTrue    the spherical-harmonic coefficients, a closed formula plus
%                  ONE planted anomaly at (n,m) = (5,3), cosine term
%     cfg.yTrue    the per-arc empirical parameters of each system
%     cfg.trueSigma / cfg.trueAR   the noise law actually used to generate the
%                  observations. The NOISE DECLARED to the toolbox (cfg.cov,
%                  stated in the task as lsqvce.Noise.arma / .white) has the
%                  right AR shape but unit sigma for every stream, so the
%                  variance components VCE estimates should come out as the
%                  planted sigmas: 1.0 / 2.0 / 0.5.
%
%   CFG IN DETAIL
%     cfg  1x1 struct  the configuration and the planted truth. Fields:
%            .nmax     1x1 double   maximum SH degree (8 -> 77 coefficients)
%            .nc       1x1 double   common-parameter count, numel(cfg.xTrue)
%            .meta     1x1 struct   .n .m (1xnc double) .cs (1xnc char), the
%                                   degree/order/term label of every column
%            .days     1xD double   the day list
%            .arcsOfDay 1xD cell    arc IDs present on each day (ragged)
%            .mOfArc   DxA double   epochs of each (day, arc); 0 where absent
%            .arcOrd   DxA double   global arc ordinal; 0 where absent. KEPT
%                                   only so the tests can check makeLoader's
%                                   plan against an independent table -- the
%                                   Builds never read it (WP-52)
%            .nArcTot  1x1 double   total arcs over all days (16)
%            .sys      1x3 struct   .id .name .obsType .nlPerArc .nl
%            .xTrue    ncx1 double  the planted coefficients
%            .yTrue    1x1 struct   .potC .potD .kbr, the planted locals
%            .trueSigma/.trueAR 1x1 struct  the planted noise law per stream
%            .cov      1x1 struct   the DECLARED noise per stream (law and
%                                   parameters: .type .ar .ma .sigma), from
%                                   which gravityTask states its lsqvce.Noise
%            .r .incl .dt .lagC2D   orbit geometry
%
%   ACCURACY
%     Not applicable: a deterministic generator (the noise drawn from one
%     RandStream per day and stream, seeded), no floating-point result of its
%     own; DEMO_GRAVITY_DESIGN carries the design's ACCURACY. The WP-161
%     migration keeps every seed, the order of every draw and the DENSE local
%     design of the v6 builder, so the rows are those of the v6 declaration by
%     construction (not measured in the migrating session, which had no MATLAB).
%
%   ERRORS
%     None of its own: the fixture is self-contained and deterministic. A
%     declaration mistake raises lsqvce:TaskDefinition:* at the add-call that
%     makes it; a mistake inside a Build raises lsqvce:DayContext:* while the
%     day is read.
%
%   EXAMPLE
%     [cfg, task] = demo_gravity_dataLoader();
%     R = lsqvce.fit(task, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, MaxPasses=4, ...
%             Verbose=false));
%     max(abs(R.x - cfg.xTrue) ./ sqrt(diag(R.sol.Cxx)))   % < 5 (sol.Cxx is sigma0^2-scaled)
%
%   See also DEMO_GRAVITY, DEMO_GRAVITY_DESIGN, DEMO_GRACE_DATALOADER,
%   LSQVCE.TASKDEFINITION, LSQVCE.SYSTEMDEFINITION, LSQVCE.NOISE, LSQVCE.DAYCONTEXT,
%   LSQVCE.FIT
%
%   Updated 2026-11-03 (CET, UTC+01:00), WP-162: the outputs are [cfg, task];
%   the compiled struct is task.declaration() wherever the engine tier needs
%   it. With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-11-02 (CET, UTC+01:00), WP-161: declared as a lsqvce.TaskDefinition, one Build per stream. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-09-09 (CEST, UTC+02:00), R6 WP-102 (v6): buildGravityDay
%   returns one cell per scheduled arc (the cells are the arcs; arcIdx /
%   arcId are gone). With assistance from Claude (Anthropic), model Claude
%   Fable 5.1.
%   Updated 2026-09-07 (CEST, UTC+02:00), R6 WP-95: the contract is named v5,
%   and the stale third output L "ready for LSQVCE.ENGINE.ADJUSTVCE" is gone --
%   the signature is [cfg, data] and that engine unit no longer exists; the
%   remaining field list is retitled CFG IN DETAIL so it no longer poses as a
%   second OUTPUT block. With assistance from Claude (Anthropic), model
%   Claude Opus 5.
%   Updated 2026-08-31 (CEST, UTC+02:00), v4 SCOPED tables (Level L11,
%   WP-77): the declaration tables carry the scope column (commonParams
%   'common'|'day', localParams 'arc'|'day'|'system'); localsPerArc left the
%   USER surface. With assistance from Claude (Anthropic), model
%   Claude Fable 5.
%   Updated 2026-08-31 (CEST, UTC+02:00), v3 data contract (Level L10,
%   WP-74): USER surface only -- returns [cfg, data], the v3 declaration
%   (commonParams = the degree/order/term table, obsTypes from the
%   covariance fieldnames, systemId -> id); the ENGINE-tier loader
%   construction moved to dev/tests/fixture_gravity.m. With assistance
%   from Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-31 (CEST, UTC+02:00), R4 Level L7 (WP-67): third output
%   DATA -- the same declaration as a 2.0 data struct for LSQVCE.FIT. With
%   assistance from Claude (Anthropic), model Claude Fable 5.
%   Created 2026-08-21 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.
%
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-124: references to the data
%   contract point at the API reference. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-29 (CET, UTC+01:00), WP-159: the declaration is the data
%   struct of the API reference (it was called "v5"; since v6 a day is one
%   cell per arc). With assistance from Claude (Anthropic), model Claude
%   Opus 5.5.

    cfg  = gravityCfg();
    task = gravityTask(cfg);
end

%==========================================================================
function task = gravityTask(cfg)
%GRAVITYTASK  The declaration (WP-161): the names declare, the counts are
%   derived. One common parameter per SH coefficient, named by term, degree
%   and order ("C2_0", "S2_1", ...) and carrying them as Info; three systems
%   on the ragged schedule, each with its own per-arc locals and ONE stream
%   whose Build reads the day.
%   INPUT   cfg   1x1 struct  see gravityCfg
%   OUTPUT  task  1x1 lsqvce.TaskDefinition
    task  = lsqvce.TaskDefinition(Days=cfg.days);
    names = arrayfun(@(k) sprintf('%s%d_%d', cfg.meta.cs(k), cfg.meta.n(k), cfg.meta.m(k)), ...
                     1:cfg.nc, 'UniformOutput', false);
    info  = table(cfg.meta.n(:), cfg.meta.m(:), cellstr(cfg.meta.cs(:)), ...
                  'VariableNames', {'degree', 'order', 'term'});
    task.addCommon(string(names), Info=info);           % the columns of A, in meta order

    % the declared noise: the right AR shape, unit sigma (cfg.cov)
    noise = struct('potC', lsqvce.Noise.arma(cfg.cov.potC.ar, cfg.cov.potC.sigma, MA=cfg.cov.potC.ma), ...
                   'potD', lsqvce.Noise.arma(cfg.cov.potD.ar, cfg.cov.potD.sigma, MA=cfg.cov.potD.ma), ...
                   'kbr',  lsqvce.Noise.white(cfg.cov.kbr.sigma));
    build = struct('potC', @(day, ctx) buildPotC(day, ctx, cfg), ...
                   'potD', @(day, ctx) buildPotD(day, ctx, cfg), ...
                   'kbr',  @(day, ctx) buildKbr(day, ctx, cfg));
    perArc = ["bias" "drift" "quadratic"];             % the per-arc empirical parameters
    for s = 1:numel(cfg.sys)
        q   = cfg.sys(s);
        sys = task.addSystem(string(q.name), Id=q.id, ...
                             Arcs=@(day, ~) cfg.arcsOfDay{cfg.days == day});   % the ragged schedule
        sys.addLocal(perArc(1:q.nlPerArc), Scope="arc");      % C 2, D 3, KBR 1 per arc
        sys.addStream(string(q.obsType), Build=build.(q.obsType), Noise=noise.(q.obsType));
    end
end

%==========================================================================
function cfg = gravityCfg()
%GRAVITYCFG  The deterministic configuration and planted truth.
%   OUTPUT cfg 1x1 struct  see the header of DEMO_GRAVITY_DATALOADER.
    cfg.nmax = 8;
    cfg.days = 1:5;
    cfg.r    = 1.07;                       % ~450 km altitude, Earth radii
    cfg.incl = deg2rad(89);                % near-polar: a real polar gap
    cfg.dt   = 5;                          % seconds between epochs
    cfg.lagC2D = 0.0035;                   % along-track lag C -> D [rad],
                                           % ~220 km at this radius (GRACE-like)

    % ---- the ragged arc schedule ---------------------------------------
    % Arc IDs present on each day; note days 3 and 5 are missing arcs by ID,
    % not merely short of them.
    cfg.arcsOfDay = {[1 2 3], [1 2 3 4], [1 3], [1 2 3 4], [2 3 4]};
    cfg.maxArc    = 4;
    % Epochs of each (day, arc): every arc has its own length.
    cfg.mOfArc = [180 150 200   0; ...
                  160 210 140 190; ...
                  220   0 170   0; ...
                  150 160 200 180; ...
                    0 190 210 150];
    % Global arc ordinal, day-then-arc: the numbering the toolbox gives the
    % local-parameter blocks itself (the Builds never read it; the tests
    % check the loader's plan against it).
    cfg.arcOrd = zeros(numel(cfg.days), cfg.maxArc);
    k = 0;
    for d = 1:numel(cfg.days)
        for a = cfg.arcsOfDay{d}
            k = k + 1;
            cfg.arcOrd(d, a) = k;
        end
    end
    cfg.nArcTot = k;

    % ---- the three systems, each with its OWN local count --------------
    cfg.sys = struct( ...
        'id',       {1, 2, 3}, ...
        'name',     {'GRACE-C', 'GRACE-D', 'KBR'}, ...
        'obsType',  {'potC', 'potD', 'kbr'}, ...
        'nlPerArc', {2, 3, 1});
    for s = 1:numel(cfg.sys)
        cfg.sys(s).nl = cfg.sys(s).nlPerArc * cfg.nArcTot;
    end

    % ---- the planted coefficient field ---------------------------------
    probe = struct('u', 0, 'lam0', 0, 'incl', cfg.incl, 'r', cfg.r);
    [~, meta] = demo_gravity_design(cfg.nmax, probe);
    cfg.meta = meta;
    cfg.nc   = numel(meta.n);
    xT = (cos(3*meta.n + meta.m) ./ meta.n.^2).';
    isS = meta.cs == 'S';
    xT(isS) = (sin(2*meta.n(isS) - meta.m(isS)) ./ meta.n(isS).^2).';
    anom = meta.n == 5 & meta.m == 3 & meta.cs == 'C';
    xT(anom) = xT(anom) + 0.15;            % THE planted anomaly
    cfg.xTrue = xT;

    % ---- the planted per-arc empirical parameters ----------------------
    cfg.yTrue = struct();
    for s = 1:numel(cfg.sys)
        ot = cfg.sys(s).obsType;
        nl = cfg.sys(s).nl;
        cfg.yTrue.(ot) = 0.02 * cos(0.7 * (1:nl).' + s);
    end

    % ---- the planted noise law, and what we DECLARE to the toolbox -----
    % Declared: the right AR shape, unit sigma. Planted: the same shape,
    % these sigmas. So the variance components should come out as 1 / 2 / 0.5.
    cfg.trueSigma = struct('potC', 1.0,  'potD', 2.0,  'kbr', 0.5);
    cfg.trueAR    = struct('potC', [1 -0.6 0.2], 'potD', [1 -0.3], ...
                           'kbr',  1);
    cfg.cov = struct( ...
        'potC', struct('type','arma', 'ar',[1 -0.6 0.2], 'ma',1, 'sigma',1), ...
        'potD', struct('type','arma', 'ar',[1 -0.3],     'ma',1, 'sigma',1), ...
        'kbr',  struct('type','diag', 'sigma',1));

    % Seeds: the noise REALISATION is random but reproducible; every other
    % number in this fixture is deterministic.
    cfg.seed = struct('potC', 510001, 'potD', 510002, 'kbr', 510003);
end

%==========================================================================
function out = buildPotC(day, ctx, cfg)
%BUILDPOTC  GRACE-C's potential on DAY: the SH design along C's orbit, and
%   per arc a bias and a drift.
%   INPUT day 1x1 double; ctx 1x1 lsqvce.DayContext; cfg 1x1 struct.
%   OUTPUT out 1x1 lsqvce.DayData.
    design   = @(a, m) demo_gravity_design(cfg.nmax, arcOrbit(cfg, day, a, m, 0));
    partials = @(tau) [ones(numel(tau), 1), tau];                   % bias, drift
    out = gravityArcs(day, ctx, cfg, design, partials);
end

function out = buildPotD(day, ctx, cfg)
%BUILDPOTD  GRACE-D's potential on DAY: the SH design along D's orbit (C's,
%   trailing by cfg.lagC2D), and per arc a bias, a drift and a quadratic.
%   INPUT day 1x1 double; ctx 1x1 lsqvce.DayContext; cfg 1x1 struct.
%   OUTPUT out 1x1 lsqvce.DayData.
    design   = @(a, m) demo_gravity_design(cfg.nmax, arcOrbit(cfg, day, a, m, cfg.lagC2D));
    partials = @(tau) [ones(numel(tau), 1), tau, tau.^2];           % + quadratic
    out = gravityArcs(day, ctx, cfg, design, partials);
end

function out = buildKbr(day, ctx, cfg)
%BUILDKBR  The KBR stream on DAY, read by energy balance: the potential
%   DIFFERENCE along the baseline (D's design minus C's), so its design
%   couples both orbits; per arc one range bias.
%   INPUT day 1x1 double; ctx 1x1 lsqvce.DayContext; cfg 1x1 struct.
%   OUTPUT out 1x1 lsqvce.DayData.
    design   = @(a, m) demo_gravity_design(cfg.nmax, arcOrbit(cfg, day, a, m, cfg.lagC2D)) ...
                       - demo_gravity_design(cfg.nmax, arcOrbit(cfg, day, a, m, 0));
    partials = @(tau) ones(numel(tau), 1);                           % range bias
    out = gravityArcs(day, ctx, cfg, design, partials);
end

function out = gravityArcs(day, ctx, cfg, design, partials)
%GRAVITYARCS  One satellite-day, shared by the three Builds: physical design,
%   planted truth, noise -- one arc per scheduled id of the day.
%   INPUT   day       1x1 double            the day (an entry of cfg.days)
%           ctx       1x1 lsqvce.DayContext  ArcIds (the day's schedule), the
%                                           columns by name (ctx.cols), NumLocal
%           cfg       1x1 struct            the configuration and planted truth
%           design    function_handle       A = design(arc, m), m x nc
%           partials  function_handle       P = partials(tau), m x nlPerArc,
%                                           tau the normalised arc time
%   OUTPUT  out       1x1 lsqvce.DayData     cell ia is arc ctx.ArcIds(ia):
%                                           A m x nc, B m x nl (dense), l, t m x 1
%   The local partials of an arc go into ONE dense B in the columns
%   ctx.cols(names, Arc=id) announces (as the v6 builder held it; Local= would
%   build it sparse); l uses the same columns to pick the planted locals.
    ot    = char(ctx.Stream);                % potC | potD | kbr (seed, truth)
    di    = find(cfg.days == day);
    names = reshape(ctx.Local.Name, 1, []);  % this system's per-arc locals, in order
    rs    = RandStream('mt19937ar', 'Seed', cfg.seed.(ot) + 97*day);
    out   = ctx.newDay();
    tOff  = 0;
    for a = ctx.ArcIds                       % == cfg.arcsOfDay{di}, the resolver's list
        m = cfg.mOfArc(di, a);

        % ---- the physical design and this arc's empirical parameters ----
        A    = design(a, m);
        tau  = (0:m-1).' / m;                % normalised arc time
        Barc = partials(tau);
        cols = ctx.cols(names, Arc=a);       % the arc's local columns, by name

        % ---- the planted observation ------------------------------------
        w = cfg.trueSigma.(ot) * randn(rs, m, 1);
        e = filter(1, cfg.trueAR.(ot), w);

        B = zeros(m, ctx.NumLocal);  B(:, cols) = Barc;
        l = A * cfg.xTrue + Barc * cfg.yTrue.(ot)(cols) + e;
        t = tOff + (0:m-1).' * cfg.dt;
        out.addArc(A, l, t, B=B, Arc=a);

        tOff = tOff + m * cfg.dt + 600;      % a gap between arcs
    end
end

function orbit = arcOrbit(cfg, day, arc, m, lag)
%ARCORBIT  Deterministic orbit samples for one arc of one satellite.
%   Circular, inclined, with an Earth-rotating node; the day and the arc ID
%   shift the ground track so that no two arcs repeat, and LAG offsets the
%   trailing satellite along-track.
%   INPUT   cfg; day 1x1; arc 1x1; m 1x1 epochs; lag 1x1 along-track [rad]
%   OUTPUT  orbit 1x1 struct  .u .lam0 (mx1) .incl .r  as DEMO_GRAVITY_DESIGN
%                             requires
    k     = (0:m-1).';
    du    = 2*pi * 1.15 / 160;                 % ~1.15 revolutions per 160 epochs
    phase = 0.7*day + 1.31*arc + lag;
    orbit = struct('u',    phase + du*k, ...
                   'lam0', 0.4*day + 0.17*arc - (du/15.3)*k, ...
                   'incl', cfg.incl, 'r', cfg.r);
end
