function [cfg, data] = demo_gravity_dataLoader()
%DEMO_GRAVITY_DATALOADER  GRACE-like gravity-recovery data for DEMO_GRAVITY / DEMO_FIT.
%
%   [cfg, data] = demo_gravity_dataLoader()
%
%   Returns a configuration CFG and the declaration as a v5 data struct DATA
%   (the form LSQVCE.FIT takes -- see the API reference), so a PHYSICALLY
%   STRUCTURED gravity recovery can be run end to end without any external
%   data:
%
%       [cfg, data] = demo_gravity_dataLoader();
%       R = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false));
%
%   OUTPUT
%     cfg   1x1 struct  the configuration (days, nc, sys, cov, arcsOfDay, truth)
%     data  1x1 struct  .commonParams (one row per SH coefficient: degree,
%                       order, term behind the 'common' token -- nc derived), .days, .systems(1x3) with
%                       buildDay, covariance (fieldname = the stream),
%                       arcs (per-day resolver), localParams ('arc' rows), id, name
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
%     System 1 'C'  obsType 'potC'  potential along GRACE-C's orbit
%     System 2 'D'  obsType 'potD'  potential along GRACE-D's orbit, which
%                                   trails C along-track by a fixed lag
%     System 3 'K'  obsType 'kbr'   the inter-satellite POTENTIAL DIFFERENCE,
%                                   A_kbr = A_D - A_C. This is the energy-
%                                   balance reading of the KBR observable and
%                                   it is what makes the demo GRACE-like: one
%                                   stream whose design couples both orbits
%                                   and whose noise is far quieter than the
%                                   two absolute streams.
%
%   THE THREE THINGS THIS FIXTURE EXERCISES THAT THE OTHERS DO NOT
%     1. A DIFFERENT NUMBER OF LOCAL UNKNOWNS PER SATELLITE. Locals are
%        per-arc empirical parameters and each system declares its own count:
%        C carries 2 per arc (bias, drift), D carries 3 (bias, drift,
%        quadratic), the KBR stream carries 1 (a range bias). nl is therefore
%        32 / 48 / 16 -- the block-arrow solver must handle three different
%        local block sizes in one adjustment.
%     2. A DIFFERENT NUMBER OF ARCS PER DAY. The schedule is ragged (3, 4, 2,
%        4, 3 arcs), and the missing arcs are missing BY ID, not by count:
%        day 3 has arcs 1 and 3, day 5 has arcs 2, 3 and 4. This is what a
%        real arc schedule looks like once manoeuvres and data gaps have
%        taken their share.
%     3. A DIFFERENT ARC LENGTH FOR EVERY ARC (140 to 220 epochs). Nothing in
%        the toolbox may assume a fixed arc length.
%        Points 1 to 3 are DECLARED, not hand-wired (WP-52): 'arcs' takes a
%        per-day resolver (required here: the 'arc' template rows make nl a
%        function of the schedule), the day returns ONE CELL PER SCHEDULED
%        ARC (v6, R6 WP-102: the cells are the arcs), and
%        the 'arc' template rows hand the toolbox the local bookkeeping -- it
%        numbers every (day, arc), DERIVES nl, and passes this day's column
%        offsets to the builder. An earlier version of this file wrapped the
%        loader's closure and carried its own arc-ordinal table; both are
%        gone, which is the honest proof that the support is native.
%
%   PLANTED TRUTH (deterministic -- the loop must recover it)
%     cfg.xTrue    the spherical-harmonic coefficients, a closed formula plus
%                  ONE planted anomaly at (n,m) = (5,3), cosine term
%     cfg.yTrue    the per-arc empirical parameters of each system
%     cfg.trueSigma / cfg.trueAR   the noise law actually used to generate the
%                  observations. The COVARIANCE DECLARED to the toolbox has
%                  the right AR shape but unit sigma for every stream, so the
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
%                                   builder no longer reads it (WP-52)
%            .nArcTot  1x1 double   total arcs over all days (16)
%            .sys      1x3 struct   .id .name .obsType .nlPerArc .nl
%            .xTrue    ncx1 double  the planted coefficients
%            .yTrue    1x1 struct   .potC .potD .kbr, the planted locals
%            .trueSigma/.trueAR 1x1 struct  the planted noise law per stream
%            .cov      1x1 struct   the DECLARED covariance per stream
%            .r .incl .dt .lagC2D   orbit geometry
%
%   ERRORS
%     None. The fixture is self-contained and deterministic.
%
%   EXAMPLE
%     [cfg, data] = demo_gravity_dataLoader();
%     R = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, MaxPasses=4, ...
%             Verbose=false));
%     max(abs(R.x - cfg.xTrue) ./ (sqrt(diag(R.sol.Cxx))*R.sigma0))   % < 5
%
%   See also DEMO_GRAVITY, DEMO_GRAVITY_DESIGN, DEMO_GRACE_DATALOADER,
%   LSQVCE.FIT
%
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

    cfg  = gravityCfg();
    data = gravityData(cfg);
end

%==========================================================================
function data = gravityData(cfg)
%GRAVITYDATA  The declaration as a v5 data struct (L11 WP-77): the tables
%   declare, the counts are derived -- commonParams is the degree/order/term
%   table of the SH coefficients (nc = rows, and the columns MEAN something),
%   the observation type is the covariance fieldname, id is the system id.
%   fit builds the Loaders from this internally.
    data = struct();
    data.commonParams = [repmat({'common'}, numel(cfg.meta.n), 1), ...
                         num2cell(cfg.meta.n(:)), num2cell(cfg.meta.m(:)), ...
                         cellstr(cfg.meta.cs(:))];        % scoped: nc x 4
    data.days = cfg.days;
    for s = 1:numel(cfg.sys)
        ot = cfg.sys(s).obsType;
        d  = struct();
        d.buildDay     = @(dd, o, cm, lo, po, ai) buildGravityDay(cfg, dd, o, ai);
        d.covariance   = struct(ot, cfg.cov.(ot));        % fieldname = the stream
        d.arcs         = @(dd, o) cfg.arcsOfDay{cfg.days == dd};
        d.localParams  = [repmat({'arc'}, cfg.sys(s).nlPerArc, 1), ...
                          num2cell((1:cfg.sys(s).nlPerArc).')];   % v5 'arc' templates
        d.id           = cfg.sys(s).id;
        d.name         = cfg.sys(s).name;
        if s == 1, data.systems = d; else, data.systems(s) = d; end
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
    % Global arc ordinal, day-then-arc. It indexes the local-parameter
    % blocks, so it must be deterministic and shared by all three systems.
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
function D = buildGravityDay(cfg, day, obsType, arcInfo)
%BUILDGRAVITYDAY  One satellite-day: physical design, planted truth, noise.
%   INPUT   cfg 1x1 struct; day 1x1 double; obsType 1xK char;
%           arcInfo 1x1 struct  the day's local-block plan from makeLoader
%                               (WP-52): .arcId .offset .count .nl .offsetOf
%   OUTPUT  D 1x1 struct  .A .B .l .t  {1xnArc} cells, ONE CELL PER ARC of
%                         the day's schedule (v6): cell ia is arc arcs(ia),
%                         A m x nc, B m x nl, l and t m x 1
    di   = find(cfg.days == day);
    arcs = cfg.arcsOfDay{di};              % == the resolver's list for this day
    s    = find(strcmp({cfg.sys.obsType}, obsType));
    nlPA = cfg.sys(s).nlPerArc;
    nl   = arcInfo.nl;                     % DERIVED by makeLoader, not restated

    nArc = numel(arcs);
    D = struct('A', {cell(1, nArc)}, 'B', {cell(1, nArc)}, 'l', {cell(1, nArc)}, 't', {cell(1, nArc)});

    rs = RandStream('mt19937ar', 'Seed', cfg.seed.(obsType) + 97*day);
    tOff = 0;
    for ia = 1:nArc
        a    = arcs(ia);
        m    = cfg.mOfArc(di, a);

        % ---- the physical design ---------------------------------------
        Ac = demo_gravity_design(cfg.nmax, arcOrbit(cfg, day, a, m, 0));
        switch obsType
            case 'potC'
                A = Ac;
            case 'potD'
                A = demo_gravity_design(cfg.nmax, ...
                        arcOrbit(cfg, day, a, m, cfg.lagC2D));
            case 'kbr'
                % energy-balance reading: the potential DIFFERENCE along the
                % baseline, so this design couples both orbits
                A = demo_gravity_design(cfg.nmax, ...
                        arcOrbit(cfg, day, a, m, cfg.lagC2D)) - Ac;
        end

        % ---- this arc's empirical parameters ---------------------------
        tau = (0:m-1).' / m;                   % normalised arc time
        switch nlPA
            case 1, Barc = ones(m, 1);                       % range bias
            case 2, Barc = [ones(m,1), tau];                 % bias, drift
            case 3, Barc = [ones(m,1), tau, tau.^2];         % + quadratic
        end
        cols = arcInfo.offsetOf(a) + (1:nlPA);   % offsets owned by makeLoader

        % ---- the planted observation ------------------------------------
        w = cfg.trueSigma.(obsType) * randn(rs, m, 1);
        e = filter(1, cfg.trueAR.(obsType), w);

        B = zeros(m, nl);  B(:, cols) = Barc;
        D.A{ia} = A;
        D.B{ia} = B;
        D.l{ia} = A * cfg.xTrue + Barc * cfg.yTrue.(obsType)(cols) + e;
        D.t{ia} = tOff + (0:m-1).' * cfg.dt;

        tOff = tOff + m * cfg.dt + 600;        % a gap between arcs
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
