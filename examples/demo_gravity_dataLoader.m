function [cfg, L] = demo_gravity_dataLoader()
%DEMO_GRAVITY_DATALOADER  GRACE-like gravity-recovery data for DEMO_GRAVITY.
%
%   [cfg, L] = demo_gravity_dataLoader()
%
%   Returns a configuration CFG and a cell array L of three Loaders, so a
%   PHYSICALLY STRUCTURED gravity recovery can be run end to end without any
%   external data:
%
%       [cfg, L] = demo_gravity_dataLoader();
%       R = lsqvce.adjustVCE(L, cfg.days).groupComponentsBy('obsType').run();
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
%        Points 2 and 3 are served by answering the accumulator's PER-DAY
%        'meta' probe with that day's actual arc list (see raggedClosure
%        below) -- the loader contract supports this directly, and it is why
%        the engines never see an arc that does not exist.
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
%   OUTPUT
%     cfg  1x1 struct  the configuration and the planted truth. Fields:
%            .nmax     1x1 double   maximum SH degree (8 -> 77 coefficients)
%            .nc       1x1 double   common-parameter count, numel(cfg.xTrue)
%            .meta     1x1 struct   .n .m (1xnc double) .cs (1xnc char), the
%                                   degree/order/term label of every column
%            .days     1xD double   the day list
%            .arcsOfDay 1xD cell    arc IDs present on each day (ragged)
%            .mOfArc   DxA double   epochs of each (day, arc); 0 where absent
%            .arcOrd   DxA double   global arc ordinal; 0 where absent
%            .nArcTot  1x1 double   total arcs over all days (16)
%            .sys      1x3 struct   .id .name .obsType .nlPerArc .nl
%            .xTrue    ncx1 double  the planted coefficients
%            .yTrue    1x1 struct   .potC .potD .kbr, the planted locals
%            .trueSigma/.trueAR 1x1 struct  the planted noise law per stream
%            .cov      1x1 struct   the DECLARED covariance per stream
%            .r .incl .dt .lagC2D   orbit geometry
%     L    1x3 cell    the Loaders, ready for LSQVCE.ADJUSTVCE
%
%   ERRORS
%     None. The fixture is self-contained and deterministic.
%
%   EXAMPLE
%     [cfg, L] = demo_gravity_dataLoader();
%     R = lsqvce.adjustVCE(L, cfg.days, 'MaxOuter', 4, 'Verbose', false) ...
%             .groupComponentsBy('obsType').run();
%     max(abs(R.x - cfg.xTrue) ./ (sqrt(diag(R.sol.Cxx))*R.sigma0))   % < 5
%
%   See also DEMO_GRAVITY, DEMO_GRAVITY_DESIGN, DEMO_GRACE_DATALOADER,
%   LSQVCE.MAKELOADER
%
%   Created 2026-08-21 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.

    cfg = gravityCfg();
    L   = gravityLoaders(cfg);
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
function L = gravityLoaders(cfg)
%GRAVITYLOADERS  Build one Loader per system, each with a per-day arc list.
%   INPUT  cfg 1x1 struct   the configuration
%   OUTPUT L   1x3 cell     the Loaders
    L = cell(1, numel(cfg.sys));
    for s = 1:numel(cfg.sys)
        ot  = cfg.sys(s).obsType;
        cov = struct(ot, cfg.cov.(ot));
        L0 = lsqvce.makeLoader(@(d, o) buildGravityDay(cfg, d, o), ...
                'observationTypes', {ot}, ...
                'covariance',       cov, ...
                'arcs',             cfg.maxArc, ...
                'systemId',         cfg.sys(s).id, ...
                'name',             cfg.sys(s).name, ...
                'numLocal',         cfg.sys(s).nl, ...
                'numCommon',        cfg.nc, ...
                'days',             cfg.days);
        % THE RAGGED-SCHEDULE HOOK. makeLoader resolves 'arcs' once, without
        % reading a day, so on its own it would declare four arcs on every
        % day. The accumulator, however, probes 'meta' PER DAY -- so we wrap
        % the closure and answer that probe with the day's ACTUAL arc list.
        % The engines then never ask for an arc that does not exist.
        L0.closure = @(day, obsType, arcID, opt) ...
                        raggedClosure(L0, cfg, day, obsType, arcID, opt);
        L{s} = L0;
    end
end

function out = raggedClosure(L0, cfg, day, obsType, arcID, opt)
%RAGGEDCLOSURE  Answer the per-day 'meta' probe with that day's arcs; delegate
%   everything else to the Loader built by LSQVCE.MAKELOADER.
%   INPUT   L0 1x1 struct the wrapped Loader; cfg; then the 4-arg callback
%           signature (day, obsType, arcID, opt), where day == 'meta' is the
%           discovery probe and obsType then carries the day number.
%   OUTPUT  out 1x1 struct  the meta struct, or one arc's batch.
    if (ischar(day) || isstring(day)) && strcmp(day, 'meta')
        d = obsType;                                 % the probe carries the day
        out.obsTypes = L0.obsTypes;
        out.arcs     = repmat({cfg.arcsOfDay{cfg.days == d}}, ...
                              1, numel(L0.obsTypes));
        return
    end
    out = L0.batch(day, obsType, arcID, opt);
end

%==========================================================================
function D = buildGravityDay(cfg, day, obsType)
%BUILDGRAVITYDAY  One satellite-day: physical design, planted truth, noise.
%   INPUT   cfg 1x1 struct; day 1x1 double; obsType 1xK char
%   OUTPUT  D 1x1 struct  .A (mxnc) .B (mxnl) .l (mx1) .t (mx1)
%                         .arcIdx {1xmaxArc} cell, EMPTY for absent arcs
%                         (never requested -- meta does not list them)
    di   = find(cfg.days == day);
    arcs = cfg.arcsOfDay{di};
    s    = find(strcmp({cfg.sys.obsType}, obsType));
    nlPA = cfg.sys(s).nlPerArc;
    nl   = cfg.sys(s).nl;

    mTot = sum(cfg.mOfArc(di, arcs));
    D.A  = zeros(mTot, cfg.nc);
    D.B  = zeros(mTot, nl);
    D.l  = zeros(mTot, 1);
    D.t  = zeros(mTot, 1);
    D.arcIdx = cell(1, cfg.maxArc);

    rs = RandStream('mt19937ar', 'Seed', cfg.seed.(obsType) + 97*day);
    r0 = 0;
    tOff = 0;
    for a = arcs
        m    = cfg.mOfArc(di, a);
        rows = r0 + (1:m);
        ord  = cfg.arcOrd(di, a);              % global arc ordinal

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
        cols = (ord-1)*nlPA + (1:nlPA);

        % ---- the planted observation ------------------------------------
        w = cfg.trueSigma.(obsType) * randn(rs, m, 1);
        e = filter(1, cfg.trueAR.(obsType), w);

        D.A(rows, :)    = A;
        D.B(rows, cols) = Barc;
        D.l(rows)       = A * cfg.xTrue + Barc * cfg.yTrue.(obsType)(cols) + e;
        D.t(rows)       = tOff + (0:m-1).' * cfg.dt;
        D.arcIdx{a}     = rows;

        r0   = r0 + m;
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
