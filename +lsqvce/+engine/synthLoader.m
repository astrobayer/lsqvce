function out = synthLoader(cfg, sat, day, obsType, arcID, opt)
%SYNTHLOADER  Synthetic ground-truth loader -- and the loader TEMPLATE.
%
%   Bind it as a closure over cfg and the satellite, so the accumulator sees
%   the four-argument signature it expects:
%
%       loader = @(d,ot,a,uo) synthLoader(cfg, satID, d, ot, a, uo);
%       acc.addSystem(satID, nl, loader, days);
%
%   COPY THIS FILE as the starting point for a real loader. The structure --
%   meta call, single-day cache, arc slicing, covariance closure -- is exactly
%   what a GRACE-FO loader needs; only buildDay changes.
%
%   TWO SUPPORTED WAYS TO USE YOUR COPY (both first-class):
%     (1) RAW CLOSURE (this file's pattern). Bind it as the 4-arg closure above
%         and hand it straight to acc.addSystem, or to LSQVCE.ENGINE.ADJUSTBLOCK's
%         struct-array form (.sysID/.nl/.loader). Simplest; no plot hooks or
%         parametrisation pipe (the engines guard those with isfield, so their
%         absence is safe -- see the Loader contract in LSQVCE.ENGINE.MAKELOADER).
%     (2) WRAPPED by LSQVCE.ENGINE.MAKELOADER. Keep your buildDay, pass it as
%         makeLoader(@(d,ot) myBuildDay(...), 'obsTypes',..., 'cov',..., ...)
%         to get a full Loader STRUCT with meta/cache plus the optional
%         'CacheDays', 'Priors', 'Scoped' tables and 'CommonParams'/'LocalParams'/
%         'ProcessingOptions' features. Use this when you want those conveniences.
%   Do NOT hand-roll a wrapper named makeLoader in your own code -- that name
%   is taken by LSQVCE.ENGINE.MAKELOADER; call it or bind a raw closure instead.
%
%   CONFIGURATION (cfg)
%     .nc, .nl        parameter counts
%     .x_c            nc x 1  true common parameters
%     .y              containers.Map: sat -> nl x 1 true local parameters
%     .obsTypes       cellstr, e.g. {'pos','kbr'}
%     .arcsPerDay     arcs per (day, obsType)
%     .mArc           rows per arc
%     .dt             sampling interval [s]
%     .cov.(obsType)  covariance spec (see LSQVCE.ENGINE.COVMODEL)
%     .seed           RNG seed
%     .collinear      make two local columns near-collinear (for tests)
%     .collinearEps   how near: 0 = exactly singular
%
%   THE TWO CALL SIGNATURES
%     meta = loader('meta', day, [], opt)   -> .obsTypes, .arcs
%     arc  = loader(day, obsType, arcID, opt)
%
%   WHY THERE IS A CACHE
%     Real design matrices are built per DAY (one integration of the
%     variational equations), but the accumulator consumes them per ARC. The
%     persistent single-day cache bridges that: buildDay runs once per
%     (day, obsType), and the arcs of that day are then sliced out of it.
%     Because the accumulator visits day -> obsType -> arc in order, all arcs
%     of a day are consumed consecutively and only ONE day is ever resident.
%
%     If even one day of A is too large, drop the cache and make buildDay
%     arc-aware: integrate only over the requested arc's epochs. The meta call
%     must then get the arc boundaries from a lightweight scan of the
%     observation timestamps, without touching the design matrix.
%
%   THE CACHE KEY IS A CORRECTNESS ISSUE
%     The key includes everything that changes the design matrix. If an option
%     in opt alters A or B (e.g. switching empirical accelerations on) it MUST
%     appear in the key, or the second pass in LSQVCE.ENGINE.RESIDUALS will silently
%     reuse a stale day and sigma0 will not reproduce. That mismatch is
%     detected and warned about, but a correct key avoids it entirely.
%
%   ARC-WISE LOCAL PARAMETERS
%     If the local parameters are per-arc (arc-wise initial state, say), then
%     nl counts ALL arcs of the system and B is mostly zeros. Return B sparse,
%     or the m x (6*nArc) dense block will dominate memory.
%
%   EXAMPLE
%     cfg = struct('nc',6, 'nl',4, 'obsTypes',{{'pos','kbr'}}, ...
%                  'arcsPerDay',2, 'mArc',400, 'dt',5, 'seed',42, ...
%                  'collinear',false, 'collinearEps',1e-7);
%     cfg.x_c = randn(6,1);
%     cfg.y   = containers.Map(1, randn(4,1));
%     cfg.cov.pos = struct('type','arma','ar',[1 -1.6 0.9 -0.25],'ma',1,'sigma',0.7);
%     cfg.cov.kbr = struct('type','diag','sigma',0.2);
%
%     acc = lsqvce.engine.BlockAccumulator(6, 'Backend','qr');
%     acc.addSystem(1, 4, @(d,ot,a,uo) synthLoader(cfg,1,d,ot,a,uo), 1:3);
%     sol = acc.solve();
%     max(abs(sol.x - cfg.x_c))          % recovers the truth
%
%   INPUT
%     cfg      1x1 struct  the synthetic configuration
%     sysID    1x1         which system to build for
%     day      1x1         day index
%     obsType  1xk  char   which stream
%     arcID    1x1         arc index within the day
%     opt      1x1 struct  runtime options
%
%   OUTPUT
%     out  1x1 struct  one batch: .A .B .l .t and the covariance spec
%
%   ACCURACY
%     A generator: its claims are determinism and the declared noise law.
%     The residual pass re-reads the same day, so sigma0 reproduces to 1e-10
%     relative (test_precision (testResidualPassReproducesSigma0)); with an
%     'arma' pos and a 'diag' kbr stream sigma0 = 1 within 0.05 and the
%     whitened residuals have std 1 within 0.08, lags 1-5 inside 2.5/sqrt(n)
%     (test_precision (testSigma0IsUnbiased, testWhitenedResidualsAreWhite)).
%     The 'psd' generator (an AR(200) fit) is not pinned. WP-139.
%
%   ERRORS
%     lsqvce:synthLoader:badArc          the arc index is out of range
%     lsqvce:synthLoader:noNoiseModel    the configuration declares no noise model
%     lsqvce:synthLoader:robustFactors   opt.rowWeights / epochWeights carry factors
%                                        for the arc: the template applies none
%                                        (makeLoader does)
%
%   See also LSQVCE.ENGINE.BLOCKACCUMULATOR, LSQVCE.ENGINE.COVMODEL
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-139: the help states its
%   accuracy, every number traced to the test that pins it. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-14 (CEST, UTC+02:00), WP-135 (review 3, lane 4): robust
%   factors for an arc are refused (lsqvce:synthLoader:robustFactors) -- the
%   template ignored them while the driver reported them applied. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-09-09 (CEST, UTC+02:00), R6 WP-102 (v6): the internal day
%   is cut into arc CELLS (the v6 day contract; the raw closure's arc slice
%   is a cell pick). With assistance from Claude (Anthropic), model Claude
%   Fable 5.1.
%   Updated 2026-08-31 (CEST, UTC+02:00), R4 G-06 (Level L1): the covariance
%   override is resolved by the package-private APPLYCOVOVERRIDE, so the
%   per-satellite shape LSQVCE.ENGINE.ADJUSTVCEJOINT emits is honoured by this template
%   (it used to read only the flat shape and silently dropped every variance
%   component under adjustVCE). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-40 r2: the 'psd' GENERATION
%   path floors cm.S before psd2acv (the one raw-S call in the toolbox).
%   With assistance from Claude (Anthropic), model Claude Fable 5.
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    persistent CACHE
    if isempty(CACHE), CACHE = struct('key','', 'D', []); end

    %--------------------------------------------------------------- meta call
    if (ischar(day) || isstring(day)) && strcmp(day, 'meta')
        out.obsTypes = cfg.obsTypes;
        out.arcs     = cell(1, numel(cfg.obsTypes));
        for k = 1:numel(cfg.obsTypes)
            out.arcs{k} = 1:cfg.arcsPerDay;
        end
        return
    end

    %---------------------------------------------------------------- arc call
    % Everything that can change A, B or l belongs in the key. That includes the
    % PROBLEM DIMENSIONS: cfg.nc sets the width of A and cfg.nl the width of B,
    % so two fixtures sharing (sat, day, obsType) but differing in nc/nl would
    % otherwise collide and the second would silently receive the first's
    % wrongly-sized arcs. cfg.seed likewise changes the generated data.
    tag  = lsqvce.util.getfielddef(opt, 'forceModelTag', '');
    seed = lsqvce.util.getfielddef(cfg, 'seed', 0);
    key  = sprintf('%d|%d|%s|%s|nc%d|nl%d|sd%g', sat, day, obsType, tag, ...
                   cfg.nc, cfg.nl, seed);

    if ~strcmp(CACHE.key, key)
        CACHE.D   = buildDay(cfg, sat, day, obsType);
        CACHE.key = key;                        % the previous day is released here
    end
    D = CACHE.D;

    if arcID < 1 || arcID > numel(D.A)
        error('lsqvce:synthLoader:badArc', ...
              'Arc %d does not exist (day %d has %d arcs).', ...
              arcID, day, numel(D.A));
    end

    % WP-135 (review 3, lane 4): the template applies no robust factors --
    % makeLoader does (epochWeights as a nugget before the law, rowWeights
    % after it). A factor meant for this arc was ignored while the driver
    % reported it applied; it is refused instead.
    rk = lsqvce.util.rowWeightKey(lsqvce.util.getfielddef(opt, 'sysID', sat), day, obsType, arcID);
    for f = {'rowWeights', 'epochWeights'}
        if isfield(opt, f{1}) && ~isempty(opt.(f{1})) && isKey(opt.(f{1}), rk)
            error('lsqvce:synthLoader:robustFactors', ...
                  ['opt.%s carries factors for arc %s, which the template does not apply. ' ...
                   'Wrap your buildDay with lsqvce.engine.makeLoader for a robust pass.'], f{1}, rk);
        end
    end

    out.A = D.A{arcID};                         % v6: the cells are the arcs
    out.B = D.B{arcID};
    out.l = D.l{arcID};
    out.t = D.t{arcID};

    % Runtime override of the covariance model, if the caller asked for one.
    % R4 G-06: resolved by the ONE package-private merge, so the per-satellite
    % shape LSQVCE.ENGINE.ADJUSTVCEJOINT emits (opt.covOverride.bySys.(satKey).(obsType))
    % reaches a template-derived loader exactly as it reaches a makeLoader
    % batch. The flat shape opt.covOverride.(obsType) keeps working.
    cm = cfg.cov.(obsType);
    if isfield(opt,'covOverride')
        cm = applyCovOverride(cm, opt.covOverride, obsType, opt, day, arcID);   % R7-13: + the per-arc factor
    end

    [out.Cfcn, out.Wfcn, info] = lsqvce.engine.covModel(cm, out.t, opt);
    out.filterMemory = info.filterMemory;       % the accumulator reads this
    if lsqvce.util.getfielddef(info, 'wholeArc', false)
        out.wholeArc = true;                    % a dense law is whitened whole (WP-125)
        out.filterMemory = numel(out.l);
    end
end

%==========================================================================
function D = buildDay(cfg, sat, day, obsType)
%BUILDDAY  Build one day of design matrix and observations.
%   In a real loader this is where the orbit integration, partial derivatives
%   and observation modelling happen. Everything else in this file is plumbing
%   that carries over unchanged.

    M  = cfg.arcsPerDay * cfg.mArc;

    % Deterministic per-(sat, day, obsType) stream: the loader MUST return the
    % same data on the second pass, or LSQVCE.ENGINE.RESIDUALS cannot reproduce sigma0.
    sd = mod(cfg.seed + 1000*sat + 7*day + 31*sum(double(obsType)), 2^31 - 1);
    rs = RandStream('mt19937ar', 'Seed', sd);

    D.t = (0:M-1).' * cfg.dt;
    D.A = randn(rs, M, cfg.nc);
    D.B = randn(rs, M, cfg.nl);

    % Obs types carry different information about the common parameters.
    switch obsType
        case 'kbr'
            k = 1:min(3, cfg.nc);
            D.A(:,k) = 0.1 * D.A(:,k);
    end

    % Optional near-collinear local columns -> ill-conditioned local block.
    if isfield(cfg,'collinear') && cfg.collinear && cfg.nl >= 2
        D.B(:,end) = D.B(:,end-1) + cfg.collinearEps * randn(rs, M, 1);
    end

    e   = synthNoise(rs, cfg.cov.(obsType), M, cfg.dt);
    y   = cfg.y(sat);
    l   = D.A*cfg.x_c + D.B*y + e;

    % v6 (R6 WP-102): THE CELLS ARE THE ARCS -- cut the day's M rows into
    % arcsPerDay cells of mArc rows each
    A = D.A;  B = D.B;  t = D.t;
    D = struct('A', {cell(1, cfg.arcsPerDay)}, 'B', {cell(1, cfg.arcsPerDay)}, ...
               'l', {cell(1, cfg.arcsPerDay)}, 't', {cell(1, cfg.arcsPerDay)});
    for a = 1:cfg.arcsPerDay
        r = (a-1)*cfg.mArc + (1:cfg.mArc);
        D.A{a} = A(r, :);  D.B{a} = B(r, :);  D.l{a} = l(r);  D.t{a} = t(r);
    end
end

%==========================================================================
function e = synthNoise(rs, cm, M, dt)
%SYNTHNOISE  Draw noise with exactly the covariance the model claims.
%   The burn-in matters: without it the first samples are not stationary and
%   sigma0 comes out slightly below 1, which looks like a bug in the estimator
%   when it is really a bug in the data generator.

    nBurn = 2000;

    switch lower(cm.type)
    case 'diag'
        s = cm.sigma;
        if isscalar(s)
            e = s * randn(rs, M, 1);
        else
            e = s(:) .* randn(rs, M, 1);
        end

    case 'arma'
        a = cm.ar(:).';
        b = lsqvce.util.getfielddef(cm, 'ma', 1);
        w = cm.sigma * randn(rs, M + nBurn, 1);
        e = filter(b, a, w);                    % COLOURING filter: b/a
        e = e(nBurn+1:end);

    case 'psd'
        % Fit a long AR to the PSD and drive it -- exact for AR-type spectra,
        % and a close approximation otherwise.
        % WP-40 r2 (latent-bug fix found while chasing F-15): this is the
        % ONE psd2acv call fed RAW cm.S -- a fed-back or estimated PSD with
        % negative ordinates crashed the GENERATION path while the modelling
        % path (covModel.psdModel) floors. Same floor here.
        Sgen = cm.S(:);
        Sgen(~isfinite(Sgen) | Sgen < 0) = 0;
        Sgen = max(Sgen, 1e-10 * max(Sgen));
        [a, s2] = lsqvce.kernel.levinsonDurbin(lsqvce.kernel.psd2acv(cm.f, Sgen, dt, 201), 200);
        w = sqrt(s2) * randn(rs, M + nBurn, 1);
        e = filter(1, a, w);                    % COLOURING filter: 1/a
        e = e(nBurn+1:end);

    otherwise
        error('lsqvce:synthLoader:noNoiseModel', ...
              'synthLoader cannot generate noise for covariance type ''%s''.', ...
              cm.type);
    end
end
