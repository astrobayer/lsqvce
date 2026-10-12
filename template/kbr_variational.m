function [data, kit] = kbr_variational(cfg)
%KBR_VARIATIONAL  GRACE-FO range-rate day builder by the variational equations.
%
%   SYNTAX
%     data        = kbr_variational(cfg)
%     [data, kit] = kbr_variational(cfg)
%
%   DESCRIPTION
%     A real-data template for the lsqvce `data` declaration (the Manual,
%     chapter 3): the GRACE-FO Level-1B range-rates (KBR1B) of the C/D pair,
%     and their reduced-dynamic positions (GNV1B) that make the initial
%     states observable, as observation equations of a gravity field in
%     fully normalized spherical harmonics -- by the VARIATIONAL EQUATIONS.
%     Per arc and satellite the reference orbit and
%         Psi(t) = d[r; v](t) / d[r0; v0 | C_nm, S_nm | ACC bias | ACC scale]
%     are integrated together,
%         d/dt [r; v] = [v; a(t, r)],
%         d/dt Psi    = [0 I; G 0] Psi + [0; da/dp],
%     with the field's acceleration a, its gradient G = da/dr and the
%     partials da/dC_nm, da/dS_nm all ANALYTIC: the potential is expanded in
%     normalized solid harmonics by the Cartesian (pole-safe) recursion, and
%     d/dx, d/dy, d/dz are exact sparse maps degree n -> n+1 between
%     coefficient vectors, so a, G and da/dp come out of ONE harmonic sweep
%     per force evaluation. The integrator is an 8th-order Adams-Bashforth-
%     Moulton (PECE, exact rational coefficients) on the observation grid,
%     started by RK4 substeps. A range-rate row is
%         d rho_dot = [ h' e' ] (Psi_D - Psi_C),
%         e = (r_D - r_C)/rho,  h = (v_D - v_C - rho_dot e)/rho,
%     formed on the fly while the pair is integrated (Psi is never stored --
%     at degree 60 one satellite-day of Psi would be 3 GB).
%
%     The parameters (scoped tables, the Manual, chapter 3 section 6):
%       common   C_nm, S_nm for n = cfg.degree(1) .. cfg.degree(2)
%       'arc'    per satellite r0, v0 (6 + 6); KBR bias, drift, cos u, sin u (4)
%       'day'    per satellite ACC bias (3) and scale (3)            (6 + 6)
%     Every column carries a declared UNIT (column 5 of commonParams, column 3
%     of localParams): the column is the partial times the unit, so the
%     estimate is in that unit -- physical value = unit .* estimate
%     (kit.unitCommon .* R.x for the field). The units bring every column to
%     O(1): unscaled, the columns span ~10 decades and the QR rank decision
%     (1e3 eps on the pivots) would see the ACC scale as rank-deficient.
%
%   INPUT (cfg, 1x1 struct; [] = default)
%     days          1xD  double  MJD of the days to process            REQUIRED
%     dataDir       char         folder of the L1B files               REQUIRED
%     gravityField  char | struct  ICGEM .gfc file, or struct .GM .R .C .S
%                                  ((N+1)x(N+1) fully normalized)     REQUIRED
%     trf2crf       function_handle [R, Rdot] = f(tGps): r_crf = R r_trf,
%                                  v_crf = R v_trf + Rdot r_trf     REQUIRED
%                   (the IERS 2010 transformation for real data, e.g. from
%                   SOFA with EOP; kit.eraRotation is the Earth-rotation-
%                   angle-only stand-in for synthetic data -- see LIMITATIONS)
%     degree        1x2  [nmin nmax] estimated degrees           ([2 60])
%     fieldDegree   1x1  degree of the reference field     (degree(2), >= it)
%     arcsPerDay    1x1  arcs per day; arc length 86400/arcsPerDay       (1)
%     arcs          fh   @(day, obsType) arc ids (a resolver, the Manual
%                        chapter 3 section 4.2) instead of 1:arcsPerDay ([])
%     step          1x1  integration and KBR step, s                    (5)
%     posSampling   1x1  GNV1B position sampling used, s (multiple of step) (30)
%     sigma         struct .kbr (m/s, 2e-7), .pos (m, 0.02): the starting
%                        LEVEL of each stream (Estimate="noise" updates it)
%     obsTypes      cellstr  subset of {'kbr','posC','posD'}  (all three)
%     files         struct .gnv .kbr .acc .sca: sprintf patterns of the file
%                        names, called with (date 'yyyy-mm-dd', satellite
%                        letter) -- default GRACE-FO v04: 'GNV1B_%s_%s_04.txt',
%                        'KBR1B_%s_%s_04.txt' (letter 'Y'), 'ACT1B_%s_%s_04.txt',
%                        'SCA1B_%s_%s_04.txt'
%     formats       struct .gnv .kbr .acc .sca: textscan formats and the
%                        field positions (see L1BFORMATS below) -- default the
%                        GRACE-FO Level-1B v04 ASCII columns
%     quatConvention 'srf2crf' | 'crf2srf'  what the SCA1B quaternion rotates
%                                                             ('srf2crf')
%     accToSrf      3x3  rotation of the ACC product's axes into the SRF (eye(3))
%     accBias0, accScale0  struct .C .D (3x1): the a-priori ACC calibration
%                        (0 and 1) -- the linearization point of the 'day' rows
%     fieldUpdate   np x 1  a-priori correction of the estimated C_nm, S_nm
%                        (physical, commonParams order) -- a second Gauss-Newton
%                        iteration passes kit.unitCommon .* R.x here      (0)
%     background    function_handle  a_bg = f(tGps, r_crf, sat): every other
%                        force (tides, third bodies, AOD1B, relativity),
%                        m/s^2 in the CRF                              (none)
%     maxAccGap     1x1  longest ACC/SCA gap bridged by interpolation, s (10)
%
%   OUTPUT
%     data  1x1 struct  the lsqvce declaration: commonParams, days, systems(1)
%           (id 1, name 'GRACE-FO', buildDay, covariance, arcs, localParams,
%           processingOptions) -- R = lsqvce.fit(data, opts)
%     kit   1x1 struct  the physics as function handles, for tests and for a
%           second iteration: unitCommon, unitLocal, coeffIndex, harmonics,
%           makeField, fieldEval, propagate, accModel, eraRotation, readGfc,
%           readL1B, buildDay
%
%   ACCURACY (dev/validation/kbr_variational_mirror.py, the numerical twin)
%     harmonics vs Legendre functions 7e-11; d/dx operators vs differences
%     2e-9; acceleration vs closed-form J2 2e-16; gradient vs differences
%     8e-10; integrator vs Kepler after one day at 5 s 6e-6 m; Psi vs
%     re-integration 1e-6 (rounding floor); range-rate partials vs
%     differences 9e-5 (rounding of rho_dot); closed loop (one day, degree 8,
%     2 % ACC scale error planted) 3.5e-5 of the planted field after one
%     linearization, 8e-6 after a second (cfg.fieldUpdate, accBias0/Scale0).
%
%   ERRORS
%     kbrvar:badConfig       a REQUIRED cfg field is missing, or a field is malformed
%     kbrvar:noFile          an L1B or .gfc file does not exist
%     kbrvar:badFile         a file has no YAML header end, or no data rows
%     kbrvar:gfcTimeVariable the .gfc carries time-variable terms (gfct/trnd/asin/acos)
%     kbrvar:noInitialState  an arc has no GNV1B epoch on the grid for a satellite
%     kbrvar:accGap          ACC1B/SCA1B leave a gap longer than maxAccGap in an arc
%     kbrvar:offGrid         a KBR1B epoch is not on the integration grid
%     kbrvar:shortArc        an arc holds fewer than 8 integration steps
%     kbrvar:emptyArc        a scheduled arc carries no data (declare cfg.arcs)
%     (+ every declaration check of lsqvce.fit, the Manual chapter 3 section 3.5)
%
%   EXAMPLE
%     cfg = struct('days', 60000:60029, 'dataDir', '/data/gracefo/L1B', ...
%                  'gravityField', '/data/models/GOCO06s.gfc', ...
%                  'trf2crf', @myIers2010, 'degree', [2 60], 'arcsPerDay', 4);
%     [data, kit] = kbr_variational(cfg);
%     R  = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Noise=struct('Grid', 1080)));
%     dC = kit.unitCommon .* R.x;                  % physical corrections
%     cfg.fieldUpdate = dC;  data = kbr_variational(cfg);   % iteration 2
%
%   LIMITATIONS
%     - The Earth-rotation stand-in kit.eraRotation omits precession-nutation
%       (0.35 deg since 2000) and polar motion; on real data use the IERS
%       2010 transformation (cfg.trf2crf). It is fine for synthetic data.
%     - The L1B column layouts, the quaternion direction and the ACC frame
%       are the GRACE-FO v04 conventions as declared in cfg.formats /
%       quatConvention / accToSrf: check them against the product handbook
%       of the release you process.
%     - Every force but the static field and the accelerometer enters through
%       cfg.background; its velocity dependence is not in Psi.
%     - GNV1B is a reduced-dynamic orbit; for a gravity solution replace it by
%       a kinematic orbit (a reader with the same output) and its covariance.
%     - KBR1B quality flags are not screened (NaN rows are dropped).
%     - Memory per day: 8 * m * nc bytes per stream (KBR at 5 s, degree 60:
%       0.5 GB); shorten the day or the degree, or thin posSampling.
%     - One linearization per call: iterate with cfg.fieldUpdate and the ACC
%       a-priori; the initial states always start from GNV1B.
%
%   See also LSQVCE.FIT, LSQVCE.OPTIONS, LSQVCE.LAYOUT, GRAVFIRE_DATA
%
%   Created 2026-10-12 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- WP-122.

arguments
    cfg (1,1) struct
end
    cfg = withDefaults(cfg);
    gf  = readGravity(cfg.gravityField);
    N   = cfg.fieldDegree;
    if size(gf.C, 1) < N + 1
        error('kbrvar:badConfig', 'the gravity field holds degree %d, fieldDegree is %d.', ...
              size(gf.C, 1) - 1, N);
    end

    % ---- the common table: C_nm, S_nm of the estimated degrees, coeffIndex order
    [nn, mm, cs] = coeffIndex(N);
    est = find(nn >= cfg.degree(1) & nn <= cfg.degree(2));
    np  = numel(est);
    csName = {'C', 'S'};
    unitSh = 1e-10;                                       % the ONE unit of the SH columns
    data.commonParams = [repmat({'common'}, np, 1), csName(cs(est) + 1).', ...
                         num2cell(nn(est)), num2cell(mm(est)), num2cell(unitSh * ones(np, 1))];
    data.days = cfg.days(:).';

    % ---- the field: reference coefficients + the a-priori correction
    c = fieldVector(gf, N);
    if ~isempty(cfg.fieldUpdate)
        if numel(cfg.fieldUpdate) ~= np
            error('kbrvar:badConfig', 'fieldUpdate has %d entries; the estimated set has %d.', ...
                  numel(cfg.fieldUpdate), np);
        end
        c(est) = c(est) + cfg.fieldUpdate(:);
    end
    F = makeField(c, N, est, gf.GM, gf.R);

    % ---- the local table: 16 'arc' rows, 12 'day' rows, unit in column 3
    sat = {'C', 'D'};
    icName = {'x0', 'y0', 'z0', 'vx0', 'vy0', 'vz0'};
    icUnit = [1 1 1 1e-3 1e-3 1e-3];
    accName = {'accBias_x', 'accBias_y', 'accBias_z', 'accScale_x', 'accScale_y', 'accScale_z'};
    accUnit = [1e-9 1e-9 1e-9 1e-3 1e-3 1e-3];
    arcRows = {};  dayRows = {};
    for k = 1:2
        arcRows = [arcRows; strcat(sat{k}, '_', icName(:)), num2cell(icUnit(:))]; %#ok<AGROW>
        dayRows = [dayRows; strcat(sat{k}, '_', accName(:)), num2cell(accUnit(:))]; %#ok<AGROW>
    end
    if any(strcmp(cfg.obsTypes, 'kbr'))                   % the KBR empiricals need the KBR stream
        arcRows = [arcRows; {'kbr_bias'; 'kbr_drift'; 'kbr_cos1cpr'; 'kbr_sin1cpr'}, num2cell(1e-7 * ones(4, 1))];
    end
    localParams = [[repmat({'arc'}, size(arcRows, 1), 1), arcRows]; ...
                   [repmat({'day'}, size(dayRows, 1), 1), dayRows]];

    % ---- covariance: the temporal law per stream (white to start); positions
    %      carry their sigmas with the day (Cep), so their level is 1
    cov = struct();
    for k = 1:numel(cfg.obsTypes)
        switch cfg.obsTypes{k}
            case 'kbr',  cov.kbr = struct('type', 'diag', 'sigma', cfg.sigma.kbr);
            otherwise,   cov.(cfg.obsTypes{k}) = struct('type', 'diag', 'sigma', 1);
        end
    end

    % ---- what buildDay needs, bound once (cache: a handle, one day resident)
    po = cfg;
    po.field = F;  po.np = np;  po.unitSh = unitSh;
    po.unitIc = icUnit;  po.unitAcc = accUnit;
    po.cache = containers.Map('KeyType', 'double', 'ValueType', 'any');

    s = struct();
    s.id = 1;
    s.name = 'GRACE-FO';
    s.buildDay = @buildKbrDay;
    s.covariance = cov;
    if isempty(cfg.arcs), s.arcs = cfg.arcsPerDay; else, s.arcs = cfg.arcs; end
    s.localParams = localParams;
    s.processingOptions = po;
    data.systems = s;

    kit = struct('unitCommon', unitSh * ones(np, 1), 'unitLocal', [cell2mat(arcRows(:, 2)); cell2mat(dayRows(:, 2))], ...
                 'coeffIndex', @coeffIndex, 'harmonics', @harmonics, 'makeField', @makeField, ...
                 'fieldEval', @fieldEval, 'propagate', @propagate, 'accModel', @accModel, ...
                 'eraRotation', @eraRotation, 'readGfc', @readGfc, 'readL1B', @readL1B, ...
                 'buildDay', @buildKbrDay, 'field', F);
end

% ==========================================================================
%  DECLARATION HELPERS
% ==========================================================================
function cfg = withDefaults(cfg)
%WITHDEFAULTS  The documented defaults; the REQUIRED fields are refused by name.
    for f = {'days', 'dataDir', 'gravityField', 'trf2crf'}
        if ~isfield(cfg, f{1}) || isempty(cfg.(f{1}))
            error('kbrvar:badConfig', 'cfg.%s is required (see help kbr_variational).', f{1});
        end
    end
    if ~isa(cfg.trf2crf, 'function_handle')
        error('kbrvar:badConfig', 'cfg.trf2crf must be a function handle [R, Rdot] = f(tGps).');
    end
    d = struct('degree', [2 60], 'fieldDegree', [], 'arcsPerDay', 1, 'arcs', [], 'step', 5, ...
               'posSampling', 30, 'sigma', struct('kbr', 2e-7, 'pos', 0.02), ...
               'obsTypes', {{'kbr', 'posC', 'posD'}}, 'files', l1bFiles(), 'formats', l1bFormats(), ...
               'quatConvention', 'srf2crf', 'accToSrf', eye(3), ...
               'accBias0', struct('C', zeros(3, 1), 'D', zeros(3, 1)), ...
               'accScale0', struct('C', ones(3, 1), 'D', ones(3, 1)), ...
               'fieldUpdate', [], 'background', [], 'maxAccGap', 10);
    for f = fieldnames(d).'
        if ~isfield(cfg, f{1}) || isempty(cfg.(f{1})), cfg.(f{1}) = d.(f{1}); end
    end
    for f = fieldnames(d.sigma).'                          % a partial sigma keeps the other default
        if ~isfield(cfg.sigma, f{1}), cfg.sigma.(f{1}) = d.sigma.(f{1}); end
    end
    cfg.obsTypes = cellstr(cfg.obsTypes);
    cfg.obsTypes = cfg.obsTypes(:).';
    if isempty(cfg.fieldDegree), cfg.fieldDegree = cfg.degree(2); end
    if numel(cfg.degree) ~= 2 || cfg.degree(1) < 2 || cfg.degree(2) < cfg.degree(1) ...
            || cfg.fieldDegree < cfg.degree(2)
        error('kbrvar:badConfig', 'degree must be [nmin nmax], 2 <= nmin <= nmax <= fieldDegree.');
    end
    if mod(cfg.posSampling, cfg.step) ~= 0 || mod(86400 / cfg.arcsPerDay, cfg.step) ~= 0
        error('kbrvar:badConfig', 'posSampling and the arc length must be multiples of step.');
    end
    bad = setdiff(cfg.obsTypes, {'kbr', 'posC', 'posD'});
    if ~isempty(bad)
        error('kbrvar:badConfig', 'unknown observation type ''%s''.', bad{1});
    end
    if ~any(strcmp(cfg.quatConvention, {'srf2crf', 'crf2srf'}))
        error('kbrvar:badConfig', 'quatConvention must be ''srf2crf'' or ''crf2srf''.');
    end
end

function f = l1bFiles()
%L1BFILES  Default GRACE-FO v04 file names: sprintf(pattern, 'yyyy-mm-dd', letter).
    f = struct('gnv', 'GNV1B_%s_%s_04.txt', 'kbr', 'KBR1B_%s_%s_04.txt', ...
               'acc', 'ACT1B_%s_%s_04.txt', 'sca', 'SCA1B_%s_%s_04.txt');
end

function f = l1bFormats()
%L1BFORMATS  The GRACE-FO Level-1B v04 ASCII columns: a textscan format and
%   the positions of the fields used, in the textscan output.
%     GNV1B  gps_time id coord_ref xpos ypos zpos (3 err) xvel yvel zvel (3 err) qualflg
%     KBR1B  gps_time biased_range range_rate range_accl iono_corr lighttime_corr
%            lighttime_rate lighttime_accl ant_centr_corr ant_centr_rate
%            ant_centr_accl (4 SNR) qualflg
%     ACT1B  gps_time id lin_accl_x y z ang_accl_x y z acl_x_res y z qualflg
%     SCA1B  gps_time id sca_id quatangle quaticoeff quatjcoeff quatkcoeff qual_rss qualflg
    f.gnv = struct('fmt', ['%f %s %s' repmat(' %f', 1, 12) ' %s'], 't', 1, 'pos', 4:6, 'vel', 10:12);
    f.kbr = struct('fmt', [repmat('%f ', 1, 15) '%s'], 't', 1, 'rr', 3, 'ltr', 7, 'acr', 10);
    f.acc = struct('fmt', ['%f %s' repmat(' %f', 1, 9) ' %s'], 't', 1, 'acc', 3:5);
    f.sca = struct('fmt', ['%f %s' repmat(' %f', 1, 6) ' %s'], 't', 1, 'q', 4:7);
end

% ==========================================================================
%  THE DAY BUILDER
% ==========================================================================
function D = buildKbrDay(day, obsType, ~, ~, po, info)
%BUILDKBRDAY  One (day, stream): A, B, l, t (+ Cep) -- one cell per arc. The
%   whole day (every stream) is built by ONE integration of each arc and kept
%   in po.cache; the other streams of the day are served from it.
    if ~isKey(po.cache, day)
        if po.cache.Count > 0, remove(po.cache, keys(po.cache)); end   % one day resident
        po.cache(day) = computeDay(day, po, info);
    end
    built = po.cache(day);
    D = built.(obsType);
end

function out = computeDay(day, po, info)
%COMPUTEDAY  Read the day's L1B files, integrate every arc of the pair with
%   its variational equations, and form the rows of every stream.
    dateStr = char(datetime(day, 'ConvertFrom', 'modifiedjuliandate', 'Format', 'yyyy-MM-dd'));
    t0day = (day - 51544.5) * 86400;                      % GPS seconds of the day start
    sats = {'C', 'D'};
    for k = 1:2
        g  = readL1B(fullfile(po.dataDir, sprintf(po.files.gnv, dateStr, sats{k})), po.formats.gnv);
        ac = readL1B(fullfile(po.dataDir, sprintf(po.files.acc, dateStr, sats{k})), po.formats.acc);
        q  = readL1B(fullfile(po.dataDir, sprintf(po.files.sca, dateStr, sats{k})), po.formats.sca);
        S.(sats{k}) = struct('gnv', g, 'acc', ac, 'sca', q);
    end
    kb = struct('t', zeros(0, 1));                        % read only when the stream is asked for
    if any(strcmp(po.obsTypes, 'kbr'))
        kb = readL1B(fullfile(po.dataDir, sprintf(po.files.kbr, dateStr, 'Y')), po.formats.kbr);
    end
    ids = info.arcId;
    nArc = numel(ids);
    arcLen = 86400 / po.arcsPerDay;
    for s = po.obsTypes(:).'
        out.(s{1}) = struct('A', {cell(1, nArc)}, 'B', {cell(1, nArc)}, 'l', {cell(1, nArc)}, ...
                            't', {cell(1, nArc)});
        if s{1}(1) == 'p', out.(s{1}).Cep = cell(1, nArc); end
    end
    for a = 1:nArc
        span = t0day + (ids(a) - 1) * arcLen + [0 arcLen];
        arc = computeArc(span, S, kb, po, info, ids(a));
        for s = po.obsTypes(:).'
            o = arc.(s{1});
            out.(s{1}).A{a} = o.A;  out.(s{1}).B{a} = o.B;
            out.(s{1}).l{a} = o.l;  out.(s{1}).t{a} = o.t;
            if s{1}(1) == 'p', out.(s{1}).Cep{a} = o.Cep; end
        end
        if all(cellfun(@(s) isempty(arc.(s).l), po.obsTypes))
            error('kbrvar:emptyArc', ['day %d, arc %d carries no data; its ''arc'' unknowns would be ' ...
                  'unobserved -- declare cfg.arcs (a resolver) without it.'], day, ids(a));
        end
    end
end

function arc = computeArc(span, S, kb, po, info, arcId)
%COMPUTEARC  One arc of the pair: the common start t0 on the grid, the two
%   initial states from GNV1B, the integration of [r v Psi] for C and D
%   together, and the rows of every stream, formed at each grid epoch.
    h = po.step;  np = po.np;  nc = info.nc;  nl = info.nl;
    sats = {'C', 'D'};
    % ---- start: the first grid epoch of the arc where BOTH satellites have a GNV1B state
    onGrid = @(t) abs(t - span(1) - h * round((t - span(1)) / h)) < 1e-3;
    tc = cell(1, 2);
    for k = 1:2
        g = S.(sats{k}).gnv;
        tc{k} = g.t(g.t >= span(1) & g.t < span(2) & onGrid(g.t));
    end
    tBoth = intersect(round(tc{1} * 1e3), round(tc{2} * 1e3)) / 1e3;
    if isempty(tBoth)
        error('kbrvar:noInitialState', 'arc [%.0f %.0f] (GPS s): no GNV1B epoch on the grid for both satellites.', span);
    end
    t0 = tBoth(1);
    nstep = round((span(2) - h - t0) / h);
    if nstep < 8
        error('kbrvar:shortArc', 'arc [%.0f %.0f] (GPS s) holds %d steps from its first state; 8 are needed.', span, nstep);
    end
    tEnd = t0 + nstep * h;
    ncol = 6 + np + 6;  nY = 6 + 6 * ncol;
    Y0 = zeros(2 * nY, 1);
    for k = 1:2
        g = S.(sats{k}).gnv;
        i = find(abs(g.t - t0) < 1e-3, 1);
        [R, Rd] = po.trf2crf(t0);
        r = g.pos(i, :).';  v = g.vel(i, :).';
        Y0((k - 1) * nY + (1:nY)) = [R * r; R * v + Rd * r; reshape([eye(6), zeros(6, ncol - 6)], [], 1)];
        a = S.(sats{k}).acc;  q = S.(sats{k}).sca;
        ctx.(sats{k}) = struct('F', po.field, 'np', np, 'ncol', ncol, 'rot', po.trf2crf, ...
            'accfun', accModel(a.t, a.acc, q.t, q.q, po.accToSrf, po.quatConvention, [t0 tEnd], po.maxAccGap), ...
            'b0', po.accBias0.(sats{k})(:), 's0', po.accScale0.(sats{k})(:), 'bg', po.background, 'sat', sats{k});
    end

    % ---- which grid epoch carries which row
    want = @(s) any(strcmp(po.obsTypes, s));
    kbrAt = zeros(nstep + 1, 1);  obsK = zeros(0, 1);  tK = zeros(0, 1);
    if want('kbr')
        sel = kb.t >= t0 & kb.t <= tEnd + 1e-3;
        tk = kb.t(sel);  rr = kb.rr(sel) + kb.ltr(sel) + kb.acr(sel);
        ok = isfinite(rr);  tk = tk(ok);  rr = rr(ok);
        kk = (tk - t0) / h;
        if any(abs(kk - round(kk)) > 1e-3)
            error('kbrvar:offGrid', 'KBR1B epochs off the %g s grid of arc [%.0f %.0f] (GPS s).', h, span);
        end
        kbrAt(round(kk) + 1) = 1:numel(tk);  obsK = rr;  tK = tk;
    end
    posAt = zeros(nstep + 1, 2);  obsP = cell(1, 2);  tP = cell(1, 2);
    for k = 1:2
        if ~want(['pos' sats{k}]), continue, end
        g = S.(sats{k}).gnv;
        sel = g.t >= t0 & g.t <= tEnd + 1e-3 & ...
              abs(g.t - span(1) - po.posSampling * round((g.t - span(1)) / po.posSampling)) < 1e-3 & ...
              all(isfinite(g.pos), 2);
        tP{k} = g.t(sel);  obsP{k} = g.pos(sel, :);
        posAt(round((tP{k} - t0) / h) + 1, k) = 1:numel(tP{k});
    end

    % ---- the rows, preallocated; B as triplets (it spans the whole adjustment)
    off = info.offsetOf(arcId);  dOff = info.dayOffset;
    colIc = {off + (1:6), off + (7:12)};  colKbr = off + (13:16);
    colAcc = {dOff + (1:6), dOff + (7:12)};
    accIdx = 6 + np + (1:6);  shIdx = 6 + (1:np);  uSh = po.unitSh;
    mK = numel(tK);
    AK = zeros(mK, nc);  lK = zeros(mK, 1);
    iK = zeros(mK * 28, 1);  jK = iK;  vK = iK;
    AP = cell(1, 2);  lP = AP;  iP = AP;  jP = AP;  vP = AP;
    for k = 1:2
        m3 = 3 * numel(tP{k});
        AP{k} = zeros(m3, nc);  lP{k} = zeros(m3, 1);
        iP{k} = zeros(m3 * 12, 1);  jP{k} = iP{k};  vP{k} = iP{k};
    end
    tMid = mean(span);

    rhs = @(t, Y) [satRhs(t, Y(1:nY), ctx.C); satRhs(t, Y(nY + 1:end), ctx.D)];
    abm8(rhs, t0, Y0, h, nstep, @emit);

    % ---- the streams of this arc
    arc.kbr = struct('A', AK, 'B', sparse(iK, jK, vK, mK, nl), 'l', lK, 't', tK(:));
    for k = 1:2
        nE = numel(tP{k});
        % B DENSE on the multi-axis streams: the toolbox's per-epoch whitener of
        % a stream with Cep reshapes B to 3-D (lsqvce.kernel.separableWhiten>tfilt),
        % which MATLAB refuses for sparse arrays (MATLAB:reshape:sparseNDmatrix,
        % WP-122 CI round 1) -- although the Manual allows a sparse B
        arc.(['pos' sats{k}]) = struct('A', AP{k}, 'B', full(sparse(iP{k}, jP{k}, vP{k}, 3 * nE, nl)), ...
            'l', lP{k}, 't', repelem(tP{k}(:), 3), 'Cep', po.sigma.pos * ones(nE, 3));
    end

    function emit(ek, eY)
    %EMIT  The rows of grid epoch ek: KBR (both satellites' Psi), positions.
    %   (Its own names only: a nested function shares every name it has in
    %   common with computeArc.)
        et = t0 + ek * h;
        eyc = eY(1:6);  ePc = reshape(eY(7:nY), 6, ncol);
        eyd = eY(nY + (1:6));  ePd = reshape(eY(nY + 7:end), 6, ncol);
        ei = kbrAt(ek + 1);
        if ei > 0
            edr = eyd(1:3) - eyc(1:3);  edv = eyd(4:6) - eyc(4:6);
            erho = norm(edr);  ee = edr / erho;  erd = ee.' * edv;
            eJ = [((edv - erd * ee) / erho).', ee.'];
            egC = eJ * ePc;  egD = eJ * ePd;
            AK(ei, 1:np) = (egD(shIdx) - egC(shIdx)) * uSh;
            eu = argLat(eyc);
            evals = [-egC(1:6) .* po.unitIc, egD(1:6) .* po.unitIc, ...
                     [1, (et - tMid) / 3600, cos(eu), sin(eu)] * 1e-7, ...
                     -egC(accIdx) .* po.unitAcc, egD(accIdx) .* po.unitAcc];
            ecols = [colIc{1}, colIc{2}, colKbr, colAcc{1}, colAcc{2}];
            ep = (ei - 1) * 28 + (1:28);
            iK(ep) = ei;  jK(ep) = ecols;  vK(ep) = evals;
            lK(ei) = obsK(ei) - erd;
        end
        for es = 1:2
            ej = posAt(ek + 1, es);
            if ej == 0, continue, end
            if es == 1, ey = eyc; eP = ePc; else, ey = eyd; eP = ePd; end
            eR = po.trf2crf(et);
            erows = 3 * (ej - 1) + (1:3);
            AP{es}(erows, 1:np) = eP(1:3, shIdx) * uSh;
            lP{es}(erows) = eR * obsP{es}(ej, :).' - ey(1:3);
            evals = [eP(1:3, 1:6) .* po.unitIc, eP(1:3, accIdx) .* po.unitAcc];
            [ejj, eii] = meshgrid([colIc{es}, colAcc{es}], erows);
            ep = (ej - 1) * 36 + (1:36);
            iP{es}(ep) = eii(:);  jP{es}(ep) = ejj(:);  vP{es}(ep) = evals(:);
        end
    end
end

function u = argLat(y)
%ARGLAT  Argument of latitude of the state y = [r; v] (rad).
    r = y(1:3);  hh = cross(r, y(4:6));  hh = hh / norm(hh);
    nod = cross([0; 0; 1], hh);  nod = nod / norm(nod);
    u = atan2(dot(cross(nod, r), hh), dot(nod, r));
end

% ==========================================================================
%  DYNAMICS: force model, variational equations, integrator
% ==========================================================================
function dY = satRhs(t, Y, C)
%SATRHS  d/dt [r; v; Psi] of one satellite in the CRF. Psi columns:
%   [r0 v0 (6) | C_nm, S_nm (np) | ACC bias (3) | ACC scale (3)].
    r = Y(1:3);  v = Y(4:6);  Psi = reshape(Y(7:end), 6, C.ncol);
    R = C.rot(t);
    [a, G, P] = fieldEval(C.F, R.' * r);
    [Q, acc] = C.accfun(t);
    a = R * a + Q * (C.b0 + C.s0 .* acc);
    if ~isempty(C.bg), a = a + C.bg(t, r, C.sat); end
    G = R * G * R.';
    Fp = zeros(3, C.ncol);
    Fp(:, 6 + (1:C.np)) = R * P;
    Fp(:, 6 + C.np + (1:3)) = Q;
    Fp(:, 9 + C.np + (1:3)) = Q .* acc.';
    dPsi = [Psi(4:6, :); G * Psi(1:3, :) + Fp];
    dY = [v; a; dPsi(:)];
end

function abm8(rhs, t0, Y0, h, nstep, emit)
%ABM8  Adams-Bashforth-Moulton, order 8, PECE, fixed step h on the grid
%   t0 + k h (k = 0..nstep); the first 7 steps by RK4 at h/20. emit(k, Y) is
%   called at every grid epoch. Coefficients: the exact rationals /120960.
    AB = [434241 -1152169 2183877 -2664477 2102243 -1041723 295767 -36799].' / 120960;
    AM = [36799 139849 -121797 123133 -88547 41499 -11351 1375].' / 120960;
    if nstep < 8
        error('kbrvar:shortArc', 'an arc needs at least 8 integration steps (got %d).', nstep);
    end
    Y = Y0;
    H = zeros(numel(Y0), 8);  H(:, 1) = rhs(t0, Y);
    emit(0, Y);
    sub = 20;  hs = h / sub;
    for k = 1:7
        for j = 0:sub - 1
            t = t0 + (k - 1) * h + j * hs;
            k1 = rhs(t, Y);                 k2 = rhs(t + hs / 2, Y + hs / 2 * k1);
            k3 = rhs(t + hs / 2, Y + hs / 2 * k2);  k4 = rhs(t + hs, Y + hs * k3);
            Y = Y + hs / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
        end
        H = [rhs(t0 + k * h, Y), H(:, 1:7)];
        emit(k, Y);
    end
    for k = 8:nstep
        t = t0 + k * h;
        Yp = Y + h * (H * AB);
        Y = Y + h * ([rhs(t, Yp), H(:, 1:7)] * AM);
        H = [rhs(t, Y), H(:, 1:7)];
        emit(k, Y);
    end
end

function [states, PsiEnd] = propagate(F, rot, accfun, b0, s0, t0, y0, h, nstep, bg)
%PROPAGATE  One satellite: the states [r v] at t0 + k h (k = 0..nstep) and
%   Psi at the last epoch -- the simulator of the tests and of a user's
%   closed loop (the same force model and integrator as buildDay).
    if nargin < 10, bg = []; end
    np = numel(F.est);  ncol = 6 + np + 6;
    C = struct('F', F, 'np', np, 'ncol', ncol, 'rot', rot, 'accfun', accfun, ...
               'b0', b0(:), 's0', s0(:), 'bg', bg, 'sat', '');
    states = zeros(nstep + 1, 6);  PsiEnd = [];
    abm8(@(t, Y) satRhs(t, Y, C), t0, [y0(:); reshape([eye(6), zeros(6, ncol - 6)], [], 1)], ...
         h, nstep, @store);
    function store(k, Y)
        states(k + 1, :) = Y(1:6).';
        if k == nstep, PsiEnd = reshape(Y(7:end), 6, ncol); end
    end
end

function fun = accModel(tAcc, acc, tSca, quat, accToSrf, convention, window, maxGap)
%ACCMODEL  [Q, a] = fun(t): the SRF->CRF rotation from SCA1B and the ACC
%   acceleration in the SRF, interpolated linearly (quaternions sign-aligned,
%   renormalized). A gap longer than maxGap inside the window is refused.
    checkCover(tAcc, window, maxGap, 'ACC1B');
    checkCover(tSca, window, maxGap, 'SCA1B');
    accS = (accToSrf * acc.').';
    for i = 2:size(quat, 1)                               % q and -q are one rotation
        if quat(i, :) * quat(i - 1, :).' < 0, quat(i, :) = -quat(i, :); end
    end
    transposeIt = strcmp(convention, 'crf2srf');
    fun = @evalAt;
    function [Q, a] = evalAt(t)
        a = interp1(tAcc, accS, t, 'linear').';
        q = interp1(tSca, quat, t, 'linear');  q = q / norm(q);
        Q = quat2rot(q);
        if transposeIt, Q = Q.'; end
    end
end

function checkCover(t, window, maxGap, name)
%CHECKCOVER  The samples t must cover window with no gap longer than maxGap.
    s = t(t >= window(1) - maxGap & t <= window(2) + maxGap);
    if isempty(s) || s(1) > window(1) || s(end) < window(2) || any(diff(s) > maxGap)
        error('kbrvar:accGap', '%s does not cover [%.0f %.0f] (GPS s) with gaps <= %g s.', ...
              name, window, maxGap);
    end
end

function R = quat2rot(q)
%QUAT2ROT  Rotation matrix of the unit quaternion q = [q0 q1 q2 q3] (scalar first).
    a = q(1);  b = q(2);  c = q(3);  d = q(4);
    R = [a*a + b*b - c*c - d*d, 2*(b*c - a*d),         2*(b*d + a*c); ...
         2*(b*c + a*d),         a*a - b*b + c*c - d*d, 2*(c*d - a*b); ...
         2*(b*d - a*c),         2*(c*d + a*b),         a*a - b*b - c*c + d*d];
end

function [R, Rd] = eraRotation(tGps)
%ERAROTATION  TRF -> CRF by the Earth rotation angle only (UT1 ~ GPS - 18 s):
%   NO precession-nutation, NO polar motion. A stand-in for synthetic data;
%   real data need the IERS 2010 transformation (see LIMITATIONS).
    du = (tGps - 18) / 86400;
    th = 2 * pi * mod(0.7790572732640 + 1.00273781191135448 * du, 1);
    w = 2 * pi * 1.00273781191135448 / 86400;
    c = cos(th);  s = sin(th);
    R = [c -s 0; s c 0; 0 0 1];
    Rd = w * [-s -c 0; c -s 0; 0 0 0];
end

% ==========================================================================
%  GRAVITY: normalized solid harmonics, exact derivative operators
% ==========================================================================
function [nn, mm, cs] = coeffIndex(N)
%COEFFINDEX  The coefficient order up to degree N: per degree C_n0..C_nn,
%   then S_n1..S_nn; (n, m, cs) with cs 0 = C, 1 = S. Index of C_nm is
%   n^2 + m + 1, of S_nm n^2 + n + m + 1.
    K = (N + 1)^2;
    nn = zeros(K, 1);  mm = nn;  cs = nn;
    for n = 0:N
        i = n^2 + (1:2 * n + 1);
        nn(i) = n;
        mm(i) = [0:n, 1:n];
        cs(i) = [zeros(1, n + 1), ones(1, n)];
    end
end

function [V, W] = harmonics(r, L, R)
%HARMONICS  Vb_nm + i Wb_nm = (R/r)^(n+1) Pbar_nm(sin phi) exp(i m lambda),
%   n, m <= L (row n+1, column m+1), by the Cartesian recursion: no
%   trigonometry, no pole singularity.
    x = r(1);  y = r(2);  z = r(3);
    r2 = x * x + y * y + z * z;  q = R / r2;
    V = zeros(L + 1);  W = zeros(L + 1);
    V(1, 1) = R / sqrt(r2);
    for m = 1:L
        if m == 1, c = sqrt(3); else, c = sqrt((2 * m + 1) / (2 * m)); end
        V(m + 1, m + 1) = c * q * (x * V(m, m) - y * W(m, m));
        W(m + 1, m + 1) = c * q * (x * W(m, m) + y * V(m, m));
    end
    for n = 1:L
        m = 0:n - 1;
        a = sqrt((2 * n - 1) * (2 * n + 1) ./ ((n - m) .* (n + m)));
        V(n + 1, m + 1) = a .* (z * q) .* V(n, m + 1);
        W(n + 1, m + 1) = a .* (z * q) .* W(n, m + 1);
        if n >= 2
            m = 0:n - 2;
            b = sqrt((2 * n + 1) * (n + m - 1) .* (n - m - 1) ./ ((2 * n - 3) * (n + m) .* (n - m)));
            V(n + 1, m + 1) = V(n + 1, m + 1) - b .* (R * q) .* V(n - 1, m + 1);
            W(n + 1, m + 1) = W(n + 1, m + 1) - b .* (R * q) .* W(n - 1, m + 1);
        end
    end
end

function D = operators(N)
%OPERATORS  {Dx, Dy, Dz}: R * d/dx_i of an expansion of degree N, written in
%   the degree-(N+1) basis -- (N+2)^2 x (N+1)^2 sparse, exact.
    [nn, mm, cs] = coeffIndex(N);
    K = numel(nn);
    iC = @(n, m) n.^2 + m + 1;  iS = @(n, m) n.^2 + n + m + 1;
    I = cell(3, 1);  Jc = cell(3, 1);  Vv = cell(3, 1);
    for j = 1:K
        n = nn(j);  m = mm(j);
        up = sqrt((1 + (m > 0)) * (2*n + 1) * (n + m + 1) * (n + m + 2) / (2 * (2*n + 3)));
        same = sqrt((2*n + 1) * (n + m + 1) / ((2*n + 3) * (n - m + 1)));
        if m >= 1
            down = sqrt(2 * (2*n + 1) / ((1 + (m > 1)) * (2*n + 3) * (n - m + 2) * (n - m + 1)));
        else
            down = 0;
        end
        k = (n - m + 2) * (n - m + 1);
        if cs(j) == 0
            if m == 0
                e = {[iC(n+1, 1), -up], [iS(n+1, 1), -up]};
            else
                e = {[iC(n+1, m+1), -up/2; iC(n+1, m-1), k*down/2], ...
                     [iS(n+1, m+1), -up/2; sIdx(n+1, m-1), -k*down/2]};
            end
            e{3} = [iC(n+1, m), -(n - m + 1) * same];
        else
            e = {[iS(n+1, m+1), -up/2; sIdx(n+1, m-1), k*down/2], ...
                 [iC(n+1, m+1), up/2; iC(n+1, m-1), k*down/2], ...
                 [iS(n+1, m), -(n - m + 1) * same]};
        end
        for d = 1:3
            ed = e{d};  ed = ed(ed(:, 1) > 0, :);           % W_{n,0} does not exist
            I{d} = [I{d}; ed(:, 1)];  Jc{d} = [Jc{d}; j * ones(size(ed, 1), 1)];  Vv{d} = [Vv{d}; ed(:, 2)];
        end
    end
    D = cell(1, 3);
    for d = 1:3
        D{d} = sparse(I{d}, Jc{d}, Vv{d}, (N + 2)^2, K);
    end
end

function i = sIdx(n, m)
%SIDX  Index of S_nm, or 0 for m = 0 (no such term).
    if m > 0, i = n^2 + n + m + 1; else, i = 0; end
end

function F = makeField(c, N, est, GM, R)
%MAKEFIELD  A field of degree N (coefficient vector c, coeffIndex order) with
%   its acceleration and gradient expansions and the partial operators of the
%   estimated coefficients est -- built once; fieldEval is one harmonic sweep.
    D1 = operators(N);  D2 = operators(N + 1);
    F.N = N;  F.GM = GM;  F.R = R;  F.est = est(:);  F.c = c(:);
    F.ac = [D1{1} * c(:), D1{2} * c(:), D1{3} * c(:)];                  % (N+2)^2 x 3
    F.gc = [D2{1} * F.ac(:, 1), D2{1} * F.ac(:, 2), D2{1} * F.ac(:, 3), ...
            D2{2} * F.ac(:, 2), D2{2} * F.ac(:, 3), D2{3} * F.ac(:, 3)];   % xx xy xz yy yz zz
    F.Dest = {D1{1}(:, est), D1{2}(:, est), D1{3}(:, est)};
    L = N + 2;
    for d = [1 2]
        [nn, mm, cs] = coeffIndex(N + d);
        lin = sub2ind([L + 1, L + 1], nn + 1, mm + 1);
        F.map{d} = struct('c', find(cs == 0), 'lc', lin(cs == 0), 's', find(cs == 1), 'ls', lin(cs == 1));
    end
end

function [a, G, P] = fieldEval(F, r)
%FIELDEVAL  Acceleration a (3x1), gradient G (3x3) and partials P (3 x np) of
%   the field F at the TRF point r -- one sweep of the harmonics to N+2.
    [V, W] = harmonics(r, F.N + 2, F.R);
    h = cell(1, 2);
    for d = [1 2]
        m = F.map{d};
        h{d} = zeros((F.N + d + 1)^2, 1);
        h{d}(m.c) = V(m.lc);  h{d}(m.s) = W(m.ls);
    end
    s1 = F.GM / F.R^2;  s2 = F.GM / F.R^3;
    a = s1 * (F.ac.' * h{1});
    g = s2 * (F.gc.' * h{2});
    G = [g(1) g(2) g(3); g(2) g(4) g(5); g(3) g(5) g(6)];
    P = s1 * full([h{1}.' * F.Dest{1}; h{1}.' * F.Dest{2}; h{1}.' * F.Dest{3}]);
end

function c = fieldVector(gf, N)
%FIELDVECTOR  The (N+1)^2 coefficient vector of gf.C / gf.S, coeffIndex order.
    [nn, mm, cs] = coeffIndex(N);
    c = zeros(numel(nn), 1);
    for j = 1:numel(nn)
        if cs(j) == 0, c(j) = gf.C(nn(j) + 1, mm(j) + 1); else, c(j) = gf.S(nn(j) + 1, mm(j) + 1); end
    end
end

% ==========================================================================
%  READERS: ICGEM .gfc, GRACE-FO Level-1B ASCII
% ==========================================================================
function gf = readGravity(g)
%READGRAVITY  A .gfc file name, or a struct with GM, R, C, S.
    if isstruct(g)
        if ~all(isfield(g, {'GM', 'R', 'C', 'S'}))
            error('kbrvar:badConfig', 'cfg.gravityField as a struct needs GM, R, C, S.');
        end
        gf = g;
    else
        gf = readGfc(char(g));
    end
end

function gf = readGfc(file)
%READGFC  An ICGEM .gfc file (static part): GM, R and the fully normalized
%   C, S as (N+1) x (N+1) matrices (row n+1, column m+1).
    if ~isfile(file), error('kbrvar:noFile', 'no such file: %s', file); end
    txt = fileread(file);
    head = regexp(txt, 'end_of_head[^\n]*\n', 'end', 'once');
    if isempty(head), error('kbrvar:badFile', '%s: no end_of_head line.', file); end
    hdr = txt(1:head);  body = txt(head + 1:end);
    gf.GM = headerValue(hdr, 'earth_gravity_constant');  gf.R = headerValue(hdr, 'radius');
    N = headerValue(hdr, 'max_degree');
    if any(isnan([gf.GM gf.R N])), error('kbrvar:badFile', '%s: header lacks GM, radius or max_degree.', file); end
    if ~isempty(regexp(body, '(?m)^\s*(gfct|trnd|asin|acos)\s', 'once'))
        error('kbrvar:gfcTimeVariable', '%s carries time-variable terms; give its static field.', file);
    end
    body = strrep(body, 'D', 'E');                         % Fortran exponents
    tok = textscan(body, '%s %f %f %f %f %*[^\n]');
    sel = strcmp(tok{1}, 'gfc');
    gf.C = zeros(N + 1);  gf.S = zeros(N + 1);
    idx = sub2ind([N + 1, N + 1], tok{2}(sel) + 1, tok{3}(sel) + 1);
    gf.C(idx) = tok{4}(sel);  gf.S(idx) = tok{5}(sel);
end

function v = headerValue(hdr, key)
%HEADERVALUE  The number after a keyword of an ICGEM header, NaN if absent.
    tok = regexp(hdr, ['(?m)^\s*' key '\s+(\S+)'], 'tokens', 'once');
    v = NaN;
    if ~isempty(tok), v = str2double(strrep(tok{1}, 'D', 'E')); end
end

function o = readL1B(file, fm)
%READL1B  One GRACE-FO Level-1B ASCII file: the YAML header is skipped, the
%   columns read by fm.fmt, and every field of fm (but fmt) returned as the
%   numeric columns it names.
    if ~isfile(file), error('kbrvar:noFile', 'no such file: %s', file); end
    txt = fileread(file);
    e = regexp(txt, 'End of YAML header[^\n]*\n', 'end', 'once');
    if isempty(e), error('kbrvar:badFile', '%s: no ''End of YAML header'' line.', file); end
    C = textscan(txt(e + 1:end), fm.fmt);
    o = struct();
    for f = setdiff(fieldnames(fm).', {'fmt'})
        o.(f{1}) = [C{fm.(f{1})}];
    end
    if isempty(o.t), error('kbrvar:badFile', '%s: no data rows.', file); end
end
