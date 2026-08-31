function [cfg, data] = demo_grace_dataLoader()
%DEMO_GRACE_DATALOADER  Synthetic GRACE-FO data generator for the VCE demo.
%
%   [cfg, data] = demo_grace_dataLoader()
%
%   Returns a ready-to-run configuration CFG and the declaration as a v3 data
%   struct DATA (DATA_CONTRACT.md), so the separable VCE loop can be
%   exercised end to end without real data:
%
%       [cfg, data] = demo_grace_dataLoader();
%       R = lsqvce.fit(data, lsqvce.Options(Model="vce", GroupBy="obsType"));
%
%   TOPOLOGY (real GRACE-FO, "Option B", plus an independent third satellite).
%     System 1  -- the C/D PAIR.  TWO satellites C, D each with a 3-D position
%       stream (posC, posD, separable covariance C = D (R_t kron I_3) D with a
%       FIXED per-epoch 3x3 D), and ONE inter-satellite range-rate stream (rr,
%       scalar) whose design couples BOTH satellites' arc states. Because rr
%       touches two satellites' local parameters and the block-arrow model has
%       one local block per system, C and D are MERGED into one pair-system with
%       a stacked local block [y_C | y_D].
%     System 2  -- satellite E, INDEPENDENT.  A third satellite with ONLY a 3-D
%       position stream (posE) and NO inter-satellite link, so it is a clean,
%       separate block-arrow system with its own local (arc) parameters y_E. It
%       observes the SAME common gravity field x_c as the pair (that is what
%       makes combining the systems meaningful), but its locals are its own.
%
%   Both systems share the common (gravity) parametrisation; VCE groups by
%   obsType, giving FOUR variance components (posC, posD, rr, posE). Each stream
%   starts from a deliberately WRONG a-priori (temporally white, sigma 1).
%
%   PARAMETRISATION PIPE (v3 contract).  The unknowns are described by two
%   cell arrays with one ROW per unknown (data.commonParams and the system's
%   localParams), piped through to the day-builder as its 3rd/4th args; nc
%   and nl are DERIVED from their row counts, and the column count is YOUR
%   convention (this demo uses 14). A 'processingOptions' struct is piped as
%   the 5th arg. The demo's column convention (yours may differ) is:
%
%     col  1  kind            'gravity' | 'arcBias' | ...   (char)
%     col  2  label           human-readable name          (char)
%     col  3  system          system id this unknown belongs to (numeric; 0 =
%                             common/shared)
%     col  4  satellite       'C'|'D'|'E'|'' ('' = shared)  (char)
%     col  5  degree  n       SH degree   (gravity) or []   (numeric)
%     col  6  order   m       SH order    (gravity) or []   (numeric)
%     col  7  cs              'C'|'S' coeff (gravity) or ''  (char)
%     col  8  arc             arc index (local) or []        (numeric)
%     col  9  component       'x'|'y'|'z'|'' axis or ''      (char)
%     col 10  unit            e.g. '-'|'m'|'m/s'             (char)
%     col 11  a-priori value  x0                             (numeric)
%     col 12  a-priori sigma  or Inf if free                 (numeric)
%     col 13  estimate flag   true/false                     (logical)
%     col 14  group tag       free grouping label            (char)
%
%   The toolbox only checks the SHAPE (a 2-D cell); the MEANINGS above are
%   this demo's convention, read inside demoBuildDay. This is the TEMPLATE:
%   copy demoBuildDay / demoEpochCov and the two paramTable builders as the
%   starting point for your own case. The epochCovariance hook each system
%   declares is what gives the VCE Shape step its exact separable-descale
%   path (DATA_CONTRACT.md section 7).
%
%   See also LSQVCE.FIT, DEMO_GRACE, DEMO_SEPARABLE
%
%   Updated 2026-08-31 (CEST, UTC+02:00), v3 data contract (Level L10,
%   WP-74): USER surface only -- returns [cfg, data], the v3 declaration
%   (obsTypes from the covariance fieldnames, localParams tables, id,
%   epochCovariance hooks); the ENGINE-tier loader construction moved to
%   dev/tests/fixture_grace.m. With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   Created 2026-07-16 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.
    cfg  = demoCfg();
    data = graceData(cfg);
end

%==========================================================================
function data = graceData(cfg)
%GRACEDATA  The declaration as a v3 data struct: two systems (the merged
%   C/D pair and independent E), each with its covariance (the fieldnames
%   ARE the streams, in processing order), its localParams table, and the
%   epochCovariance hook for the exact VCE Shape path.
    data = struct();
    data.commonParams = cfg.common;                   % nc = rows (14 cols, OUR convention)
    data.days = cfg.days;

    s = struct();
    s.buildDay          = @(d,ot,common,local,po) demoBuildDay(cfg, d, ot, common, local, po);
    s.covariance        = whiteCov(cfg.obsTypes);     % posC, posD, rr -- wrong on purpose
    s.arcs              = cfg.arcsPerDay;
    % NOTE: id must survive [store.sysID] concatenation as ONE element --
    % numeric or a single char, never a multi-char string like 'CD'.
    s.id                = 1;
    s.name              = 'demo-CD-pair';
    s.localParams       = cfg.localPair;
    s.processingOptions = cfg.procOpt;
    s.epochCovariance   = @(d,ot) demoEpochCov(cfg, d, ot);
    data.systems = s;

    s.covariance  = whiteCov(cfg.obsTypesE);          % posE only
    s.id          = 2;
    s.name        = 'demo-E';
    s.localParams = cfg.localE;
    data.systems(2) = s;
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

    % ---- parametrisation tables (14 columns, one row per unknown) ----------
    cfg.common = commonParamTable(cfg);      % nc rows -> nc derived
    cfg.localPair = localParamTable(cfg, {'C','D'}, [cfg.y_C, cfg.y_D]); % 2*nlSat
    cfg.localE    = localParamTable(cfg, {'E'},     cfg.y_E);            % nlSat
    cfg.procOpt   = struct('mEpoch',cfg.mEpoch, 'dt',cfg.dt, ...
                           'arcsPerDay',cfg.arcsPerDay, 'seed',cfg.seed);
end

function cov = whiteCov(obsTypes)
%WHITECOV  A-priori white covariance per obsType (separable for positions).
    cov = struct();
    for k = 1:numel(obsTypes)
        ot = obsTypes{k};
        if demoAxes(ot) > 1
            cov.(ot) = struct('type','separable', 'ar',1, 'sigma',1, 'd',3);
        else
            cov.(ot) = struct('type','diag', 'sigma', 1);
        end
    end
end

function D = demoBuildDay(cfg, day, obsType, common, local, procOpt)  %#ok<INUSD>
%DEMOBUILDDAY  One synthetic day, parameter-driven (5-arg buildDay form).
%   common/local are the 14-column parametrisation tables; nc = rows(common),
%   nlThis = rows(local). This demo reads the COUNTS from the tables (proving
%   the pipe works) and could equally read col-5/6 (degree/order) etc. to build
%   the design; here the design is synthetic so we only use the counts and the
%   a-priori truth carried in cfg.
    nc     = size(common, 1);           % DERIVED common count (== cfg.nc)
    nlThis = size(local,  1);           % this system's local count
    nlSat  = cfg.nlSat;

    nArc  = cfg.arcsPerDay;
    nAx   = demoAxes(obsType);
    mE    = cfg.mEpoch;                 % epochs per arc
    m     = nArc * mE * nAx;            % rows per day (epoch-major)

    s = RandStream('mt19937ar', 'Seed', ...
            mod(cfg.seed*1e5 + day*997 + sum(double(obsType))*31, 2^31-1));

    A = A_design(s, m, nc);
    B = zeros(m, nlThis);               % this system's local block
    switch obsType
        case 'posC'                      % pair system: fill C sub-columns
            B(:, 1:nlSat) = B_design(s, m, nlSat, nArc, mE*nAx);
            yContrib = B * [cfg.y_C; zeros(nlSat,1)];
        case 'posD'                      % pair system: fill D sub-columns
            B(:, nlSat+(1:nlSat)) = B_design(s, m, nlSat, nArc, mE*nAx);
            yContrib = B * [zeros(nlSat,1); cfg.y_D];
        case 'rr'                        % pair system: couples BOTH, weakly
            B(:, 1:nlSat)         = 0.1 * B_design(s, m, nlSat, nArc, mE*nAx);
            B(:, nlSat+(1:nlSat)) = 0.1 * B_design(s, m, nlSat, nArc, mE*nAx);
            yContrib = B * [cfg.y_C; cfg.y_D];
        case 'posE'                      % independent system: its own locals
            B(:, 1:nlSat) = B_design(s, m, nlSat, nArc, mE*nAx);
            yContrib = B * cfg.y_E;
        otherwise
            error('demo_grace_dataLoader:obsType', ...
                  'Unknown demo obsType ''%s''.', obsType);
    end

    l0 = A * cfg.x_c + yContrib;
    if nAx > 1
        sigT = cfg.trueSigma.(obsType);
        Ld   = chol(cfg.Dep, 'lower');       % 3x3
        nEd  = mE * nArc;                    % epochs in the day
        e    = zeros(m, 1);
        for a = 1:nArc
            base = (a-1) * mE * nAx;
            W = zeros(nAx, mE);
            for ax = 1:nAx
                w = randn(s, mE, 1);
                W(ax,:) = filter(1, cfg.trueAR.(obsType), w).';  % time colour
            end
            W = sigT * (Ld * W);             % cross-axis mix + scale
            e(base + (1:mE*nAx)) = W(:);      % epoch-major [x1 y1 z1 ...]
        end
        D.Cep = repmat(reshape(cfg.Dep, 1, nAx, nAx), nEd, 1, 1);
    else
        e = colouredNoise(s, m, cfg.trueAR.(obsType), ...
                          cfg.trueSigma.(obsType), nAx, nArc);
    end
    D.A = A;  D.B = B;  D.l = l0 + e;
    D.t = reshape(repmat((0:mE*nArc-1)*cfg.dt, nAx, 1), [], 1);  % epoch-major
    D.arcIdx = arcIndex(nArc, mE*nAx);
    if nAx > 1, D.nAxis = nAx; end       % tell the accumulator the axis count
end

function n = demoAxes(obsType)
%DEMOAXES  Axis count for a demo obsType (posC/posD/posE -> 3, rr -> 1).
    if numel(obsType) >= 3 && strcmpi(obsType(1:3), 'pos'), n = 3; else, n = 1; end
end

function dd = demoEpochCov(cfg, day, obsType)  %#ok<INUSD>
%DEMOEPOCHCOV  Per-day epoch-covariance in the LSQVCE.KERNEL.RESIDUALACV
%   contract: .arc(i).epochs (day-epoch indices) and .Cep (nEpochDay x d x d).
%   Only meaningful for the multi-axis position streams.
    nAx  = demoAxes(obsType);
    mE   = cfg.mEpoch;  nArc = cfg.arcsPerDay;
    nEd  = mE * nArc;                          % epochs in the day
    dd.Cep = repmat(reshape(cfg.Dep, 1, nAx, nAx), nEd, 1, 1);
    dd.arc = struct('epochs', {}, 'rows', {});
    for a = 1:nArc
        ep = (a-1)*mE + (1:mE);                % day-epoch indices
        dd.arc(a).epochs = ep;
        dd.arc(a).rows   = (a-1)*mE*nAx + (1:mE*nAx);
    end
end

%% ===================================================== parametrisation tables
function C = commonParamTable(cfg)
%COMMONPARAMTABLE  nc rows x 14 cols describing the shared gravity unknowns.
%   Synthetic degree/order just enumerate the nc coefficients; a real table
%   would carry the true SH (n,m,C/S). Shared -> system 0, satellite ''.
    nc = cfg.nc;
    C  = cell(nc, 14);
    for i = 1:nc
        n = floor(sqrt(i-1));  m = (i-1) - n*n;      % toy (n,m) enumeration
        C(i,:) = { 'gravity', sprintf('Cnm_%d_%d', n, m), 0, '', ...
                   n, m, 'C', [], '', '-', cfg.x_c(i), Inf, true, 'gravity' };
    end
end

function T = localParamTable(cfg, sats, yTrue)
%LOCALPARAMTABLE  (nlSat*numel(sats)) rows x 14 cols: arc-bias locals per sat.
%   yTrue is nlSat x numel(sats) (or nlSat x 1) of a-priori truth values.
    nlSat = cfg.nlSat;  ns = numel(sats);
    T = cell(nlSat*ns, 14);
    r = 0;
    for is = 1:ns
        sat = sats{is};
        for k = 1:nlSat
            r = r + 1;
            T(r,:) = { 'arcBias', sprintf('bias_%s_%d', sat, k), NaN, sat, ...
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

function idx = arcIndex(nArc, mArc)
    idx = cell(1, nArc);
    for a = 1:nArc
        idx{a} = (a-1)*mArc + (1:mArc);
    end
end
