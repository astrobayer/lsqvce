function uo = jointModelOpt(model, streams, opt)
%LSQVCE.ENGINE.JOINTMODELOPT  Turn a joint-mode model struct into loader UserOpt.
%
%   uo = LSQVCE.ENGINE.JOINTMODELOPT(model, streams)
%   uo = LSQVCE.ENGINE.JOINTMODELOPT(model, streams, 'Name', value, ...)
%
%   The missing piece of the LSQVCE.ENGINE.ADJUSTVCEJOINT contract. That driver hands
%   the caller's accFcn a MODEL and expects the caller to weight the data with
%   it; this function does that translation, so a caller writes
%
%       accFcn = @(model) myBuild(cfg, lsqvce.engine.jointModelOpt(model, streams));
%
%   instead of hand-rolling the override and getting the conventions wrong.
%
%   THREE CONVENTIONS THAT ARE EASY TO GET WRONG (all encoded here)
%     1. VARIANCE vs STANDARD DEVIATION. model.sigma2 holds VARIANCE factors,
%        but a covariance spec's 'sigma' field is a STANDARD DEVIATION. The
%        scale therefore enters as sqrt(sigma2), not sigma2. Getting this wrong
%        does not error -- the loop simply diverges, which is exactly how it was
%        first found (a test fixture scaled by sigma2 and the VCE ran away).
%     2. WHERE THE OVERRIDE LIVES. Loaders read it from
%        opt.covOverride.bySys.(satKey).(obsType), with satKey formed by
%        matlab.lang.makeValidName on the system id -- NOT from the covariance
%        spec directly, and NOT from a flat per-obsType field alone.
%     3. psdMode IS READ FROM THE UserOpt, not from the covariance spec. Putting
%        'psdMode' inside the spec has NO effect. It matters: the default is
%        'ar', which fits an AR to the model and can produce a non-stationary
%        polynomial; 'toeplitz' builds the covariance directly (Ellmer eqs.
%        6.5.6/6.5.17) and is what joint mode wants.
%
%   INPUT
%     model    1x1 struct   as handed to accFcn by LSQVCE.ENGINE.ADJUSTVCEJOINT:
%                             .C       nLag x 1  current covariance function
%                             .sigma2  M x 1     current arc variance factors
%                             .iter    1x1       iteration number
%                           On iter == 1 .C may be the starting guess and
%                           .sigma2 may be empty; both are handled.
%     streams  1xS cell     the (sysID, obsType) pairs the override must cover,
%                           each either 'satKey|obsType' (as LSQVCE.ENGINE.ADJUSTVCE
%                           names them) or a 1x2 cell {sysID, obsType}.
%
%   OPTIONS (Name/Value)
%     'PsdMode'  char    'toeplitz' (default) | 'ar' | 'arma'. See convention 3.
%     'Sigma'    1x1     an extra overall standard deviation applied to every
%                        stream (default 1), for callers that carry a separate
%                        global scale.
%     'Base'     1x1 struct  extra fields merged into every stream's spec, e.g.
%                        struct('type','psd') if the loader needs a type hint.
%
%   OUTPUT
%     uo   1x1 struct   a UserOpt ready to hand to LSQVCE.ENGINE.BLOCKACCUMULATOR:
%                         .covOverride.bySys.(satKey).(obsType) = spec
%                         .psdMode
%
%   EXAMPLE
%     streams = {'x1|posC', 'x1|posD', 'x1|rr', 'x2|posE'};
%     accFcn  = @(model) buildAcc(cfg, lsqvce.engine.jointModelOpt(model, streams));
%     out     = lsqvce.engine.adjustVCEjoint(accFcn, 'NumLags', 60);
%
%     function acc = buildAcc(cfg, uo)
%         acc = lsqvce.engine.BlockAccumulator(cfg.nc, 'KeepLocal', true, 'UserOpt', uo);
%         for s = cfg.sats
%             acc.addSystem(s, cfg.nl, cfg.loader(s), cfg.days);
%         end
%     end
%
%   NOTE: the 'vceMode','joint' path on LSQVCE.ENGINE.ADJUSTVCE does this internally
%   (see its covOverrideOpt), so this helper is for callers driving
%   LSQVCE.ENGINE.ADJUSTVCEJOINT themselves.
%
%   ERRORS
%     lsqvce:jointModelOpt:badModel   the model spec is unrecognised
%     lsqvce:jointModelOpt:badStream  the stream key is unknown
%     lsqvce:jointModelOpt:badC       the cross-covariance has the wrong shape
%
%   See also LSQVCE.ENGINE.ADJUSTVCEJOINT, LSQVCE.ENGINE.ADJUSTVCE, LSQVCE.ENGINE.COVMODEL
%
%   Updated 2026-09-01 (CEST, UTC+02:00), R4 G-08 (Level L3): the acf lookup
%   converts covModel's TIME lags to SAMPLE lags with dt taken from the
%   batch's own lag matrix; before, dt ~= 1 s gave a wrongly sampled, then
%   truncated (indefinite) covariance. With assistance from Claude
%   (Anthropic), model Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (b) b2 (F-07):
%   typed `arguments` signature, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    % WP-39 stage (b) b2: typed `arguments` signature (canonical
    % UpperCamelCase; legacy casings keep matching -- MATLAB NV
    % semantics; enum VALUES case-sensitive lowercase).
arguments
    model
    streams
    opt.PsdMode (1,1) string {mustBeMember(opt.PsdMode, ["toeplitz","ar","arma"])} = "toeplitz"
    opt.Sigma (1,1) double {mustBePositive} = 1
    opt.Base (1,1) struct = struct()
    opt.Separable cell = {}                 % R5 M3 (H-15): stream keys whose base law is
                                            % 'separable'; they get an AR(p) TEMPORAL law
                                            % fitted to C (Levinson) and KEEP their type,
                                            % axis count and per-epoch Cep
    opt.ArOrder (1,1) double {mustBeInteger, mustBePositive} = 20
end

    o = struct('PsdMode', char(opt.PsdMode), 'Sigma', opt.Sigma, 'Base', opt.Base);

    if ~isstruct(model) || ~isfield(model, 'C')
        error('lsqvce:jointModelOpt:badModel', ...
             ['model must be the struct LSQVCE.ENGINE.ADJUSTVCEJOINT passes to accFcn ' ...
              '(fields .C, .sigma2, .iter).']);
    end
    C = model.C(:);
    if isempty(C) || ~all(isfinite(C))
        error('lsqvce:jointModelOpt:badC', ...
              'model.C must be a finite covariance function.');
    end

    % ---- the overall scale: VARIANCE -> STANDARD DEVIATION (convention 1) ---
    %   model.sigma2 holds per-arc VARIANCE factors. A covariance spec wants a
    %   standard deviation, so the scale is sqrt of the (mean) variance. The
    %   arc-to-arc detail cannot be expressed in a per-stream override, so the
    %   mean is used and the arc factors remain the driver's business.
    s2 = 1;
    if isfield(model, 'sigma2') && ~isempty(model.sigma2)
        v = model.sigma2(:);
        v = v(isfinite(v) & v > 0);
        if ~isempty(v), s2 = mean(v); end
    end
    sig = o.Sigma * sqrt(max(s2, eps));

    % ---- the covariance spec, shared by every stream -----------------------
    %   'acf' takes a handle over a LAG MATRIX, which is how covModel calls it.
    spec = o.Base;
    spec.type  = 'acf';
    spec.sigma = sig;
    spec.acf   = makeAcf(C);
    % R5 M3 (H-15): the SEPARABLE form of the same law -- an AR(p) fitted to C
    % by Levinson, innovation std scaled by sig; the base spec's type, axis
    % count d and per-epoch Cep survive the merge because they are not set here
    specSep = [];
    if ~isempty(opt.Separable)
        p = min(opt.ArOrder, numel(C) - 1);
        [ar, s2i] = lsqvce.kernel.levinsonDurbin(C(:) / max(C(1), eps), p);
        specSep = struct('ar', reshape(ar, 1, []), 'sigma', sig * sqrt(max(s2i, eps)));
    end
    sepKeys = cellfun(@(s) char(string(s)), opt.Separable, 'UniformOutput', false);

    % ---- pack into opt.covOverride.bySys.(satKey).(obsType) (convention 2) --
    bySys = struct();
    for k = 1:numel(streams)
        [satKey, ot] = splitStreamKey(streams{k});
        if ~isfield(bySys, satKey), bySys.(satKey) = struct(); end
        if ~isempty(specSep) && any(strcmp(sepKeys, sprintf('%s|%s', satKey, ot)))
            bySys.(satKey).(ot) = specSep;
        else
            bySys.(satKey).(ot) = spec;
        end
    end

    uo = struct('covOverride', struct('bySys', bySys), ...
                'psdMode', lower(char(o.PsdMode)));      % convention 3
end

%==========================================================================
function fh = makeAcf(c)
%MAKEACF  Wrap a discrete autocovariance vector as an acf(lagMatrix, opt) handle.
%   Lags beyond the modelled range return zero -- the model simply says "no
%   correlation out there" rather than extrapolating.
    c = c(:);
    L = numel(c);
    fh = @(lag, ~) acvLookup(c, L, lag);
end

function C = acvLookup(c, L, lag)
%ACVLOOKUP  c at the batch's lags in SAMPLE units, zero beyond the modelled range.
%   covModel's 'acf' contract hands over TIME differences (t(idx) - t(idx)',
%   seconds); the driver's C is per SAMPLE lag (the DCT of a PSD over sample
%   lags). The sampling interval is the smallest non-zero |lag| of the batch
%   itself (regular sampling), so no caller has to know dt (R4 G-08: the old
%   round(|lag|)+1 read every dt-th entry for dt ~= 1 s and returned 0 past
%   index L -- an indefinite Toeplitz matrix, every arc skipped as nonPD).
    a  = abs(lag);
    nz = a(a > 0);
    if isempty(nz), dt = 1; else, dt = min(nz(:)); end
    k  = round(a / dt) + 1;            % lag 0 -> index 1
    C  = zeros(size(k));
    in = k >= 1 & k <= L;
    C(in) = c(k(in));
end

%==========================================================================
function [satKey, ot] = splitStreamKey(s)
%SPLITSTREAMKEY  Accept 'satKey|obsType' or {sysID, obsType}.
    if iscell(s) && numel(s) == 2
        satKey = matlab.lang.makeValidName(char(string(s{1})));
        ot     = char(string(s{2}));
        return
    end
    parts = split(string(s), '|');
    if numel(parts) ~= 2
        error('lsqvce:jointModelOpt:badStream', ...
             ['A stream must be ''satKey|obsType'' or {sysID, obsType} ' ...
              '(got ''%s'').'], char(string(s)));
    end
    satKey = matlab.lang.makeValidName(char(parts(1)));
    ot     = char(parts(2));
end
