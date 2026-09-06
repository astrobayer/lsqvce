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
%                           each either 'satKey|obsType' (as LSQVCE.FIT
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
%   NOTE: LSQVCE.FIT drives LSQVCE.ENGINE.ADJUSTVCEJOINT through this helper;
%   it is public for callers driving the joint estimator themselves. The
%   override carries the temporal law and a level for every stream; a
%   stream's per-epoch error information (the day's Cep) composes with it in
%   the loader, so multi-axis and scalar streams need no special case here
%   (R6 WP-90: the former 'Separable' option is gone).
%
%   ERRORS
%     lsqvce:jointModelOpt:badModel   the model spec is unrecognised
%     lsqvce:jointModelOpt:badStream  the stream key is unknown
%     lsqvce:jointModelOpt:badC       the cross-covariance has the wrong shape
%
%   See also LSQVCE.ENGINE.ADJUSTVCEJOINT, LSQVCE.ENGINE.COVMODEL, LSQVCE.ENGINE.MAKELOADER
%
%   Updated 2026-09-06 (CEST, UTC+02:00), R6 WP-90: the 'Separable' option
%   and the AR form of the law retired -- the loader composes the PSD law
%   with the arc's Cep. With assistance from Claude (Anthropic), model
%   Claude Fable 5.1.
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

    % ---- the covariance spec, shared by every stream -----------------------
    % R5 M3 3b (finding H-17): the law leaves here in the PSD DOMAIN, on
    % normalised frequency (the lags are sample lags): S(fNorm) =
    % 2 [C_0 + 2 sum_n C_n cos(2 pi fNorm n)] on L points 0 .. 0.5, floored
    % at eps*max -- positive by construction; covModel places it on the arc's
    % physical axis with the arc's dt and resamples it to whatever grid the
    % arc's length needs. Before 3b the override was a truncated ACF, which is
    % a non-PD Toeplitz for every arc longer than NumLags (the G-08 class).
    % The pooled C is a SHAPE (normalised to C(1) = 1); every stream keeps its
    % declared level times its own factor (sigmaScale, applied by
    % applyCovOverride) -- one global level would erase the scale difference
    % between, say, range-rate and positions (found on the demo: RCOND 1e-50).
    spec = lawSpecs(C, o.Base);
    % R6 WP-90: the override carries the PSD law for EVERY 'law' stream; the
    % LOADER composes it with the per-epoch Cep the arc carries (multi-axis
    % or scalar) -- the former 'Separable' AR form is gone, the decision sits
    % where the data is.
    perStream = [];
    if isfield(model, 'streamSigma2'), perStream = model.streamSigma2; end

    % ---- pack into opt.covOverride.bySys.(satKey).(obsType) (convention 2) --
    % R5 M3 3d: per-stream POLICY (model.Estimate: 'none' -> no override at
    % all; 'scale' -> the level only, on the declared law; 'law' -> the
    % PSD shape of the stream's pool group (model.Cgroups / model.Pool)
    % plus the level).
    E = [];  Cg = [];  pool = 'all';
    if isfield(model, 'Estimate'), E = model.Estimate; end
    if isfield(model, 'Cgroups'),  Cg = model.Cgroups; end
    if isfield(model, 'Pool'),     pool = lower(char(model.Pool)); end
    bySys = struct();
    for k = 1:numel(streams)
        [satKey, ot] = splitStreamKey(streams{k});
        key = sprintf('%s|%s', satKey, ot);
        pol = 'law';
        if ~isempty(E) && isKey(E, key), pol = lower(char(E(key))); end
        if strcmp(pol, 'none'), continue, end            % untouched: no override
        if ~isfield(bySys, satKey), bySys.(satKey) = struct(); end
        lvl = o.Sigma;                                   % relative to the declared sigma
        if ~isempty(perStream) && isKey(perStream, key)
            lvl = lvl * sqrt(max(perStream(key), eps));
        elseif s2 ~= 1
            lvl = lvl * sqrt(max(s2, eps));              % no per-stream factor: the mean
        end
        if strcmp(pol, 'scale')
            bySys.(satKey).(ot) = struct('sigmaScale', lvl);
            continue
        end
        % 'law': the group's C when a per-group estimate exists, else the pooled C
        switch pool
            case 'obstype', gk = ot;
            case 'stream',  gk = key;
            otherwise,      gk = 'all';
        end
        sp = spec;
        if ~isempty(Cg) && isKey(Cg, gk) && ~strcmp(gk, 'all')
            sp = lawSpecs(Cg(gk), o.Base);
        end
        sp.sigmaScale = lvl;
        bySys.(satKey).(ot) = sp;
    end

    uo = struct('covOverride', struct('bySys', bySys), ...
                'psdMode', lower(char(o.PsdMode)));      % convention 3
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

function spec = lawSpecs(C, base)
%LAWSPECS  The PSD-shape spec of one covariance function C (sample lags):
%   S(fNorm) = 2 [C_0 + 2 sum_n C_n cos(2 pi fNorm n)] on L points 0..0.5,
%   normalised to C(1) = 1, floored at eps*max.
    C = C(:);
    L = numel(C);
    Cn = C / max(C(1), eps);
    fNorm = (0:L-1).' / (2*(L-1));
    n = 1:L-1;
    S = 2 * (Cn(1) + 2 * (cos(2*pi*fNorm*n) * Cn(2:end)));
    S = max(S, eps * max(S));
    spec = base;
    spec.type  = 'psd';
    spec.fNorm = fNorm;
    spec.S     = S;
end
