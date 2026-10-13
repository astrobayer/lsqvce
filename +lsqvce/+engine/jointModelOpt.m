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
%                             .S       L x 1     the current law: the bin
%                                      spectrum (>= 0) on f_j = j / (2(L-1))
%                                      cycles per sample, LSQVCE.KERNEL.HATACV;
%                                      .Sgroups the same per pool group
%                             .sigma2  M x 1     current arc variance factors
%                             .iter    1x1       iteration number -> uo.pass
%                                      (the local priors' 'from' reads it)
%                             .rowWeights containers.Map (optional) the robust
%                                      pass's per-epoch factors -> uo.epochWeights
%                                      (the loader's nugget, before the law)
%                             .arcFactor containers.Map (optional)
%                                      LSQVCE.UTIL.ROWWEIGHTKEY -> the arc's
%                                      variance factor RELATIVE to its stream's
%                                      level (R7-13); emitted as
%                                      covOverride.byArc (sqrt: a std)
%                           On iter == 1 .S may be the starting guess and
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
%                         .covOverride.byArc  (when model.arcFactor is given)
%                         .psdMode
%                         .pass          model.iter (1 without it)
%                         .epochWeights  model.rowWeights, when it holds any
%                       WP-134: .pass and .epochWeights were added by LSQVCE.FIT
%                       alone, so an engine caller following the EXAMPLE below
%                       never activated a local prior with from > 1 and never
%                       applied the robust factors (while the loop still divided
%                       the levels by their truncation factor)
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
%     lsqvce:jointModelOpt:badS       the law is not a finite non-negative spectrum
%
%   See also LSQVCE.ENGINE.ADJUSTVCEJOINT, LSQVCE.ENGINE.COVMODEL, LSQVCE.ENGINE.MAKELOADER
%
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-134: the pass number and the
%   robust factors ride on the UserOpt here (uo.pass, uo.epochWeights), the one
%   place for every caller; LSQVCE.FIT's withPass is gone. With assistance from
%   Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-111 part 3 (R7-02): the law
%   arrives as the bin spectrum model.S / .Sgroups and leaves as the 'psd'
%   spec with ordinates 2 S on fNorm = f_j -- exactly the hat-basis law, with
%   no cosine transform of a truncated ACF and no floor (badC -> badS). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-116 part 2 (R7-13): the
%   per-arc factors reach the loader (covOverride.byArc, a relative standard
%   deviation per arc on top of the stream's level), so they weight the
%   solve. With assistance from Claude (Anthropic), model Claude Opus 5.5.
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

    if ~isstruct(model) || ~isfield(model, 'S')
        error('lsqvce:jointModelOpt:badModel', ...
             ['model must be the struct LSQVCE.ENGINE.ADJUSTVCEJOINT passes to accFcn ' ...
              '(fields .S, .sigma2, .iter).']);
    end
    S = checkedLaw(model.S);

    % ---- the overall scale: VARIANCE -> STANDARD DEVIATION (convention 1) ---
    %   model.sigma2 holds per-arc VARIANCE factors. A covariance spec wants a
    %   standard deviation, so the scale is sqrt of the (mean) variance. The
    %   arc-to-arc detail rides separately (model.arcFactor -> byArc, below).
    s2 = 1;
    if isfield(model, 'sigma2') && ~isempty(model.sigma2)
        v = model.sigma2(:);
        v = v(isfinite(v) & v > 0);
        if ~isempty(v), s2 = mean(v); end
    end

    % ---- the covariance spec, shared by every stream -----------------------
    % R5 M3 3b (finding H-17) / R7 WP-111 part 3: the law leaves here in the
    % PSD DOMAIN, on normalised frequency (the lags are sample lags): the bin
    % spectrum S itself, as covModel's ONE-sided ordinates 2 S on fNorm = f_j;
    % covModel places it on the arc's physical axis with the arc's dt and
    % LSQVCE.KERNEL.PSD2ACV integrates it exactly -- r = K S at every lag,
    % positive semi-definite for every arc length. (3b sent the cosine
    % transform of a truncated ACF, floored at eps*max: the floor is what made
    % the declared covariance near-singular, review R7-02.)
    % The pooled law is a SHAPE (normalised to r(0) = 1); every stream keeps its
    % declared level times its own factor (sigmaScale, applied by
    % applyCovOverride) -- one global level would erase the scale difference
    % between, say, range-rate and positions (found on the demo: RCOND 1e-50).
    spec = lawSpecs(S, o.Base);
    % R6 WP-90: the override carries the PSD law for EVERY 'law' stream; the
    % LOADER composes it with the per-epoch Cep the arc carries (multi-axis
    % or scalar) -- the former 'Separable' AR form is gone, the decision sits
    % where the data is.
    perStream = [];
    if isfield(model, 'streamSigma2'), perStream = model.streamSigma2; end

    % ---- pack into opt.covOverride.bySys.(satKey).(obsType) (convention 2) --
    % R5 M3 3d: per-stream POLICY (model.Estimate: 'none' -> no override at
    % all; 'scale' -> the level only, on the declared law; 'law' -> the
    % PSD shape of the stream's pool group (model.Sgroups / model.Pool)
    % plus the level).
    E = [];  Sg = [];  pool = 'all';
    if isfield(model, 'Estimate'), E = model.Estimate; end
    if isfield(model, 'Sgroups'),  Sg = model.Sgroups; end
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
        % 'law': the group's S when a per-group estimate exists, else the pooled S
        switch pool
            case 'obstype', gk = ot;
            case 'stream',  gk = key;
            otherwise,      gk = 'all';
        end
        sp = spec;
        if ~isempty(Sg) && isKey(Sg, gk) && ~strcmp(gk, 'all')
            sp = lawSpecs(checkedLaw(Sg(gk)), o.Base);
        end
        sp.sigmaScale = lvl;
        bySys.(satKey).(ot) = sp;
    end

    ovr = struct('bySys', bySys);
    % R7-13: the per-arc factors WEIGHT the solve (Ellmer's convention). Each
    % is relative to its stream's level above, so the arc's variance is level
    % x factor -- the number eq. 6.5.28 multiplies in adjustVCEjoint. Before,
    % only the stream mean reached the loader and a bad arc kept its stream's
    % weight (review R7-13: one arc 10x in sigma, formal/true sd 2.6).
    if isfield(model, 'arcFactor') && ~isempty(model.arcFactor) && model.arcFactor.Count > 0
        ak = keys(model.arcFactor);
        av = cellfun(@(k) sqrt(max(model.arcFactor(k), eps)), ak);
        ovr.byArc = containers.Map(ak, num2cell(av));
    end
    uo = struct('covOverride', ovr, ...
                'psdMode', lower(char(o.PsdMode)));      % convention 3
    % WP-134: the pass rides on the UserOpt (the local priors' 'from',
    % makeLoader>priorArc), and so do the robust factors of the previous pass
    % (R6 WP-105: model.rowWeights -> UserOpt.epochWeights, the additive nugget
    % makeLoader>batchOf composes BEFORE the temporal law -- the Koch convention)
    uo.pass = 1;
    if isfield(model, 'iter') && ~isempty(model.iter), uo.pass = model.iter; end
    if isfield(model, 'rowWeights') && ~isempty(model.rowWeights) && model.rowWeights.Count > 0
        uo.epochWeights = model.rowWeights;
    end
end

%==========================================================================
function [satKey, ot] = splitStreamKey(s)
%SPLITSTREAMKEY  Accept 'satKey|obsType' or {sysID, obsType}.
    if iscell(s) && numel(s) == 2
        [~, satKey] = lsqvce.util.streamKey(s{1});
        ot     = char(string(s{2}));
        return
    end
    parts = split(string(s), '|');
    if numel(parts) ~= 2
        error('lsqvce:jointModelOpt:badStream', ...
             ['A stream must be ''satKey|obsType'' or {sysID, obsType} ' ...
              '(got ''%s'').'], char(string(s)));
    end
    [~, satKey] = lsqvce.util.streamKey(parts(1));
    ot     = char(parts(2));
end

function S = checkedLaw(S)
%CHECKEDLAW  The bin spectrum as a column, refused unless finite, >= 0, not all zero.
    S = S(:);
    if numel(S) < 2 || ~all(isfinite(S)) || any(S < 0) || ~any(S > 0)
        error('lsqvce:jointModelOpt:badS', ...
              'model.S must be a finite non-negative bin spectrum (>= 2 bins, not all zero).');
    end
end

function spec = lawSpecs(S, base)
%LAWSPECS  The 'psd' spec of the bin spectrum S: covModel's ONE-sided
%   ordinates 2 S on fNorm = j / (2(L-1)) (LSQVCE.KERNEL.HATACV's convention),
%   scaled to r(0) = 1 (the law is a shape; the level is sigmaScale), merged
%   into base. r(0) = d sum_j c_j S_j, d = 1/(2(L-1)), c_j = 2 (1 at the ends).
    L = numel(S);
    c = 2 * ones(L, 1);  c([1 L]) = 1;
    r0 = sum(c .* S) / (2*(L-1));
    spec = base;
    spec.type  = 'psd';
    spec.fNorm = (0:L-1).' / (2*(L-1));
    spec.S     = 2 * S / r0;
end
