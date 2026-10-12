function [Cfcn, Wfcn, info] = covModel(cm, t, opt)
%LSQVCE.ENGINE.COVMODEL  Build covariance / whitening operators for one arc.
%
%   [Cfcn, Wfcn, info] = LSQVCE.ENGINE.COVMODEL(cm, t, opt)
%
%   Turns a covariance SPECIFICATION into the two operators the accumulator
%   can use. It returns closures over the arc's epochs, so the caller never
%   handles the model parameters again.
%
%   INPUT
%     cm   struct, .type selects the model:
%       'diag'    .sigma   scalar or m x 1 standard deviation
%       'acf'     .acf     @(dt, opt) -> correlation at lag dt
%                 .sigma   scalar
%       'psd'     .f       nf x 1 frequency grid [Hz], one-sided
%                 .S       nf x 1 PSD [unit^2/Hz]
%       'arma'    .ar      [1 a1 ... ap]
%                 .ma      [1 b1 ... bq]   (optional, default 1)
%                 .sigma   innovation std
%       'var'     .ar      d x d x p coefficient matrices (vector AR)
%                 .sigma   d x d innovation covariance
%                 .acv     d x d x >=p Gamma(0..), optional (else computed
%                          from ar/sigma via LSQVCE.KERNEL.VARACV). Rows epoch-major;
%                          the arc must set nAxis = d. The 'needs-var' case
%                          of LSQVCE.KERNEL.CROSSCORR; fit with LSQVCE.KERNEL.WHITTLE.
%                          Takes NO Cep (it carries its own cross-axis structure).
%       'custom'  .fcn     @(t_batch, opt) -> dense covariance
%
%       Every type above is a TEMPORAL law. ONE optional field composes it
%       with the per-epoch ERROR INFORMATION of the arc (R6 WP-90):
%         .Cep     nEpoch x 1 | nEpoch x d | nEpoch x d x d   the per-epoch
%                  standard deviation of a scalar stream, the per-axis
%                  standard deviations (a diagonal D), or the full per-epoch
%                  covariance (a PPP solution). Then C = D (R_t kron I_d) D:
%                  every epoch is descaled by D^-1/2 and the temporal whitener
%                  runs down the epochs of each axis (LSQVCE.KERNEL.SEPARABLEWHITEN's
%                  order); Cfcn is the TEMPORAL covariance over epochs; the
%                  law is built on the EPOCH grid t(1:d:end); filterMemory is
%                  p*d rows. The loader sets it from the day buildDay returned
%                  (the API reference, chapter buildDay, Cep); a user never writes it
%                  on the spec. d = size(Cep, 2), and a .d on the spec must
%                  agree. The former 'separable' type was 'arma' plus this field.
%         .nugget  nEpoch x 1, >= 0   (R7 WP-115 part 3b; with Cep only) the
%                  joint estimator's robust factors f as an ADDITIVE nugget,
%                  q = 1/f - 1: every axis of epoch e gains q_e * Sigma_ee on
%                  the diagonal of its temporal covariance, so the epoch's
%                  marginal variance is Sigma_ee / f_e (as the multiplicative
%                  fold gave) while its covariances with the neighbours stay
%                  those of the law. The fold scaled them by 1/sqrt(f) too,
%                  and under a correlated law a zero-weight epoch then became
%                  a phantom observation for its neighbours (formal/true sd
%                  0.31 at phi 0.9, HANDOVER 12.43). The whitener is EXACT
%                  and needs the whole arc (info.wholeArc): the law whitener
%                  F, then a rank-k correction G = I - Q V diag(g) V' Q' over
%                  the k nugget epochs (Q R = qr(F U), Woodbury), O(m k + k^3).
%                  Cfcn stays the LAW; the nugget travels as info.nugget for
%                  the joint sweep, which composes it itself.
%
%     t    m x 1  epochs of the current arc [s] (epoch-major with Cep: every
%                 d-th value is one epoch)
%     opt  struct runtime options (see below)
%
%   OPTIONS (opt)
%     psdMode      'ar' (default) | 'arma' | 'toeplitz'   -- for cm.type 'psd'
%     arOrder      AR order p                             (default 30)
%     maOrder      MA order q, only if psdMode = 'arma'   (default 0)
%     psdFloor     relative floor applied to the PSD      (default 1e-10)
%     whitenMode   'auto' (default) | 'filter' | 'chol'
%     BatchSize, BatchOverlap                             (used by 'toeplitz')
%
%   OUTPUT
%     Cfcn  @(idx, opt) -> dense numel(idx) x numel(idx) covariance,
%           or [] if the model has no dense form worth building.
%     Wfcn  @(idx, A, B, l, opt) -> whitened A, B, l. Applies an operator F
%           with F'*F = inv(C). [] if no filter form exists.
%     info  .filterMemory  required BatchOverlap (0 if none)
%           .type, .arOrder, .sigma, .acv, .cost
%           .nugget, .wholeArc   (with cm.nugget) the relative nugget per
%                         epoch and true when any is > 0: the arc must be
%                         whitened in ONE batch (arcBatchWhiten honours it)
%
%   THE ACCUMULATOR PREFERS Wfcn
%     If both are returned, LSQVCE.ENGINE.BLOCKACCUMULATOR uses Wfcn and never builds a
%     covariance matrix at all. Cfcn is kept as a reference path: run one arc
%     both ways and check the normals agree (they do, to ~1e-11).
%
%   COST TABLE (b = batch size, p = filter memory)
%     type      path       cost/batch   memory    BatchSize constraint
%     diag      filter     O(b)         O(1)      none
%     arma      filter     O(b*p)       O(p)      none (overlap p)
%     psd/ar    filter     O(b*p)       O(p)      none (overlap p)
%     psd/toep  chol       O(b^3)       O(b^2)    b >> correlation length
%     acf       chol       O(b^3)       O(b^2)    b >> correlation length
%     var       filter     O(b*p*d^2)   O(p*d^2)  epoch-aligned, overlap p*d
%     custom    chol       O(b^3)       O(b^2)    b >> correlation length
%
%     At b = 5000 the dense path needs a 200 MB matrix and ~40 GFLOP. That,
%     not RAM, is the real ceiling. The filter path has no ceiling: set
%     BatchSize to the whole arc.
%
%   PSD -> AR IS THE RECOMMENDED PATH, AND IT IS EXACT
%     A residual PSD defines a stationary process. Fitting an AR(p) to it and
%     whitening with the FIR filter is not an approximation of the dense
%     Cholesky -- the filter IS the Cholesky inverse factor. Both are lower
%     triangular square roots of inv(C), and that factor is unique up to sign,
%     so the two agree to ~1e-14 row for row. You pay O(b*p) instead of O(b^3)
%     and give up nothing.
%
%     VALIDATE ONCE per observation type with LSQVCE.ENGINE.CHECKPSD before trusting it:
%         chk = lsqvce.engine.checkPSD(cm, dt, opt, true);
%         % aim for chk.acvRelErr < 1e-3; raise arOrder, or use 'toeplitz'
%
%   BATCH OVERLAP
%     info.filterMemory tells the accumulator how much to overlap batches.
%     Leave BatchOverlap unset and it is picked up automatically. Setting it
%     to 0 with a filter model costs ~12 orders of magnitude of accuracy.
%
%   EXAMPLE
%     % coloured 'pos', white 'kbr'
%     dt = 5;
%     f  = linspace(0, 0.5/dt, 4096).';
%     H  = 1 ./ (exp(-1i*2*pi*f*dt*(0:3)) * [1 -1.6 0.9 -0.25].');
%     cfg.cov.pos = struct('type','psd', 'f',f, 'S', 0.7^2*dt*2*abs(H).^2);
%     cfg.cov.kbr = struct('type','diag','sigma',0.2);
%
%     opt = struct('psdMode','ar', 'arOrder',20);
%     [Cf, Wf, info] = lsqvce.engine.covModel(cfg.cov.pos, t, opt);
%     info.filterMemory      % = 20 -> the accumulator overlaps by 20 rows
%
%   ERRORS
%     The spec is not a recognisable model
%       lsqvce:covModel:noType         cm has no .type field
%       lsqvce:covModel:unknownType    .type is not one of the eight above
%       lsqvce:covModel:noCovForm      the requested form cannot be built from
%                                      this spec
%       lsqvce:covModel:noFilterForm   ... likewise for the whitening filter
%
%     A required field for the CHOSEN type is missing or malformed
%       lsqvce:covModel:noSigma        .sigma missing        (most types)
%       lsqvce:covModel:badSigma       .sigma not a positive scalar
%       lsqvce:covModel:sigmaLength    per-observation .sigma has the wrong length
%       lsqvce:covModel:noAR           .ar missing           ('arma')
%       lsqvce:covModel:noACF          .acv missing          ('acf')
%       lsqvce:covModel:noPSD          .psd missing          ('psd')
%       lsqvce:covModel:noVAR          .ar/.sigma missing    ('var')
%       lsqvce:covModel:noFcn          .fcn missing          ('custom')
%       lsqvce:covModel:maNotInvertible the MA part ('arma', or a psd law in
%                                      'arma' mode) has a root on or outside the
%                                      unit circle: its whitener is not stable
%
%     The per-epoch error information does not fit the arc (R6 WP-90)
%       lsqvce:covModel:badCep         .Cep is neither nEpoch x d nor nEpoch x d x d,
%                                      does not cover the arc's rows, or disagrees
%                                      with a .d on the spec
%       lsqvce:covModel:cepOnVar       a 'var' law with a .Cep
%       lsqvce:covModel:badNugget      .nugget without a Cep, of the wrong length,
%                                      negative or not finite
%       lsqvce:covModel:nuggetNeedsWholeArc  a nugget arc's whitener called on a
%                                      batch that is not the whole arc
%       lsqvce:covModel:notWholeEpochs a batch (or a dense VAR request) does not
%                                      hold whole epochs; rows must be EPOCH-MAJOR
%                                      and the arc must carry nAxis = d
%
%     The model is not usable as given
%       lsqvce:covModel:degeneratePSD  the PSD has a zero or negative ordinate,
%                                      so its inverse does not exist
%       lsqvce:covModel:badPsdMode     unknown PSD interpolation mode
%       lsqvce:covModel:badWhitenMode  unknown whitening mode
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:covModel:negativePSD                [warning] The PSD has X negative values (min X)
%     lsqvce:covModel:psdSanitationEscape        PSD sanitation escape
%
%   ACCURACY
%     The whiteners are exact factors of the covariance they are built from
%     (the filter IS the Cholesky inverse factor, see above): diag / arma /
%     psd are pinned in test_precision and test_signalCore; the composition
%     with Cep is pinned against the dense reference D (R_t kron I_d) D at
%     1e-8 relative on the whitened normals (test_streamUtils, through the
%     loader and directly), and the PSD law on the epoch grid against
%     psd2acv at the epoch dt (1e-7, the psdRidge) -- R6 WP-90.
%
%   See also LSQVCE.ENGINE.CHECKPSD, LSQVCE.KERNEL.ARMAWHITEN, LSQVCE.KERNEL.SEPARABLEWHITEN,
%            LSQVCE.KERNEL.DECORRFACTOR, LSQVCE.ENGINE.BLOCKACCUMULATOR, LSQVCE.ENGINE.MAKELOADER
%
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-125: a dense law is whitened over
%   the whole arc (info.wholeArc); an MA part extends filterMemory by the decay
%   of 1/B(z); lsqvce:covModel:maNotInvertible. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-111 part 2 (R7-01): the
%   'toeplitz' PSD law is no longer resampled onto floor(m/2)+1 points; its
%   ACV is psd2acv's exact integral (no circulant). With assistance from
%   Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-115 part 3b (R7-10, user
%   decision (b)): the robust factors compose as an additive nugget
%   (cm.nugget) instead of the multiplicative fold through Cep. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-09-09 (CEST, UTC+02:00), R6 WP-106: the retired 'separable'
%   type no longer listed under noAR (found by audit check 28 on its first
%   run). With assistance from Claude (Anthropic), model Claude Fable 5.1.
%   Updated 2026-09-06 (CEST, UTC+02:00), R6 WP-90: ONE composition of any
%   temporal law with the per-epoch Cep replaces the 'separable' type; the
%   law is built on the epoch grid. With assistance from Claude (Anthropic),
%   model Claude Fable 5.1.
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-124: Cep reference -> API reference.
%   With assistance from Claude (Anthropic), model Claude Opus 5.5.

    narginchk(2,3);
    if nargin < 3 || isempty(opt), opt = struct(); end

    if ~isstruct(cm) || ~isfield(cm,'type')
        error('lsqvce:covModel:noType', ...
              'The covariance model must be a struct with a ''type'' field.');
    end
    validateattributes(t, {'double'}, {'vector','finite','nonempty'}, mfilename, 't', 2);

    t = t(:);
    % R6 WP-90: the per-epoch ERROR INFORMATION D (cm.Cep, set by the loader
    % from the day buildDay returned) composes with the TEMPORAL law that
    % cm.type declares: C = D (R_t kron I_d) D. The temporal law is built on
    % the EPOCH grid (every d-th row) -- a law placed on the row grid of a
    % 3-axis stream would see dt / 3 (the round-3b class of bug).
    Cep = lsqvce.util.getfielddef(cm, 'Cep', []);
    hasCep = ~isempty(Cep);
    d = 1;
    if hasCep
        [nEp, d] = cepShape(Cep);
        if numel(t) ~= nEp * d
            error('lsqvce:covModel:badCep', ...
                  ['cm.Cep has %d epochs of %d axes (%d rows) but the arc has %d rows; rows ' ...
                   'must be epoch-major and Cep must cover every epoch of the arc.'], ...
                   nEp, d, nEp * d, numel(t));
        end
        dDecl = lsqvce.util.getfielddef(cm, 'd', []);
        if ~isempty(dDecl) && dDecl ~= d
            error('lsqvce:covModel:badCep', ...
                  'cm.d = %d but cm.Cep has %d axes; the axis count is stated once, by Cep.', dDecl, d);
        end
        if strcmpi(cm.type, 'var')
            error('lsqvce:covModel:cepOnVar', ...
                  ['''var'' carries its own cross-axis structure; a per-epoch Cep cannot be ' ...
                   'composed with it. A stream that returns Cep declares a scalar temporal ' ...
                   'law (diag / arma / psd / acf).']);
        end
        t = t(1:d:end);                              % the epoch grid
    end
    m = numel(t);
    nug = nuggetOf(cm, hasCep, m);                   % R7 WP-115 part 3b

    dt = median(diff(t));
    if m < 2 || ~isfinite(dt) || dt <= 0
        dt = 1;
    end

    Cfcn = [];
    Wfcn = [];
    info = struct('filterMemory',0, 'type',lower(cm.type), 'arOrder',0, ...
                  'sigma',NaN, 'acv',[], 'cost','', 'nAxis',1);

    mode = lower(lsqvce.util.getfielddef(opt, 'whitenMode', 'auto'));

    switch lower(cm.type)

    %------------------------------------------------------------------ diag
    case 'diag'
        s = lsqvce.util.getfielddef(cm, 'sigma', []);
        if isempty(s)
            error('lsqvce:covModel:noSigma', '''diag'' requires cm.sigma.');
        end
        if isscalar(s), s = s*ones(m,1); else, s = s(:); end
        if numel(s) ~= m
            error('lsqvce:covModel:sigmaLength', ...
                  'cm.sigma has %d entries but the arc has %d epochs.', numel(s), m);
        end
        if any(s <= 0) || any(~isfinite(s))
            error('lsqvce:covModel:badSigma', ...
                  'cm.sigma must be finite and strictly positive.');
        end

        Cfcn = @(idx,o) diag(s(idx).^2);
        Wfcn = @(idx,A,B,l,o) diagWhiten(s(idx), A, B, l);
        info.sigma = mean(s);
        info.cost  = 'O(b)';

    %------------------------------------------------------------------- acf
    case 'acf'
        if ~isfield(cm,'acf') || ~isa(cm.acf,'function_handle')
            error('lsqvce:covModel:noACF', '''acf'' requires cm.acf as a function handle.');
        end
        s = lsqvce.util.getfielddef(cm, 'sigma', 1);
        if ~isscalar(s), s = mean(s(:)); end
        if s <= 0
            error('lsqvce:covModel:badSigma', 'cm.sigma must be positive.');
        end

        Cfcn = @(idx,o) (s^2) * cm.acf( t(idx) - t(idx).', o );
        info.sigma = s;
        info.cost  = 'O(b^3)';

    %------------------------------------------------------------------- psd
    case 'psd'
        [Cfcn, Wfcn, info] = psdModel(cm, t, dt, opt, info);

    %------------------------------------------------------------------ arma
    case 'arma'
        a = lsqvce.util.getfielddef(cm, 'ar', []);
        if isempty(a)
            error('lsqvce:covModel:noAR', '''arma'' requires cm.ar.');
        end
        a = a(:).';
        b = lsqvce.util.getfielddef(cm, 'ma', 1);
        b = b(:).';
        s = lsqvce.util.getfielddef(cm, 'sigma', []);
        if isempty(s) || s <= 0
            error('lsqvce:covModel:badSigma', '''arma'' requires a positive cm.sigma.');
        end

        p = numel(a) - 1;
        r = s^2 * lsqvce.kernel.armaACV(a, b, max(p, 2));      % errors if non-stationary

        Wfcn = @(idx,A,B,l,o) lsqvce.kernel.armaWhiten(a, b, s, r, A, B, l, idx(1)==1);
        Cfcn = @(idx,o) s^2 * toeplitz(lsqvce.kernel.armaACV(a, b, numel(idx)));

        info.filterMemory = filterMemoryOf(p, b);  % WP-125: an MA part makes it IIR
        info.arOrder      = p;
        info.sigma        = s;
        info.acv          = r;
        info.cost         = sprintf('O(b*%d)', p);

    %------------------------------------------------------------------- var
    % Vector AR: temporal cross-axis correlation, the 'needs-var' verdict of
    % LSQVCE.KERNEL.CROSSCORR. Rows must be epoch-major; the arc must carry
    % nAxis = d so the accumulator batches on whole epochs.
    case 'var'
        Av = lsqvce.util.getfielddef(cm, 'ar', []);
        if isempty(Av) || ndims(Av) ~= 3 || size(Av,1) ~= size(Av,2)
            error('lsqvce:covModel:noVAR', ...
                  '''var'' requires cm.ar as a d x d x p array.');
        end
        d  = size(Av, 1);
        p  = size(Av, 3);
        Sg = lsqvce.util.getfielddef(cm, 'sigma', []);
        if isempty(Sg) || ~isequal(size(Sg), [d d])
            error('lsqvce:covModel:badSigma', ...
                  '''var'' requires cm.sigma as the d x d innovation covariance.');
        end
        Gv = lsqvce.util.getfielddef(cm, 'acv', []);
        if isempty(Gv)
            Gv = lsqvce.kernel.varACV(Av, Sg, max(p, 1));   % startup block; stationarity check
        end

        Wfcn = @(idx,A,B,l,o) lsqvce.kernel.varWhiten(Av, Sg, Gv, d, A, B, l, idx(1)==1);
        Cfcn = @(idx,o) varDense(Av, Sg, d, numel(idx));

        info.filterMemory = p * d;                % ROWS of overlap
        info.nAxis        = d;
        info.arOrder      = p;
        info.cost         = sprintf('O(b*%d*%d^2)', p, d);

    %---------------------------------------------------------------- custom
    case 'custom'
        if ~isfield(cm,'fcn') || ~isa(cm.fcn,'function_handle')
            error('lsqvce:covModel:noFcn', '''custom'' requires cm.fcn as a function handle.');
        end
        Cfcn      = @(idx,o) cm.fcn(t(idx), o);
        info.cost = 'O(b^3)';

    otherwise
        error('lsqvce:covModel:unknownType', ...
            ['Unknown covariance type ''%s''. Valid: diag, acf, psd, arma, var, custom. ' ...
             'The former ''separable'' type is a temporal law (arma / diag / psd) plus the ' ...
             'per-epoch Cep that buildDay returns with the day (WP-90).'], cm.type);
    end

    %---------------------------------------------------- explicit user override
    if strcmp(mode, 'chol')
        if isempty(Cfcn)
            error('lsqvce:covModel:noCovForm', ...
                  'whitenMode=''chol'' but model ''%s'' has no dense form.', cm.type);
        end
        Wfcn = [];
        info.filterMemory = 0;
        info.cost = 'O(b^3)';
    elseif strcmp(mode, 'filter')
        if isempty(Wfcn)
            error('lsqvce:covModel:noFilterForm', ...
                ['whitenMode=''filter'' but model ''%s'' has no filter form. ' ...
                 'Only diag, arma, var and psd (psdMode ar/arma) can be filtered.'], ...
                 cm.type);
        end
    elseif ~strcmp(mode, 'auto')
        error('lsqvce:covModel:badWhitenMode', ...
              'whitenMode must be ''auto'', ''filter'' or ''chol'' (got ''%s'').', mode);
    end

    %------------------------------------------- a dense law is whitened whole
    % WP-125 (user decision 2026-10-12): a law without a filter form (acf,
    % custom, psd 'toeplitz', or whitenMode 'chol') is factorised per batch by
    % the iterator; batches of a dense law cannot be joined by an overlap (the
    % conditional on a finite past is not the conditional on the whole past),
    % so an arc longer than one batch lost the correlation between batches.
    % The arc is now factorised whole -- exact, at O(m^3) and m^2 memory.
    if isempty(Wfcn) && ~isempty(Cfcn)
        info.wholeArc = true;
    end

    %------------------------------------------------ the composition (WP-90)
    if hasCep
        [Cfcn, Wfcn, info] = composeWithCep(Cfcn, Wfcn, info, Cep, d, nug);
    end
end

%==========================================================================
function [nEp, d] = cepShape(Cep)
%CEPSHAPE  [nEpoch, d] of a per-epoch scaling: nEp x d (sigmas) or nEp x d x d.
    if ~isnumeric(Cep) || isempty(Cep)
        error('lsqvce:covModel:badCep', ...
              'cm.Cep must be a numeric nEpoch x d (sigmas) or nEpoch x d x d (covariances) array.');
    end
    sz = size(Cep);
    if numel(sz) == 2
        nEp = sz(1);  d = sz(2);
    elseif numel(sz) == 3 && sz(2) == sz(3)
        nEp = sz(1);  d = sz(2);
    else
        error('lsqvce:covModel:badCep', ...
              'cm.Cep must be nEpoch x d or nEpoch x d x d (got %s).', mat2str(sz));
    end
end

%==========================================================================
function [Cfcn, Wfcn, info] = composeWithCep(Ct, Wt, info, Cep, d, nug)
%COMPOSEWITHCEP  C = D (R_t kron I_d) D from the temporal law and the per-epoch D.
%   Ct / Wt are the temporal covariance / whitener over EPOCHS (built above
%   on the epoch grid). The row whitener descales every epoch by D^-1/2 and
%   then applies the temporal whitener PER AXIS (rows k:d:end are one series
%   down the epochs -- the order LSQVCE.KERNEL.SEPARABLEWHITEN proves). Cfcn
%   stays the TEMPORAL covariance over epochs: what the joint estimator forms
%   per axis after the same descale. A dense-only law (psd 'toeplitz', acf,
%   custom) is composed through the dense decorrelation factor of Ct.
    if isempty(Wt) && isempty(Ct)
        error('lsqvce:covModel:noFilterForm', ...
              'the temporal law ''%s'' has neither a filter nor a dense form to compose with Cep.', info.type);
    end
    % (the runtime opt is passed through when the caller gives it -- every
    % engine does -- and omitted otherwise, as the temporal closures allow)
    if ~isempty(Wt)
        Wfcn = @(idx, A, B, l, varargin) composeRows(@(ep, A, B, l) Wt(ep, A, B, l, varargin{:}), ...
                                                     Cep, d, idx, A, B, l);
    else
        Wfcn = @(idx, A, B, l, varargin) composeRows(@(ep, A, B, l) denseWhiten(Ct(ep, varargin{:}), A, B, l), ...
                                                     Cep, d, idx, A, B, l);
    end
    Cfcn = Ct;
    info.filterMemory = info.filterMemory * d;       % ROWS of overlap (epoch-major)
    info.nAxis = d;
    info.Cep   = Cep;                                % the fixed per-epoch part travels with the arc
    if ~isempty(nug) && any(nug > 0)
        % R7 WP-115 part 3b: the nugget. The temporal whitener per axis becomes
        % G F (nuggetRows); it is not causal, so the arc is whitened whole.
        if isempty(Ct)
            error('lsqvce:covModel:badNugget', ...
                  'the temporal law ''%s'' has no dense form: a nugget cannot be scaled to its variance.', info.type);
        end
        if ~isempty(Wt)
            Wlaw = @(ep, A, B, l, varargin) Wt(ep, A, B, l, varargin{:});
        else
            Wlaw = @(ep, A, B, l, varargin) denseWhiten(Ct(ep, varargin{:}), A, B, l);
        end
        Wep0 = Wfcn;
        Wfcn = @(idx, A, B, l, varargin) nuggetRows(Wep0, Wlaw, Ct, nug, d, idx, A, B, l, varargin{:});
        info.nugget   = nug;
        info.wholeArc = true;
    end
end

function nug = nuggetOf(cm, hasCep, m)
%NUGGETOF  The validated cm.nugget (m x 1, >= 0, finite; with Cep only), or [].
    nug = lsqvce.util.getfielddef(cm, 'nugget', []);
    if isempty(nug), return, end
    nug = nug(:);
    if ~hasCep || numel(nug) ~= m || any(~isfinite(nug)) || any(nug < 0)
        error('lsqvce:covModel:badNugget', ...
              ['cm.nugget must be a finite, non-negative nEpoch x 1 vector (%d epochs) and ' ...
               'composes with the per-epoch Cep only (got %d values, Cep given: %d).'], ...
              m, numel(nug), hasCep);
    end
end

function [A, B, l] = nuggetRows(Wc, Wlaw, Ct, q, d, idx, A, B, l, varargin)
%NUGGETROWS  Whiten a WHOLE arc under C = D ((R_t + diag(c)) kron I_d) D,
%   c_e = q_e Sigma_ee (R7 WP-115 part 3b, mirrors/wp115/nugget.py).
%   Wc is the composed whitener (descale by D, then the law's whitener F
%   down the epochs of each axis); Wlaw is F alone, over epochs. The nugget
%   is a rank-k correction over the k epochs with q > 0: with U their unit
%   columns, F U = Q R (thin QR) and K = diag(1/c) + R'R,
%   inv(R_t + U diag(c) U') = F' (I - Q B Q') F with B = R K^-1 R'
%   (Woodbury), and G = I - Q V diag(1 - sqrt(1 - beta)) V' Q' (B = V
%   diag(beta) V') satisfies G'G = I - Q B Q'. G is applied to every axis
%   series after Wc. K holds 1/c, never c, so a zero-weight epoch (c ~ 1e12
%   Sigma_ee) costs no accuracy: the rows reproduce A' C^-1 A to rounding
%   (1e-15 in the mirror).
    mE = numel(q);
    if numel(idx) ~= mE * d || idx(1) ~= 1
        error('lsqvce:covModel:nuggetNeedsWholeArc', ...
              ['A nugget arc is whitened whole (%d rows); this batch holds rows %d..%d. ' ...
               'The batch iterator honours the arc''s wholeArc flag; a caller that batches ' ...
               'by hand must pass the whole arc.'], mE * d, idx(1), idx(end));
    end
    [A, B, l] = Wc(idx, A, B, l, varargin{:});
    z  = find(q > 0);  k = numel(z);
    Cz = Ct(z, varargin{:});                          % the law over the nugget epochs
    if isvector(Cz), cz = Cz(:); else, cz = diag(Cz); end
    c  = q(z) .* cz;
    E  = zeros(mE, k);
    E(sub2ind([mE, k], z(:), (1:k).')) = 1;
    [FU, ~, ~] = Wlaw(1:mE, E, zeros(mE, 0), zeros(mE, 1), varargin{:});
    [Q, R] = qr(FU, 0);
    Km = diag(1 ./ c) + R.' * R;
    Bm = R * (Km \ R.');
    [V, beta] = eig((Bm + Bm.') / 2, 'vector');
    g = 1 - sqrt(max(1 - beta, 0));
    P = Q * V;
    for j = 1:d
        rows = j:d:mE * d;
        A(rows, :) = A(rows, :) - P * (g .* (P.' * A(rows, :)));
        B(rows, :) = B(rows, :) - P * (g .* (P.' * B(rows, :)));
        l(rows)    = l(rows)    - P * (g .* (P.' * l(rows)));
    end
end

function [A, B, l] = composeRows(Wep, Cep, d, idx, A, B, l)
%COMPOSEROWS  One batch: descale whole epochs, then whiten each axis down the epochs.
    idx = idx(:).';
    if mod(numel(idx), d) ~= 0 || mod(idx(1) - 1, d) ~= 0
        error('lsqvce:covModel:notWholeEpochs', ...
              ['A batch of %d rows starting at row %d does not hold whole epochs of %d axes; ' ...
               'rows must be epoch-major and a batch must not split an epoch.'], numel(idx), idx(1), d);
    end
    ep = (idx(1) - 1) / d + (1:numel(idx) / d);      % the batch's epochs
    De = epochSlice(Cep, idx, d);
    [A, B, l] = lsqvce.kernel.separableWhiten(1, 1, 1, De, d, A, B, l, false);   % D^-1/2 only (AR(0))
    for k = 1:d
        rows = k:d:numel(idx);
        [A(rows, :), B(rows, :), l(rows)] = Wep(ep, A(rows, :), B(rows, :), l(rows));
    end
end

function [A, B, l] = denseWhiten(C, A, B, l)
%DENSEWHITEN  Whiten one axis series with the dense factor of its temporal covariance.
    W = lsqvce.kernel.decorrFactor(C);
    A = W * A;  B = W * B;  l = W * l;
end

%==========================================================================
function [Aw, Bw, lw] = diagWhiten(s, A, B, l)
%DIAGWHITEN  Scale rows by 1/sigma. Shape-safe for empty A or B.
    Aw = A;  Bw = B;
    if ~isempty(A), Aw = A ./ s;  end
    if ~isempty(B), Bw = B ./ s;  end
    lw = l ./ s;
end

%==========================================================================
function De = epochSlice(Cep, idx, d)
%EPOCHSLICE  Select the per-epoch scaling rows for a batch.
%   idx are epoch-major ROW indices (length multiple of d); the corresponding
%   epochs are unique(ceil(idx/d)). Cep is nEpoch x d (diagonal) or
%   nEpoch x d x d (full). Returns the sub-array SEPARABLEWHITEN expects.
    ep = unique(ceil(idx(:)/d), 'stable');
    if ismatrix(Cep) && size(Cep,2) == d
        De = Cep(ep, :);                 % nBatchEpoch x d diagonal sigma
    elseif ndims(Cep) == 3 && size(Cep,2) == d && size(Cep,3) == d
        De = Cep(ep, :, :);              % nBatchEpoch x d x d full covariance
    else
        error('lsqvce:covModel:badCep', ...
              'cm.Cep must be nEpoch x %d or nEpoch x %d x %d.', d, d, d);
    end
end

%==========================================================================
function [Cfcn, Wfcn, info] = psdModel(cm, t, dt, opt, info)
%PSDMODEL  PSD -> dense Toeplitz, or AR/ARMA whitening filter.

    % R5 M3 3b: a PSD may be given on NORMALISED frequency (cm.fNorm, cycles
    % per sample, 0 .. 0.5) -- the form the estimators emit, since their lags
    % are sample lags; it is placed on the physical axis with THIS arc's dt.
    if isfield(cm, 'fNorm') && ~isfield(cm, 'f')
        cm.f = cm.fNorm(:) / dt;
        cm.S = cm.S(:) * dt;      % r(0) = int S df: the axis shrinks by dt, S grows by it
    end
    if ~isfield(cm,'f') || ~isfield(cm,'S')
        error('lsqvce:covModel:noPSD', '''psd'' requires cm.f and cm.S (or cm.fNorm and cm.S).');
    end
    % R5 M4 4b: a 'psd' law may carry a LEVEL, cm.sigma (a standard deviation
    % multiplying the process; the PSD is then a SHAPE) -- the form the joint
    % estimator emits (shape + per-stream level via sigmaScale). Ignoring it
    % made every pass re-estimate the level from scratch and the product ran
    % away (found by the gravity demo: 1.5, 5, 22, 416 on a planted 2).
    sLvl = lsqvce.util.getfielddef(cm, 'sigma', []);
    if ~isempty(sLvl)
        if ~isscalar(sLvl) || ~isfinite(sLvl) || sLvl <= 0
            error('lsqvce:covModel:badSigma', '''psd'' cm.sigma must be a positive scalar level.');
        end
        cm.S = cm.S(:) * sLvl^2;
    end

    md = lower(lsqvce.util.getfielddef(opt, 'psdMode',  'ar'));
    p  =       lsqvce.util.getfielddef(opt, 'arOrder',  30);
    q  =       lsqvce.util.getfielddef(opt, 'maOrder',  0);
    fl =       lsqvce.util.getfielddef(opt, 'psdFloor', 1e-10);

    f = cm.f(:);
    S = cm.S(:);

    if any(S < 0)
        lsqvce.util.log.warn('lsqvce:covModel:negativePSD', ...
            ['The PSD has %d negative values (min %.3e). They are being ' ...
             'floored, but a negative PSD is not physical -- check how it was ' ...
             'estimated (a raw periodogram difference can go negative).'], ...
             sum(S < 0), min(S));
    end
    Smax = max(S);
    if ~isfinite(Smax) || Smax <= 0
        error('lsqvce:covModel:degeneratePSD', 'The PSD is zero or non-finite everywhere.');
    end
    % WP-40 r2 hardening: the elementwise max() floor SHOULD already heal
    % every negative/NaN entry, yet CI observed psd2acv rejecting S as
    % non-nonnegative directly below (stack pinned in the WP-40 r2 round
    % note). Sanitize by explicit indexing -- beyond doubt for every value
    % class -- and assert with diagnostics so any residual escape names its
    % value instead of failing five lines downstream.
    S(~isfinite(S) | S < fl*Smax) = fl*Smax;
    if any(S <= 0) || any(~isfinite(S))
        error('lsqvce:covModel:psdSanitationEscape', ...
             ['PSD sanitation escape: min(S) = %.6e, #nonfinite = %d, ' ...
              'fl = %.3e, Smax = %.6e. Report this state.'], ...
              min(S), sum(~isfinite(S)), fl, Smax);
    end

    Cfcn = [];
    Wfcn = [];

    switch md
    case 'toeplitz'
        % WP-125: a dense law is factorised over the WHOLE arc (info.wholeArc
        % below), so its ACV must reach every lag of the arc, not one batch
        nmax = numel(t);
        % RESOLUTION (R7 WP-111 part 2, finding R7-01). The ACV is the EXACT
        % integral of the piecewise-linear spectrum through the given
        % ordinates (LSQVCE.KERNEL.PSD2ACV): no lag period, so a coarse law
        % (the joint estimator's Grid of 60 bins) yields its own ACV at every
        % lag of any arc, and a positive S a positive-definite Toeplitz. The
        % law used to be resampled onto floor(m/2)+1 points and integrated by
        % trapezoids, which is periodic in the lag: r(m-1) = r(1) on even
        % arcs, = r(0) on odd ones (a circulant, Toeplitz cond 3.3e8 at
        % phi 0.5, measured on the runner -- HANDOVER 12.44).
        r    = lsqvce.kernel.psd2acv(f, S, dt, nmax);      % once per arc, reused per batch
        rg   = lsqvce.util.getfielddef(opt, 'psdRidge', 1e-8);
        if rg > 0 && ~isempty(r), r(1) = r(1) * (1 + rg); end
        Cfcn = @(idx,o) toeplitz(r(1:numel(idx)));
        info.filterMemory = 0;
        info.acv          = r;
        info.cost         = 'O(b^3)';

    case {'ar','arma'}
        validateattributes(p, {'numeric'}, {'scalar','integer','positive'}, ...
                           mfilename, 'opt.arOrder');
        r0      = lsqvce.kernel.psd2acv(f, S, dt, p+1);
        [a, s2] = lsqvce.kernel.levinsonDurbin(r0, p);

        if strcmp(md,'arma') && q > 0
            [a, b, s2] = lsqvce.kernel.ar2arma(a, s2, q, opt);
        else
            b = 1;
        end

        s   = sqrt(s2);
        pEf = numel(a) - 1;                       % may be < p if truncated
        rr  = lsqvce.kernel.psd2acv(f, S, dt, max(pEf, 2)); % exact stationary startup

        Wfcn = @(idx,A,B,l,o) lsqvce.kernel.armaWhiten(a, b, s, rr, A, B, l, idx(1)==1);
        Cfcn = @(idx,o) s2 * toeplitz(lsqvce.kernel.armaACV(a, b, numel(idx)));

        info.filterMemory = filterMemoryOf(pEf, b);   % WP-125: an MA part makes it IIR
        info.arOrder      = pEf;
        info.sigma        = s;
        info.acv          = rr;
        info.cost         = sprintf('O(b*%d)', pEf);

    otherwise
        error('lsqvce:covModel:badPsdMode', ...
              'psdMode must be ''ar'', ''arma'' or ''toeplitz'' (got ''%s'').', md);
    end
end

%==========================================================================
function n = filterMemoryOf(p, b)
%FILTERMEMORYOF  The batch overlap that makes the whitening filter A(z)/B(z) exact.
%   The whitener filter(a, b, x) restarts from zero state at the start of
%   every overlap. With B = 1 (pure AR) it is FIR of length p+1, and p rows of
%   overlap reproduce the whole-arc output exactly. An MA part B(z) puts
%   1/B(z) -- an IIR filter -- behind it, whose memory is infinite: the
%   overlap is p plus the number of samples after which the impulse response
%   of 1/B(z) has decayed below eps of its peak (WP-125; before, the overlap
%   counted p alone and the restart transient of the MA part leaked into
%   every batch after the first).
%
%   INPUT   p  1x1  AR order;  b  1x(q+1)  MA polynomial [1 b1 ... bq]
%   OUTPUT  n  1x1  rows of overlap (p when q = 0)
%   ERRORS  lsqvce:covModel:maNotInvertible  a root of B(z) on or outside the
%           unit circle: the whitening filter 1/B(z) is not stable
    n = p;
    b = b(:).';
    if numel(b) < 2 || all(b(2:end) == 0), return, end
    rho = max(abs(roots(b)));
    if ~(rho < 1)
        error('lsqvce:covModel:maNotInvertible', ...
              ['The MA polynomial has a root of modulus %.6g >= 1: the whitening ' ...
               'filter 1/B(z) is not stable. Give an invertible MA part.'], rho);
    end
    N = max(16, ceil(log(eps) / log(rho)));
    for it = 1:12                                   % doubling; 16 * 2^12 > 65 000
        h = filter(1, b, [1; zeros(N - 1, 1)]);
        if all(abs(h(end - numel(b) + 2:end)) <= eps * max(abs(h))), break, end
        N = 2 * N;
    end
    n = p + N;
end

%==========================================================================
function C = varDense(Av, Sg, d, mrows)
%VARDENSE  Block-Toeplitz covariance of a VAR over whole epochs.
    if mod(mrows, d) ~= 0
        error('lsqvce:covModel:notWholeEpochs', ...
            ['A dense VAR covariance was requested for %d rows, not a ' ...
             'multiple of d = %d. Batches must hold whole epochs.'], mrows, d);
    end
    nE = mrows / d;
    G  = lsqvce.kernel.varACV(Av, Sg, nE-1);
    C  = zeros(mrows);
    for i = 1:nE
        for j = 1:nE
            k = i - j;
            if k >= 0, Bk = G(:,:,k+1);
            else,      Bk = G(:,:,-k+1).'; end
            C((i-1)*d + (1:d), (j-1)*d + (1:d)) = Bk;
        end
    end
end
