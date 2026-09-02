function [Aw, Bw, lw] = separableWhiten(a, s, r, sigEpoch, d, A, B, l, isFirstBatch)
%LSQVCE.KERNEL.SEPARABLEWHITEN  Whiten multi-axis observations with a separable covariance.
%
%   [Aw,Bw,lw] = LSQVCE.KERNEL.SEPARABLEWHITEN(a, s, r, sigEpoch, d, A, B, l, isFirstBatch)
%
%   THE MODEL
%       C = D * (R_t kron I_d) * D
%
%   D        per-epoch scaling (diagonal, or a full d x d covariance per epoch)
%   R_t      temporal correlation, the SAME for every axis
%   kron I_d NO cross-axis correlation
%
%   This is the realistic case for 3-D kinematic positions: you have an epochwise
%   covariance from the PPP solution, and a PSD describing the temporal colour.
%
%   THE WHITENER FACTORISES
%       W = (F_t kron I_d) * D^-1
%
%   so the operation is:
%     1. DE-SCALE each epoch by its Cholesky factor  (D^-1)
%     2. FILTER each axis independently down epochs  (F_t, the scalar AR filter)
%
%   ORDER MATTERS: de-scale FIRST, then filter. The reverse is wrong unless the
%   per-epoch sigma is constant.
%
%   COST
%     O(nEpoch*d^3) for the de-scale + O(m*p) for the filter. Multi-axis is
%     therefore NO MORE EXPENSIVE per row than the scalar case -- the d^3 term is
%     trivial (27 flops per epoch for d = 3). That is the payoff of separability.
%
%   ROW ORDERING -- EPOCH-MAJOR, NOT AXIS-MAJOR
%       [x1 y1 z1 x2 y2 z2 ...]     correct
%       [x1 x2 ... y1 y2 ... ]      WRONG
%   The block-Toeplitz structure is contiguous only in epoch-major order, and the
%   filter reshape below depends on it. m must be a whole multiple of d, and a
%   batch must never split an epoch (LSQVCE.ENGINE.BLOCKACCUMULATOR rounds BatchSize down to a
%   multiple of d).
%
%   INPUT
%     a          1 x (p+1)  scalar temporal AR polynomial, shared by all axes
%     s          scalar     innovation std of the temporal model
%     r          >= p x 1   temporal autocovariance (already scaled by s^2)
%     sigEpoch   nEpoch x d          per-epoch, per-axis sigma      (diagonal D)
%             or nEpoch x d x d      per-epoch full covariance      (D = chol)
%     d          number of axes
%     A, B, l    the arc's design blocks and observations, epoch-major
%     isFirstBatch  true for the first batch of the arc (stationary startup)
%
%   USE THE FULL d x d IF YOU HAVE IT
%     A kinematic PPP solution gives you a full 3x3 per-epoch covariance: the
%     vertical is 2-3x worse than horizontal AND correlated with it. Passing only
%     a diagonal sigma throws that away and mis-weights everything. The diagonal
%     form is supported, but it is the special case, not the default you should
%     want. LSQVCE.KERNEL.CROSSCORR will tell you if you got this wrong: it reports
%     'cross-axis-zero-lag'.
%
%   WHEN SEPARABILITY FAILS
%     If the axes are cross-correlated IN TIME (axis i at epoch k with axis j at
%     epoch k-tau), this model is wrong and you need a vector-AR. LSQVCE.KERNEL.CROSSCORR
%     detects that case and reports 'needs-var'. Do not assume; test.
%
%   EXAMPLE
%     % 3-D positions, 5 s sampling, full per-epoch covariance from PPP
%     d   = 3;
%     Cep = ppp.covariance;                     % nEpoch x 3 x 3
%     [a, s2] = lsqvce.kernel.levinsonDurbin(lsqvce.kernel.psd2acv(f, S, dt, p+1), p);
%     r   = lsqvce.kernel.psd2acv(f, S, dt, max(p,2));
%
%     [Aw,Bw,lw] = lsqvce.kernel.separableWhiten(a, sqrt(s2), r, Cep, d, A, B, l, true);
%
%   INPUT
%     a            1x(p+1) double  temporal AR polynomial [1 a1 ... ap]
%     s            1x1     double  innovation std of the temporal law
%     r            1xL     double  temporal ACV for the stationary startup
%     sigEpoch     nEp x d | nEp x d x d  the per-epoch scaling D
%     d            1x1     double  axes per epoch
%     A            m x nc  double  common design rows (EPOCH-MAJOR)
%     B            m x nl  double  local design rows  (may be [])
%     l            m x 1   double  observations
%     isFirstBatch 1x1     logical apply the stationary startup
%
%   OUTPUT
%     Aw, Bw, lw   the same arrays, whitened: solving with them is equivalent
%                  to a GLS solve against C = D (R_t kron I_d) D'
%
%   ERRORS
%     lsqvce:separable:notWholeEpochs  m is not a multiple of d, i.e. the
%                  batch does not contain whole epochs. Rows MUST be epoch-major
%                  and the arc must set nAxis = d.
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:separable:acvTooShort               The autocovariance has X lags but the startup needs X.lower
%     lsqvce:separable:badSigma                  Innovation std s = X must be positive
%     lsqvce:separable:badSigmaShape             sigEpoch must be nEpoch x X (diagonal) or nEpoch x X x X (full per-ep...
%     lsqvce:separable:epochNotPD                The X x X covariance of epoch X is not positive definite
%     lsqvce:separable:sigmaLength               sigEpoch has X epochs but the block has X., nE*d, 1); Xpoch-major if...
%     lsqvce:separable:startupNotPD              [warning] The stationary startup block is not PD; keeping the transient
%
%   See also LSQVCE.KERNEL.CROSSCORR, LSQVCE.KERNEL.ARMAWHITEN, LSQVCE.ENGINE.COVMODEL
%
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+kernel (R4 decision R4-D1,
%   WP-64a): the KERNEL tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    narginchk(9,9);

    a = a(:).';
    r = r(:);
    p = numel(a) - 1;
    m = size(l,1);
    d = double(d);

    if mod(m, d) ~= 0
        error('lsqvce:separable:notWholeEpochs', ...
            ['The block has %d rows, which is not a multiple of d = %d. Rows must ' ...
             'be epoch-major and a batch must not split an epoch.'], m, d);
    end
    nE = m / d;

    if s <= 0 || ~isfinite(s)
        error('lsqvce:separable:badSigma', 'Innovation std s = %g must be positive.', s);
    end

    % ---- 1. de-scale each epoch -------------------------------------------
    Ad = descale(A, sigEpoch, nE, d);
    Bd = descale(B, sigEpoch, nE, d);
    ld = descale(l, sigEpoch, nE, d);

    % ---- 2. filter each axis independently, down epochs -------------------
    Aw = tfilt(Ad, a, s, r, nE, d, p, isFirstBatch);
    Bw = tfilt(Bd, a, s, r, nE, d, p, isFirstBatch);
    lw = tfilt(ld, a, s, r, nE, d, p, isFirstBatch);
end

%==========================================================================
function X = descale(X, sigEpoch, nE, d)
%DESCALE  Apply D^-1: divide each epoch by its (Cholesky of the) covariance.
    if isempty(X), return, end
    nCol = size(X,2);

    if ismatrix(sigEpoch) && size(sigEpoch,2) == d
        % --- diagonal case: sigEpoch is nEpoch x d -------------------------
        if size(sigEpoch,1) ~= nE
            error('lsqvce:separable:sigmaLength', ...
                  'sigEpoch has %d epochs but the block has %d.', ...
                  size(sigEpoch,1), nE);
        end
        sc = reshape(sigEpoch.', nE*d, 1);         % epoch-major
        if any(sc <= 0) || any(~isfinite(sc))
            error('lsqvce:separable:badSigma', 'sigEpoch must be finite and positive.');
        end
        X = X ./ sc;

    elseif ndims(sigEpoch) == 3 && size(sigEpoch,2) == d && size(sigEpoch,3) == d
        % --- full d x d per epoch ------------------------------------------
        if size(sigEpoch,1) ~= nE
            error('lsqvce:separable:sigmaLength', ...
                  'sigEpoch has %d epochs but the block has %d.', ...
                  size(sigEpoch,1), nE);
        end
        Xr = reshape(X, d, nE, nCol);              % axis, epoch, col
        for k = 1:nE
            Ck = squeeze(sigEpoch(k,:,:));
            Ck = (Ck + Ck.')/2;
            [Lk, fail] = chol(Ck, 'lower');
            if fail > 0
                error('lsqvce:separable:epochNotPD', ...
                    ['The %d x %d covariance of epoch %d is not positive ' ...
                     'definite.'], d, d, k);
            end
            Xr(:,k,:) = Lk \ squeeze(Xr(:,k,:));
        end
        X = reshape(Xr, nE*d, nCol);

    else
        error('lsqvce:separable:badSigmaShape', ...
            ['sigEpoch must be nEpoch x %d (diagonal) or nEpoch x %d x %d ' ...
             '(full per-epoch covariance).'], d, d, d);
    end
end

%==========================================================================
function Y = tfilt(X, a, s, r, nE, d, p, isFirstBatch)
%TFILT  Apply the scalar temporal filter independently down each axis.
%   Rows are epoch-major (row = (epoch-1)*d + axis). To filter each axis (and
%   each design column) independently DOWN epochs, we must put epochs on
%   dimension 1 with every (axis,col) pair as its own column. A plain
%   reshape(X, d*nCol, nE).' interleaves columns for nCol > 1 -- it groups the
%   first d*nCol linear elements, which spans only ONE epoch's worth of the
%   first few columns, so consecutive "rows" of the reshape mix design columns
%   and the filter convolves across the wrong series. Reshape to
%   (d, nEpoch, nCol) then permute epochs to the front instead.
    if isempty(X), Y = X; return, end
    nCol = size(X,2);

    Xr = reshape(X, d, nE, nCol);                  % axis, epoch, col
    Xr = permute(Xr, [2 1 3]);                      % epoch, axis, col
    Xr = reshape(Xr, nE, d*nCol);                   % nEpoch x (d*nCol), one series/col
    Yr = filter(a, 1, Xr, [], 1) / s;

    if p > 0 && isFirstBatch
        nS = min(p, nE);
        if numel(r) < nS
            error('lsqvce:separable:acvTooShort', ...
                  'The autocovariance has %d lags but the startup needs %d.', ...
                  numel(r), nS);
        end
        G = toeplitz(r(1:nS));
        [Lc, fail] = chol(G, 'lower');
        if fail == 0
            Wp = Lc \ eye(nS);
            Yr(1:nS,:) = Wp * Xr(1:nS,:);
        else
            lsqvce.util.log.warn('lsqvce:separable:startupNotPD', ...
                'The stationary startup block is not PD; keeping the transient.');
        end
    end

    Yr = reshape(Yr, nE, d, nCol);                  % epoch, axis, col
    Yr = permute(Yr, [2 1 3]);                      % axis, epoch, col
    Y  = reshape(Yr, nE*d, nCol);                   % back to epoch-major rows
end
