function [Aw, Bw, lw] = varWhiten(Avar, Sigma, G, d, A, B, l, isFirstBatch)
%LSQVCE.KERNEL.VARWHITEN  Whiten multi-axis observations against a vector-AR model.
%
%   [Aw,Bw,lw] = LSQVCE.KERNEL.VARWHITEN(Avar, Sigma, G, d, A, B, l, isFirstBatch)
%
%   The case LSQVCE.KERNEL.CROSSCORR flags as 'needs-var': axes cross-correlated IN
%   TIME (axis i at epoch k with axis j at epoch k-tau), which the separable
%   model cannot represent. The whitener is a MATRIX FIR:
%
%       w_t = L^-1 * ( z_t - sum_{k=1}^p A_k z_{t-k} ),   L*L' = Sigma
%
%   applied down epochs, plus an exact stationary startup for the first p
%   epochs (the block analogue of LSQVCE.KERNEL.ARMAWHITEN's startup rows). Cost
%   O(nEpoch * p * d^2) -- for d = 3 essentially the scalar cost times 9.
%
%   EXACTNESS
%     With the startup block, the operator equals the dense block-Cholesky
%     inverse chol(C)^-1 ROW FOR ROW (C the block-Toeplitz covariance):
%     relative difference ~1e-10, and the GLS solution through it matches the
%     dense solve to 1e-9 (validation/v6verify.py tests 5-6). Not an
%     approximation.
%
%   ROW ORDERING -- EPOCH-MAJOR, exactly as LSQVCE.KERNEL.SEPARABLEWHITEN
%       [x1 y1 z1 x2 y2 z2 ...]   and m must be a whole multiple of d.
%
%   BATCHING
%     On any batch after the first the FIR restarts from zero state; the
%     first p EPOCHS (p*d rows) are transient and must be overlapped and
%     discarded -- LSQVCE.ENGINE.BLOCKACCUMULATOR does this from arc.filterMemory = p*d and
%     rounds batches to whole epochs from arc.nAxis = d. Same contract as the
%     scalar filter path; same failure mode if the overlap is forced to 0.
%
%   INPUT
%     Avar    d x d x p     VAR coefficient matrices (LSQVCE.KERNEL.WHITTLE output)
%     Sigma   d x d         innovation covariance
%     G       d x d x >=p   Gamma(0..) for the startup block (LSQVCE.KERNEL.VARACV, or
%                           the empirical cross-ACV the model was fitted on)
%     d       scalar        number of axes
%     A, B    m x *         design blocks, epoch-major rows (may be m x 0)
%     l       m x 1         observations, epoch-major
%     isFirstBatch          true -> exact stationary startup for epochs 1..p
%
%   INPUT
%     Avar         d x d x p  vector-AR coefficient matrices
%     Sigma        d x d      innovation covariance
%     G            d x d x >=p  stationary cross-ACV (from LSQVCE.KERNEL.VARACV)
%     d            1x1        axes per epoch
%     A            m x nc     common design rows (EPOCH-MAJOR)
%     B            m x nl     local design rows (may be [] or sparse; sparse in, sparse out)
%     l            m x 1      observations
%     isFirstBatch 1x1 logical  apply the exact stationary startup
%
%   OUTPUT
%     Aw, Bw, lw   the same arrays, whitened against the vector-AR model
%
%   ERRORS
%     lsqvce:varWhiten:badA         the VAR coefficient array has the wrong shape
%     lsqvce:varWhiten:badSigma     the innovation covariance is not d x d
%     lsqvce:varWhiten:badStartup   the stationary startup block is malformed
%     lsqvce:varWhiten:notWholeEpochs  the batch does not contain whole epochs:
%                                   rows must be EPOCH-MAJOR and nAxis = d
%
%   See also LSQVCE.KERNEL.WHITTLE, LSQVCE.KERNEL.VARACV, LSQVCE.KERNEL.SEPARABLEWHITEN, LSQVCE.KERNEL.ARMAWHITEN
%
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+kernel (R4 decision R4-D1,
%   WP-64a): the KERNEL tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.
%
%   Sparse blocks: 2026-10-12 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- WP-123.

    narginchk(8,8);
    % WP-123: sparse blocks (makeLoader allows a sparse B) run the dense
    % whitener on their non-zero columns -- N-D reshape and filter refuse
    % sparse input. The dense path below is unchanged.
    if issparse(A) || issparse(B) || issparse(l)
        dense = @(A, B, l) lsqvce.kernel.varWhiten(Avar, Sigma, G, d, A, B, l, isFirstBatch);
        [Aw, Bw, lw] = activeColumns(dense, A, B, l);
        return
    end
    [dd, d2, p] = size(Avar);
    m = size(l,1);
    d = double(d);
    if dd ~= d || d2 ~= d
        error('lsqvce:varWhiten:badA', 'Avar must be %d x %d x p.', d, d);
    end
    if mod(m, d) ~= 0
        error('lsqvce:varWhiten:notWholeEpochs', ...
            ['The block has %d rows, not a multiple of d = %d. Rows must be ' ...
             'epoch-major and a batch must never split an epoch.'], m, d);
    end
    nE = m / d;

    % WP-43 (T-03a family): the help text promises badSigma for "not d x d",
    % but only the PD half was enforced -- a wrong-sized PD Sigma fell through
    % to a raw dimension error deeper in. The documented id is the contract.
    if ~isequal(size(Sigma), [d d])
        error('lsqvce:varWhiten:badSigma', ...
              'Sigma must be %d x %d; got %d x %d.', ...
              d, d, size(Sigma,1), size(Sigma,2));
    end

    [Lf, fail] = chol(Sigma, 'lower');
    if fail > 0
        error('lsqvce:varWhiten:badSigma', ...
              'Sigma is not positive definite.');
    end

    % One pass over the stacked matrix: whitening is a row operator.
    J  = [A, B, l(:)];
    q  = size(J, 2);
    Je = permute(reshape(J.', q, d, nE), [2 1 3]);   % d x q x nE epoch slabs
    W  = zeros(d, q, nE);

    t0 = 1;
    if isFirstBatch
        k0 = min(p, nE);
        Cs = blockToeplitz(G, k0);
        [Rs, f2] = chol(Cs);
        if f2 > 0
            error('lsqvce:varWhiten:badStartup', ...
                ['The p-epoch startup covariance is not positive definite. ' ...
                 'G is not a valid ACV sequence up to lag %d.'], k0-1);
        end
        Z0 = zeros(d*k0, q);
        for t = 1:k0
            Z0((t-1)*d + (1:d), :) = Je(:,:,t);
        end
        W0 = Rs.' \ Z0;                               % chol(Cs)^-1 * Z0
        for t = 1:k0
            W(:,:,t) = W0((t-1)*d + (1:d), :);
        end
        t0 = k0 + 1;
    end

    for t = t0:nE
        acc = Je(:,:,t);
        for k = 1:min(p, t-1)
            acc = acc - Avar(:,:,k) * Je(:,:,t-k);
        end
        W(:,:,t) = Lf \ acc;
    end

    Jw = zeros(m, q);
    for t = 1:nE
        Jw((t-1)*d + (1:d), :) = W(:,:,t);
    end

    na = size(A,2);  nb = size(B,2);
    Aw = Jw(:, 1:na);
    Bw = Jw(:, na + (1:nb));
    lw = Jw(:, end);
end

%==========================================================================
function C = blockToeplitz(G, nE)
%BLOCKTOEPLITZ  Dense symmetric block-Toeplitz covariance from a cross-ACV.
%   C = blockToeplitz(G, nE) assembles the (nE*d) x (nE*d) covariance of nE
%   consecutive epochs of a d-variate stationary process whose cross-
%   autocovariance is Gamma(k) = G(:,:,k+1), using Gamma(-k) = Gamma(k)'.
%
%   INPUT
%     G   d x d x (>=nE)  cross-ACV slabs; G(:,:,k+1) = Gamma(k), k = 0..nE-1
%     nE  1 x 1           number of epochs (block rows/cols)
%   OUTPUT
%     C   (nE*d) x (nE*d) symmetric block-Toeplitz covariance, epoch-major
    d = size(G,1);
    C = zeros(nE*d);
    for i = 1:nE
        for j = 1:nE
            k = i - j;
            if k >= 0, Bk = G(:,:,k+1);
            else,      Bk = G(:,:,-k+1).'; end
            C((i-1)*d + (1:d), (j-1)*d + (1:d)) = Bk;
        end
    end
end
