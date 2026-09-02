function [A, Sigma, info] = whittle(G, p)
%LSQVCE.KERNEL.WHITTLE  Multichannel Levinson-Durbin: cross-ACV -> VAR(p).
%
%   [A, Sigma] = LSQVCE.KERNEL.WHITTLE(G, p)
%
%   Fits  z_t = sum_{k=1}^p A_k z_{t-k} + e_t,  C{e} = Sigma  to the given
%   cross-autocovariance sequence by the Whittle recursion -- the exact
%   multichannel generalisation of LSQVCE.KERNEL.LEVINSONDURBIN. On the EXACT ACV of a
%   VAR(p) it returns the generating coefficients to ~1e-9, and fitting a
%   higher order makes the extra matrices vanish (validation/v6verify.py
%   tests 3-4) -- the same order-sufficiency signature the scalar
%   LSQVCE.ENGINE.SELECTARORDER exploits.
%
%   INPUT
%     G   d x d x (>= p+1)   Gamma(0..), Gamma(j) = E[z_t z_{t-j}']
%     p   scalar             VAR order
%
%   OUTPUT
%     A       d x d x p      forward coefficient matrices
%     Sigma   d x d          forward innovation covariance
%     info    struct         .Sigmab backward innovation covariance
%                            .trace  trace(Sigma) per order 0..p (the
%                                    multichannel analogue of the scalar
%                                    prediction-variance curve)
%
%   ASYMMETRY IS REAL
%     Unlike the scalar case, Gamma(j)' ~= Gamma(j) in general, so forward
%     and backward coefficient sets differ and BOTH are carried through the
%     recursion. Feeding a symmetrised ACV of a genuinely cross-correlated
%     process silently destroys the lead/lag information -- exactly the part
%     that made a VAR necessary.
%
%   ERRORS
%     lsqvce:whittle:badG    the periodogram argument has the wrong shape
%     lsqvce:whittle:notPD   the fitted innovation covariance is not positive
%                            definite
%
%   See also LSQVCE.KERNEL.VARACV, LSQVCE.KERNEL.VARWHITEN, LSQVCE.KERNEL.LEVINSONDURBIN
%
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    narginchk(2,2);
    d = size(G,1);
    if size(G,2) ~= d || size(G,3) < p+1
        error('lsqvce:whittle:badG', ...
              'G must be d x d x (>= p+1); got %d x %d x %d for p = %d.', ...
              size(G,1), size(G,2), size(G,3), p);
    end

    A  = zeros(d, d, 0);
    B  = zeros(d, d, 0);
    Sf = (G(:,:,1) + G(:,:,1).')/2;
    Sb = Sf;
    tr = zeros(1, p+1);  tr(1) = trace(Sf);

    for n = 0:p-1
        D = G(:,:,n+2);
        for k = 1:n
            D = D - A(:,:,k) * G(:,:,n-k+2);
        end
        Kf = D / Sb;                       % D * inv(Sb)
        Kb = D.' / Sf;

        SfOld = Sf;
        An = zeros(d, d, n+1);  Bn = zeros(d, d, n+1);
        An(:,:,n+1) = Kf;       Bn(:,:,n+1) = Kb;
        for k = 1:n
            An(:,:,k) = A(:,:,k) - Kf * B(:,:,n-k+1);
            Bn(:,:,k) = B(:,:,k) - Kb * A(:,:,n-k+1);
        end
        Sf = Sf - Kf * Sb    * Kf.';
        Sb = Sb - Kb * SfOld * Kb.';
        Sf = (Sf + Sf.')/2;  Sb = (Sb + Sb.')/2;
        A = An;  B = Bn;
        tr(n+2) = trace(Sf);

        [~, cf] = chol(Sf);
        if any(~isfinite(Sf(:))) || cf > 0
            error('lsqvce:whittle:notPD', ...
                ['The forward innovation covariance lost positive ' ...
                 'definiteness at order %d. The ACV sequence is not a valid ' ...
                 'covariance (estimation noise, or the biased/unbiased ' ...
                 'normalisation was mixed).'], n+1);
        end
    end

    Sigma = Sf;
    if nargout > 2
        info = struct('Sigmab', Sb, 'trace', tr);
    end
end
