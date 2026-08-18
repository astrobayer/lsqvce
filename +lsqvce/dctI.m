function [y, X] = dctI(x, varargin)
%LSQVCE.DCTI  Type-I discrete cosine transform (covariance <-> PSD pairing).
%
%   y      = LSQVCE.DCTI(x)
%   y      = LSQVCE.DCTI(x, 'Name', value, ...)
%   [y, X] = LSQVCE.DCTI(...)
%
%   The discrete transform Ellmer (2018, eqs. 6.5.12-6.5.14) uses to pair the
%   discretised autocovariance function C^n_xx with the power spectral density
%   S^j_xx of a wide-sense-stationary noise process:  S = X C  and  C = X^-1 S.
%   For the ORTHONORMAL normalisation (the default) X is self-inverse (X^-1 = X),
%   so the same call inverts itself.
%
%   THE MATRIX (Rao & Yip 1990; Ellmer eq. 6.5.12)
%     X_mn = sqrt(2/(N-1)) * k_m k_n cos( m n pi / (N-1) ),   m,n in [0,N),
%   with the endpoint weights k_i selected by 'Norm':
%     'ortho'  (default)  k_0 = k_{N-1} = 1/sqrt(2), else 1  -> X orthonormal,
%                         self-inverse (X X = I). Numerically sound; standard
%                         type-I DCT (== scipy dct type=1, norm='ortho').
%     'ellmer'            k_0 = 1/sqrt(2), else 1 (the thesis as literally
%                         written). NOT orthonormal but well-conditioned
%                         (cond ~ 1.5); inverse is via 'Inverse',true (uses X^-1).
%                         Use for bit-exact ITSG-Grace2016 reproduction. The two
%                         normalisations differ only in the Nyquist-bin scaling
%                         (a sqrt(2) factor on the highest-frequency amplitude).
%
%   INPUT
%     x    N x K double   K columns of length-N sequences to transform (each
%                         column independently). N >= 2. A row vector is treated
%                         as one column. For the covariance->PSD direction pass
%                         the autocovariance C^n_xx (n = 0..N-1) as x.
%
%   OPTIONS (Name/Value)
%     'Norm'     char     'ortho' (default) | 'ellmer'  (see above)
%     'Inverse'  logical  false (default) forward S = X C; true applies the
%                         inverse C = X^-1 S. For 'ortho' the inverse equals the
%                         forward (X is self-inverse) so this is a no-op flag
%                         there; for 'ellmer' it forms and applies X^-1.
%
%   OUTPUT
%     y    N x K double   the transform of each column of x
%     X    N x N double   the transform matrix actually applied (X for forward,
%                         X^-1 for inverse). Handy to reuse across many columns
%                         or to cache; forming it is O(N^2).
%
%   INVARIANTS (tested)
%     'ortho':  X = X.',  X*X = I  (self-inverse, symmetric)
%     round-trip: dctI(dctI(x)) = x for 'ortho'; and for 'ellmer'
%                 dctI(dctI(x,'Norm','ellmer'),'Norm','ellmer','Inverse',true)=x
%
%   EXAMPLE
%     C = [1.0 0.5 0.2 0.05 0.0].';        % an autocovariance sequence
%     S = lsqvce.dctI(C);                  % its type-I DCT (a PSD proxy)
%     C2 = lsqvce.dctI(S);                 % ortho is self-inverse -> C2 == C
%     max(abs(C - C2))                     % ~ 1e-15
%
%   ERRORS
%     lsqvce:dctI:tooShort          the series is shorter than the transform needs
%
%   See also LSQVCE.PSD2ACV, LSQVCE.TOEPLITZCOFACTOR, LSQVCE.COVCOMPONENT
%
%   Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.

    p = inputParser;  p.FunctionName = 'lsqvce.dctI';
    p.addParameter('Norm'   , 'ortho', @(s) any(strcmpi(s, {'ortho','ellmer'})));
    p.addParameter('Inverse', false,   @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
    p.parse(varargin{:});
    o = p.Results;
    inv = logical(o.Inverse);

    wasRow = isrow(x);
    if wasRow, x = x(:); end
    N = size(x, 1);
    if N < 2
        error('lsqvce:dctI:tooShort', ...
              'The transform needs at least N = 2 samples (got %d).', N);
    end

    X = dctIMatrix(N, o.Norm);
    if inv
        if strcmpi(o.Norm, 'ortho')
            A = X;                       % orthonormal: X^-1 = X
        else
            A = invWellCond(X);          % ellmer: explicit inverse
        end
    else
        A = X;
    end

    y = A * x;
    if wasRow, y = y.'; end
    if nargout > 1, X = A; end
end

%==========================================================================
function X = dctIMatrix(N, nrm)
%DCTIMATRIX  The N x N type-I DCT matrix with the chosen endpoint weighting.
    m = (0:N-1).';
    n = (0:N-1);
    k = ones(N, 1);
    k(1) = 1/sqrt(2);                    % index 0
    if strcmpi(nrm, 'ortho')
        k(N) = 1/sqrt(2);               % index N-1 too -> orthonormal
    end
    X = sqrt(2/(N-1)) * (k * k.') .* cos(m * n * pi / (N-1));
    X = (X + X.')/2;                     % enforce exact symmetry (rounding)
end

%==========================================================================
function Xi = invWellCond(X)
%INVWELLCOND  Inverse of the (non-orthonormal but well-conditioned) ellmer X.
%   X is cond ~ 1.5, so a direct solve is safe; use \ against the identity.
    Xi = X \ eye(size(X, 1));
end
