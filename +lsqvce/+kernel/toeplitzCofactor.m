function out = toeplitzCofactor(c, opt)
%LSQVCE.KERNEL.TOEPLITZCOFACTOR  Symmetric Toeplitz cofactor matrix from a lag vector.
%
%   Q     = LSQVCE.KERNEL.TOEPLITZCOFACTOR(c)                 (dense N x N)
%   apply = LSQVCE.KERNEL.TOEPLITZCOFACTOR(c, 'Form', 'apply')(function handle)
%   out   = LSQVCE.KERNEL.TOEPLITZCOFACTOR(c, 'Name', value, ...)
%
%   Builds the cofactor matrix Ellmer (2018, eq. 6.5.6) assigns to a short arc of
%   a wide-sense-stationary noise process: the symmetric Toeplitz matrix whose
%   first row/column is the discretised autocovariance function c = [C^0 ... C^{N-1}],
%
%       Q_ll = sum_{n=0}^{N-1} c(n) * V_n              (eq. 6.5.17)
%
%   with V_n the lag components (LSQVCE.ENGINE.COVCOMPONENT). All arcs of the same
%   observation type and length share this Q; shorter arcs use the leading
%   principal submatrix. Q is a COFACTOR (unit-variance shape); the arc's full
%   covariance is Sigma = sigma_m^2 * Q (eq. 6.5.7).
%
%   TWO FORMS
%     'dense' (default)  returns the explicit N x N matrix. Use for short arcs or
%                        when a Cholesky whitener chol(Q) is wanted next
%                        (LSQVCE.KERNEL.DECORRFACTOR).
%     'apply'            returns a handle q = @(x) Q*x that never forms Q, using
%                        the banded V_n applies (LSQVCE.ENGINE.COVCOMPONENT). Use for
%                        long arcs where N x N is too large.
%
%   INPUT
%     c    N x 1 double   autocovariance sequence C^n_xx, n = 0..N-1; c(1) = C^0
%                         is the variance (must be > 0). N >= 1.
%
%   OPTIONS (Name/Value)
%     'Form'  char        'dense' (default) | 'apply'
%     'Check' logical     if true (default), warn when c(1) <= 0 or |c(n)| > c(1)
%                         (a valid autocovariance peaks at lag 0). Set false to
%                         silence in inner loops.
%
%   OUTPUT (depends on 'Form')
%     Q      N x N double     the dense cofactor         ('dense')
%     apply  function_handle  q = @(x) Q*x, x is N x K   ('apply')
%
%   INVARIANTS (tested)
%     'dense' == toeplitz(c)   and is symmetric
%     apply(x) == toeplitz(c)*x   for any x
%     Q positive definite iff the sequence is a valid (PD) autocovariance;
%       LSQVCE.KERNEL.DECORRFACTOR floors non-PD sequences and reports it.
%
%   EXAMPLE
%     c = [1 0.5 0.2 0.05].';
%     Q = lsqvce.kernel.toeplitzCofactor(c);            % 4x4 symmetric Toeplitz
%     W = lsqvce.kernel.decorrFactor(Q);                % whitener chol(Q)^-1
%     % long arc, matrix-free:
%     q = lsqvce.kernel.toeplitzCofactor(c, 'Form','apply');
%     y = q(randn(4, 8));                        % = Q * X, never forms Q
%
%   ERRORS
%     lsqvce:toeplitzCofactor:empty  the autocovariance is empty
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:toeplitzCofactor:badVariance        [warning] C^0 = c(1) = X is not positive; a cofactor matrix needs a positive la...
%     lsqvce:toeplitzCofactor:notACV             [warning] |c(n)| exceeds c(1) at some lag; a valid autocovariance peaks at lag 0
%
%   See also LSQVCE.ENGINE.COVCOMPONENT, LSQVCE.KERNEL.DCTI, LSQVCE.KERNEL.DECORRFACTOR,
%            LSQVCE.KERNEL.PSD2ACV
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (c) c1 (F-07):
%   typed `arguments` signature, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+kernel (R4 decision R4-D1,
%   WP-64a): the KERNEL tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

    % WP-39 stage (c) c1: typed `arguments` signature (canonical
    % UpperCamelCase; legacy casings keep matching -- MATLAB NV
    % semantics; enum VALUES case-sensitive lowercase).
arguments
    c
    opt.Form (1,1) string {mustBeMember(opt.Form, ["dense","apply"])} = "dense"
    opt.Check (1,1) logical = true
end

    o = struct('Form', char(opt.Form), 'Check', opt.Check);

    c = c(:);
    N = numel(c);
    if N < 1
        error('lsqvce:toeplitzCofactor:empty', 'The lag vector c is empty.');
    end
    if logical(o.Check)
        if c(1) <= 0
            lsqvce.util.log.warn('lsqvce:toeplitzCofactor:badVariance', ...
                ['C^0 = c(1) = %.3g is not positive; a cofactor matrix needs a ' ...
                 'positive lag-0 variance.'], c(1));
        elseif any(abs(c) > c(1) * (1 + 1e-10))
            lsqvce.util.log.warn('lsqvce:toeplitzCofactor:notACV', ...
                ['|c(n)| exceeds c(1) at some lag; a valid autocovariance peaks ' ...
                 'at lag 0. The Toeplitz matrix may not be positive definite.']);
        end
    end

    if strcmpi(o.Form, 'apply')
        out = @(x) toeplitzApply(c, x);
    else
        out = toeplitz(c);               % symmetric: first row = first col = c
    end
end

%==========================================================================
function y = toeplitzApply(c, x)
%TOEPLITZAPPLY  y = toeplitz(c) * x via banded lag components (no dense matrix).
    N = numel(c);
    if isrow(x), x = x(:); end
    y = zeros(size(x));
    for n = 0:N-1
        if c(n+1) ~= 0
            y = y + c(n+1) * lsqvce.engine.covComponent(x, n);
        end
    end
end
