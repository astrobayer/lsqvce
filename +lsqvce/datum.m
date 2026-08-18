function [B, info] = datum(N, varargin)
%LSQVCE.DATUM  Detect a datum defect and construct the datum matrix B.
%
%   [B, info] = LSQVCE.DATUM(N) inspects a (possibly rank-deficient) normal
%   matrix and returns a basis B for its null space -- the datum. The columns
%   of B are the transformations the observations CANNOT see: for a levelling
%   net, an overall height shift; for a 3-D geodetic net, 3 translations (+3
%   rotations, + scale if unconstrained).
%
%   [B, info] = LSQVCE.DATUM(N, 'datumB', Buser) validates a datum matrix you
%   supply. ALWAYS PREFER THIS. You know what your datum is; the eigenvalue
%   threshold used by the automatic detector does not, and can be fooled --
%   in testing, a regularisation of only 1e-10*||N|| was enough to hide a
%   genuine defect from it. Automatic detection is a fallback that announces
%   itself as a guess.
%
%   NAME-VALUE
%     'datumB'  u x d  datum matrix. If given, it is CHECKED (||NB||/||N|| must
%                      be small) rather than re-derived.
%     'tol'     relative eigenvalue threshold for auto-detection (default 1e-10)
%     'verbose' announce what was found (default true)
%
%   OUTPUT info
%     .defect     d, the dimension of the null space (0 = full rank)
%     .rank       u - d
%     .method     'user' | 'auto' | 'none'
%     .evRatio    eigenvalues of N normalised by the largest
%     .residual   ||N*B|| / ||N||, i.e. how well B really spans the null space
%
%   WHY THIS MATTERS
%   ----------------
%   On a free network the redundancy is n - rank(A) = n - (u - d), NOT n - u.
%   Getting this wrong does not merely mis-state a degree of freedom: it makes
%   the identity sum(r_i) = n - u fail, and it biases the variance factor by a
%   factor of (n-u+d)/(n-u). On a small levelling net (n=14, u=6, d=1) that is
%   a 12.5% error in s0^2 -- easily enough to fail a global test that should
%   have passed, and to send you hunting for outliers that do not exist.
%
%   INPUT
%     N         u x u  numeric  the normal matrix (may be singular)
%     varargin  Name/Value: 'tol' (rank tolerance), 'type' ('auto' | explicit)
%
%   OUTPUT
%     B     u x k  numeric  basis of the null space -- the datum defect
%     info  1x1    struct   .rank, .defect, .tol, .detected (true when the
%                           defect was found automatically rather than supplied)
%
%   EXAMPLE
%       % A levelling network measures only HEIGHT DIFFERENCES, so the overall
%       % level is unobservable: exactly one datum defect.
%       u = 6;  edges = [1 2;2 3;3 4;4 5;5 6;1 3;2 5;4 6];
%       A = zeros(size(edges,1), u);
%       for i = 1:size(edges,1)
%           A(i, edges(i,1)) = -1;  A(i, edges(i,2)) = 1;
%       end
%
%       [B, info] = lsqvce.datum(A'*A);
%       % info.defect -> 1        B -> spans ones(u,1), the constant shift
%
%       [sol, D] = lsqvce.adjust(u, src, nB, 'datumB', B);
%
%       % B must span null(N) -- an INNER constraint. "Fix station 3"
%       % (B = e_3) is a HARD constraint and is REFUSED: that is a
%       % reparameterisation (drop the column), not a datum.
%
%   ERRORS
%     lsqvce:datum:badB             B is not a numeric matrix of the right width
%     lsqvce:datum:rankDeficientB   B does not have full column rank, so it is
%                                   not a basis
%     lsqvce:datum:zeroN            the normal matrix is empty or all zero
%
%   See also LSQVCE.ACCUMULATOR, LSQVCE.ADJUST.
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.

ip = inputParser;
ip.FunctionName = 'lsqvce.datum';
addRequired (ip, 'N', @(x) isnumeric(x) && ismatrix(x) && size(x,1) == size(x,2));
addParameter(ip, 'datumB', [],    @(x) isempty(x) || isnumeric(x));
addParameter(ip, 'tol', 1e-10, @(x) isscalar(x) && x > 0);
addParameter(ip, 'verbose', true, @(x) islogical(x) && isscalar(x));
parse(ip, N, varargin{:});
o = ip.Results;

u = size(N, 1);
N = (N + N') / 2;                    % symmetrise
nrm = norm(N, 'fro');
if nrm == 0
    error('lsqvce:datum:zeroN', 'The normal matrix is identically zero.');
end

info = struct('defect', 0, 'rank', u, 'method', 'none', ...
              'evRatio', [], 'residual', 0);

% ---------------------------------------------------------------- user datum
if ~isempty(o.datumB)
    B = o.datumB;
    if size(B, 1) ~= u
        error('lsqvce:datum:badB', ...
              'B has %d rows but there are %d parameters.', size(B,1), u);
    end
    d = size(B, 2);

    % Orthonormalise. The caller's B need not be orthonormal (a levelling
    % datum is usually just ones(u,1)), but every downstream formula assumes
    % a well-conditioned basis, and the pseudoinverse correction
    % B*inv(B'B)^2*B' is numerically nasty if B'B is ill-conditioned.
    [Q, Rb] = qr(B, 0);
    if min(abs(diag(Rb))) < eps(max(abs(diag(Rb)))) * max(d, 1) * 10
        error('lsqvce:datum:rankDeficientB', ...
              ['The supplied datum matrix B is itself rank deficient: its ' ...
               '%d columns are not independent. Each column must be a ' ...
               'distinct datum transformation.'], d);
    end
    B = Q;

    resid = norm(N * B, 'fro') / nrm;
    info.residual = resid;
    info.defect = d;
    info.rank   = u - d;
    info.method = 'user';

    if resid > 1e-6
        lsqvce.util.log.warn('lsqvce:datum:BnotNull', ...
            ['The supplied datum matrix does NOT span the null space of N: ' ...
             '||N*B||/||N|| = %.2e, which should be ~0.\n' ...
             'Either B is wrong, or the network is not actually free in those ' ...
             'directions (perhaps you already applied constraints). ' ...
             'Constraining a direction the observations DO see will bias the ' ...
             'solution, not just fix its origin.'], resid);
    end

    if lsqvce.util.verbosity.progress(o.verbose)
        lsqvce.util.log.msg(['Datum: %d supplied constraint(s), residual ||NB||/||N|| = ' ...
                 '%.2e. Redundancy will use n - %d.\n'], d, resid, u - d);
    end
    return;
end

% ---------------------------------------------------------------- auto-detect
[V, Dg] = eig(N);
ev = real(diag(Dg));
[ev, ix] = sort(ev, 'ascend');
V = V(:, ix);

evr = ev / max(abs(ev));
info.evRatio = evr;

d = sum(evr < o.tol);
if d == 0
    B = zeros(u, 0);
    info.method = 'auto';
    if lsqvce.util.verbosity.progress(o.verbose)
        lsqvce.util.log.msg('Datum: N is full rank (smallest eigenvalue ratio %.2e).\n', ...
                evr(1));
    end
    return;
end

B = V(:, 1:d);
info.defect   = d;
info.rank     = u - d;
info.method   = 'auto';
info.residual = norm(N * B, 'fro') / nrm;

% Is the gap between the "null" and "non-null" eigenvalues actually clear? If
% the spectrum decays smoothly there is no honest cut point, and reporting a
% definite defect would be false precision.
if d < u
    gap = evr(d+1) / max(evr(d), realmin);
else
    gap = Inf;
end

lsqvce.util.log.warn('lsqvce:datum:autoDetected', ...
    ['Datum defect was DETECTED AUTOMATICALLY, not supplied. Found d = %d ' ...
     'null direction(s) using a relative eigenvalue threshold of %.0e.\n' ...
     'This is a guess. The eigen-gap is %.1e (large = confident, ~1 = the ' ...
     'spectrum decays smoothly and there is no honest cut point).\n' ...
     'If you know your datum -- and you usually do, it is the set of ' ...
     'transformations your observations cannot see -- pass it explicitly as ' ...
     '''datumB'' and this guess will not be used.'], d, o.tol, gap);

if lsqvce.util.verbosity.progress(o.verbose)
    lsqvce.util.log.msg(['Datum: auto-detected %d defect(s). Eigenvalue ratios: %s\n' ...
             '       Redundancy will use n - %d.\n'], ...
             d, mat2str(evr(1:min(d+2, u)).', 3), u - d);
end
end
