function E = ellipse(D, pair, varargin)
%LSQVCE.ELLIPSE  Confidence ellipse for a pair of correlated parameters.
%
%   E = LSQVCE.ELLIPSE(D, [i j]) returns the confidence ellipse of parameters
%   i and j from the 2x2 submatrix of Qxx.
%
%   WHEN THE REPORT SAYS "look at the error ellipse instead", THIS IS IT.
%
%   Two parameters correlated at 0.99 may each have a small standard deviation
%   while being almost completely inseparable. The individual sigmas describe
%   the two SIDES of a bounding box; the ellipse inside it can be a needle. The
%   sigmas are then true but deeply misleading: they suggest you know each
%   parameter to that precision, when in fact you only know one combination of
%   them well and the orthogonal combination hardly at all.
%
%   The ellipse says what you actually know. Its semi-major axis is the
%   direction you know WORST, and its length is how badly.
%
%   INPUT
%     D     1x1 struct  diagnostics from lsqvce.adjust / lsqvce.diagnose; uses
%                       the parameter estimate and Qxx (cofactor) plus sigma0/s0.
%     pair  1x2 double  the two parameter indices [i j] to plot the ellipse for.
%
%   NAME-VALUE
%     'conf'     confidence level (default 0.95)
%     'scale'    'apriori' (use sigma0) | 'aposteriori' (use s0)  [default
%                'aposteriori', matching what the report tabulates]
%     'npts'     points on the returned outline (default 200)
%
%   OUTPUT E
%     .semiMajor, .semiMinor    axis half-lengths at the requested confidence
%     .angle                    orientation of the major axis, radians CCW from
%                               the +i axis
%     .ratio                    semiMajor/semiMinor -- the elongation. Above
%                               ~10 the parameters are effectively one.
%     .xy                       2 x npts outline, centred on the estimate
%     .centre                   [x_i; x_j]
%     .corr                     the correlation between them
%     .sigInd                   [sigma_i; sigma_j] -- the misleading numbers
%     .worstDir                 unit vector of the semi-major axis: the
%                               combination of the two parameters you know least
%
%   EXAMPLE
%       [~, D] = lsqvce.adjust(u, src, 4, 'verbose', false);
%       E = lsqvce.ellipse(D, [1 2]);            % parameters 1 and 2
%
%       plot(E.xy(:,1), E.xy(:,2), 'LineWidth', 1.5);  hold on
%       plot(E.centre(1), E.centre(2), 'r+');
%       axis equal; grid on
%
%       % E.xy         npts x 2, ready for plot()
%       % E.semiMajor  / E.semiMinor / E.angle
%       % E.ratio      elongation -- large means the pair is poorly separated
%
%       % Draw this whenever |rho| > 0.95: the individual sigmas are then
%       % meaningless in isolation, however small they look.
%
%   ERRORS
%     lsqvce:ellipse:badPair        the parameter pair is not two valid indices
%     lsqvce:ellipse:samePair       both indices are the same parameter
%
%   See also LSQVCE.DIAGNOSE.
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.

ip = inputParser;
addParameter(ip, 'conf',  0.95, @(x) isscalar(x) && x > 0 && x < 1);
addParameter(ip, 'scale', 'aposteriori', ...
             @(s) any(strcmpi(s, {'apriori', 'aposteriori'})));
addParameter(ip, 'npts',  200,  @(x) isscalar(x) && x >= 8);
parse(ip, varargin{:});
o = ip.Results;

if numel(pair) ~= 2 || any(pair < 1) || any(pair > D.u) || any(pair ~= fix(pair))
    error('lsqvce:ellipse:badPair', ...
          'pair must be two integer parameter indices in 1..%d.', D.u);
end
i = pair(1);  j = pair(2);
if i == j
    error('lsqvce:ellipse:samePair', 'The two indices must differ.');
end

Q = D.Qxx([i j], [i j]);
Q = (Q + Q') / 2;

if strcmpi(o.scale, 'aposteriori')
    fac = D.s0sq;
    if ~isfinite(fac) || fac <= 0, fac = D.sigma0^2; end
else
    fac = D.sigma0^2;
end
C = fac * Q;

[V, Lam] = eig(C);
lam = diag(Lam);
[lam, ix] = sort(lam, 'descend');
V = V(:, ix);

if any(lam < 0)
    lam = max(lam, 0);
end

% Chi-square with 2 dof gives the scale factor for a joint confidence region.
% (Using 1 dof would give the marginal interval on each axis, which is a
% DIFFERENT and smaller region -- the classic factor-of-2.45-vs-1.96 confusion.)
k2 = lsqvce.stat.chi2Inv(o.conf, 2);

a = sqrt(k2 * lam(1));      % semi-major
b = sqrt(k2 * lam(2));      % semi-minor
ang = atan2(V(2,1), V(1,1));

th = linspace(0, 2*pi, o.npts);
% npts x 2, NOT 2 x npts: this is the orientation plot() expects, so the caller
% can write plot(E.xy(:,1), E.xy(:,2)) directly.
%
% It WAS 2 x npts, which silently broke demo_workflow: E.xy(:,1) plotted a
% single POINT instead of the ellipse, and no error was raised. Caught only by
% running in MATLAB.
xy = (V * [a*cos(th); b*sin(th)]).';

cx = [D.x(i); D.x(j)];
si = sqrt(fac * [Q(1,1); Q(2,2)]);
rho = Q(1,2) / sqrt(max(Q(1,1)*Q(2,2), realmin));

E = struct('semiMajor', a, 'semiMinor', b, 'angle', ang, ...
           'ratio', a / max(b, realmin), ...
           'xy', xy + cx.', 'centre', cx, 'corr', rho, ...
           'sigInd', si, 'conf', o.conf, 'pair', [i j], ...
           'worstDir', V(:,1), 'bestDir', V(:,2), ...
           'lambda', lam, 'scale', o.scale);
end
