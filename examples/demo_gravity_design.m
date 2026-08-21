function [A, meta] = demo_gravity_design(nmax, orbit)
%DEMO_GRAVITY_DESIGN  Spherical-harmonic potential design along an orbit.
%   [A, meta] = demo_gravity_design(nmax, orbit)
%
%   The ONE physically structured design generator in this tree: the
%   partials of the (unitless, GM = R = 1) gravitational potential
%       V(r,th,lam) = (1/r) * sum_n (1/r)^n * sum_m Pbar_nm(cos th)
%                     * (C_nm cos(m lam) + S_nm sin(m lam))
%   with respect to every fully-normalised coefficient C_nm / S_nm of
%   degrees n = 2..nmax, evaluated at points sampled along a circular
%   inclined orbit with Earth rotation. Unlike the Gaussian fixtures, the
%   columns carry the REAL structure of a gravity design: order-block
%   correlation, near-polar-gap conditioning, and (1/r)^n degree
%   attenuation. Everything here is DETERMINISTIC -- no random numbers.
%
%   INPUTS
%     nmax    (1,1) double   maximum degree, integer >= 2.
%     orbit   (1,1) struct   .u      (M,1) double  argument of latitude [rad]
%                            .lam0   (M,1) double  inertial longitude of the
%                                                  ascending node minus Earth
%                                                  rotation angle [rad]
%                            .incl   (1,1) double  inclination [rad]
%                            .r      (1,1) double  orbit radius in Earth
%                                                  radii, > 1
%
%   OUTPUTS
%     A       (M,K) double   design matrix; K = (nmax+1)^2 - 4 columns in
%                            the order given by meta.n/meta.m/meta.cs.
%     meta    (1,1) struct   .n  (1,K) double  degree per column
%                            .m  (1,K) double  order per column
%                            .cs (1,K) char    'C' or 'S' per column
%
%   ACCURACY
%     Measured against the analytic SH algebra (HANDOVER 6.1.4, authority
%     row "Analytic SH algebra"): the fully-normalised Legendre recursion
%     satisfies the addition-theorem identity
%     Pbar_n0^2 + sum_{m=1..n} Pbar_nm^2 = 2n+1 to 4.3e-14 absolute at
%     nmax = 12 across near-pole-to-equator latitudes (MEASURED 2026-08-20,
%     user machine; the test asserts <= 5e-13, a decade of headroom on an
%     O(nmax)*eps recursion-growth bound), and the closed forms
%     Pbar_10 = sqrt(3) t, Pbar_11 = sqrt(3) sqrt(1-t^2),
%     Pbar_20 = sqrt(5) (3t^2-1)/2 to <= 1e-14 (test_physOracle asserts
%     both; the tolerances are eps-level growth of an O(nmax) recursion).
%
%   ERRORS
%     Input:
%       lsqvce:gravityDesign:badNmax   nmax is not an integer scalar >= 2
%       lsqvce:gravityDesign:badOrbit  orbit misses a field, u/lam0 differ
%                                      in length, or r <= 1
%
%   EXAMPLE
%     M = 200;  k = (0:M-1).';
%     orbit = struct('u', 0.05*k, 'lam0', -0.003*k, ...
%                    'incl', deg2rad(89), 'r', 1.07);
%     [A, meta] = demo_gravity_design(8, orbit);   % 200 x 77
%
%   LIMITATIONS
%     Potential values only (no gradients); degrees 0/1 excluded by design
%     (datum); unitless scaling GM = R = 1 -- the STRUCTURE is the point,
%     not geophysical units.
%
%   CONSUMERS
%     DEMO_GRAVITY_DATALOADER  the GRACE-like demo (three streams, per-system
%                              local counts, ragged arcs)
%     TEST_PHYSORACLE          the closed-loop oracle suite
%   Both call THIS function -- the Legendre recursion and its ACCURACY claim
%   exist once (WP-51; the test reaches examples/ the same way
%   TEST_ADJUSTVCE reaches DEMO_GRACE_DATALOADER).
%
%   See also DEMO_GRAVITY, DEMO_GRAVITY_DATALOADER, TEST_PHYSORACLE
%
%   Updated 2026-08-21 (CEST, UTC+02:00), WP-51: moved from tests/private to
%   examples/ and renamed, so the demo and the oracle suite share ONE
%   implementation instead of duplicating the recursion. With assistance
%   from Claude (Anthropic), model Claude Fable 5.
%   Created 2026-08-20 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.

    narginchk(2, 2);
    if ~isnumeric(nmax) || ~isscalar(nmax) || nmax < 2 || nmax ~= round(nmax)
        error('lsqvce:gravityDesign:badNmax', ...
              'nmax must be an integer scalar >= 2 (got %g).', nmax);
    end
    need = {'u', 'lam0', 'incl', 'r'};
    if ~isstruct(orbit) || ~all(isfield(orbit, need))
        error('lsqvce:gravityDesign:badOrbit', ...
              'orbit must carry the fields u, lam0, incl, r.');
    end
    u    = orbit.u(:);
    lam0 = orbit.lam0(:);
    if numel(u) ~= numel(lam0)
        error('lsqvce:gravityDesign:badOrbit', ...
              'orbit.u and orbit.lam0 must have equal length (%d vs %d).', ...
              numel(u), numel(lam0));
    end
    if ~isscalar(orbit.r) || orbit.r <= 1
        error('lsqvce:gravityDesign:badOrbit', ...
              'orbit.r must be a scalar > 1 Earth radius (got %g).', orbit.r);
    end

    M = numel(u);
    % geometry of a circular inclined orbit (spherical Earth):
    sinphi = sin(orbit.incl) .* sin(u);                 % sin(latitude)
    lam    = lam0 + atan2(cos(orbit.incl) .* sin(u), cos(u));
    t      = sinphi;                                    % Pbar argument

    K  = (nmax + 1)^2 - 4;                              % n = 2..nmax, C and S
    A  = zeros(M, K);
    mn = zeros(1, K);  mm = zeros(1, K);  mc = repmat('C', 1, K);

    P = legendreBar(nmax, t);                           % M x (n,m) table
    att = (1 / orbit.r) * (1 / orbit.r) .^ (0:nmax);    % (1/r)^(n+1)

    k = 0;
    for n = 2:nmax
        for m = 0:n
            Pnm = P{n + 1}(:, m + 1) * att(n + 1);
            k = k + 1;
            A(:, k) = Pnm .* cos(m * lam);
            mn(k) = n;  mm(k) = m;  mc(k) = 'C';
            if m > 0
                k = k + 1;
                A(:, k) = Pnm .* sin(m * lam);
                mn(k) = n;  mm(k) = m;  mc(k) = 'S';
            end
        end
    end
    meta = struct('n', mn, 'm', mm, 'cs', mc);
end

function P = legendreBar(nmax, t)
%LEGENDREBAR  Fully-normalised associated Legendre functions, all degrees.
%   P{n+1} is M x (n+1) holding Pbar_nm(t) for m = 0..n. Standard stable
%   forward-column recursion (sectorial seed, then upward in degree).
    M = numel(t);
    s = sqrt(1 - t.^2);
    P = cell(1, nmax + 1);
    P{1} = ones(M, 1);                                  % Pbar_00
    for n = 1:nmax
        P{n + 1} = zeros(M, n + 1);
        % sectorial Pbar_nn
        if n == 1
            P{2}(:, 2) = sqrt(3) * s;
        else
            P{n + 1}(:, n + 1) = sqrt((2*n + 1) / (2*n)) .* s ...
                                 .* P{n}(:, n);
        end
        % Pbar_n,n-1
        P{n + 1}(:, n) = sqrt(2*n + 1) .* t .* P{n}(:, n);
        % upward recursion for m = 0..n-2
        for m = 0:(n - 2)
            a = sqrt((2*n - 1) * (2*n + 1) / ((n - m) * (n + m)));
            b = sqrt((2*n + 1) * (n + m - 1) * (n - m - 1) / ...
                     ((n - m) * (n + m) * (2*n - 3)));
            P{n + 1}(:, m + 1) = a .* t .* P{n}(:, m + 1) ...
                                 - b .* P{n - 1}(:, m + 1);
        end
    end
end
