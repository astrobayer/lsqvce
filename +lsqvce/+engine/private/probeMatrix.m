function [raw, seedUsed] = probeMatrix(n, nProbes, dist, seed)
%PROBEMATRIX  The ONE raw probe generator (Rademacher / Gaussian).
%   [raw, seedUsed] = probeMatrix(n, nProbes, dist, seed)
%
%   Every stochastic probe block in the toolbox comes from THIS function
%   (WP-36, fixes F-06; rule 3): LSQVCE.ENGINE.TRACEESTIMATE,
%   LSQVCE.ENGINE.STREAMCONTRIBUTION and LSQVCE.ENGINE.SHAREDPROBE previously carried
%   three private generators with TWO different Rademacher recipes and,
%   with an empty seed, drew from -- and advanced -- the GLOBAL stream
%   while their help promised a private one. The one true contract:
%
%     seed given  -> RandStream('twister', 'Seed', seed): bit-reproducible,
%                    IDENTICAL probes across all callers for the same
%                    (n, nProbes, dist, seed).
%     seed empty  -> a concrete seed is DRAWN ONCE from the package-private
%                    session stream (PROBESTREAM; the global RNG is never
%                    touched or advanced) and then used exactly as above;
%                    it is returned in seedUsed so the caller can PIN it
%                    (reproducibility of "unseeded" runs, and the vce
%                    final-pass == loop-pass guarantee).
%
%   The Rademacher recipe is sign(randn) with sign(0) -> +1 -- the recipe
%   two of the three callers already used, kept as THE one so their seeded
%   realisations are unchanged.
%
%   INPUT
%     n        1x1 double   rows of the probe block (operator dimension)
%     nProbes  1x1 double   columns Z (number of probes)
%     dist     1x1 char     'rademacher' | 'gaussian'
%     seed     1x1 double   RNG seed, or [] for a private draw (see above)
%
%   OUTPUT
%     raw      n x nProbes double   UN-normalised probes (+-1, or N(0,1));
%                                   normalisation (e.g. eq. 6.5.43's
%                                   1/sqrt(Z)) stays with the caller
%     seedUsed 1x1 double           the seed actually used (== seed when
%                                   given; the private draw otherwise)
%
%   ERRORS
%     lsqvce:probeMatrix:badDist   dist is neither 'rademacher' nor
%                                  'gaussian'.
%
%   EXAMPLE
%     [Z1, s] = probeMatrix(60, 100, 'rademacher', []);   % private draw
%     Z2      = probeMatrix(60, 100, 'rademacher', s);    % pinned replay
%     % isequal(Z1, Z2) is true.
%
%   See also PROBESTREAM, LSQVCE.ENGINE.TRACEESTIMATE, LSQVCE.ENGINE.SHAREDPROBE,
%            LSQVCE.ENGINE.STREAMCONTRIBUTION
%
%   Created 2026-08-19 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.

    if isempty(seed)
        seedUsed = randi(probeStream(), [0, 2^31 - 2], 1, 1);
    else
        seedUsed = seed;
    end
    rs = RandStream('twister', 'Seed', seedUsed);

    switch lower(dist)
        case 'rademacher'
            raw = sign(randn(rs, n, nProbes));
            raw(raw == 0) = 1;                       % sign(0) -> +1
        case 'gaussian'
            raw = randn(rs, n, nProbes);
        otherwise
            error('lsqvce:probeMatrix:badDist', ...
                  'dist must be ''rademacher'' or ''gaussian'' (got ''%s'').', ...
                  char(string(dist)));
    end
end
