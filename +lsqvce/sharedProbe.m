function out = sharedProbe(U, varargin)
%LSQVCE.SHAREDPROBE  The once-per-iteration shared probe block UiZ = U^-1*Zbar.
%
%   out = LSQVCE.SHAREDPROBE(U)
%   out = LSQVCE.SHAREDPROBE(U, 'Name', value, ...)
%
%   Builds the single most-reused object in Ellmer's (2018) stochastic-model
%   estimation: the triangular solve U^-1 * Zbar of a Rademacher probe block
%   against the factored normal matrix N = U'*U.
%
%   WHY THIS IS ITS OWN FUNCTION (Ellmer p. 69, eq. 6.5.48)
%     U comes from the ACCUMULATED normal equations, so it is IDENTICAL for every
%     short arc in the iteration. Ellmer's central optimisation is to solve
%     U^-1 Zbar ONCE and reuse it for all arcs, all lags, and (eq. 6.5.61) the
%     geophysical redundancy too. Recomputing it per arc would dominate the cost.
%     Build it once at the top of each outer iteration and pass it down.
%
%   NORMALISATION
%     Zbar = Z / sqrt(Z) with Z the RAW +1/-1 Rademacher matrix and Z the probe
%     count (eq. 6.5.43), the same convention as LSQVCE.TRACEESTIMATE, so that
%     trace(X) ~ trace(Zbar' X Zbar) is unbiased.
%
%   INPUT
%     U    P x P double   upper-triangular factor of the normal matrix, N = U'*U
%                         (e.g. chol(N), or the QR triangular factor R of the
%                         whitened [A l], which IS this U -- see the note below).
%
%   OPTIONS (Name/Value)
%     'numProbes' 1x1 int      Z, number of probes (default 100, Ellmer's value)
%     'seed'      1x1 int|[]   RNG seed for reproducible probes (default []).
%                              Empty: a seed is drawn ONCE from the package's
%                              private session stream (probeMatrix contract,
%                              WP-36). Either way a private RandStream is used
%                              and the global RNG is
%                              never disturbed.
%     'probes'    P x Z double a FIXED raw +1/-1 probe block to use instead of
%                              drawing one (normalisation is applied here).
%                              Overrides 'numProbes'/'seed'.
%     'distribution' char      'rademacher' (default) | 'gaussian'
%
%   OUTPUT
%     out   1x1 struct
%       .UiZ     P x Z double   U^-1 * Zbar, the shared block to pass around
%       .Zbar    P x Z double   the normalised probe block itself
%       .probes  1x1            Z actually used
%       .P       1x1            parameter count
%
%   QR NOTE
%     With the QR backend the triangular factor R of the whitened [A l] already
%     satisfies N = R'*R, so R IS Ellmer's U -- pass it directly and skip forming
%     and re-factoring N. That is strictly better conditioned than chol(A'A).
%
%   EXAMPLE
%     sp  = lsqvce.sharedProbe(U, 'numProbes', 100, 'seed', 42);
%     for each arc
%         r = lsqvce.redundancyLags(Abar, W, sp.UiZ, ebar, 60);   % reuse sp.UiZ
%     end
%
%   ERRORS
%     lsqvce:sharedProbe:notSquare  the operator is not square
%     lsqvce:sharedProbe:badProbes  the probe count is not a positive integer
%
%   See also LSQVCE.REDUNDANCYLAGS, LSQVCE.TRACEESTIMATE,
%            LSQVCE.STREAMCONTRIBUTION
%
%   Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-36 (fixes F-06): probes come
%   from the ONE private generator probeMatrix (empty seed = private
%   session-stream draw, global RNG untouched; seeded realisations
%   unchanged). With assistance from Claude (Anthropic), model Claude
%   Fable 5.

    p = inputParser;  p.FunctionName = 'lsqvce.sharedProbe';
    p.addParameter('numProbes'   , 100, @(v) isscalar(v) && v >= 1 && v == floor(v));
    p.addParameter('seed'        , []);
    p.addParameter('probes'      , [], @(v) isnumeric(v));
    p.addParameter('distribution', 'rademacher', ...
                   @(s) any(strcmpi(s, {'rademacher','gaussian'})));
    p.parse(varargin{:});
    o = p.Results;

    P = size(U, 1);
    if size(U, 2) ~= P
        error('lsqvce:sharedProbe:notSquare', ...
              'U must be square (got %d x %d).', P, size(U, 2));
    end

    if ~isempty(o.probes)
        raw = o.probes;
        if size(raw, 1) ~= P
            error('lsqvce:sharedProbe:badProbes', ...
                  'The supplied probe block has %d rows but U is %d x %d.', ...
                  size(raw, 1), P, P);
        end
        Z = size(raw, 2);
    else
        Z = o.numProbes;
        % WP-36 (fixes F-06): the ONE generator; empty seed = a private
        % session-stream draw, the global RNG is never touched.
        raw = probeMatrix(P, Z, o.distribution, o.seed);
    end

    Zbar = raw / sqrt(Z);                    % eq. 6.5.43 normalisation
    UiZ  = U \ Zbar;                         % the shared triangular solve

    out = struct('UiZ', UiZ, 'Zbar', Zbar, 'probes', Z, 'P', P);
end

