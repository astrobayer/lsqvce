function con = contribution(acc, sol, varargin)
%LSQVCE.CONTRIBUTION  Per-group information share of the solution (contribution).
%
%   con = LSQVCE.CONTRIBUTION(acc, sol)
%   con = LSQVCE.CONTRIBUTION(acc, sol, 'Name', value, ...)
%
%   Quantifies how much each system / observation-type / group determines the
%   common (gravity) solution -- the "contribution matrix" C_g = N^-1 N_g^cc, with
%   the defining property sum_g C_g = I. Two levels of detail:
%
%   LEVEL 1 -- the scalar SHARE (default, exact, cheap).
%     share_g = tr(C_g) / u_red = r_g / sum_g r_g, the fraction of the effectively
%     determined parameters attributable to group g. r_g is the EXACT per-group
%     redundancy from LSQVCE.REDUNDANCY (r_g = sum_{i in g}(1 - h_i)); no probing,
%     no approximation. The shares sum to 1.
%
%   LEVEL 2 -- the per-parameter DIAGONAL (opt-in, streamed, approximate).
%     diag(C_g) = diag(N^-1 N_g^cc), one value per common parameter: the fraction
%     of information on THAT coefficient coming from group g (the interpretable
%     GRACE view -- e.g. which system pins the low degrees vs the high). Computed
%     matrix-free by LSQVCE.STREAMCONTRIBUTION (Hutchinson probes over one data
%     sweep), so it never forms N_g^cc. The diagonals sum to ones(nc,1) across
%     groups (exactly, up to probe noise, because sum_g N_g^cc = N). This is a
%     stochastic estimate: check con.group(k).relSE.
%
%   WEIGHT DEPENDENCE (read this)
%     N and every N_g^cc depend on the whitening, so the contribution is
%     WEIGHT-DEPENDENT and is only meaningful at the CONVERGED variance
%     components. Call this AFTER the VCE loop has settled (on the final acc/sol),
%     not per outer pass. In LSQVCE.ADJUSTVCE the Level-1 share is already tracked
%     per pass in R.history.contribution; this function is the standalone /
%     Level-2 door.
%
%   INPUT
%     acc   1x1 lsqvce.BlockAccumulator   built with 'KeepLocal', true
%     sol   1x1 struct                    solve output carrying sol.x (nc x 1)
%
%   OPTIONS (Name/Value)
%     'GroupBy'      char | fh   'system' (default) | 'obsType' | 'day' | 'arc' |
%                                @(sat,obsType,arc,day)->key. Same spec as
%                                LSQVCE.REDUNDANCY / LSQVCE.VCE.
%     'Level'        1x1 int     1 (share only, default) or 2 (also the diagonal).
%     'Probes'       1x1 int     Hutchinson probes for Level 2 (default 100).
%     'Seed'         1x1 int|[]  RNG seed for reproducible Level-2 probes ([]).
%     'MinRedundancy' 1x1        passed to LSQVCE.REDUNDANCY (default 50).
%     'Verbose'      1x1 logical progress banner (default false).
%
%   OUTPUT
%     con   1x1 struct with fields
%       .group   1xG struct, in GroupBy order:
%                  .key    the group key
%                  .r      1x1     exact redundancy r_g = tr(N^-1 N_g^cc)
%                  .share  1x1     r_g / sum_g r_g  (Level-1 scalar share, sums 1)
%                  .diag   nc x 1  diag(N^-1 N_g^cc) (Level 2 only, else [])
%                  .relSE  1x1     relative MC standard error of the streamed
%                                  trace (Level 2 only, else NaN)
%       .u        1x1        total parameters (m - redundancy = sum of r plus ...)
%       .rTotal   1x1        sum_g r_g = m - u_full (the exact denominator)
%       .level    1x1        1 or 2 (what was computed)
%       .nc       1x1        common-parameter count (Level 2; else [])
%       .note     char       provenance / weight-dependence reminder
%
%   EXAMPLE
%     R   = P.run();                              % adjustVCE, converged
%     acc = R.history ... ;                       % or a direct BlockAccumulator
%     con = lsqvce.contribution(acc, sol, 'GroupBy','obsType', 'Level',2);
%     [con.group.share]                           % scalar share per stream
%     bar(con.group(1).diag);                     % per-coefficient share, group 1
%
%   See also LSQVCE.STREAMCONTRIBUTION, LSQVCE.REDUNDANCY, LSQVCE.VCE,
%            LSQVCE.ADJUSTVCE
%
%   Created 2026-07-17 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.

    p = inputParser;
    p.addParameter('GroupBy'      , 'system');
    p.addParameter('Level'        , 1, @(v) isscalar(v) && any(v == [1 2]));
    p.addParameter('Probes'       , 100, @(v) isscalar(v) && v >= 1 && v == floor(v));
    p.addParameter('Seed'         , []);
    p.addParameter('MinRedundancy', 50);
    p.addParameter('Verbose'      , false, @(v) isscalar(v));
    p.parse(varargin{:});
    o = p.Results;

    % ---- Level 1: exact per-group redundancy share -----------------------
    red = lsqvce.redundancy(acc, sol, 'GroupBy', o.GroupBy, ...
                            'MinRedundancy', o.MinRedundancy, ...
                            'Collect', false, 'Verbose', logical(o.Verbose));
    rKeys = {red.group.key};
    rVals = [red.group.r];
    rTot  = red.rTotal;
    if rTot > 0
        shares = rVals / rTot;
    else
        shares = nan(size(rVals));
    end

    G = struct('key', {}, 'r', {}, 'share', {}, 'diag', {}, 'relSE', {});
    for k = 1:numel(rKeys)
        G(k) = struct('key', rKeys{k}, 'r', rVals(k), 'share', shares(k), ...
                      'diag', [], 'relSE', NaN); %#ok<AGROW>
    end

    ncOut = [];
    if o.Level >= 2
        % ---- Level 2: streamed per-parameter diagonal --------------------
        % Pin the group order to the Level-1 keys so diag aligns with share.
        sc = lsqvce.streamContribution(acc, sol, 'GroupBy', o.GroupBy, ...
                    'Groups', {red.group.key}, 'Probes', o.Probes, ...
                    'Seed', o.Seed, 'WantDiag', true, 'Verbose', logical(o.Verbose));
        ncOut = sc.nc;
        scKeys = cellfun(@(k) char(string(k)), {sc.group.key}, ...
                         'UniformOutput', false);
        for k = 1:numel(G)
            kk = char(string(G(k).key));
            [tf, j] = ismember(kk, scKeys);
            if tf
                G(k).diag  = sc.group(j).diag;
                G(k).relSE = sc.group(j).relSE;
            end
        end
    end

    con = struct('group', G, 'u', red.u, 'rTotal', rTot, ...
                 'level', o.Level, 'nc', ncOut, ...
                 'note', ['contribution is weight-dependent -- meaningful only ' ...
                          'at the converged variance components (compute once, ' ...
                          'at the end). Level 1 share is exact; Level 2 diagonal ' ...
                          'is a Hutchinson estimate (check relSE).']);

    if o.Verbose
        lsqvce.util.log.msg(['%s%d group(s), share sums to %.4f' ...
            '%s\n'], lsqvce.util.verbosity.banner('Contribution'), numel(G), sum(shares(isfinite(shares))), ...
            ternaryStr(o.Level >= 2, sprintf(', diag over %d probe(s)', o.Probes), ''));
    end
end

function s = ternaryStr(c, a, b)
%TERNARYSTR  Inline string conditional: a if c else b.
%   INPUT  c 1x1 logical; a,b char   OUTPUT s char
    if c, s = a; else, s = b; end
end
