function con = contribution(acc, sol, opt)
%LSQVCE.ENGINE.CONTRIBUTION  Per-group information share of the solution (contribution).
%
%   con = LSQVCE.ENGINE.CONTRIBUTION(acc, sol)
%   con = LSQVCE.ENGINE.CONTRIBUTION(acc, sol, 'Name', value, ...)
%
%   Quantifies how much each system / observation-type / group determines the
%   common (gravity) solution -- the "contribution matrix" C_g = N^-1 N_g^cc, with
%   the defining property sum_g C_g = I. Two levels of detail:
%
%   LEVEL 1 -- the scalar SHARE (default, exact, cheap).
%     share_g = t_g / sum_g t_g with t_g = tr(C_g) = sum_{i in g} z_i' N^-1 z_i,
%     the common part of the leverage of the group's rows (z_i the row reduced
%     by the elimination of the locals): the fraction of the common solution
%     attributable to group g. Exact, from the same sweep as the redundancy
%     (LSQVCE.ENGINE.REDUNDANCY, field tcc); no probing. Over the data the
%     t_g sum to nc less what a common prior absorbs, tr(N^-1 K'K) (WP-131:
%     the help said nc), and the shares, relative to the DATA's part, to 1.
%
%   LEVEL 2 -- the per-parameter DIAGONAL (opt-in, streamed, approximate).
%     diag(C_g) = diag(N^-1 N_g^cc), one value per common parameter: the fraction
%     of information on THAT coefficient coming from group g (the interpretable
%     GRACE view -- e.g. which system pins the low degrees vs the high). Computed
%     matrix-free by LSQVCE.ENGINE.STREAMCONTRIBUTION (Hutchinson probes over one data
%     sweep), so it never forms N_g^cc. The diagonals sum to ones(nc,1) across
%     groups (exactly, up to probe noise, because sum_g N_g^cc = N) -- less
%     diag(N^-1 K'K) when a common prior carries a share. This is a
%     stochastic estimate: check con.group(k).relSE.
%
%   WEIGHT DEPENDENCE (read this)
%     N and every N_g^cc depend on the whitening, so the contribution is
%     WEIGHT-DEPENDENT and is only meaningful at the CONVERGED variance
%     components. Call this AFTER the VCE loop has settled (on the final acc/sol),
%     not per outer pass.
%
%   INPUT
%     acc   1x1 lsqvce.engine.BlockAccumulator   built with 'KeepLocal', true
%     sol   1x1 struct                    solve output carrying sol.x (nc x 1)
%
%   OPTIONS (Name/Value)
%   (WP-39 stage (b): canonical UpperCamelCase names below; legacy
%   casings still match -- MATLAB NV semantics -- and enum VALUES are
%   case-sensitive lowercase.)
%     'GroupBy'      char | fh   'system' (default) | 'obsType' | 'day' | 'arc' |
%                                @(sat,obsType,arc,day)->key. Same spec as
%                                LSQVCE.ENGINE.REDUNDANCY / LSQVCE.ENGINE.VCE.
%     'Level'        1x1 int     1 (share only, default) or 2 (also the diagonal).
%     'Probes'       1x1 int     Hutchinson probes for Level 2 (default 100).
%     'Seed'         1x1 int|[]  RNG seed for reproducible Level-2 probes ([]).
%     'MinRedundancy' 1x1        passed to LSQVCE.ENGINE.REDUNDANCY (default 50).
%     'Verbose'      1x1 logical progress banner (default false).
%
%   OUTPUT
%     con   1x1 struct with fields
%       .group   1xG struct, in GroupBy order:
%                  .key    the group key
%                  .t      1x1     exact common trace t_g = tr(N^-1 N_g^cc),
%                                  N_g^cc = Z_g'Z_g, Z_g = A_g - B_g G (the
%                                  rows reduced by the local elimination)
%                  .r      1x1     the group's redundancy r_g = sum(1 - h_i)
%                  .share  1x1     t_g / sum_g t_g (sums to 1); the Level-2
%                                  diagonal of the same group sums to t_g up to
%                                  the probe noise (.relSE) -- both levels use
%                                  this ONE split of N (WP-131)
%                  .diag   nc x 1  diag(N^-1 N_g^cc), each entry >= 0 in
%                                  expectation (Level 2 only, else [])
%                  .relSE  1x1     relative MC standard error of the streamed
%                                  trace (Level 2 only, else NaN)
%       .u        1x1        total parameters (m - redundancy = sum of r plus ...)
%       .rTotal   1x1        sum_g r_g: m - u_full without a common prior, m + nk -
%                            u_full - r_K with one (r_K the prior rows'
%                            redundancy, LSQVCE.ENGINE.REDUNDANCY; WP-134)
%       .level    1x1        1 or 2 (what was computed)
%       .nc       1x1        common-parameter count (Level 2; else [])
%       .seedUsed 1x1        the Level-2 probe seed (replays the run; [] at
%                            Level 1) -- WP-134
%       .note     char       provenance / weight-dependence reminder
%
%   EXAMPLE
%     R   = lsqvce.fit(data, Estimate="noise");   % converged; R.acc, R.sol
%     con = lsqvce.engine.contribution(R.acc, R.sol, 'GroupBy','obsType', 'Level',2);
%     [con.group.share]                           % scalar share per stream
%     bar(con.group(1).diag);                     % per-coefficient share, group 1
%     % (R.contribution(...) is the same call on the Results)
%
%   See also LSQVCE.ENGINE.STREAMCONTRIBUTION, LSQVCE.ENGINE.REDUNDANCY, LSQVCE.ENGINE.VCE,
%            LSQVCE.RESULTS
%
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-134 (review 3 lows): .seedUsed
%   (Level 2) and .rTotal's value with a common prior documented. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-131: Level 2 streams the same
%   (leverage) split as Level 1 -- before, its half-Schur split differed per
%   group whenever a system's locals spanned several groups. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-125: the Level-1 share is the
%   common trace t_g = tr(inv(N) N_g^cc), not the redundancy. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (b) (F-07): typed
%   `arguments` signature, canonical UpperCamelCase option names (legacy
%   casings keep matching; enum VALUES case-sensitive lowercase). With
%   assistance from Claude (Anthropic), model Claude Fable 5.
%   Created 2026-07-17 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.
arguments
    acc
    sol
    opt.GroupBy = 'system'
    opt.Level (1,1) double {mustBeMember(opt.Level, [1 2])} = 1
    opt.Probes (1,1) double {mustBeInteger, mustBePositive} = 100
    opt.Seed double {mustBeScalarOrEmpty} = []
    opt.MinRedundancy (1,1) double {mustBePositive} = 50
    opt.Verbose (1,1) logical = false
end

    % WP-39 stage (b): typed `arguments` signature (names already
    % UpperCamelCase; NAME matching case-insensitive per MATLAB semantics).
    o = opt;

    % ---- Level 1: exact per-group common-trace share ----------------------
    red = lsqvce.engine.redundancy(acc, sol, 'GroupBy', o.GroupBy, ...
                            'MinRedundancy', o.MinRedundancy, ...
                            'Collect', false, 'Verbose', logical(o.Verbose));
    rKeys = {red.group.key};
    rVals = [red.group.r];
    rTot  = red.rTotal;
    % WP-125: the share is the COMMON TRACE t_g (what the group determines),
    % not the redundancy r_g (what it leaves undetermined) -- the two are
    % nearly complementary, so the redundancy share ranked the groups upside
    % down on the information they carry
    tVals = [red.group.tcc];
    tTot  = sum(tVals);
    if tTot > 0
        shares = tVals / tTot;
    else
        shares = nan(size(tVals));
    end

    G = struct('key', {}, 't', {}, 'r', {}, 'share', {}, 'diag', {}, 'relSE', {});
    for k = 1:numel(rKeys)
        G(k) = struct('key', rKeys{k}, 't', tVals(k), 'r', rVals(k), ...
                      'share', shares(k), 'diag', [], 'relSE', NaN); %#ok<AGROW>
    end

    ncOut = [];  seedUsed = [];
    if o.Level >= 2
        % ---- Level 2: streamed per-parameter diagonal --------------------
        % Pin the group order to the Level-1 keys so diag aligns with share.
        sc = lsqvce.engine.streamContribution(acc, sol, 'GroupBy', o.GroupBy, ...
                    'Groups', {red.group.key}, 'Probes', o.Probes, ...
                    'Seed', o.Seed, 'WantDiag', true, 'Verbose', logical(o.Verbose));
        ncOut = sc.nc;  seedUsed = sc.seedUsed;
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
                 'level', o.Level, 'nc', ncOut, 'seedUsed', seedUsed, ...
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
