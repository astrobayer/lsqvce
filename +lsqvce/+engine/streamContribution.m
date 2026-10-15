function out = streamContribution(acc, sol, opt)
%LSQVCE.ENGINE.STREAMCONTRIBUTION  Streaming, matrix-free per-group N^-1 N_g^cc probes.
%
%   out = LSQVCE.ENGINE.STREAMCONTRIBUTION(acc, sol)
%   out = LSQVCE.ENGINE.STREAMCONTRIBUTION(acc, sol, 'Name', value, ...)
%
%   The shared primitive behind BOTH the streaming Hutchinson VCE trace and the
%   Level-2 contribution diagonal. For every group g it estimates, by Hutchinson
%   probing over a SINGLE data sweep and WITHOUT ever forming N_g^cc,
%
%       t_g       = tr( N^-1 N_g^cc )              scalar contribution / redundancy
%       diag(C_g) = diag( N^-1 N_g^cc )   (nc x 1) per-parameter information share
%
%   where N = Ncc_red = sum_g N_g^cc is the TOTAL reduced common normal (already
%   factored once as chol) and N_g^cc = Z_g'Z_g is group g's contribution to it,
%   Z_g = A_g - B_g G the group's rows reduced by the local elimination -- the
%   LEVERAGE split, the one Level 1 uses (WP-131), so t_g is in expectation the
%   group's common leverage sum_{i in g} z_i' N^-1 z_i and the diagonal sums to
%   it. The group-summed invariants hold EXACTLY (up to probe noise) without a
%   common prior: sum_g t_g = nc and sum_g diag(C_g) = ones(nc,1), because
%   sum_g N_g^cc = N. N includes a common prior's normals K'K / sigma^2, which
%   belong to no data group: then sum_g t_g = nc - tr(N^-1 K'K) / sigma^2
%   (WP-134: the help claimed nc regardless).
%
%   MEMORY-HONEST DESIGN (why this exists)
%     N_g^cc is nc x nc and a group can span many arcs, so holding all G of them
%     is O(G*nc^2) and destroys the streaming property that makes adjustBlock
%     worth using. We never store a group normal. What Hutchinson needs is only
%     the ACTION N_g^cc*Z on an nc x Z probe block, and that action is exactly
%     group g's contribution to a mat-vec, re-accumulated by STREAMING that
%     group's arcs -- mirroring the exact leverage path in LSQVCE.ENGINE.REDUNDANCY:
%
%       N_g^cc Z = sum_{arc in g} A_arc'(A_arc Z) - Nlc_g'(G_s Z) - G_s'(Nlc_g Z)
%                  + G_s' Mll_g (G_s Z)
%
%     with A_arc the whitened common design of the arc (m x nc), Nlc_g / Mll_g
%     the group's cross block B'A and local Gram B'B (streamed with it), and
%     G_s = inv(Nll_s) Nlc_s the per-SYSTEM Schur factor (built once, identical
%     to the G of LSQVCE.ENGINE.BLOCKCOFACTOR). All groups are
%     accumulated in ONE sweep -- each arc routes its A'(A Z) into its own group's
%     nc x Z accumulator, so the cost is a handful of data passes independent of G.
%
%     N^-1 is then applied to the assembled N_g^cc Z through the shared Cholesky
%     Lc of N (two triangular solves), never by forming N^-1.
%
%   RELATION TO THE EXACT PATH (read this before trusting the numbers)
%     The exact per-arc leverage path (LSQVCE.ENGINE.REDUNDANCY, LSQVCE.ENGINE.BLOCKCOFACTOR) already
%     gives t_g EXACTLY as r_g = sum_{i in g}(1 - h_i), cheaply, and is the
%     DEFAULT and the authority. This estimator is an OPT-IN accelerator for very
%     large nc where the exact per-group DIAGONAL diag(N^-1 N_g^cc) would cost
%     O(G*nc^2) to form. LSQVCE.ENGINE.BLOCKCOFACTOR documents that Hutchinson on the JOINT
%     leverage operator N^-1 N (trace u, energy in one mode) has 20-50% error and
%     was rejected there -- that caveat is real. Here the PER-GROUP operators
%     N^-1 N_g^cc have trace r_g (typically << u) and the group-sum constraint is
%     exact, which tempers but does not remove the noise: always check .relSE and
%     prefer the exact path whenever it is affordable.
%
%   INPUT
%     acc   1x1 lsqvce.engine.BlockAccumulator   built with 'KeepLocal', true (needs the
%                                         per-system Schur handles to re-stream)
%     sol   1x1 struct                    the solve output; must carry sol.x
%                                         (nc x 1, the common estimate). The
%                                         probed quadratic form A'A depends only
%                                         on the whitening (the variance
%                                         components), not on x, so sol is used
%                                         as a contract/consistency check that
%                                         the accumulator has been solved.
%
%   OPTIONS (Name/Value)
%   (WP-39 stage (b): canonical UpperCamelCase names below; legacy
%   casings still match -- MATLAB NV semantics -- and enum VALUES are
%   case-sensitive lowercase.)
%     'GroupBy'   char | function_handle  'arc' | 'obsType' | 'system' | 'day', or
%                                         a predicate @(sat,obsType,arc,day)->key
%                                         (default 'system'). Same spec as
%                                         LSQVCE.ENGINE.REDUNDANCY / LSQVCE.ENGINE.VCE.
%     'Groups'    1xG cell | []           pin the group order (keys as returned by
%                                         GroupBy). [] discovers them in sweep
%                                         order (default []).
%     'Probes'    1x1 int                 Z, number of Hutchinson probes (default
%                                         100, Ellmer's value).
%     'Seed'      1x1 int | []            RNG seed for reproducible probes
%                                         (default []: private stream, no global
%                                         side effect). Uses one FIXED probe block
%                                         across all groups so the group-sum
%                                         invariants hold on the SAME probes.
%     'WantDiag'  1x1 logical             also return the per-parameter diagonal
%                                         diag(N^-1 N_g^cc) per group (default
%                                         true). Set false for trace-only (§1.1).
%     'Verbose'   1x1 logical             progress banner (default false).
%
%   OUTPUT
%     out   1x1 struct with fields
%       .group    1xG struct   per group, in the pinned/discovered order:
%                   .key    the group key (as GroupBy returns it)
%                   .trace  1x1     estimate of tr(N^-1 N_g^cc)
%                   .diag   nc x 1  estimate of diag(N^-1 N_g^cc), or [] if
%                                   'WantDiag' is false
%                   .se     1x1     Monte-Carlo standard error of .trace (NaN
%                                   for one probe -- WP-134, it read 0)
%                   .relSE  1x1     .se / |.trace|  (NaN if .trace == 0)
%                   .n      1x1     observation rows in the group (exact count
%                                   from the same sweep)
%                   .localAbsorbed 1x1  EXACT local leverage of the group,
%                                   local_g = tr(Qll B_g'B_g), Qll = inv(Nll)
%                                   (its pseudo-inverse for a rank-deficient
%                                   block), the second term of the block
%                                   leverage h = z'Sxx z + b'Qll b
%                                   (LSQVCE.ENGINE.BLOCKCOFACTOR, WP-131); a
%                                   whole-system group gives exactly the local
%                                   rank. No probes enter this term.
%                   .r      1x1     STREAMED REDUNDANCY (WP-32, fixes F-01):
%                                   r_g = .n - .localAbsorbed - .trace
%                                   == n_g - sum_{i in g} h_i in expectation,
%                                   with sum_g .r = m - u without a common
%                                   prior, and sum_g .r + r_K = m + nk - u with
%                                   one (r_K its nk rows' redundancy,
%                                   LSQVCE.ENGINE.REDUNDANCY). THIS is the Helmert
%                                   denominator; .trace alone is only the
%                                   cc-attribution (sums to nc, NOT m-u).
%                   .rSE    1x1     Monte-Carlo se of .r (== .se: the local
%                                   term and the count are exact)
%       .nc       1x1        common-parameter count (size of the probe space)
%       .probes   1x1        Z actually used
%       .traceSum 1x1        sum_g .trace  (~ nc without a common prior; the
%                            sanity check)
%       .seedUsed 1x1        the probe seed (the given Seed, or the one the
%                            private stream drew -- passing it back replays the
%                            run; WP-134)
%       .note     char       provenance string
%
%   COST
%     One (batched) data sweep, plus Z triangular solves of size nc against the
%     shared Lc. Total O(passes*nnz + Z*nc^2). Memory is the probe block (nc x Z)
%     plus ONE nc x Z accumulator PER GROUP for the whole sweep: O(G*nc*Z), not
%     O(G*nc^2) -- 24 GB for 10^4 arcs, nc = 3000, Z = 100 (WP-134: it said one).
%
%   EXAMPLE
%     acc = lsqvce.engine.BlockAccumulator(nc, 'KeepLocal', true);
%     ... addSystem(...) ...
%     sol = acc.solve();
%     c   = lsqvce.engine.streamContribution(acc, sol, 'GroupBy','system', 'Probes',200);
%     share = [c.group.trace] ./ c.traceSum;      % relative info share per system
%     bar(c.group(1).diag);                       % per-coeff share of system 1
%
%   ACCURACY
%     .trace (with its Monte-Carlo .se) and .diag are Hutchinson estimates;
%     .n and .localAbsorbed are exact: with the exact t_g they give the exact
%     r_g to 1e-8 n_g, and the streamed trace sits within max(6 se, 5 %) of
%     that t_g at 2000 probes (test_contribution
%     (test_contribution_level2_is_level1_split_by_obsType)). Without a common
%     prior the group diagonals sum to ones(nc,1) to 1e-8 on one probe block
%     (test_contribution (test_stream_diag_sums_to_ones)). The exact path
%     stays the authority. WP-139.
%
%   ERRORS
%     lsqvce:streamContribution:badSol      the solution struct is not from this
%                                           accumulator
%     lsqvce:makeKeyFcn:badGroupBy  the grouping key is unknown (WP-38:
%                   the one shared key resolver)
%     lsqvce:streamContribution:noStore     nothing was retained; re-run with
%                                           'store','all'
%     lsqvce:streamContribution:noArcs      a system declares no arcs
%     lsqvce:streamContribution:cNotPD      a contribution block is not positive
%                                           definite
%
%   See also LSQVCE.ENGINE.CONTRIBUTION, LSQVCE.ENGINE.REDUNDANCY, LSQVCE.ENGINE.BLOCKCOFACTOR,
%            LSQVCE.ENGINE.TRACEESTIMATE, LSQVCE.ENGINE.VCE
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-139: the help states its
%   accuracy, every number traced to the test that pins it. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-134 (review 3 lows): se is NaN
%   for one probe (it read 0, the estimate as if exact); out.seedUsed; the
%   help states the sums with a common prior and the per-group memory. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-133: n_g counts observations,
%   not the rows of weight 0 (arc.excluded). With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-131: the leverage split
%   N_g = Z_g'Z_g (Level 1's) replaces the symmetric half-Schur split, so a
%   group's Level-2 diagonal sums to its Level-1 share and is non-negative;
%   localAbsorbed = tr(Qll Mll_g) on the solve's local cofactor. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-86 part 3 (R7-24): the
%   retired leverage / leverageBlocks references point at
%   lsqvce.engine.BlockCofactor. With assistance from Claude (Anthropic),
%   model Claude Opus 5.5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (b) (F-07): typed
%   `arguments` signature, canonical UpperCamelCase option names (legacy
%   casings keep matching; enum VALUES case-sensitive lowercase). With
%   assistance from Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-38 (fixes F-08): key presets
%   and the batch-and-whiten skeleton live in the ONE private makeKeyFcn /
%   arcBatchWhiten (byte-identical numerics; GroupBy error id now
%   lsqvce:makeKeyFcn:badGroupBy). With assistance from Claude
%   (Anthropic), model Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-36 (fixes F-06): probes from
%   the ONE private generator probeMatrix; empty seed = private
%   session-stream draw, global RNG untouched; seeded realisations
%   unchanged. With assistance from Claude (Anthropic), model Claude
%   Fable 5.
%   Updated 2026-08-18 (CEST, UTC+02:00), WP-32: every group additionally
%   carries .n, .localAbsorbed and the streamed redundancy .r/.rSE (see the
%   OUTPUT block for the derivation); arcCommonProbe streams B_g'B_g and the
%   row count in the same pass; per-system leverage blocks X = Sxx G' and
%   G Sxx G' are built once per system with one nc-solve. Numerical evidence:
%   validation/v10review.py T2d/T2e. With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   Created 2026-07-17 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.
arguments
    acc
    sol
    opt.GroupBy = 'system'
    opt.Groups = []
    opt.Probes (1,1) double {mustBeInteger, mustBePositive} = 100
    opt.Seed double {mustBeScalarOrEmpty} = []
    opt.WantDiag (1,1) logical = true
    opt.Verbose (1,1) logical = false
end

    % WP-39 stage (b): typed `arguments` signature (names already
    % UpperCamelCase; NAME matching case-insensitive per MATLAB semantics).
    o = opt;
    wantDiag = logical(o.WantDiag);

    ids = acc.systemIDs();
    if isempty(ids)
        error('lsqvce:streamContribution:noStore', ...
              ['The accumulator has no kept systems. Construct it with ' ...
               '''KeepLocal'', true so the per-system Schur handles survive.']);
    end

    if ~isstruct(sol) || ~isfield(sol, 'x')
        error('lsqvce:streamContribution:badSol', ...
              'sol must be the solve struct carrying the common estimate sol.x.');
    end
    Ncc_red = acc.reducedNormalMatrix();            % TOTAL N, summed over systems
    nc      = size(Ncc_red, 1);
    keyFcn  = makeKeyFcn(o.GroupBy);

    % ---- shared Cholesky of N (the ONLY nc^3 factorisation) --------------
    [Lc, pc] = chol((Ncc_red + Ncc_red.')/2, 'lower');
    if pc > 0
        error('lsqvce:streamContribution:cNotPD', ...
            ['The reduced common normal matrix is not positive definite (datum ' ...
             'defect). Add constraints or a common-part prior before probing ' ...
             'the contribution.']);
    end

    % ---- the fixed probe block (shared by ALL groups) --------------------
    %   Rademacher +-1 from the ONE generator (WP-36, fixes F-06): a given
    %   seed reproduces bit-identically across every caller; an empty seed
    %   draws once from the package-private session stream -- the global RNG
    %   is untouched (the promise this comment made before the fix while the
    %   code advanced the global stream).
    [Z, seedUsed] = probeMatrix(nc, o.Probes, 'rademacher', o.Seed);   % WP-134: the seed returned

    % ---- per-system broadcast state (serial, cheap) ----------------------
    %   G_s = inv(Nll_s) Nlc_s (nl x nc), the per-system Schur factor, identical
    %   to the G of LSQVCE.ENGINE.BLOCKCOFACTOR. Precompute G_s*Z (nl x Z) once per system.
    nSys   = numel(ids);
    stores = cell(nSys, 1);
    for k = 1:nSys
        st = acc.getStore(ids(k));
        [~, G] = lsqvce.engine.BlockAccumulator.localSolveOf(st, []);   % nl x nc (WP-128)
        % WP-131: the exact local leverage of a group is tr(Qll Mll_g), Qll the
        % solve's local cofactor (on its factor, or its pseudo-inverse: the
        % accumulator's localCovariance with Cxx = 0, s0^2 = 1)
        Qll = acc.localCovariance(st.sysID, zeros(nc), 1);
        stores{k} = struct('sysID', st.sysID, 'dataFcn', st.dataFcn, ...
                           'uo', st.uo, 'dayList', st.dayList, ...
                           'G', G, 'GZ', G * Z, 'Qll', Qll);
    end

    % ---- group bookkeeping (pinned order if 'Groups' given) --------------
    if isempty(o.Groups)
        keys    = {};
        keyVals = {};
    else
        keyVals = o.Groups(:).';
        keys    = cellfun(@(k) char(string(k)), keyVals, 'UniformOutput', false);
    end
    nG0   = numel(keys);
    AccAZ = repmat({zeros(nc, o.Probes)}, 1, nG0);    % N_g Z (the leverage split), nc x Z
    nObs  = zeros(1, nG0);                            % n_g, rows per group
    locA  = zeros(1, nG0);                            % exact local absorbed

    % ---- ONE data sweep: the LEVERAGE split of N (WP-131) ----------------
    %   N_g = Z_g' Z_g,  Z_g = A_g - B_g G_s  (the rows reduced by the local
    %   elimination) -- the split Level 1 uses (BlockCofactor.leverage's hc,
    %   t_g = sum_{i in g} z_i' Sxx z_i). It sums to N over the groups of a
    %   system (G_s' Nll_s G_s = G_s' Nlc_s), is positive semidefinite per
    %   group, and its matvec is assembled from the per-arc pieces:
    %       N_g Z = A_g'A_g Z - Nlc_g' (G Z) - G' (Nlc_g Z) + G' Mll_g (G Z).
    %   The symmetric half-Schur split used before also summed to N but not
    %   group by group to Level 1's t_g (off by tr(Sxx G' B_g'(A_g - B_g G)),
    %   zero only for a whole system): a group's Level-2 diagonal did not sum
    %   to its share, and could be negative
    %   (dev/validation/wp131_level2_split_mirror.py).
    for k = 1:nSys
        st = stores{k};
        uo = st.uo;
        % Accumulate this system's per-group A'(AZ) and cross block Nlc_g^(s).
        gAZ  = {};                                    % local: this system's groups
        gNlc = {};                                    % Nlc_g^(s), nl x nc each
        gMll = {};                                    % B_g'B_g,   nl x nl each
        gN   = [];                                    % rows per local group
        gIdx = [];                                    % global group index per local
        for d = st.dayList(:).'
            meta = st.dataFcn('meta', d, [], uo);
            for ot = 1:numel(meta.obsTypes)
                obsType = meta.obsTypes{ot};
                for a = meta.arcs{ot}(:).'
                    key  = keyFcn(st.sysID, obsType, a, d);
                    ckey = char(string(key));
                    [tf, gi] = ismember(ckey, keys);
                    if ~tf
                        keys{end+1}    = ckey;        %#ok<AGROW>
                        keyVals{end+1} = key;         %#ok<AGROW>
                        gi = numel(keys);
                        AccAZ{gi} = zeros(nc, o.Probes); %#ok<AGROW>
                        nObs(gi)  = 0;                   %#ok<AGROW>
                        locA(gi)  = 0;                   %#ok<AGROW>
                    end
                    lg = find(gIdx == gi, 1);
                    if isempty(lg)
                        gIdx(end+1) = gi;            %#ok<AGROW>
                        lg          = numel(gIdx);
                        gAZ{lg}     = zeros(nc, o.Probes);   %#ok<AGROW>
                        gNlc{lg}    = zeros(size(st.G));      %#ok<AGROW> % nl x nc
                        gMll{lg}    = zeros(size(st.G,1));    %#ok<AGROW> % nl x nl
                        gN(lg)      = 0;                      %#ok<AGROW>
                    end
                    [AZ, NlcArc, MllArc, nArc] = arcCommonProbe(st, obsType, a, d, uo, Z);
                    if ~isempty(AZ)
                        gAZ{lg}  = gAZ{lg}  + AZ;
                        gNlc{lg} = gNlc{lg} + NlcArc;
                        gMll{lg} = gMll{lg} + MllArc;
                        gN(lg)   = gN(lg)   + nArc;
                    end
                end
            end
        end
        % Fold this system's per-group pieces into the global accumulators:
        % the leverage-split matvec and the exact local leverage tr(Qll Mll_g)
        for lg = 1:numel(gIdx)
            gi   = gIdx(lg);
            NlcG = gNlc{lg};                          % Nlc_g^(s), nl x nc
            MllG = gMll{lg};                          % B_g'B_g,   nl x nl
            AccAZ{gi} = AccAZ{gi} + gAZ{lg} - NlcG.' * st.GZ - st.G.' * (NlcG * Z) ...
                                  + st.G.' * (MllG * st.GZ);
            locA(gi)  = locA(gi) + sum(sum(st.Qll .* MllG));   % tr(Qll Mll_g), both symmetric
            nObs(gi) = nObs(gi) + gN(lg);
        end
        if o.Verbose
            lsqvce.util.log.msg('%sstreamed system %s\n', lsqvce.util.verbosity.banner('Contribution'), ...
                string(st.sysID));
        end
    end

    nG = numel(keys);
    if nG == 0
        error('lsqvce:streamContribution:noArcs', ...
              'No arcs were enumerated -- nothing to probe.');
    end

    % ---- apply N^-1 and reduce to trace / diagonal per group -------------
    % WP-32: alongside the cc-attribution trace, every group now carries the
    % STREAMED REDUNDANCY  r_g = n_g - local_g - tr(N^-1 N_g^cc). Only the
    % trace term is stochastic, so r inherits the trace's Monte-Carlo se.
    grp = struct('key', {}, 'trace', {}, 'diag', {}, 'se', {}, 'relSE', {}, ...
                 'n', {}, 'localAbsorbed', {}, 'r', {}, 'rSE', {});
    traceSum = 0;
    for gi = 1:nG
        NgZ  = AccAZ{gi};                             % N_g Z (leverage split), nc x Z
        W    = Lc.' \ (Lc \ NgZ);                     % N^-1 N_g^cc Z, nc x Z
        % Hutchinson trace: mean_z z'(N^-1 N_g Z);  per-probe values for the se.
        perProbe = sum(Z .* W, 1);                    % 1 x Z
        tG       = mean(perProbe);
        se       = NaN;                               % WP-134: one draw has no spread
        if o.Probes > 1, se = std(perProbe) / sqrt(o.Probes); end
        if tG ~= 0, rse = se / abs(tG); else, rse = NaN; end
        if wantDiag
            dG = mean(Z .* W, 2);                     % nc x 1
        else
            dG = [];
        end
        grp(gi) = struct('key', keyVals{gi}, 'trace', tG, 'diag', dG, ...
                         'se', se, 'relSE', rse, ...
                         'n', nObs(gi), 'localAbsorbed', locA(gi), ...
                         'r', nObs(gi) - locA(gi) - tG, 'rSE', se); %#ok<AGROW>
        traceSum = traceSum + tG;
    end

    out = struct('group', grp, 'nc', nc, 'probes', o.Probes, ...
                 'traceSum', traceSum, 'seedUsed', seedUsed, ...
                 'note', sprintf(['streaming Hutchinson N^-1 N_g^cc over %d ' ...
                          'probe(s); sum_g trace = %.6g (nc = %d less a common ' ...
                          'prior''s share)'], o.Probes, traceSum, nc));

    if o.Verbose
        lsqvce.util.log.msg(['%s%d group(s), sum trace = %.4f ' ...
            '(nc = %d), max relSE = %.2e\n'], lsqvce.util.verbosity.banner('Contribution'), nG, traceSum, nc, ...
            max([grp.relSE]));
    end
end

%==========================================================================
function [AZ, Nlc, Mll, nRows] = arcCommonProbe(st, obsType, arc, day, uo, Z)
%ARCCOMMONPROBE  Stream one arc: A_w'(A_w Z), cross block B_w'A_w, local Gram.
%   AZ    (nc x Z)  the arc's contribution to sum_arc A'(A Z)
%   Nlc   (nl x nc) the arc's whitened cross block B_w' A_w, summed over batches
%   Mll   (nl x nl) the arc's whitened local Gram B_w' B_w (WP-32: feeds the
%                   exact local-absorbed term of the streamed redundancy)
%   nRows (1 x 1)   observation rows this arc contributes (n_g bookkeeping)
%   Whitening is applied EXACTLY as the exact leverage path (LSQVCE.ENGINE.REDUNDANCY),
%   so grouping/scaling stay consistent. Returns [] for a data gap.
    uo.obsType = obsType;
    uo.day     = day;
    uo.arc     = arc;
    a = st.dataFcn(day, obsType, arc, uo);
    if isempty(a) || ~isfield(a, 'l') || isempty(a.l)
        AZ = [];  Nlc = [];  Mll = [];  nRows = 0;  return
    end

    nc  = size(Z, 1);
    nl  = size(a.B, 2);
    m   = size(a.A, 1);
    % WP-38 (fixes F-08): the batch/overlap/whitening skeleton lives in the
    % ONE private arcBatchWhiten; this callback keeps the payload and the
    % accumulation order (byte-identical outputs). NOTE the shared Cfcn
    % branch also whitens lb (the leverage path's behaviour) -- this
    % callback never touches lb, so the outputs are unchanged.
    AZ    = zeros(nc, size(Z, 2));
    Nlc   = zeros(nl, nc);
    Mll   = zeros(nl, nl);
    nRows = m - nnz(lsqvce.util.getfielddef(a, 'excluded', []));   % WP-133: weight-0 rows are not observations
    arcBatchWhiten(a, uo, @perBatch);

    function perBatch(~, keep, Ab, Bb, ~)
        Aw  = Ab(keep, :);                   % rows that belong to this batch
        Bw  = Bb(keep, :);
        AZ  = AZ  + Aw.' * (Aw * Z);         % nc x Z, streamed
        Nlc = Nlc + Bw.' * Aw;               % nl x nc, this arc's cross block
        Mll = Mll + Bw.' * Bw;               % nl x nl, this arc's local Gram
    end
end

%==========================================================================
