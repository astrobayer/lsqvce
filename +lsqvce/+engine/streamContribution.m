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
%   factored once as chol) and N_g^cc is group g's contribution to it. The
%   group-summed invariants hold EXACTLY (up to probe noise): sum_g t_g = nc and
%   sum_g diag(C_g) = ones(nc,1), because sum_g N_g^cc = N.
%
%   MEMORY-HONEST DESIGN (why this exists)
%     N_g^cc is nc x nc and a group can span many arcs, so holding all G of them
%     is O(G*nc^2) and destroys the streaming property that makes adjustBlock
%     worth using. We never store a group normal. What Hutchinson needs is only
%     the ACTION N_g^cc*Z on an nc x Z probe block, and that action is exactly
%     group g's contribution to a mat-vec, re-accumulated by STREAMING that
%     group's arcs -- mirroring the exact leverage path in LSQVCE.ENGINE.REDUNDANCY:
%
%       N_g^cc Z = sum_{arc in g} A_arc' (A_arc Z)   -   Nlc_s' ( G_s Z )
%
%     with A_arc the whitened common design of the arc (m x nc), G_s = inv(Nll_s)
%     Nlc_s the per-SYSTEM Schur factor (built once, shared by all that system's
%     arcs, identical to the G of LSQVCE.ENGINE.BLOCKCOFACTOR). The first term streams per arc;
%     the second is a single per-system correction applied once. All groups are
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
%                   .se     1x1     Monte-Carlo standard error of .trace
%                   .relSE  1x1     .se / |.trace|  (NaN if .trace == 0)
%                   .n      1x1     observation rows in the group (exact count
%                                   from the same sweep)
%                   .localAbsorbed 1x1  EXACT local dof absorbed by the group:
%                                   local_g = tr((Nll^-1 + G Sxx G') B_g'B_g)
%                                           - tr(B_g'A_g Sxx G')
%                                   from the leverage-block identities
%                                   inv(N)_ll = Nll^-1 + G Sxx G' and
%                                   inv(N)_lc = -G Sxx (LSQVCE.ENGINE.BLOCKCOFACTOR); a
%                                   whole-system group gives exactly nl_s. No
%                                   probes enter this term.
%                   .r      1x1     STREAMED REDUNDANCY (WP-32, fixes F-01):
%                                   r_g = .n - .localAbsorbed - .trace
%                                   == n_g - sum_{i in g} h_i in expectation,
%                                   with sum_g .r = m - u. THIS is the Helmert
%                                   denominator; .trace alone is only the
%                                   cc-attribution (sums to nc, NOT m-u).
%                   .rSE    1x1     Monte-Carlo se of .r (== .se: the local
%                                   term and the count are exact)
%       .nc       1x1        common-parameter count (size of the probe space)
%       .probes   1x1        Z actually used
%       .traceSum 1x1        sum_g .trace  (should be ~ nc; the sanity check)
%       .note     char       provenance string
%
%   COST
%     One (batched) data sweep, plus Z triangular solves of size nc against the
%     shared Lc. Total O(passes*nnz + Z*nc^2). Memory is the probe block (nc x Z)
%     plus, transiently, one nc x Z group accumulator -- NOT O(G*nc^2).
%
%   EXAMPLE
%     acc = lsqvce.engine.BlockAccumulator(nc, 'KeepLocal', true);
%     ... addSystem(...) ...
%     sol = acc.solve();
%     c   = lsqvce.engine.streamContribution(acc, sol, 'GroupBy','system', 'Probes',200);
%     share = [c.group.trace] ./ c.traceSum;      % relative info share per system
%     bar(c.group(1).diag);                       % per-coeff share of system 1
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
    Z = probeMatrix(nc, o.Probes, 'rademacher', o.Seed);

    % ---- per-system broadcast state (serial, cheap) ----------------------
    %   G_s = inv(Nll_s) Nlc_s (nl x nc), the per-system Schur factor, identical
    %   to the G of LSQVCE.ENGINE.BLOCKCOFACTOR. Precompute G_s*Z (nl x Z) once per system.
    nSys   = numel(ids);
    stores = cell(nSys, 1);
    for k = 1:nSys
        st = acc.getStore(ids(k));
        [~, G] = lsqvce.engine.BlockAccumulator.localSolveOf(st, []);   % nl x nc (WP-128)
        % WP-32: the exact local-absorbed term needs the two leverage-block
        % identities inv(N)_ll = Nll^-1 + G Sxx G' and inv(N)_lc = -G Sxx
        % (LSQVCE.ENGINE.BLOCKCOFACTOR). One nc-solve per system yields X = Sxx G';
        % everything else is nl x nl.
        X   = Lc.' \ (Lc \ G.');                    % Sxx G', nc x nl
        GSG = G * X;                                 % G Sxx G', nl x nl
        stores{k} = struct('sysID', st.sysID, 'dataFcn', st.dataFcn, ...
                           'uo', st.uo, 'dayList', st.dayList, ...
                           'Nlc', st.Nlc, 'Nllr', st.Nll_r, 'GZ', G * Z, ...
                           'X', X, 'GSG', GSG);
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
    AccAZ = repmat({zeros(nc, o.Probes)}, 1, nG0);    % sum_arc A'(A Z), nc x Z
    Corr  = repmat({zeros(nc, o.Probes)}, 1, nG0);    % -1/2 Schur term, nc x Z
    nObs  = zeros(1, nG0);                            % n_g, rows per group
    locA  = zeros(1, nG0);                            % exact local absorbed

    % ---- ONE data sweep: symmetric Schur decomposition -------------------
    %   N_g^cc = Ncc_g - 1/2 ( Nlc_g' inv(Nll_s) Nlc_s + Nlc_s' inv(Nll_s) Nlc_g )
    %   which partitions EXACTLY (sum_g N_g^cc = N) and is INDEPENDENT of arc
    %   ordering -- the only correct split when a system's arcs span several
    %   groups (e.g. obsType grouping of a merged pair-system). Its matvec is
    %       N_g^cc Z = sum_{arc in g} A'(A Z)
    %                  - 1/2 ( Nlc_g' (G_s Z) + Nlc_s' (G_g Z) ),   G_g=inv(Nll)Nlc_g
    %   accumulated per (group, system): within a system's pass we build that
    %   group's cross block Nlc_g^(s) = sum_{arc in g, sys s} B'A and apply the
    %   two half-terms, so systems never cross-contaminate.
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
                        Corr{gi}  = zeros(nc, o.Probes); %#ok<AGROW>
                        nObs(gi)  = 0;                   %#ok<AGROW>
                        locA(gi)  = 0;                   %#ok<AGROW>
                    end
                    lg = find(gIdx == gi, 1);
                    if isempty(lg)
                        gIdx(end+1) = gi;            %#ok<AGROW>
                        lg          = numel(gIdx);
                        gAZ{lg}     = zeros(nc, o.Probes);   %#ok<AGROW>
                        gNlc{lg}    = zeros(size(st.Nlc));    %#ok<AGROW> % nl x nc
                        gMll{lg}    = zeros(size(st.Nlc,1));  %#ok<AGROW> % nl x nl
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
        % Fold this system's per-group pieces into the global accumulators,
        % applying the symmetric half-Schur correction now (G_s Z is st.GZ).
        for lg = 1:numel(gIdx)
            gi   = gIdx(lg);
            NlcG = gNlc{lg};                          % Nlc_g^(s), nl x nc
            GgZ  = st.Nllr \ (NlcG * Z);              % inv(Nll_s) Nlc_g Z, nl x Z
            half = 0.5 * ( NlcG.' * st.GZ + st.Nlc.' * GgZ );   % nc x Z
            AccAZ{gi} = AccAZ{gi} + gAZ{lg};
            Corr{gi}  = Corr{gi}  + half;
            % WP-32: the EXACT local dof absorbed by this group (no probes):
            %   local_g = tr((Nll^-1 + G Sxx G') Mll_g) - tr(Nlc_g Sxx G')
            % from the leverage-block identities. For a whole-system group it
            % collapses to nl_s exactly; per group it completes the streamed
            % cc-attribution to the true absorbed dof sum_{i in g} h_i.
            MllG     = gMll{lg};                      % B_g'B_g, nl x nl
            locA(gi) = locA(gi) + trace(st.Nllr \ MllG) ...
                                + trace(st.GSG * MllG) ...
                                - trace(NlcG * st.X);
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
        NgZ  = AccAZ{gi} - Corr{gi};                  % N_g^cc Z, nc x Z
        W    = Lc.' \ (Lc \ NgZ);                     % N^-1 N_g^cc Z, nc x Z
        % Hutchinson trace: mean_z z'(N^-1 N_g Z);  per-probe values for the se.
        perProbe = sum(Z .* W, 1);                    % 1 x Z
        tG       = mean(perProbe);
        se       = std(perProbe) / sqrt(o.Probes);
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
                 'traceSum', traceSum, ...
                 'note', sprintf(['streaming Hutchinson N^-1 N_g^cc over %d ' ...
                          'probe(s); sum_g trace = %.6g (target nc = %d)'], ...
                          o.Probes, traceSum, nc));

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
    nRows = m;
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
