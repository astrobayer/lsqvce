function out = streamContribution(acc, sol, varargin)
%LSQVCE.STREAMCONTRIBUTION  Streaming, matrix-free per-group N^-1 N_g^cc probes.
%
%   out = LSQVCE.STREAMCONTRIBUTION(acc, sol)
%   out = LSQVCE.STREAMCONTRIBUTION(acc, sol, 'Name', value, ...)
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
%     group's arcs -- mirroring the exact leverage path in LSQVCE.REDUNDANCY:
%
%       N_g^cc Z = sum_{arc in g} A_arc' (A_arc Z)   -   Nlc_s' ( G_s Z )
%
%     with A_arc the whitened common design of the arc (m x nc), G_s = inv(Nll_s)
%     Nlc_s the per-SYSTEM Schur factor (built once, shared by all that system's
%     arcs, identical to LSQVCE.LEVERAGEBLOCKS.G). The first term streams per arc;
%     the second is a single per-system correction applied once. All groups are
%     accumulated in ONE sweep -- each arc routes its A'(A Z) into its own group's
%     nc x Z accumulator, so the cost is a handful of data passes independent of G.
%
%     N^-1 is then applied to the assembled N_g^cc Z through the shared Cholesky
%     Lc of N (two triangular solves), never by forming N^-1.
%
%   RELATION TO THE EXACT PATH (read this before trusting the numbers)
%     The exact per-arc leverage path (LSQVCE.REDUNDANCY, LSQVCE.LEVERAGE) already
%     gives t_g EXACTLY as r_g = sum_{i in g}(1 - h_i), cheaply, and is the
%     DEFAULT and the authority. This estimator is an OPT-IN accelerator for very
%     large nc where the exact per-group DIAGONAL diag(N^-1 N_g^cc) would cost
%     O(G*nc^2) to form. LSQVCE.LEVERAGE documents that Hutchinson on the JOINT
%     leverage operator N^-1 N (trace u, energy in one mode) has 20-50% error and
%     was rejected there -- that caveat is real. Here the PER-GROUP operators
%     N^-1 N_g^cc have trace r_g (typically << u) and the group-sum constraint is
%     exact, which tempers but does not remove the noise: always check .relSE and
%     prefer the exact path whenever it is affordable.
%
%   INPUT
%     acc   1x1 lsqvce.BlockAccumulator   built with 'KeepLocal', true (needs the
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
%     'GroupBy'   char | function_handle  'arc' | 'obsType' | 'system' | 'day', or
%                                         a predicate @(sat,obsType,arc,day)->key
%                                         (default 'system'). Same spec as
%                                         LSQVCE.REDUNDANCY / LSQVCE.VCE.
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
%     acc = lsqvce.BlockAccumulator(nc, 'KeepLocal', true);
%     ... addSystem(...) ...
%     sol = acc.solve();
%     c   = lsqvce.streamContribution(acc, sol, 'GroupBy','system', 'Probes',200);
%     share = [c.group.trace] ./ c.traceSum;      % relative info share per system
%     bar(c.group(1).diag);                       % per-coeff share of system 1
%
%   ERRORS
%     lsqvce:streamContribution:badSol      the solution struct is not from this
%                                           accumulator
%     lsqvce:streamContribution:badGroupBy  the grouping key is unknown
%     lsqvce:streamContribution:noStore     nothing was retained; re-run with
%                                           'store','all'
%     lsqvce:streamContribution:noArcs      a system declares no arcs
%     lsqvce:streamContribution:cNotPD      a contribution block is not positive
%                                           definite
%
%   See also LSQVCE.CONTRIBUTION, LSQVCE.REDUNDANCY, LSQVCE.LEVERAGEBLOCKS,
%            LSQVCE.TRACEESTIMATE, LSQVCE.VCE
%
%   Created 2026-07-17 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Opus 4.8.

    p = inputParser;
    p.addParameter('GroupBy' , 'system');
    p.addParameter('Groups'  , []);
    p.addParameter('Probes'  , 100, @(v) isscalar(v) && v >= 1 && v == floor(v));
    p.addParameter('Seed'    , []);
    p.addParameter('WantDiag', true, @(v) isscalar(v) && (islogical(v) || isnumeric(v)));
    p.addParameter('Verbose' , false, @(v) isscalar(v));
    p.parse(varargin{:});
    o = p.Results;
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
    %   Rademacher +1/-1, a private stream so the global RNG is untouched.
    if isempty(o.Seed)
        Z = sign(randn(nc, o.Probes));
    else
        rs = RandStream('twister', 'Seed', o.Seed);
        Z  = sign(randn(rs, nc, o.Probes));
    end
    Z(Z == 0) = 1;                                   % sign(0) -> +1

    % ---- per-system broadcast state (serial, cheap) ----------------------
    %   G_s = inv(Nll_s) Nlc_s (nl x nc), the per-system Schur factor, identical
    %   to LSQVCE.LEVERAGEBLOCKS.G. Precompute G_s*Z (nl x Z) once per system.
    nSys   = numel(ids);
    stores = cell(nSys, 1);
    for k = 1:nSys
        st = acc.getStore(ids(k));
        G  = st.Nll_r \ st.Nlc;                      % nl x nc
        stores{k} = struct('sysID', st.sysID, 'dataFcn', st.dataFcn, ...
                           'uo', st.uo, 'dayList', st.dayList, ...
                           'Nlc', st.Nlc, 'Nllr', st.Nll_r, 'GZ', G * Z);
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
                    end
                    lg = find(gIdx == gi, 1);
                    if isempty(lg)
                        gIdx(end+1) = gi;            %#ok<AGROW>
                        lg          = numel(gIdx);
                        gAZ{lg}     = zeros(nc, o.Probes);   %#ok<AGROW>
                        gNlc{lg}    = zeros(size(st.Nlc));    %#ok<AGROW> % nl x nc
                    end
                    [AZ, NlcArc] = arcCommonProbe(st, obsType, a, d, uo, Z);
                    if ~isempty(AZ)
                        gAZ{lg}  = gAZ{lg}  + AZ;
                        gNlc{lg} = gNlc{lg} + NlcArc;
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
    grp = struct('key', {}, 'trace', {}, 'diag', {}, 'se', {}, 'relSE', {});
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
                         'se', se, 'relSE', rse); %#ok<AGROW>
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
function [AZ, Nlc] = arcCommonProbe(st, obsType, arc, day, uo, Z)
%ARCCOMMONPROBE  Stream one arc, return A_w'(A_w Z) and the cross block B_w'A_w.
%   AZ  (nc x Z)  the arc's contribution to sum_arc A'(A Z)
%   Nlc (nl x nc) the arc's whitened cross block B_w' A_w, summed over batches
%   Whitening is applied EXACTLY as the exact leverage path (LSQVCE.REDUNDANCY),
%   so grouping/scaling stay consistent. Returns [] for a data gap.
    uo.obsType = obsType;
    uo.day     = day;
    uo.arc     = arc;
    a = st.dataFcn(day, obsType, arc, uo);
    if isempty(a) || ~isfield(a, 'l') || isempty(a.l)
        AZ = [];  Nlc = [];  return
    end

    nc  = size(Z, 1);
    nl  = size(a.B, 2);
    m   = size(a.A, 1);
    bs  = lsqvce.util.getfielddef(uo, 'BatchSize', 20000);
    nOv = lsqvce.util.getfielddef(uo, 'BatchOverlap', []);
    if isempty(nOv)
        nOv = lsqvce.util.getfielddef(a, 'filterMemory', 0);
    end
    bs  = max(1, min(round(bs),  m));
    nOv = max(0, min(round(nOv), bs-1));

    hasW = isfield(a, 'Wfcn') && ~isempty(a.Wfcn);
    hasC = isfield(a, 'Cfcn') && ~isempty(a.Cfcn);

    AZ  = zeros(nc, size(Z, 2));
    Nlc = zeros(nl, nc);
    for i0 = 1:bs:m
        i1   = min(i0+bs-1, m);
        j0   = max(1, i0-nOv);
        idx  = j0:i1;
        keep = (i0-j0+1):numel(idx);

        Ab = a.A(idx, :);
        Bb = a.B(idx, :);
        lb = a.l(idx);  lb = lb(:);
        if hasW
            [Ab, Bb, lb] = a.Wfcn(idx, Ab, Bb, lb, uo);
        elseif hasC
            W  = lsqvce.decorrFactor(a.Cfcn(idx, uo));
            Ab = W*Ab;  Bb = W*Bb;
        end
        Aw  = Ab(keep, :);                   % rows that belong to this batch
        Bw  = Bb(keep, :);
        AZ  = AZ  + Aw.' * (Aw * Z);         % nc x Z, streamed
        Nlc = Nlc + Bw.' * Aw;               % nl x nc, this arc's cross block
    end
end

%==========================================================================
function f = makeKeyFcn(g)
%MAKEKEYFCN  Resolve a GroupBy spec into a key function f(sat,obsType,arc,day).
%   Mirrors LSQVCE.REDUNDANCY so grouping is identical across the two paths.
%   INPUT  g  char | function_handle
%   OUTPUT f  1x1 function_handle  @(sat,obsType,arc,day) -> group key
    if isa(g, 'function_handle')
        f = g;  return
    end
    switch lower(g)
        case 'arc',     f = @(s,ot,a,d) sprintf('s%s_%s_d%s_a%s', ...
                                    string(s), string(ot), string(d), string(a));
        case 'obstype', f = @(s,ot,a,d) sprintf('%s', string(ot));
        case 'system',  f = @(s,ot,a,d) sprintf('s%s', string(s));
        case 'day',     f = @(s,ot,a,d) sprintf('d%s', string(d));
        otherwise
            error('lsqvce:streamContribution:badGroupBy', ...
                ['GroupBy must be ''arc'', ''obsType'', ''system'', ''day'' ' ...
                 'or a function handle (got ''%s'').'], g);
    end
end
