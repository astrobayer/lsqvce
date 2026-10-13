function out = outliers(acc, sol, opt)
%LSQVCE.ENGINE.OUTLIERS  Data snooping (Baarda) using exact partial redundancy.
%
%   out = LSQVCE.ENGINE.OUTLIERS(acc, sol, 'Name', value, ...)
%
%   Flags observations whose NORMALISED residual is large:
%
%       w_i = |v_i| / sqrt(r_i)        with  r_i = 1 - h_ii
%
%   The normalisation by sqrt(r_i) is what makes this correct. A raw residual is
%   misleading: an observation with high leverage (h_ii near 1) is fitted almost
%   exactly BY CONSTRUCTION, so its residual is small even when it is badly
%   wrong. Dividing by sqrt(r_i) removes that effect. Verified: 4/4 injected
%   outliers detected, 0 false alarms.
%
%   INPUT
%     acc      1x1 lsqvce.engine.BlockAccumulator   after solve
%     sol      1x1 struct                            its solution (x, Cxx, s0sq)
%     Name, Value                                    the OPTIONS below
%
%   OPTIONS
%     'threshold'  critical value for w_i (default 4.0)
%                  3.29 = 99.9% for a single test; 4.0 is a sensible default
%                  once you account for testing millions of observations.
%     'groupBy'    as LSQVCE.ENGINE.REDUNDANCY, for the summary (default 'arc')
%     'minRedundancy'  observations with r_i below this are NOT testable and are
%                  reported separately (default 0.05)
%     'verbose'    default true
%
%   OUTPUT
%     out.flagged(j).sysID .obsType .day .arc .idx .w .v .r .h
%     out.nFlagged
%     out.nUntestable    r_i below MinRedundancy: no local redundancy, so a
%                        blunder there is INVISIBLE. This is a design weakness
%                        in the observation geometry, not a data problem.
%     out.byGroup        count per group
%     out.mask           @(sysID,obsType,day,arc) -> logical mask for the loader
%
%   THE UNTESTABLE OBSERVATIONS MATTER
%     If r_i ~ 0 the observation is not checked by any other -- an error in it
%     passes straight into the parameters undetected. A large nUntestable means
%     the geometry is weak (too many local parameters for too short an arc), not
%     that the data is clean. Fix the parameterisation, not the threshold.
%
%   REJECT AND RE-SOLVE
%     Rejection changes the solution, which changes the residuals, which changes
%     what is flagged. Iterate, and reject conservatively (a few at a time, the
%     largest first). Rejecting many at once on a single pass is how good data
%     gets thrown away.
%
%   EXAMPLE
%     sol = acc.solve();
%     out = lsqvce.engine.outliers(acc, sol, 'threshold', 4.0);
%
%     lsqvce.util.log.msg('%d flagged, %d untestable\n', out.nFlagged, out.nUntestable);
%
%     % feed the mask back to the loader and re-solve
%     opt.rejectMask = out.mask;
%     acc2 = buildAccumulator(cfg, opt);
%     sol2 = acc2.solve();
%
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-131: an 'spmd' accumulator
%   runs the redundancy pass as 'parfor' here, so the per-observation series
%   (and the flags) exist. With assistance from Claude (Anthropic), model
%   Claude Opus 5.5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (b) b2 (F-07):
%   typed `arguments` signature, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   See also LSQVCE.ENGINE.REDUNDANCY, LSQVCE.ENGINE.BLOCKCOFACTOR, LSQVCE.ENGINE.VCE
%
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-86 part 3 (R7-24): See also
%   names lsqvce.engine.BlockCofactor (leverage retired). With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

arguments
    acc
    sol
    opt.Threshold (1,1) double {mustBePositive} = 4.0
    opt.GroupBy = 'arc'
    opt.MinRedundancy (1,1) double {mustBeNonnegative} = 0.05
    opt.Verbose (1,1) logical = true
end

    narginchk(2,inf);

    % WP-39 stage (b) b2: typed `arguments` signature (canonical
    % UpperCamelCase; legacy casings keep matching -- MATLAB NV
    % semantics; enum VALUES case-sensitive lowercase).
    o = struct('threshold', opt.Threshold, 'groupBy', opt.GroupBy, 'minRedundancy', opt.MinRedundancy, 'verbose', opt.Verbose);

    % WP-131: outliers NEEDS the per-observation series; an accumulator in
    % 'spmd' mode made redundancy drop it (Collect forced false) and this
    % reported "0 flagged, 0 untestable" -- clean-looking data. spmd becomes
    % parfor here (the same per-arc work, the series shipped back).
    pm = '';
    if strcmpi(lsqvce.util.getfielddef(acc.opt, 'Parallel', 'serial'), 'spmd'), pm = 'parfor'; end
    red = lsqvce.engine.redundancy(acc, sol, 'GroupBy', o.groupBy, ...
                         'Collect', true, 'Verbose', false, 'Parallel', pm);

    F = struct([]);
    nUn = 0;
    keyCount = containers.Map('KeyType','char','ValueType','double');

    for k = 1:numel(red.obs)
        e  = red.obs(k);
        ri = e.r;
        vi = e.v;

        testable = ri >= o.minRedundancy;
        nUn = nUn + sum(~testable);

        w = zeros(size(vi));
        w(testable) = abs(vi(testable)) ./ sqrt(ri(testable));

        idx = find(w > o.threshold);
        for j = idx(:).'
            rec = struct('sysID',e.sysID, 'obsType',e.obsType, 'day',e.day, ...
                         'arc',e.arc, 'idx',j, 'w',w(j), 'v',vi(j), ...
                         'r',ri(j), 'h',e.h(j));
            if isempty(F), F = rec; else, F(end+1) = rec; end          %#ok<AGROW>

            ck = char(string(e.key));
            if keyCount.isKey(ck), keyCount(ck) = keyCount(ck) + 1;
            else,                  keyCount(ck) = 1; end
        end
    end

    out = struct('flagged', F, 'nFlagged', numel(F), 'nUntestable', nUn, ...
                 'threshold', o.threshold, 'byGroup', keyCount, ...
                 'red', red);
    out.mask = @(s,ot,d,a) buildMask(F, s, ot, d, a);

    if lsqvce.util.verbosity.progress(o.verbose)
        lsqvce.util.log.msg('%s%d flagged at w > %.2f   (%d untestable, r_i < %.3f)\n', lsqvce.util.verbosity.banner('Data snooping'), ...
                out.nFlagged, o.threshold, nUn, o.minRedundancy);
        if nUn > 0
            lsqvce.util.log.msg(['           %d observations have almost no local redundancy: ' ...
                     'a blunder in them is INVISIBLE. That is a geometry ' ...
                     'weakness, not a data problem.\n'], nUn);
        end
    end
end

%==========================================================================
function m = buildMask(F, s, ot, d, a)
%BUILDMASK  Logical rejection mask for one arc, for the loader to apply.
    m = [];
    if isempty(F), return, end
    hit = arrayfun(@(f) isequal(f.sysID,s) && strcmp(char(string(f.obsType)), ...
                       char(string(ot))) && isequal(f.day,d) && isequal(f.arc,a), F);
    m = [F(hit).idx];
end
