function res = residuals(acc, sol, opt)
%LSQVCE.ENGINE.RESIDUALS  Second pass over the data to compute residuals.
%
%   res = LSQVCE.ENGINE.RESIDUALS(acc, sol, 'Name', value, ...)
%
%   Nothing is retained after accumulation, so residuals require a SECOND
%   SWEEP through the loader. The local parameters are recovered by back-
%   substitution from the stored Schur blocks, then v = l - A*x - B*y is formed
%   arc by arc. The loader is called with exactly the same (day, obsType, arc)
%   sequence as during accumulation and MUST return the same data.
%
%   INPUT
%     acc   the LSQVCE.ENGINE.BLOCKACCUMULATOR, constructed with 'KeepLocal', true
%     sol   the struct returned by acc.solve()
%
%   OPTIONS
%     'whitened'   true (default)  also return w = W*v, whose sum of squares is
%                                  the quadratic form v'*inv(C)*v
%     'collect'    true (default)  return the per-observation residual series.
%                                  Set false for a memory-light sigma0 check.
%     'verbose'    true (default)
%
%   OUTPUT
%     res.sys(k).sysID
%     res.sys(k).y            local parameter estimate
%     res.sys(k).Cyy          its covariance (commons' uncertainty propagated)
%     res.sys(k).vtv          weighted residual square sum of this system
%     res.sys(k).m            number of observations
%     res.sys(k).rms          sqrt(vtv/m), in units of sigma0
%     res.sys(k).obs(j)       .obsType .day .arc .t .v .w .nAxis .Cep
%                             (.nAxis, .Cep: the arc's axis count and its
%                             per-epoch error information as the loader
%                             returned them -- the joint estimator's robust
%                             pass descales .v by them, R6 WP-105)
%     res.vtv, res.dof, res.s0sq, res.sigma0    recomputed globally -- WITH
%                             a common prior's pseudo-observations: its
%                             residual |(K x - x0)/sigma|^2 in vtv and its nk
%                             rows in dof, as the solve counts them (R7 WP-114,
%                             R7-15; without them every fit with data.priors
%                             warned sigma0Mismatch and blamed the loader)
%     res.prior               [] | .nk .vtv  the common prior's share
%     res.supersedesSol       1x1 logical  true when THESE numbers supersede
%                             sol.vtv/sol.sigma0 (WP-49): either the backend
%                             NaN-poisoned a cancellation-destroyed vtv, or
%                             the two routes disagree beyond 1e-6
%     res.obsType(k)          .obsType .m .vtv .rms
%                             per-obsType QUICK-LOOK: rms = sqrt(vtv/m) in
%                             units of sigma0. Uses m, not the exact
%                             redundancy, so it is biased low by
%                             1 - sqrt(r/m) (negligible when m >> u). Use it
%                             to see "kbr at 1.0, pos at 3.4" before paying
%                             for a full LSQVCE.ENGINE.VCE, not as a variance component.
%
%   WARNINGS
%     lsqvce:residuals:supersedesSol   the accumulated sigma0 is non-finite
%                                      (backend NaN-poisoned it); quote
%                                      res.sigma0/res.vtv instead
%     lsqvce:residuals:sigma0Mismatch  finite but disagreeing routes: the
%                                      loader changed between passes
%
%   ACCURACY
%     res.vtv is the DIRECT sum of squared (whitened) residuals, immune to
%     the l'Pl - n'x cancellation of the accumulated route; the two agree
%     to 1e-10 relative on healthy fixtures (regression-locked,
%     test_precision) and res is the authority when they do not (WP-49).
%
%   THE CONSISTENCY CHECK
%     res.sigma0 must equal sol.sigma0. They are computed by completely
%     different routes -- one from accumulated scalars during the first pass,
%     one from explicit residuals during the second -- so agreement is a strong
%     check that the whole chain is correct. A mismatch is reported as a
%     warning and means one of:
%
%       * the loader is not deterministic (random data selection, a changing
%         option that alters the design matrix between passes)
%       * a cache key in the loader is missing a field that opt varies over
%       * the covariance model was changed between the two calls
%
%     Do not ignore that warning: if the two passes disagree, at least one of
%     them is wrong, and you do not know which.
%
%   READING THE RESIDUALS
%     v  is in the original metric (metres, m/s, ...) and is still COLOURED --
%        plot it to see the physics: unmodelled signal, gaps, manoeuvres.
%     w  is whitened and should be white noise with unit variance. It is the
%        diagnostic that tests the STOCHASTIC model:
%
%          std(w) ~ 1               the noise level is right
%          acf(w) ~ delta           the noise COLOUR is right
%
%        If std(w) is fine but the ACF of w still shows structure, the AR order
%        is too low -- go back to LSQVCE.ENGINE.CHECKPSD and raise it. This is the single
%        most useful plot in the whole toolbox.
%
%   EXAMPLE
%     sol = acc.solve();
%     res = acc.residuals(sol);
%
%     w = vertcat(res.sys(1).obs.w);
%     lsqvce.util.log.msg('std(w) = %.4f (want 1)\n', std(w));
%
%     [a, lags] = xcorr(w - mean(w), 100, 'coeff');
%     stem(lags(lags>=0), a(lags>=0)); yline([-1 1]*1.96/sqrt(numel(w)), 'r--');
%     title('ACF of whitened residuals -- should be white');
%
%   ERRORS
%     lsqvce:residuals:badSol       the solution is not from this accumulator
%     lsqvce:residuals:noStore      nothing was retained; re-run with 'store','all'
%
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-130: dof counts a rank-deficient
%   local block by its rank (BlockAccumulator.localRankOf). With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (b) b2 (F-07):
%   typed `arguments` signature, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   See also LSQVCE.ENGINE.BLOCKACCUMULATOR, LSQVCE.ENGINE.CHECKPSD
%
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-114 (R7-15): the common prior
%   enters vtv and dof (res.prior), one estimator, one sigma0. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%
%   Updated 2026-08-20 (CEST, UTC+02:00), WP-49: NaN-aware supersede
%   logic + res.supersedesSol; this pass is the documented authority for
%   a cancellation-poisoned accumulated vtv. With assistance from Claude
%   (Anthropic), model Claude Fable 5.
%   Generated with assistance from Claude (Anthropic), model Claude
%   Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

arguments
    acc
    sol
    opt.Whitened (1,1) logical = true
    opt.Collect (1,1) logical = true
    opt.Verbose (1,1) logical = true
end

    narginchk(2,inf);

    % WP-39 stage (b) b2: typed `arguments` signature (canonical
    % UpperCamelCase; legacy casings keep matching -- MATLAB NV
    % semantics; enum VALUES case-sensitive lowercase).
    o = struct('whitened', opt.Whitened, 'collect', opt.Collect, 'verbose', opt.Verbose);

    if ~isstruct(sol) || ~isfield(sol,'x')
        error('lsqvce:residuals:badSol', 'sol must be the struct returned by acc.solve().');
    end

    x   = sol.x(:);
    ids = acc.systemIDs();
    if isempty(ids)
        error('lsqvce:residuals:noStore', ...
            ['No local blocks were kept. Construct the accumulator with ' ...
             '''KeepLocal'', true.']);
    end

    res = struct('sys', struct([]), 'vtv', 0, 'dof', 0, 's0sq', NaN, ...
                 'sigma0', NaN, 'obsType', struct([]));

    vtvTot  = 0;
    mTot    = 0;
    nParTot = numel(x);
    otMap   = containers.Map('KeyType','char','ValueType','any');
    otKeys  = {};

    stageDetail = lsqvce.util.verbosity.detail(o.verbose);
    if stageDetail
        % Mirrors the accumulation block: same heading -> day -> line
        % shape, so the SECOND pass over the data reads like the first.
        % Two passes reporting in two different layouts make a log harder
        % to follow than the work actually is.
        lsqvce.util.log.msg('%s', lsqvce.util.verbosity.heading('Compute the residuals'));
    end
    for k = 1:numel(ids)
        st = acc.getStore(ids(k));

        % ---- back-substitution of the local block ------------------------
        y   = lsqvce.engine.BlockAccumulator.localSolveOf(st, x);   % WP-128: on the factor if kept
        Cyy = acc.localCovariance(ids(k), sol.Cxx, sol.s0sq);   % R7-40: the ONE formula

        S = struct('sysID', ids(k), 'y', y, 'Cyy', Cyy, ...
                   'vtv', 0, 'm', 0, 'rms', NaN, 'obs', struct([]));

        uo = st.uo;
        % The readable label travels on UserOpt, published by
        % BlockAccumulator.processSystem, so the SECOND pass over the data
        % labels its lines exactly as the first did. Falls back to the id.
        sysName = char(string(ids(k)));
        if isfield(uo, 'systemName') && ~isempty(uo.systemName)
            sysName = char(uo.systemName);
        end
        if isnumeric(ids(k)) && isscalar(ids(k))
            sysIdx = double(ids(k));
        else
            sysIdx = k;
        end
        for d = st.dayList(:).'
            tDay = tic;
            mDay0 = S.m;
            meta = st.dataFcn('meta', d, [], uo);

            dayType = struct('obsType', {}, 'm', {}, 'seconds', {});
            for ot = 1:numel(meta.obsTypes)
                uo.obsType = meta.obsTypes{ot};
                tType = tic;  mType0 = S.m;

                for a = meta.arcs{ot}(:).'
                    uo.day = d;  uo.arc = a;
                    arc = st.dataFcn(d, uo.obsType, a, uo);
                    if isempty(arc) || ~isfield(arc,'l') || isempty(arc.l)
                        continue
                    end

                    [v, w] = arcResiduals(arc, x, y, uo, o);

                    S.vtv = S.vtv + w.'*w;
                    S.m   = S.m   + numel(v);

                    if ~otMap.isKey(uo.obsType)
                        otMap(uo.obsType) = struct('m', 0, 'vtv', 0);
                        otKeys{end+1} = uo.obsType;            %#ok<AGROW>
                    end
                    g = otMap(uo.obsType);
                    g.m   = g.m + numel(v);
                    g.vtv = g.vtv + w.'*w;
                    otMap(uo.obsType) = g;

                    if o.collect
                        e = struct('obsType', uo.obsType, 'day', d, 'arc', a, ...
                                   't', lsqvce.util.getfielddef(arc,'t',[]), ...
                                   'nAxis', lsqvce.util.getfielddef(arc, 'nAxis', 1), ...   % R6 WP-105
                                   'Cep', lsqvce.util.getfielddef(arc, 'Cep', []), ...
                                   'v', v, 'w', []);
                        if o.whitened, e.w = w; end
                        if isempty(S.obs), S.obs = e; else, S.obs(end+1) = e; end
                    end
                end
                dayType(end+1) = struct('obsType', char(uo.obsType), ...
                                        'm', S.m - mType0, ...
                                        'seconds', toc(tType));   %#ok<AGROW>
            end
            if stageDetail
                % One line per observation type, identical in shape to the
                % accumulation block. The previous form collapsed all types into
                % ONE line per day, so a system with three streams appeared to
                % have one -- and it passed the count into the obsType slot, so
                % the line read: type "4200 residual(s)" for "1".
                lsqvce.util.log.msg('%s', lsqvce.util.verbosity.dayLine(d));
                for t = 1:numel(dayType)
                    lsqvce.util.log.msg('%s', lsqvce.util.verbosity.typeLine( ...
                        sysIdx, t, dayType(t).obsType, sysName, ...
                        sprintf('%d obs', dayType(t).m), dayType(t).seconds));
                end
            elseif lsqvce.util.verbosity.progress(o.verbose)
                lsqvce.util.log.msg('%ssystem %s, day %s : %d residual(s)   %.2f s\n', ...
                    lsqvce.util.verbosity.banner('Residuals'), string(ids(k)), string(d), ...
                    S.m - mDay0, toc(tDay));
            end
        end

        S.rms   = sqrt(S.vtv / max(S.m,1));
        vtvTot  = vtvTot + S.vtv;
        mTot    = mTot   + S.m;
        nParTot = nParTot + lsqvce.engine.BlockAccumulator.localRankOf(st);   % WP-130: the rank

        if isempty(res.sys), res.sys = S; else, res.sys(end+1) = S; end

        if lsqvce.util.verbosity.progress(o.verbose)
            % A SUMMARY, and it now says so. It used to open with a bare system
            % id followed by 'm = ... rms(w) = ... |y| = ...', which read like
            % one more per-day line rather than the total for the system.
            lsqvce.util.log.msg('%s', lsqvce.util.verbosity.fieldBlock( ...
                sprintf('Summary for system "%s"', sysName), ...
                {'observations', sprintf('%d', S.m); ...
                 'rms(w)',       sprintf('%.4g', S.rms); ...
                 '|y|',          sprintf('%.4g', norm(y))}));
        end
    end

    res.vtv    = vtvTot;
    % R7 WP-114 (R7-15): the common prior is part of the estimate -- the solve
    % counts its rows in m and its residual in vtv, so this pass does too
    res.prior = [];
    nk = 0;
    if ~isempty(acc.prior)
        res.prior = struct('nk', acc.prior.nk, 'vtv', acc.prior.vtv(x));
        vtvTot = vtvTot + res.prior.vtv;  res.vtv = vtvTot;  nk = acc.prior.nk;
    end
    res.dof    = mTot + nk - nParTot;
    res.s0sq   = vtvTot / max(res.dof, 1);
    res.sigma0 = sqrt(res.s0sq);

    % ---- per-obsType sigma0 QUICK-LOOK -------------------------------------
    % rms = sqrt(vtv/m) uses m, not the exact redundancy r, so it is biased
    % LOW by exactly 1 - sqrt(r/m) -- negligible at high redundancy (measured
    % < 2% at u/m ~ 2e-3, v4verify.py) but real for small groups. It costs
    % nothing and tells you "kbr at 1.02, pos at 3.4" BEFORE you decide to pay
    % for a full LSQVCE.ENGINE.VCE. It is a triage number, not a variance component.
    for k = 1:numel(otKeys)
        g = otMap(otKeys{k});
        e = struct('obsType', otKeys{k}, 'm', g.m, 'vtv', g.vtv, ...
                   'rms', sqrt(g.vtv / max(g.m, 1)));
        if isempty(res.obsType), res.obsType = e;
        else,                    res.obsType(end+1) = e; end       %#ok<AGROW>
    end

    if lsqvce.util.verbosity.progress(o.verbose)
        % A TITLED BLOCK, like every other summary. It used to put the title in
        % the LABEL COLUMN with 'vtv = ...' on the same line, and then indent the
        % per-type rows by that column's width to line up under it -- so the
        % section had no title line of its own and everything sat 21 characters
        % to the right of every other block on screen.
        lsqvce.util.log.msg('%s', lsqvce.util.verbosity.heading('Residuals'));
        for k = 1:numel(res.obsType)
            e = res.obsType(k);
            lsqvce.util.log.msg('%s', lsqvce.util.verbosity.labelledRow(e.obsType, ...
                sprintf('%7d obs   rms(w) = %.4f', e.m, e.rms)));
        end
        % v'Pv spelled out: 'vtv' is the code's name for it and means nothing to
        % a reader. It is the WEIGHTED sum of squares -- w is the whitened
        % residual, so w'w = v'Pv -- and it is the numerator of sigma0^2.
        lsqvce.util.log.msg('%s', lsqvce.util.verbosity.labelledRow( ...
            'weighted sum of squares v''Pv', sprintf('%.8g', res.vtv)));
        lsqvce.util.log.msg('%s', lsqvce.util.verbosity.labelledRow( ...
            'degrees of freedom', sprintf('%d', res.dof)));
        lsqvce.util.log.msg('%s', lsqvce.util.verbosity.labelledRow( ...
            'sigma0', sprintf('%.8g', res.sigma0)));
    end

    % --- cross-check against the accumulation ------------------------------
    % WP-49: three outcomes, each named. (1) The accumulated sigma0 is
    % NON-FINITE -- the backend NaN-poisoned a cancellation-destroyed vtv;
    % this pass's number SUPERSEDES it (that is the documented repair path,
    % and the old relative-difference test was BLIND to NaN: |x - NaN| > tol
    % is false). (2) Finite but mismatched beyond 1e-6 -- the loader did not
    % return the same data. (3) Agreement -- nothing to say.
    res.supersedesSol = false;
    if ~isfinite(sol.sigma0)
        res.supersedesSol = true;
        lsqvce.util.log.warn('lsqvce:residuals:supersedesSol', ...
            ['The accumulated sigma0 is %g -- the backend refused a ' ...
             'cancellation-destroyed v''Pv rather than return a wrong ' ...
             'plausible number. THIS pass recomputed it honestly: ' ...
             'sigma0 = %.8g, vtv = %.8g. QUOTE res.sigma0/res.vtv, not ' ...
             'sol.sigma0/sol.vtv. To avoid the problem entirely, use ' ...
             '''Backend'',''qr''.'], sol.sigma0, res.sigma0, res.vtv);
    else
        rel = abs(res.sigma0 - sol.sigma0) / max(sol.sigma0, realmin);
        if rel > 1e-6
            res.supersedesSol = true;
            lsqvce.util.log.warn('lsqvce:residuals:sigma0Mismatch', ...
                ['sigma0 from the residual pass (%.8g) disagrees with the ' ...
                 'accumulation (%.8g), relative difference %.2e. The two are ' ...
                 'computed by independent routes, so a mismatch means the loader ' ...
                 'did NOT return the same data on the second pass. Check for ' ...
                 'non-deterministic data selection, or a cache key that is missing ' ...
                 'a field that opt varies over.'], res.sigma0, sol.sigma0, rel);
        end
    end
end

%==========================================================================
function [v, w] = arcResiduals(arc, x, y, uo, o)
%ARCRESIDUALS  Residuals of one arc, batched exactly as in accumulation.

    m = size(arc.A,1);
    v = arc.l(:) - arc.A*x - arc.B*y;      % raw residuals, original metric
    w = v;                                 % overwritten below if whitened

    if ~o.whitened
        return
    end

    hasW = isfield(arc,'Wfcn') && ~isempty(arc.Wfcn);
    hasC = isfield(arc,'Cfcn') && ~isempty(arc.Cfcn);
    hasM = isfield(arc,'W')    && ~isempty(arc.W);

    if ~(hasW || hasC || hasM)
        return                              % unit weights: w = v
    end

    % R6 WP-91 (E-09): the residual is whitened by the ONE iterator
    % arcBatchWhiten (engine-private) -- the same batch bounds, overlap and
    % whole-epoch rounding as accumulation -- on an arc whose design blocks
    % are empty and whose "observation" is v (what the old loop passed to
    % Wfcn as ze, ze, v(idx): byte-identical).
    av = struct('A', zeros(m, 0), 'B', zeros(m, 0), 'l', v);
    for f = {'Wfcn', 'Cfcn', 'W', 'filterMemory', 'nAxis', 'wholeArc'}
        if isfield(arc, f{1}), av.(f{1}) = arc.(f{1}); end
    end
    w = zeros(m,1);
    arcBatchWhiten(av, uo, @keepWhitened);

    function keepWhitened(rows, keep, ~, ~, lb)
        w(rows) = lb(keep);
    end
end
