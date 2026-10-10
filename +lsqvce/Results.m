classdef Results < handle
%LSQVCE.RESULTS  The outcome of a recovery: solution, diagnostics, provenance.
%   R = lsqvce.Results(sol, D, acc, diagSol, batchFcn, nBatch, systems)
%   R = lsqvce.Results(sol, D, acc, diagSol, batchFcn, nBatch, systems, meta)
%
%   Returned by every driver LSQVCE.FIT dispatches to -- LSQVCE.ENGINE.ADJUSTBLOCK,
%   ADJUSTVCEJOINT, ADJUSTROBUST, ADJUSTNONLINEAR. It gathers
%   everything a user needs to read and judge an adjustment in one object, and
%   keeps the lower layers reachable for experts. You rarely construct it
%   yourself; you receive it from adjustBlock and call R.report().
%
%   PROPERTIES (read-only)
%     x          nc x 1   double   Common-parameter estimate (the field).
%     sigma0     1x1      double   A-posteriori variance factor (global).
%     stopReason 1x1      string   Why the loop stopped (R5 M1): "single" |
%                                  "converged" | "maxPasses". Convergence
%                                  failure is a RESULT: the last pass is
%                                  returned with its covariance and one warning.
%     converged  1x1      logical  Did the pipeline converge? (true for a single
%                                  solve; set by the iterative engines.)
%     diag       1x1      struct   Per-observation diagnostics (LSQVCE.ENGINE.DIAGNOSE):
%                                  redundancy, leverage, w-tests, MDB, ...
%     history    1xI      struct   Per-iteration record (empty for a single
%                                  solve; filled by LSQVCE.ENGINE.ADJUSTVCEJOINT:
%                                  .maxDev .Shat .streams .sigma .white
%                                  .whitePval per pass, .skips, .nIter, and
%                                  under Robust the per-pass .robustScale (one
%                                  column per stream) .nDown .nZero -- R6 WP-105;
%                                  by LSQVCE.ENGINE.ADJUSTROBUST: .scale .dx
%                                  .nDown .nZero and, on the last pass, .factor
%                                  .keys).
%     noise      1x1      struct   Fitted noise models per stream (empty until
%                                  the variance-component loop populates it):
%                                  .stream (keys) .model{k} (.sigma the level,
%                                  .C the law) .Cgroups (one C per pool group)
%                                  .priors .priorHistory, and .robust
%                                  (R6 WP-105: .influence, .scale per stream,
%                                  .streams, .factor ONE PER EPOCH in sweep
%                                  order, .keys; empty fields under "none").
%     sol        1x1      struct   Raw output of acc.solve() (.x, .sigma0, ...).
%     acc        1x1      lsqvce.engine.BlockAccumulator   The accumulator (drop here
%                                  to re-solve, inspect blocks, get residuals).
%     source     1x1      struct   Bridged solution from LSQVCE.ENGINE.DIAGNOSEACCUM.
%     batchFcn   1x1      function_handle   Re-reader used by LSQVCE.ENGINE.DIAGNOSE.
%     nBatch     1x1      double   Batch count from the bridge.
%     systems    1xM      struct   The resolved system specs (sysID/nl/dayList).
%     meta       1x1      struct   Run provenance (WP-44). Field list is
%                                  ADD-ONLY across versions, never renamed:
%                                  .toolboxVersion 1xK char  from Contents.m,
%                                                  the ONE version source
%                                  .matlabRelease  1xK char  version('-release')
%                                  .createdUTC     1xK char  ISO 8601, trailing Z
%                                  .seeds          1x1 struct every RNG seed the
%                                                  run resolved (struct() when
%                                                  none; a resolved-but-unset
%                                                  seed is stored as [])
%                                  .options        1x1 struct the RESOLVED
%                                                  options the run actually used
%                                  A bare construction (7 arguments) gets a
%                                  meta with empty seeds/options -- the field
%                                  is never missing.
%
%   METHODS
%     V = R.report(Name,Value,...)   Print/return the level-0..4 verdict
%                                    (LSQVCE.REPORT); options pass straight through.
%     r = R.redundancy()             Total redundancy sum(r_i) = n - u.
%     res = R.residuals(Name,Value)  Whitened/raw residuals (acc.residuals).
%     con = R.contribution(...)      Per-group information share.
%     T = R.convergence(Name,Value)  Per-ITERATION trajectory: when the loop
%                                    settled (needs .history).
%     F = R.compare(Rred)            Nested-model F test against a reduced fit (R5 M7).
%     fh = R.plot(kind)              "residuals" | "contribution" | "psd" | "convergence"
%                                    (R5 M5: one figure per question, nothing exported).
%     [C, L, map] = R.toTable(data)  The expanded parameter tables with estimate and
%                                    standard deviation written in (R5-D8).
%     T = R.summary()                Per-STREAM final state: sigma, AR order,
%                                    whiteness (needs .history).
%     T = R.coefficients(Alpha=a)    The estimates as a table with formal
%                                    errors and (1-a) confidence bounds
%                                    (R4 WP-57; default a = 0.05).
%     disp(R)                        One-screen summary.
%
%   THE TWO ITERATION READERS ARE COMPLEMENTARY. convergence() reduces across
%   STREAMS (one row per pass); summary() reduces across PASSES (one row per
%   stream). A converged loop can still hide a single stream that never
%   whitened -- the per-pass view reports the worst p-value without naming it.
%
%   EXAMPLE
%     R = lsqvce.engine.adjustBlock(systems, days);
%     R.report('Level', 2);
%     if ~R.converged, lsqvce.util.log.warn('did not converge'); end
%     x = R.x;   s0 = R.sigma0;
%
%   ERRORS
%     lsqvce:Results:noVPv          compare() on a Results without a finite v'Pv
%     lsqvce:Results:noHistory      an iteration reader was called on a
%                                   single-solve Results (adjustBlock does not
%                                   iterate); use the result of an iterating
%                                   fit (Estimate = "noise", Robust, Nonlinear)
%     lsqvce:Results:noCovariance   coefficients() on a Results whose sol
%                                   carries no Cxx
%     lsqvce:Results:noAccumulator  residuals()/contribution() on a Results
%                                   without an accumulator (fit kinds B / C)
%
%   ACCURACY
%     Results adds no numerics: x and sigma0 are the engine's (sigma0 is the
%     A-POSTERIORI std of unit weight; sol.Cxx the covariance scaled by it --
%     the streaming contract every producer must establish, R4 §11.8 L4 r1),
%     report() is lsqvce.report on .diag, coefficients() is arithmetic on x
%     and diag(Cxx) with z(0.05) = 1.959963984540054 in closed form.
%
%   See also LSQVCE.FIT, LSQVCE.ENGINE.ADJUSTBLOCK, LSQVCE.ENGINE.ADJUSTVCEJOINT, LSQVCE.REPORT,
%            LSQVCE.ENGINE.DIAGNOSE, LSQVCE.ENGINE.BLOCKACCUMULATOR
%
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-118 part 1 (R7-34): the
%   noHistory remedy no longer names the deleted adjustVCE. With assistance
%   from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-09-08 (CEST, UTC+02:00), R6 WP-99: R.report() forwards the
%   loop state, so the verdict says whether the stochastic model settled.
%   With assistance from Claude (Anthropic), model Claude Opus 5.
%   Updated 2026-09-01 (CEST, UTC+02:00), R4 WP-61 (Level L3): ACCURACY block
%   (the sigma0 / Cxx contract stated once). With assistance from Claude
%   (Anthropic), model Claude Fable 5.
%   Updated 2026-08-31 (CEST, UTC+02:00), R4 WP-56b (Level L4): residuals()
%   and contribution() refuse with lsqvce:Results:noAccumulator on a Results
%   without an accumulator (the in-memory / batched kinds of fit). With
%   assistance from Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-31 (CEST, UTC+02:00), R4 WP-57 (Level L2): coefficients()
%   -- estimates, formal errors and confidence bounds as ONE table (the
%   fitlm Coefficients/coefCI pair folded into one method; z from erfcinv,
%   no Statistics Toolbox). With assistance from Claude (Anthropic), model
%   Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-44 (T-05): run provenance --
%   immutable .meta composed by the ONE private builder runMeta, constructor
%   widened to an optional 8th argument, report() forwards the provenance
%   line. With assistance from Claude (Anthropic), model Claude Fable 5.
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (c) c2 (F-07):
%   typed `arguments` signatures, canonical UpperCamelCase option names
%   (legacy casings keep matching -- MATLAB NV semantics; enum VALUES
%   case-sensitive lowercase). With assistance from Claude (Anthropic),
%   model Claude Fable 5.
%   Created 2026-07-16 12:20 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.

    properties (SetAccess = immutable)
        x
        sigma0
        sol
        diag
        acc
        source
        batchFcn
        nBatch
        systems
        meta       % WP-44: run provenance (see the PROPERTIES help block)
    end
    properties (SetAccess = private)
        converged  = true;     % single solve is, by definition, done
        stopReason = "single"; % R5 M1: WHY the loop stopped -- "single" (one
                               % solve, nothing iterated), "converged" (every
                               % criterion met), "maxPasses" (budget exhausted;
                               % the last pass is returned WITH a warning),
                               % "userStop" (Options.OnPass returned false, R5 M5)
        history    = struct([]);
        noise      = struct([]);
    end

    methods
        function obj = Results(sol, D, acc, diagSol, batchFcn, nBatch, systems, meta)
        %RESULTS  Construct a Results wrapper around a solve + diagnostics (see header).
            narginchk(7, 8);
            if nargin < 8
                % A bare construction still carries provenance -- composed by
                % the ONE builder with empty seeds/options (WP-44 robustness
                % promise: never a missing meta, never a missing field).
                meta = lsqvce.engine.runMeta(struct(), struct());
            end
            obj.sol      = sol;
            obj.x        = sol.x(:);
            obj.sigma0   = sol.sigma0;
            obj.diag     = D;
            obj.acc      = acc;
            obj.source   = diagSol;
            obj.batchFcn = batchFcn;
            obj.nBatch   = nBatch;
            obj.systems  = systems;
            obj.meta     = meta;
        end

        function V = report(obj, varargin)
        %REPORT  Judge the adjustment in words (LSQVCE.REPORT), levels 0-4.
        %   Forwards the run provenance so the printed level-0 block carries
        %   the one provenance line (WP-44), and the LOOP STATE so level 1 can
        %   say whether the stochastic model settled (R6 WP-99: every number
        %   from level 1 down rests on those weights). A single solve reports
        %   stopReason "single" and report() then says nothing about it. User
        %   options pass through and a user-supplied 'Provenance' / 'Iteration'
        %   would win (last-match NV semantics). A robust run also forwards its
        %   influence (Iteration.influence), so the verdict says the sigma0 and
        %   covariance are those of the reweighted problem (R7-18; updated
        %   2026-10-10 with assistance from Claude (Anthropic), model Claude
        %   Opus 5.5).
            nPass = 1;                       % a single solve is one pass
            if ~isempty(obj.history) && isfield(obj.history, 'nIter')
                nPass = obj.history.nIter;
            end
            it = struct('converged', obj.converged, 'stopReason', obj.stopReason, ...
                        'nIter', nPass, 'influence', '');
            % R7-18: a robust fit says so -- adjustRobust records the influence
            % in R.noise.influence, the joint loop's robust pass in R.noise.robust
            if isstruct(obj.noise) && isfield(obj.noise, 'influence')
                it.influence = char(obj.noise.influence);
            elseif isstruct(obj.noise) && isfield(obj.noise, 'robust') && ...
                    isstruct(obj.noise.robust) && isfield(obj.noise.robust, 'influence')
                it.influence = char(obj.noise.robust.influence);
            end
            V = lsqvce.report(obj.diag, 'Provenance', obj.meta, 'Iteration', it, ...
                              varargin{:});
        end

        function r = redundancy(obj)
        %REDUNDANCY  Total redundancy sum(r_i) = n - u, from the diagnostics.
            if isfield(obj.diag, 'r'),        r = sum(obj.diag.r(:));
            elseif isfield(obj.diag, 'redundancy'), r = sum(obj.diag.redundancy(:));
            else, r = NaN;
            end
        end

        function fh = plot(obj, kind, opt)
        %PLOT  One figure per question, from the Results alone (R5 M5 5c).
        %   fh = R.plot("residuals")     whitened residual per stream vs time
        %   fh = R.plot("contribution")  information share per stream
        %   fh = R.plot("psd")           the estimated spectral shape(s), and the
        %                                per-stream level (Estimate "noise" only)
        %   fh = R.plot("convergence")   the pass record (max|alpha2-1|, whiteness)
        %   fh = R.plot(kind, Stream=key, Visible=false)
        %   Replaces the v4 plot hooks (plotPostSolve / plotPostModel,
        %   ShowSeconds / ExportDir): nothing is exported by the toolbox --
        %   the figure handle is returned; saveas / exportgraphics are yours.
        %
        %   INPUT   kind 1x1 string as above;  Stream '' | 'satKey|obsType'
        %           (residuals only);  Visible logical (default true)
        %   OUTPUT  fh 1x1 matlab.ui.Figure
        %   ERRORS  lsqvce:Results:badPlotKind   an unknown kind
        %           lsqvce:Results:noHistory     "psd"/"convergence" on a single solve
        %           lsqvce:Results:noAccumulator "residuals"/"contribution" on a Results
        %                                        without an accumulator
        arguments
            obj
            kind (1,1) string
            opt.Stream (1,1) string = ""
            opt.Visible (1,1) logical = true
        end
            switch lower(char(kind))
                case 'residuals'
                    res = obj.residuals('Verbose', false);
                    fh = figure('Visible', opt.Visible, 'Name', 'lsqvce: whitened residuals');
                    n = 0;  panels = {};
                    for k = 1:numel(res.sys)
                        S = res.sys(k);
                        for j = 1:numel(S.obs)
                            key = sprintf('%s|%s', matlab.lang.makeValidName(char(string(S.sysID))), S.obs(j).obsType);
                            if strlength(opt.Stream) > 0 && ~strcmp(key, char(opt.Stream)), continue, end
                            panels{end+1} = struct('key', key, 't', S.obs(j).t(:), 'w', S.obs(j).w(:)); %#ok<AGROW>
                        end
                    end
                    streamKeys = unique(cellfun(@(p) p.key, panels, 'UniformOutput', false), 'stable');
                    for k = 1:numel(streamKeys)
                        ax = subplot(numel(streamKeys), 1, k, 'Parent', fh);
                        sel = panels(strcmp(cellfun(@(p) p.key, panels, 'UniformOutput', false), streamKeys{k}));
                        hold(ax, 'on');
                        for j = 1:numel(sel), plot(ax, sel{j}.t, sel{j}.w, '.', 'MarkerSize', 3); end
                        hold(ax, 'off');  grid(ax, 'on');
                        ylabel(ax, 'w (whitened)');  title(ax, streamKeys{k}, 'Interpreter', 'none');
                        n = n + 1;
                    end
                    if n == 0, error('lsqvce:Results:badPlotKind', 'No stream matches Stream="%s".', opt.Stream); end
                    xlabel(ax, 't');
                case 'contribution'
                    con = obj.contribution();
                    fh = figure('Visible', opt.Visible, 'Name', 'lsqvce: contribution');
                    names = arrayfun(@(g) char(string(g.key)), con.group, 'UniformOutput', false);
                    ax = axes('Parent', fh);
                    bar(ax, [con.group.share]);
                    set(ax, 'XTick', 1:numel(names), 'XTickLabel', names, 'TickLabelInterpreter', 'none');
                    ylabel(ax, 'information share');  ylim(ax, [0 1]);  grid(ax, 'on');
                    title(ax, 'share of the common parameters per stream');
                case 'psd'
                    h = obj.requireHistory('plot("psd")');
                    fh = figure('Visible', opt.Visible, 'Name', 'lsqvce: estimated spectra');
                    ax = axes('Parent', fh);
                    if isfield(obj.noise, 'Cgroups') && ~isempty(obj.noise.Cgroups) && obj.noise.Cgroups.Count > 0
                        ks = obj.noise.Cgroups.keys;
                        hold(ax, 'on');
                        for k = 1:numel(ks)
                            C = obj.noise.Cgroups(ks{k});  L = numel(C);
                            fN = (0:L-1).' / (2*(L-1));
                            S = 2 * (C(1) + 2 * (cos(2*pi*fN*(1:L-1)) * C(2:end)));
                            semilogy(ax, fN, max(S, eps*max(S)), 'DisplayName', ks{k});
                        end
                        hold(ax, 'off');  legend(ax, 'Interpreter', 'none');
                    elseif isfield(h, 'Shat') && ~isempty(h.Shat)
                        semilogy(ax, linspace(0, 0.5, size(h.Shat, 2)), max(h.Shat(end, :), realmin));
                    else
                        error('lsqvce:Results:noHistory', 'No estimated spectrum on this Results.');
                    end
                    xlabel(ax, 'normalised frequency (cycles / sample)');  ylabel(ax, 'S (shape)');
                    grid(ax, 'on');  title(ax, 'estimated spectral shape per pool group');
                case 'convergence'
                    h = obj.requireHistory('plot("convergence")');
                    fh = lsqvce.Results.plotConvergence(h, char(opt.Stream));
                    set(fh, 'Visible', opt.Visible);
                otherwise
                    error('lsqvce:Results:badPlotKind', ...
                          'plot kind must be "residuals", "contribution", "psd" or "convergence" (got "%s").', kind);
            end
        end

        function res = residuals(obj, varargin)
        %RESIDUALS  Collected residuals from the accumulator (acc.residuals).
            obj.requireAccumulator('residuals');
            res = obj.acc.residuals(obj.sol, varargin{:});
        end

        function [C, L, map] = toTable(obj, data, opt)
        %TOTABLE  Write the estimates back into the parameter tables (R5-D8 round 4).
        %   [C, L, map] = R.toTable(data)
        %   [C, L, map] = R.toTable(data, Estimate=9, Std=10)
        %   The EXPANDED tables -- one row per instantiated parameter, in the
        %   toolbox's column order -- with the estimate and its standard
        %   deviation written into the columns named (defaults 9 and 10 of the
        %   15-column table: your columns 8 and 9 after the scope token). Every
        %   other column is copied from the template row. A 'day'/'arc' row
        %   appears once per instance; map says which.
        %   INPUT   data  the declaration this Results came from (its tables are
        %                 re-expanded with lsqvce.engine.declaration -- same plans)
        %           Estimate, Std   1-based columns to write (padded if the table
        %                 is narrower)
        %   OUTPUT  C    nc x W cell   the common table, one row per column of A
        %           L    1xK cell of (nl_k x W) cells, one per system, block order
        %                (per day: arc blocks, day block; then the system block)
        %           map  struct: .common (nc x 1 struct: row, day), .local{k}
        %                (nl_k x 1 struct: row, day, arc) -- the instance of each line
        %   ERRORS  lsqvce:Results:noAccumulator   no accumulator on this Results
        %           lsqvce:Results:tableMismatch    the declaration's nc is not this Results' u
        %           lsqvce:fit:*                    the declaration does not validate
        arguments
            obj
            data (1,1) struct
            opt.Estimate (1,1) double {mustBeInteger, mustBePositive} = 9
            opt.Std (1,1) double {mustBeInteger, mustBePositive} = 10
        end
            if isempty(obj.acc)
                error('lsqvce:Results:noAccumulator', 'toTable needs a streaming (kind A) Results.');
            end
            [spec, systems] = lsqvce.engine.declaration(data);
            cp = spec.commonPlan;  T = spec.commonParams;
            if cp.nc ~= numel(obj.x)
                error('lsqvce:Results:tableMismatch', ...
                      ['The declaration expands to %d common parameters but this Results has %d: ' ...
                       'toTable needs the declaration the Results came from.'], cp.nc, numel(obj.x));
            end
            sx = sqrt(max(real(diag(obj.sol.Cxx)), 0));
            W = max([size(T, 2), opt.Estimate, opt.Std]);
            C = cell(cp.nc, W);  mc = struct('row', cell(cp.nc, 1), 'day', cell(cp.nc, 1));
            for j = 1:cp.commonCount
                r = cp.rows.common(j);
                C(j, 1:size(T, 2)) = T(r, :);  C{j, opt.Estimate} = obj.x(j);  C{j, opt.Std} = sx(j);
                mc(j).row = r;  mc(j).day = [];
            end
            for i = 1:numel(cp.days)
                if cp.dayCount(i) == 0, continue, end
                rows = cp.rows.day;
                if isfield(cp, 'dayRowsOf'), rows = cp.dayRowsOf{i}; end
                for q = 1:cp.dayCount(i)
                    j = cp.dayOffset(i) + q;  r = rows(q);
                    C(j, 1:size(T, 2)) = T(r, :);  C{j, opt.Estimate} = obj.x(j);  C{j, opt.Std} = sx(j);
                    mc(j).row = r;  mc(j).day = cp.days(i);
                end
            end
            res = obj.residuals('Verbose', false);          % the local solutions y, Cyy
            K = numel(systems);  L = cell(1, K);  ml = cell(1, K);
            for k = 1:K
                Lk = systems{k};  plan = Lk.plan;
                Tl = data.systems(k).localParams;
                if isempty(Tl), L{k} = cell(0, W);  ml{k} = struct('row', {}, 'day', {}, 'arc', {});  continue, end
                s = res.sys(k);                                 % the accumulator's system order = data.systems order
                y = s.y(:);  sy = sqrt(max(real(diag(s.Cyy)), 0));
                Wl = max([size(Tl, 2), opt.Estimate, opt.Std]);
                Lt = cell(plan.nl, Wl);  m = struct('row', cell(plan.nl, 1), 'day', cell(plan.nl, 1), 'arc', cell(plan.nl, 1));
                for i = 1:numel(plan.days)
                    d = plan.days(i);
                    arcRows = plan.rows.arc;  dayRows = plan.rows.day;
                    if isfield(plan, 'arcRowsOf'), arcRows = plan.arcRowsOf{i};  dayRows = plan.dayRowsOf{i}; end
                    for a = 1:numel(plan.arcId{i})
                        for q = 1:numel(arcRows)
                            j = plan.offset{i}(a) + q;  r = arcRows(q);
                            Lt(j, 1:size(Tl, 2)) = Tl(r, :);  Lt{j, opt.Estimate} = y(j);  Lt{j, opt.Std} = sy(j);
                            m(j).row = r;  m(j).day = d;  m(j).arc = plan.arcId{i}(a);
                        end
                    end
                    for q = 1:plan.dayCount(i)
                        j = plan.dayOffset(i) + q;  r = dayRows(q);
                        Lt(j, 1:size(Tl, 2)) = Tl(r, :);  Lt{j, opt.Estimate} = y(j);  Lt{j, opt.Std} = sy(j);
                        m(j).row = r;  m(j).day = d;  m(j).arc = [];
                    end
                end
                for q = 1:plan.sysCount
                    j = plan.sysOffset + q;  r = plan.rows.system(q);
                    Lt(j, 1:size(Tl, 2)) = Tl(r, :);  Lt{j, opt.Estimate} = y(j);  Lt{j, opt.Std} = sy(j);
                    m(j).row = r;  m(j).day = [];  m(j).arc = [];
                end
                L{k} = Lt;  ml{k} = m;
            end
            map = struct('common', mc, 'local', {ml});
        end

        function F = compare(obj, Rred, opt)
        %COMPARE  Is the extra parameter earning its place? (nested-model F test)
        %   F = R.compare(Rred)            R = the FULL model, Rred = the REDUCED (nested) one
        %   F = R.compare(Rred, Alpha=0.01, Verbose=false)
        %   R5 M7 7c: lsqvce.engine.ftest re-homed onto the Results -- the two
        %   fits must share observations and weights, the reduced one a SUBSET
        %   of the parameters; the test is F = [(vPv_red - vPv_full)/(u_full -
        %   u_red)] / s0sq_full ~ F(du, r_full) under "the extra parameters are
        %   zero". Both Results must carry a honest v'Pv (NaN-honest folds refuse
        %   with lsqvce:Results:noVPv). The residual-sum PREMISE (same rows) is
        %   checked by the engine when both Results have diagnostics.
        %   INPUT   Rred 1x1 lsqvce.Results; Alpha (0.05); Verbose (false)
        %   OUTPUT  F    1x1 struct  the engine's record (.F .p .df1 .df2 .significant
        %                            .alpha .valid) -- see LSQVCE.ENGINE.FTEST
        %   ERRORS  lsqvce:Results:noVPv    a v'Pv is NaN or missing
        %           lsqvce:ftest:*          the engine's own refusals (not nested, ...)
        arguments
            obj
            Rred (1,1) lsqvce.Results
            opt.Alpha (1,1) double {mustBePositive, mustBeLessThan(opt.Alpha, 1)} = 0.05
            opt.Verbose (1,1) logical = false
        end
            sF = solOf(obj);  sR = solOf(Rred);
            F = lsqvce.engine.ftest(sF, sR, 'Alpha', opt.Alpha, 'Verbose', opt.Verbose);

            function s = solOf(Rk)
                if isempty(Rk.sol) || ~isfield(Rk.sol, 'vtv') || ~isfinite(Rk.sol.vtv)
                    error('lsqvce:Results:noVPv', ...
                          'compare() needs a finite v''Pv on both Results (this one is NaN or missing).');
                end
                u = numel(Rk.x);  n = Rk.diag.n;
                s = struct('vPv', Rk.sol.vtv, 'n', n, 'u', u, 'rank', u, 'r', n - u);
            end
        end

        function con = contribution(obj, varargin)
        %CONTRIBUTION  Per-group information share of this solution.
        %   con = R.contribution()                 Level 1: the scalar share
        %   con = R.contribution('Level', 2)       + diag(N^-1 N_g^cc) per parameter
        %   con = R.contribution(Name, Value, ...) see LSQVCE.ENGINE.CONTRIBUTION
        %
        %   WHY THIS EXISTS. The Level-1 share is already computed every outer
        %   pass and stored in R.history.contribution, because it is free (it
        %   reuses the redundancy numbers VCE returns). LEVEL 2 IS NOT FREE: it
        %   costs one extra streamed data sweep with Hutchinson probes, so it is
        %   deliberately NOT run automatically at convergence -- most runs do not
        %   want to pay for it. This method is how you ask for it, without having
        %   to reach into R.acc and R.sol yourself.
        %
        %   INPUT   varargin  Name/Value pairs passed through to
        %                     LSQVCE.ENGINE.CONTRIBUTION ('Probes',
        %                     'Seed', ...). 'GroupBy' defaults there.
        %   OUTPUT  con  1x1 struct  .group(k).key/.r/.share/.diag/.relSE, .nc
        %
        %   EXAMPLE
        %     R   = lsqvce.fit(data, Estimate="noise");
        %     con = R.contribution('Level', 2);
        %     bar(con.group(1).diag);      % per-coefficient share of group 1
        %
        %   See also LSQVCE.ENGINE.CONTRIBUTION, LSQVCE.ENGINE.STREAMCONTRIBUTION
        %
        %   Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude
        %   (Anthropic), model Claude Opus 4.8.
            obj.requireAccumulator('contribution');
            con = lsqvce.engine.contribution(obj.acc, obj.sol, varargin{:});
        end

        function requireAccumulator(obj, caller)
        %REQUIREACCUMULATOR  Refuse cleanly when this Results has no accumulator
        %   (the in-memory / batched kinds of LSQVCE.FIT, R4 WP-56b): the second
        %   data pass and the information shares need the streaming engine; the
        %   per-observation residuals of those kinds are already in .diag.v.
            if isempty(obj.acc)
                error('lsqvce:Results:noAccumulator', ...
                    ['%s() needs the streaming engine (R.acc); this Results came ' ...
                     'from the in-memory / batched kind -- read R.diag.v instead.'], caller);
            end
        end

        % ---- iteration diagnostics: two complementary views of R.history ----
        % convergence() reduces ACROSS STREAMS to one row per pass -- "when did
        % it settle". summary() reduces ACROSS PASSES to one row per stream --
        % "what did each stream end up as". Both read only R.history, which is
        % why they live here and not on the builder: a finished Results must be
        % readable without the object that produced it.

        function [T, varargout] = convergence(obj, opt)
        %CONVERGENCE  Tabulate (and optionally plot) why the loop stopped.
        %   T = R.convergence() returns a per-iteration table of the
        %   variance-factor and stochastic-model trajectory recorded in
        %   R.history, so you can SEE which stream never whitened and which
        %   variance component oscillated -- the two reasons the loop fails to
        %   converge.
        %
        %   T = R.convergence('Plot', true) also draws two stacked panels:
        %   max|sigma_g^2 - 1| per pass (the variance-factor criterion) and the
        %   per-stream Ljung-Box whiteness p-value per pass (the whiteness
        %   criterion), with the WhiteAlpha threshold marked. Returns [T, fh]
        %   with the figure handle.
        %
        %   NAME-VALUE
        %     'Plot'   1x1 logical  draw the panels (default false).
        %     'Stream' char/string  restrict the printed/plotted streams to those
        %                           whose key CONTAINS this text (default '': all).
        %
        %   OUTPUT
        %     T   table  one row per iteration, with columns:
        %                pass, maxSigmaOffset (max_g |sigma_g^2-1|), nWhite
        %                (streams declared white), minWhitePval (worst stream this
        %                pass), converged (logical, only the last row can be true).
        %
        %   ERRORS
        %     lsqvce:Results:noHistory   this Results carries no iteration
        %                                history (it came from a single solve).
        %
        %   EXAMPLE
        %     R = lsqvce.fit(data, Estimate="noise");
        %     T = R.convergence('Plot', true);   % table + diagnostic figure
        %     % a stream stuck at minWhitePval < 0.05 for every pass needs a
        %     % finer spectral grid (Options.Noise.Grid) or more passes; an
        %     % oscillating maxSigmaOffset means too little redundancy per
        %     % pool group -- coarsen Options.Pool ("obsType" or "all").
        %
        %   See also LSQVCE.RESULTS/summary, LSQVCE.ENGINE.ADJUSTVCEJOINT
        %
        %   Created 2026-07-19 (CEST, UTC+02:00) with assistance from Claude
        %   (Anthropic), model Claude Opus 4.8. Moved here from
        %   LSQVCE.ENGINE.ADJUSTVCEJOINT: it is a POST-RUN READER of R.history and needed
        %   the builder object only to draw one threshold line.
            % WP-39 c2: typed `arguments` signature (canonical Plot/Stream).
        arguments
            obj
            opt.Plot (1,1) logical = false
            opt.Stream (1,1) string = ""
        end
            o = struct('plot', opt.Plot, 'stream', char(opt.Stream));

            h = obj.requireHistory('convergence');
            nPass = size(h.sigma, 1);

            % per-pass variance-factor offset and whiteness summaries
            maxOff = max(abs(h.sigma - 1), [], 2);            % nPass x 1
            if isfield(h, 'white') && ~isempty(h.white)
                nWhite = sum(h.white ~= 0, 2);
            else
                nWhite = nan(nPass, 1);
            end
            if isfield(h, 'whitePval') && ~isempty(h.whitePval)
                minP = min(h.whitePval, [], 2, 'omitnan');
            else
                minP = nan(nPass, 1);
            end
            conv = false(nPass, 1);
            if isfield(h, 'converged'), conv(end) = logical(h.converged); end

            T = table((1:nPass).', maxOff(:), nWhite(:), minP(:), conv(:), ...
                'VariableNames', {'pass','maxSigmaOffset','nWhite', ...
                                  'minWhitePval','converged'});

            if o.plot
                fh = lsqvce.Results.plotConvergence(h, char(o.stream));
                if nargout > 1, varargout{1} = fh; end
            elseif nargout > 1
                varargout{1} = gobjects(0);
            end
        end

        function T = summary(obj)
        %SUMMARY  What each stream ended up as, one row per stream.
        %   T = R.summary() reads the FINAL iteration of R.history and reports,
        %   per stream, the variance factor it settled on, the temporal model it
        %   fitted, and whether its residuals came out white.
        %
        %   THIS IS THE COMPLEMENT OF R.CONVERGENCE, not a duplicate of it.
        %   convergence() gives one row per PASS, reduced across streams (when
        %   did the loop settle); summary() gives one row per STREAM, reduced
        %   across passes (what did each stream become). A loop that converges
        %   can still hide one stream that never whitened, and the per-pass view
        %   reports only the WORST p-value without naming its owner.
        %
        %   OUTPUT
        %     T   table  one row per stream, with columns:
        %                stream   string   'satKey|obsType'
        %                sigma    double   final variance factor (want ~1)
        %                arOrder  double   fitted AR order (0 = white)
        %                white    logical  residuals declared white
        %                pValue   double   final Ljung-Box p (NaN if not recorded)
        %
        %   ERRORS
        %     lsqvce:Results:noHistory   this Results carries no iteration
        %                                history (it came from a single solve).
        %
        %   EXAMPLE
        %     R = lsqvce.fit(data, Estimate="noise");
        %     T = R.summary();
        %     T(~T.white, :)            % the streams that never whitened
        %     T(abs(T.sigma-1) > 0.1, :)  % the streams still badly weighted
        %
        %   See also LSQVCE.RESULTS/convergence, LSQVCE.FIT
        %
        %   Created 2026-07-19 (CEST, UTC+02:00) with assistance from Claude
        %   (Anthropic), model Claude Opus 4.8.
            h = obj.requireHistory('summary');
            streams = string(h.streams(:));
            nS = numel(streams);

            sigma   = finalRow(h, 'sigma',     nS);
            arOrder = finalRow(h, 'arOrder',   nS);
            pValue  = finalRow(h, 'whitePval', nS);
            wRaw    = finalRow(h, 'white',     nS);
            white   = ~isnan(wRaw) & wRaw ~= 0;

            T = table(streams, sigma(:), arOrder(:), white(:), pValue(:), ...
                'VariableNames', {'stream','sigma','arOrder','white','pValue'});
        end

        function setIterative(obj, history, noise, converged, reason)
        %SETITERATIVE  Populate iteration fields (used by the iterative engines).
        %   stopReason is DERIVED (R5 M1): converged -> "converged"; not
        %   converged with a history -> "maxPasses"; no history -> "single" --
        %   unless the driver NAMES it (R5 M5 5c': "userStop" when the OnPass
        %   hook returned false).
            obj.history   = history;
            obj.noise     = noise;
            obj.converged = converged;
            if nargin > 4 && strlength(string(reason)) > 0
                obj.stopReason = string(reason);
            elseif converged
                obj.stopReason = "converged";
            elseif ~isempty(history)
                obj.stopReason = "maxPasses";
            else
                obj.stopReason = "single";
            end
        end

        function T = coefficients(obj, opt)
        %COEFFICIENTS  Estimates with formal errors and confidence bounds, as a table.
        %   T = R.coefficients()            95 % bounds
        %   T = R.coefficients(Alpha=0.01)  99 % bounds
        %
        %   The fitlm pair Coefficients (table) + coefCI (matrix) folded into ONE
        %   method (rule: one way to do each thing). The bounds are
        %   x +- z * se with se = sqrt(diag(Cxx)) -- Cxx already carries the
        %   a-posteriori sigma0^2 -- and z = norminv(1 - Alpha/2), computed as
        %   sqrt(2)*erfcinv(Alpha) so the Statistics Toolbox is not needed.
        %   Normal, not Student t: the redundancy of a streamed adjustment is in
        %   the thousands, where the two agree to < 1e-3 in z.
        %
        %   INPUT
        %     Alpha   1x1 double  two-sided significance level, 0 < Alpha < 1
        %                         (default 0.05)
        %   OUTPUT
        %     T   table, nc rows  index (double) | estimate | se | lower | upper
        %
        %   ACCURACY
        %     z(0.05) = 1.959963984540054 (closed form, 16 digits); the table is
        %     arithmetic on x and diag(Cxx), no further approximation.
        %
        %   ERRORS
        %     lsqvce:Results:noCovariance   sol carries no Cxx
        arguments
            obj
            opt.Alpha (1,1) double {mustBePositive, mustBeLessThan(opt.Alpha, 1)} = 0.05
        end
            if ~isstruct(obj.sol) || ~isfield(obj.sol, 'Cxx') || isempty(obj.sol.Cxx)
                error('lsqvce:Results:noCovariance', ...
                    'This Results carries no covariance (sol.Cxx); coefficients() needs one.');
            end
            se = sqrt(max(diag(obj.sol.Cxx), 0));
            se = se(:);
            z  = sqrt(2) * erfcinv(opt.Alpha);         % norminv(1 - Alpha/2)
            x  = obj.x(:);
            T  = table((1:numel(x)).', x, se, x - z*se, x + z*se, ...
                       'VariableNames', {'index', 'estimate', 'se', 'lower', 'upper'});
        end

        function disp(obj)
        %DISP  One-line custom display of the Results object.
            u = numel(obj.x);
            lsqvce.util.log.msg('  lsqvce.Results\n');
            lsqvce.util.log.msg('    common parameters : %d\n', u);
            lsqvce.util.log.msg('    systems           : %d\n', numel(obj.systems));
            lsqvce.util.log.msg('    sigma0            : %.4f%s\n', obj.sigma0, ...
                    tern(abs(obj.sigma0-1) < 0.05, '   (~1, good)', ''));
            r = obj.redundancy();
            if ~isnan(r), lsqvce.util.log.msg('    redundancy        : %.1f\n', r); end
            lsqvce.util.log.msg('    converged         : %s\n', tern(obj.converged, 'yes', 'no'));
            lsqvce.util.log.msg('    -> R.report()  for the level-0..4 verdict; R.coefficients() for x +- bounds; R.acc for the accumulator\n');
        end
    end

    methods (Access = private)
        function h = requireHistory(obj, caller)
        %REQUIREHISTORY  Return R.history, or explain clearly why there is none.
        %   A Results from a SINGLE solve legitimately has no iteration history,
        %   so this is a normal condition to hit, not a corruption -- the message
        %   says which engine populates it rather than reporting a missing field.
        %
        %   INPUT
        %     caller  1xk  char    the method name, for the message
        %   OUTPUT
        %     h       1x1  struct  the iteration history
        %   ERRORS
        %     lsqvce:Results:noHistory
            h = obj.history;
            if isempty(h) || ~isfield(h, 'sigma') || isempty(h.sigma)
                error('lsqvce:Results:noHistory', ...
                    ['R.%s needs an iteration history, and this Results has ' ...
                     'none (a single solve does not iterate). Use the RESULTS ' ...
                     'from a run with Estimate = "noise".'], ...
                     caller);
            end
        end
    end

    methods (Access = private, Static)
        function fh = plotConvergence(h, streamFilter)
        %PLOTCONVERGENCE  Two stacked panels: the variance-factor criterion
        %   (max|sigma^2-1| per pass) and the whiteness criterion (per-stream
        %   Ljung-Box p-value per pass, with the WhiteAlpha threshold).
        %
        %   INPUT
        %     h             1x1  struct  the iteration history
        %     streamFilter  1xk  char    show only streams whose key contains
        %                                this ('' = all; no match falls back to
        %                                all rather than drawing an empty axes)
        %   OUTPUT
        %     fh            1x1  figure handle
            nPass = size(h.sigma, 1);
            passes = 1:nPass;

            % choose streams to show
            streams = h.streams;
            keep = true(1, numel(streams));
            if ~isempty(streamFilter)
                keep = contains(streams, streamFilter);
                if ~any(keep), keep = true(1, numel(streams)); end
            end

            fh = figure('Name', 'lsqvce convergence');
            ax1 = subplot(2, 1, 1, 'Parent', fh);
            semilogy(ax1, passes, max(abs(h.sigma - 1), [], 2), '-o', ...
                     'LineWidth', 1.2);
            grid(ax1, 'on');
            xlabel(ax1, 'iteration');
            ylabel(ax1, 'max|\sigma_g^2 - 1|');
            title(ax1, 'Variance components (the loop stops when this is small)');

            ax2 = subplot(2, 1, 2, 'Parent', fh);
            hold(ax2, 'on');
            if isfield(h, 'whitePval') && ~isempty(h.whitePval)
                P = h.whitePval(:, keep);
                lbl = streams(keep);
                for k = 1:size(P, 2)
                    plot(ax2, passes, P(:, k), '-o', 'LineWidth', 1.2, ...
                         'DisplayName', lbl{k});
                end
                % R7-34: the threshold is the joint estimator's fixed Ljung-Box
                % level (adjustVCEjoint>whitenessRow, alpha 0.05). The former
                % reader of hist.whiteAlpha drew nothing: no producer wrote it
                % since the separable driver was retired.
                yline(ax2, 0.05, '--', 'alpha = 0.05', 'Color', [0.2 0.5 0.2]);
                legend(ax2, 'Interpreter', 'none', 'Location', 'best');
            else
                text(0.5, 0.5, 'no whitePval recorded', 'Parent', ax2, ...
                     'HorizontalAlignment', 'center', 'Units', 'normalized');
            end
            grid(ax2, 'on');  hold(ax2, 'off');
            xlabel(ax2, 'iteration');
            ylabel(ax2, 'Ljung-Box p-value');
            title(ax2, 'Stochastic model (streams whiten when p > WhiteAlpha)');
        end
    end
end

function v = finalRow(h, field, nS)
%FINALROW  The LAST row of a per-pass history field, as a 1xnS double.
%   Returns NaN(1,nS) when the field is absent, empty, or the wrong width, so a
%   partially-recorded history yields a readable table instead of an error --
%   summary() is a diagnostic, and refusing to print because one optional field
%   is missing would defeat it.
%
%   INPUT
%     h      1x1  struct  the iteration history
%     field  1xk  char    field name (nPass x nStream)
%     nS     1x1  double  expected stream count
%   OUTPUT
%     v      1xnS double  the final pass's values (NaN where unavailable)
    v = nan(1, nS);
    if ~isfield(h, field) || isempty(h.(field)), return; end
    M = h.(field);
    if size(M, 2) ~= nS, return; end
    v = double(M(end, :));
end

function s = tern(c, a, b)
%TERN  Inline conditional: a if c else b.
    if c, s = a; else, s = b; end
end
