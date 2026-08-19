classdef Results < handle
%LSQVCE.RESULTS  The outcome of a recovery: solution, diagnostics, provenance.
%   R = lsqvce.Results(sol, D, acc, diagSol, batchFcn, nBatch, systems)
%
%   Returned by LSQVCE.ADJUSTBLOCK (and, later, LSQVCE.ADJUSTVCE). It gathers
%   everything a user needs to read and judge an adjustment in one object, and
%   keeps the lower layers reachable for experts. You rarely construct it
%   yourself; you receive it from adjustBlock and call R.report().
%
%   PROPERTIES (read-only)
%     x          nc x 1   double   Common-parameter estimate (the field).
%     sigma0     1x1      double   A-posteriori variance factor (global).
%     converged  1x1      logical  Did the pipeline converge? (true for a single
%                                  solve; set by the iterative engines.)
%     diag       1x1      struct   Per-observation diagnostics (LSQVCE.DIAGNOSE):
%                                  redundancy, leverage, w-tests, MDB, ...
%     history    1xI      struct   Per-iteration record (empty for a single
%                                  solve; filled by LSQVCE.ADJUSTVCE).
%     noise      1x1      struct   Fitted noise models per stream (empty until
%                                  the separable/VCE loop populates it).
%     sol        1x1      struct   Raw output of acc.solve() (.x, .sigma0, ...).
%     acc        1x1      lsqvce.BlockAccumulator   The accumulator (drop here
%                                  to re-solve, inspect blocks, get residuals).
%     source     1x1      struct   Bridged solution from LSQVCE.DIAGNOSEACCUM.
%     batchFcn   1x1      function_handle   Re-reader used by LSQVCE.DIAGNOSE.
%     nBatch     1x1      double   Batch count from the bridge.
%     systems    1xM      struct   The resolved system specs (sysID/nl/dayList).
%
%   METHODS
%     V = R.report(Name,Value,...)   Print/return the level-0..4 verdict
%                                    (LSQVCE.REPORT); options pass straight through.
%     r = R.redundancy()             Total redundancy sum(r_i) = n - u.
%     res = R.residuals(Name,Value)  Whitened/raw residuals (acc.residuals).
%     con = R.contribution(...)      Per-group information share.
%     T = R.convergence(Name,Value)  Per-ITERATION trajectory: when the loop
%                                    settled (needs .history).
%     T = R.summary()                Per-STREAM final state: sigma, AR order,
%                                    whiteness (needs .history).
%     disp(R)                        One-screen summary.
%
%   THE TWO ITERATION READERS ARE COMPLEMENTARY. convergence() reduces across
%   STREAMS (one row per pass); summary() reduces across PASSES (one row per
%   stream). A converged loop can still hide a single stream that never
%   whitened -- the per-pass view reports the worst p-value without naming it.
%
%   EXAMPLE
%     R = lsqvce.adjustBlock(systems, days);
%     R.report('Level', 2);
%     if ~R.converged, lsqvce.util.log.warn('did not converge'); end
%     x = R.x;   s0 = R.sigma0;
%
%   ERRORS
%     lsqvce:Results:noHistory      an iteration reader was called on a
%                                   single-solve Results (adjustBlock does not
%                                   iterate); use the result of adjustVCE
%
%   See also LSQVCE.ADJUSTBLOCK, LSQVCE.ADJUSTVCE, LSQVCE.REPORT,
%            LSQVCE.DIAGNOSE, LSQVCE.BLOCKACCUMULATOR
%
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
    end
    properties (SetAccess = private)
        converged = true;      % single solve is, by definition, done
        history   = struct([]);
        noise     = struct([]);
    end

    methods
        function obj = Results(sol, D, acc, diagSol, batchFcn, nBatch, systems)
        %RESULTS  Construct a Results wrapper around a solve + diagnostics (see header).
            narginchk(7, 7);
            obj.sol      = sol;
            obj.x        = sol.x(:);
            obj.sigma0   = sol.sigma0;
            obj.diag     = D;
            obj.acc      = acc;
            obj.source   = diagSol;
            obj.batchFcn = batchFcn;
            obj.nBatch   = nBatch;
            obj.systems  = systems;
        end

        function V = report(obj, varargin)
        %REPORT  Judge the adjustment in words (LSQVCE.REPORT), levels 0-4.
            V = lsqvce.report(obj.diag, varargin{:});
        end

        function r = redundancy(obj)
        %REDUNDANCY  Total redundancy sum(r_i) = n - u, from the diagnostics.
            if isfield(obj.diag, 'r'),        r = sum(obj.diag.r(:));
            elseif isfield(obj.diag, 'redundancy'), r = sum(obj.diag.redundancy(:));
            else, r = NaN;
            end
        end

        function res = residuals(obj, varargin)
        %RESIDUALS  Collected residuals from the accumulator (acc.residuals).
            res = obj.acc.residuals(obj.sol, varargin{:});
        end

        function con = contribution(obj, varargin)
        %CONTRIBUTION  Per-group information share of this solution.
        %   con = R.contribution()                 Level 1: the scalar share
        %   con = R.contribution('Level', 2)       + diag(N^-1 N_g^cc) per parameter
        %   con = R.contribution(Name, Value, ...) see LSQVCE.CONTRIBUTION
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
        %                     LSQVCE.CONTRIBUTION ('GroupBy', 'Level', 'Probes',
        %                     'Seed', ...). 'GroupBy' defaults there.
        %   OUTPUT  con  1x1 struct  .group(k).key/.r/.share/.diag/.relSE, .nc
        %
        %   EXAMPLE
        %     R   = lsqvce.adjustVCE(L, days).groupComponentsBy('obsType').run();
        %     con = R.contribution('GroupBy', 'obsType', 'Level', 2);
        %     bar(con.group(1).diag);      % per-coefficient share of group 1
        %
        %   See also LSQVCE.CONTRIBUTION, LSQVCE.STREAMCONTRIBUTION
        %
        %   Created 2026-07-18 (CEST, UTC+02:00) with assistance from Claude
        %   (Anthropic), model Claude Opus 4.8.
            con = lsqvce.contribution(obj.acc, obj.sol, varargin{:});
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
        %     R = lsqvce.adjustVCE(L, days).groupComponentsBy('obsType').run();
        %     T = R.convergence('Plot', true);   % table + diagnostic figure
        %     % a stream stuck at minWhitePval < WhiteAlpha for every pass needs a
        %     % higher stochasticModel MaxOrder or ARMA; an oscillating
        %     % maxSigmaOffset means too little group redundancy -- coarsen
        %     % groupComponentsBy.
        %
        %   See also LSQVCE.RESULTS/summary, LSQVCE.ADJUSTVCE/run
        %
        %   Created 2026-07-19 (CEST, UTC+02:00) with assistance from Claude
        %   (Anthropic), model Claude Opus 4.8. Moved here from
        %   LSQVCE.ADJUSTVCE: it is a POST-RUN READER of R.history and needed
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
        %     R = lsqvce.adjustVCE(L, days).groupComponentsBy('obsType').run();
        %     T = R.summary();
        %     T(~T.white, :)            % the streams that never whitened
        %     T(abs(T.sigma-1) > 0.1, :)  % the streams still badly weighted
        %
        %   See also LSQVCE.RESULTS/convergence, LSQVCE.ADJUSTVCE/run
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

        function setIterative(obj, history, noise, converged)
        %SETITERATIVE  Populate iteration fields (used by the separable engine).
            obj.history   = history;
            obj.noise     = noise;
            obj.converged = converged;
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
            lsqvce.util.log.msg('    -> R.report()  for the level-0..4 verdict; R.acc for the accumulator\n');
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
                     'from LSQVCE.ADJUSTVCE run() or empiricalStochastic().'], ...
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
                % The threshold travels WITH the history (adjustVCE stores
                % hist.whiteAlpha), so this reads nothing from the builder.
                % Older Results predating that field simply get no line.
                if isfield(h, 'whiteAlpha') && ~isempty(h.whiteAlpha)
                    alpha = h.whiteAlpha;
                    yline(ax2, alpha, '--', sprintf('WhiteAlpha = %.3g', alpha), ...
                          'Color', [0.2 0.5 0.2]);
                end
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
