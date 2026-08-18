function R = demo_grace(varargin)
%DEMO_GRACE  End-to-end walk-through of the separable VCE loop.
%
%   R = DEMO_SEPARABLE()
%   R = DEMO_SEPARABLE('maxOuter', 8, 'verbose', true)
%   R = DEMO_SEPARABLE('Plots', true, 'exportDir', '/tmp/figs')
%
%   Builds a GRACE-FO-shaped synthetic problem with TWO observation systems:
%     * System 1 -- satellites C and D, each with a 3-D position stream (posC,
%       posD) whose covariance is the SEPARABLE model C = D (R_t kron I_3) D --
%       a FIXED per-epoch 3x3 D (as from a PPP solution) times a shared temporal
%       AR law -- plus ONE inter-satellite range-rate stream (rr, scalar) whose
%       design couples BOTH satellites' arc states. Because rr touches two
%       satellites' local parameters and the block-arrow model has one local
%       block per system, C and D are merged into a single pair-system with a
%       stacked local block [y_C | y_D].
%     * System 2 -- satellite E, INDEPENDENT: 3-D positions only (posE), no
%       inter-satellite link, its own local block y_E, observing the SAME common
%       gravity field. A clean separate block-arrow system.
%   The loop starts every stream from a deliberately WRONG a-priori (temporally
%   white, sigma 1) and lets LSQVCE.ADJUSTVCE discover the covariance:
%
%     stochasticModel  auto-AR whitening per stream       (selectArAcv on descaled ACV)
%     varianceComponents  one per group            (lsqvce.vce, Helmert)
%
%   The parametrisation is passed via makeLoader's 'common'/'local' 14-column
%   tables (nc/nl derived) and piped to the day-builder -- see demo_grace_dataLoader.
%
%   It prints the sigma trajectory and the final verdict, and returns the
%   lsqvce.Results object (with R.history attached) so you can inspect it.
%
%   WHAT TO LOOK FOR
%     * sigma moves AWAY from 1 on the first pass (the misweighting is detected)
%       and settles NEAR 1 by the last pass (the model now matches the data).
%     * the 'rr' streams are declared white (AR order 0); 'pos' picks up a low
%       AR order (the injected colour). Four streams now: posC, posD, rr, posE.
%
%   This is the copy-paste starting point. To run on YOUR data, replace
%   demo_grace_dataLoader (in this folder) with your own lsqvce.makeLoader (see the
%   "USER-SPECIFIC ADAPTATION" block in `help lsqvce.adjustVCE`).
%
%   See also LSQVCE.ADJUSTVCE, LSQVCE.VCE, LSQVCE.ADJUSTBLOCK,
%            LSQVCE.MAKELOADER
%
%   Created 2026-07-16 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.

    p = inputParser;
    p.addParameter('maxOuter', 6);
    p.addParameter('verbose',  true);
    p.addParameter('Plots',    false);   % attach example per-obsType plot hooks
    p.addParameter('exportDir', '');     % where the hook PNGs go ('' = temp)
    p.parse(varargin{:});
    o = p.Results;

    [cfg, L] = demo_grace_dataLoader();   % synthetic GRACE-FO data (this folder)

    % Optional: attach EXAMPLE per-(system,obsType) plot hooks to the loaders,
    % to demonstrate the mechanism. postSolve -> a residual histogram for the
    % stream; postModel -> the fitted AR order + VCE weight as a tiny bar. Each
    % hook just BUILDS and RETURNS a figure; the engine shows/export/closes it.
    if o.Plots
        expDir = o.ExportDir;
        if isempty(expDir), expDir = fullfile(tempdir, 'lsqvce_demo_figs'); end
        for i = 1:numel(L)
            % NOTE the Map. makeLoader broadcasts a BARE HANDLE to every obsType
            % ('plotPostSolve', @fcn) -- but that normalisation happens at
            % CONSTRUCTION. Here the loaders already exist, and the engine reads
            % this field expecting a containers.Map (adjustVCE.plotHooks skips
            % anything else SILENTLY), so assigning a bare handle after the fact
            % would disable the hooks with no error. Build the Map, or pass the
            % handle to makeLoader in the first place.
            L{i}.plotPostSolve = broadcastHooks(L{i}.obsTypes, @demoPostSolvePlot);
            L{i}.plotPostModel = broadcastHooks(L{i}.obsTypes, @demoPostModelPlot);
        end
        fprintf('plot hooks ON: figures export to %s\n', expDir);
    end

    fprintf(['\n=== separable VCE demo: GRACE-FO topology (+ independent sat E) ===\n' ...
             '    System 1: satellites C,D with a 3-D position stream each\n' ...
             '    (posC, posD), covariance C = D (R_t kron I_3) D with a FIXED\n' ...
             '    per-epoch 3x3 D, plus ONE inter-satellite range-rate stream\n' ...
             '    (rr) coupling both -- C,D merged into one pair-system.\n' ...
             '    System 2: satellite E, INDEPENDENT, 3-D positions only (posE),\n' ...
             '    no link -- its own local block y_E, same common gravity field.\n' ...
             '    VCE by obsType -> FOUR variance components (posC,posD,rr,posE).\n\n']);
    fprintf('true per-stream sigma (posC, posD, rr, posE):\n');
    disp(cfg.trueSigma);
    fprintf('a-priori start: sigma = 1, white (deliberately wrong)\n');

    % Two systems (pair + E) -> group the variance components by obsType.
    % THE THREE OVERRIDES BELOW ARE DEMO-SPEED CHOICES, NOT RECOMMENDATIONS.
    % The toolbox defaults are MaxOrder 60, MaxIter 10, Tol 5e-3; this demo
    % tightens them so a synthetic run finishes quickly (MaxOrder 40, MaxIter 8)
    % and settles visibly (Tol 1e-3). For your own case, start from the defaults
    % -- just omit the two configurators -- and change one only when a
    % diagnostic tells you to (see `help lsqvce.adjustVCE`).
    P = lsqvce.adjustVCE(L, cfg.days, ...
            'maxOuter', o.MaxOuter, 'verbose', o.Verbose) ...
          .groupComponentsBy('obsType') ...
          .stochasticModel('MaxOrder', 40) ...
          .varianceComponents('MaxIter', 8, 'Tol', 1e-3);
    if o.Plots
        expDir = o.ExportDir;
        if isempty(expDir), expDir = fullfile(tempdir, 'lsqvce_demo_figs'); end
        P = P.plots('showSeconds', 1, 'exportDir', expDir);   % show 1s, export PNG
    end

    % Run empirical stochastic modelling end to end: the step-0 homogenisation
    % preflight, the SOLVE<->stochastic model <-> variance components loop, then the formal-vs-empirical
    % acceptance test (attached at R.history.empiricalStochastic).
    R = P.empiricalStochastic('report', o.Verbose);

    fprintf('\n--- sigma trajectory (rows = iterations, cols = streams) ---\n');
    fprintf('streams: %s\n', strjoin(R.history.streams, '  '));
    disp(R.history.sigma);

    % The two iteration readers, both on the Results. summary() is PER STREAM
    % (what each ended up as); convergence() is PER ITERATION (when it settled).
    fprintf('\n--- per-stream final state (R.summary) ---\n');
    disp(R.summary());

    fprintf('\n--- per-iteration trajectory (R.convergence) ---\n');
    disp(R.convergence());

    fprintf('\n--- verdict ---\n');
    R.report();

    fprintf('\nconverged: %d   iterations: %d\n', R.converged, R.history.nOuter);

    % ---------------------------------------- ABOUT THE HUTCHINSON ALTERNATIVE
    % The variance-component step of this loop runs Helmert VCE through lsqvce.vce on the
    % block-arrow accumulator. The matrix-free Hutchinson VCE (the large-u
    % method, O(u^3 + G*Z*u^2) vs O(G*u^3)) is wired as a REAL substitution on
    % the FLAT path only -- lsqvce.vceHutchinson, or lsqvce.varComp(...,
    % 'backend','hutchinson'); see examples/demo_lsqvce.m and demo_workflow.m
    % for a side-by-side run. Bringing it into the block-arrow variance-component step needs
    % per-group N_g^cc exposed from the accumulator core (roadmap item 0, still
    % open). For a real high-degree gravity field this is exactly where it would
    % pay off, and it is the next core change on the roadmap.
    fprintf(['\nnote: at gravity-field scale the variance-component step would use the ' ...
             'matrix-free\n      Hutchinson VCE; see demo_lsqvce / ' ...
             'demo_workflow for the flat-path comparison.\n']);
end

%% ===================================================== example plot hooks
function m = broadcastHooks(obsTypes, fcn)
%BROADCASTHOOKS  A containers.Map obsType -> fcn for every obsType of a loader.
    m = containers.Map('KeyType','char','ValueType','any');
    for k = 1:numel(obsTypes), m(obsTypes{k}) = fcn; end
end

function fh = demoPostSolvePlot(ctx)
%DEMOPOSTSOLVEPLOT  Example postSolve hook: histogram of this stream's residuals.
%   ctx carries the estimate (ctx.R) and identifies the stream (ctx.stream).
%   Returns a figure handle; the engine shows/exports/closes it.
    res = ctx.R.residuals('Whitened', true);
    try
        v = lsqvce.collectResiduals(res, ctx.sysID, ctx.obsType, 'v');
    catch
        v = [];
    end
    fh = figure('Visible','off');
    if isempty(v)
        title(sprintf('%s  (no residuals)', ctx.stream), 'Interpreter','none');
    else
        histogram(v, 40);
        title(sprintf('%s  postSolve residuals (pass %d)', ctx.stream, ctx.pass), ...
              'Interpreter','none');
        xlabel('whitened residual');  ylabel('count');  grid on;
    end
end

function fh = demoPostModelPlot(ctx)
%DEMOPOSTMODELPLOT  Example postModel hook: the fitted AR order and VCE weight.
%   ctx.model is the fitted noise model; ctx.sigmaG is the stream's VCE sigma.
    arOrd = 0;
    if isfield(ctx.model, 'ar') && ~isscalar(ctx.model.ar)
        arOrd = numel(ctx.model.ar) - 1;
    end
    fh = figure('Visible','off');
    bar([arOrd, ctx.sigmaG]);
    set(gca, 'XTickLabel', {'AR order', 'VCE \sigma_g'});
    title(sprintf('%s  postModel summary', ctx.stream), 'Interpreter','none');
    grid on;
end
