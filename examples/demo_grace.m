function R = demo_grace(varargin)
%DEMO_GRACE  End-to-end walk-through of the VCE loop on multi-axis and scalar streams.
%
%   R = DEMO_GRACE()
%   R = DEMO_GRACE('maxPasses', 8, 'verbose', true)
%
%   Runs a GRACE-FO-shaped synthetic problem through the USER surface:
%   the declaration comes from DEMO_GRACE_DATALOADER as a v5 data struct
%   (the API reference), the run is ONE call to LSQVCE.FIT, and everything
%   you read afterwards is on the LSQVCE.RESULTS. Two observation systems:
%     * System 1 -- satellites C and D, each with a 3-D position stream (posC,
%       posD) whose covariance is C = D (R_t kron I_3) D -- a FIXED per-epoch
%       3x3 D (as from a PPP solution, returned by buildDay with the day as
%       D.Cep) times a shared temporal law the toolbox estimates -- plus ONE
%       inter-satellite range-rate stream (rr, scalar) whose design couples
%       BOTH satellites' arc states. Because rr touches two satellites' local
%       parameters and the block-arrow model has one local block per system,
%       C and D are merged into a single pair-system with a stacked local
%       block [y_C | y_D].
%     * System 2 -- satellite E, INDEPENDENT: 3-D positions only (posE), no
%       inter-satellite link, its own local block y_E, observing the SAME
%       common gravity field. A clean separate block-arrow system.
%   The loop starts every stream from a deliberately WRONG a-priori
%   (temporally white, sigma 1) and lets Estimate="noise", Pool="stream",
%   ArcFactors=false discover the covariance: one estimated spectral shape and
%   one level per stream -- four levels, one each for posC, posD, rr, posE.
%
%   INPUT  (Name-Value, all optional)
%     maxPasses 1x1 double   passes of the loop (default 6)
%     verbose   1x1 logical  print the engine log (default true)
%   OUTPUT
%     R         1x1 lsqvce.Results   with R.history / R.noise attached
%
%   WHAT TO LOOK FOR
%     * sigma moves AWAY from 1 on the first pass (the misweighting is
%       detected) and settles NEAR the planted values by the last pass.
%     * the 'rr' stream's true noise IS white; the 'pos' streams carry
%       injected colour, which the loop recovers as the estimated spectral
%       shape (R.plot("psd") / R.noise.S) and reports in the white and
%       pValue columns of R.summary(). The joint estimator fits a spectrum,
%       not an AR order, so summary()'s arOrder column stays NaN here. Four
%       streams: posC, posD, rr, posE.
%
%   This is the copy-paste starting point. To run on YOUR data, replace
%   demo_grace_dataLoader (in this folder) with your own buildDay and data
%   struct -- the contract is the API reference (doc lsqvce).
%
%   ERRORS  none of its own; declaration errors carry lsqvce:fit ids.
%
%   See also LSQVCE.FIT, LSQVCE.OPTIONS, LSQVCE.RESULTS, DEMO_GRACE_DATALOADER
%
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-111 part 3: R.noise.S names
%   the estimated spectrum. With assistance from Claude (Anthropic), model
%   Claude Opus 5.5.
%   Updated 2026-09-07 (CEST, UTC+02:00), R6 WP-95: v5 vocabulary -- the
%   declaration is the v5 data struct (the API reference), the demo's own
%   knob is 'maxPasses' (it sets Options.MaxPasses, the ONE budget), and the
%   noise fit is described as one estimated spectral shape and one level per
%   stream instead of an AR order per stream and a component per obsType.
%   With assistance from Claude (Anthropic), model Claude Opus 5.
%   Updated 2026-08-31 (CEST, UTC+02:00), v4 SCOPED tables (Level L11,
%   WP-77): the declaration tables carry the scope column (commonParams
%   'common'|'day', localParams 'arc'|'day'|'system'); localsPerArc left the
%   USER surface. With assistance from Claude (Anthropic), model
%   Claude Fable 5.
%   Updated 2026-08-31 (CEST, UTC+02:00), v3 data contract (Level L10,
%   WP-74): USER surface only -- the run is lsqvce.fit over the v4 data
%   struct; the engine builder chain, the plot-hook wiring and the
%   empiricalStochastic walkthrough (ENGINE tier) are gone from the
%   example. With assistance from Claude (Anthropic), model Claude Fable 5.
%   Created 2026-07-16 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.
%
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-124: references to the data
%   contract point at the API reference. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.

    p = inputParser;
    p.addParameter('maxPasses', 6);
    p.addParameter('verbose',  true);
    p.parse(varargin{:});
    o = p.Results;

    [cfg, data] = demo_grace_dataLoader();   % the v5 declaration (this folder)

    fprintf(['\n=== VCE demo: GRACE-FO topology (+ independent sat E) ===\n' ...
             '    System 1: satellites C,D with a 3-D position stream each\n' ...
             '    (posC, posD), covariance C = D (R_t kron I_3) D with a FIXED\n' ...
             '    per-epoch 3x3 D, plus ONE inter-satellite range-rate stream\n' ...
             '    (rr) coupling both -- C,D merged into one pair-system.\n' ...
             '    System 2: satellite E, INDEPENDENT, 3-D positions only (posE),\n' ...
             '    no link -- its own local block y_E, same common gravity field.\n' ...
             '    VCE per stream -> FOUR variance components (posC,posD,rr,posE).\n\n']);
    fprintf('true per-stream sigma (posC, posD, rr, posE):\n');
    disp(cfg.trueSigma);
    fprintf('a-priori start: sigma = 1, white (deliberately wrong)\n');

    % THE OVERRIDES BELOW ARE DEMO-SPEED CHOICES, NOT RECOMMENDATIONS. The
    % toolbox defaults are MaxPasses 10, Tolerance 1e-2 (v5: one budget, one
    % tolerance for every iterative run); this demo runs a few passes and
    % settles visibly (Tolerance 1e-3). For your own case, start from the
    % defaults -- just omit them -- and change one only when a diagnostic
    % tells you to. Estimate="noise" with Pool="stream" and ArcFactors=false is
    % v4's Model "vce": one spectral shape and one level per stream.
    R = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, ...
            MaxPasses=o.maxPasses, Verbose=o.verbose, ...
            Tolerance=1e-3));

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

    fprintf('\nconverged: %d   iterations: %d\n', R.converged, R.history.nIter);
end
