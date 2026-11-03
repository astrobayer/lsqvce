%% DEMO_WORKFLOW  Build your OWN least-squares case, from scratch.
%
%   This is the from-scratch recipe. demo_lsqvce.m shows six FAILURE MODES on a
%   ready-made problem; demo_tutorialMatrix.m, demo_tutorialBatch.m and
%   demo_tutorialStreaming.m are the Manual's three tutorials, line by line. This script shows how to WIRE YOUR OWN DATA into a
%   lsqvce.TaskDefinition, which constructor makes which KIND of task, how to
%   choose what to ESTIMATE with lsqvce.Options, and how to read what comes back.
%
%   THE USER PATH is one declaration and one verb (the API reference is the
%   binding spec; examples/demo_fit.m shows it on every kind of task):
%       task = lsqvce.TaskDefinition.fromBatches(...);   % or fromMatrices, or (Days=...)
%       R = lsqvce.fit(task, lsqvce.Options(...));       R.report
%
%   THE MAP (five questions, in order)
%     0.  What does lsqvce need from me?          -> a TaskDefinition (Section 0)
%     1.  How do I turn my data into that?        -> Section 1
%     2.  Which KIND, what to ESTIMATE?           -> Section 2 (decision tree)
%     3.  Common vs local parameters?             -> Section 3
%     4.  Wrong weights? Estimate them.           -> Section 4
%     5.  How do I read the result?               -> Section 5
%
%   Run it section by section (Ctrl+Enter in the Editor). A script: no input,
%   it leaves its variables (task, R, vtask, Rv, ...) in the workspace.
%
%   ACCURACY  Not pinned: a walk-through. Its numbers come from seeded draws
%             (mt19937ar, seed 2026) and are printed, not asserted; the gate
%             runs it as a smoke test (test_examples).
%   ERRORS    none of its own; a declaration error is raised by the add-call
%             that makes it (lsqvce:TaskDefinition:*), a mistake inside a Build
%             while the day is read (lsqvce:DayContext:*).
%
%   See also LSQVCE.TASKDEFINITION, LSQVCE.FIT, LSQVCE.OPTIONS, LSQVCE.RESULTS,
%   DEMO_TUTORIAL, DEMO_FIT, DEMO_LSQVCE
%
%   Updated 2026-11-03 (CET, UTC+01:00), WP-163: points at the three tutorials
%   (demo_tutorialMatrix / Batch / Streaming). With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-11-02 (CET, UTC+01:00), WP-161: declared as lsqvce.TaskDefinition, one Build per group. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-09-07 (CEST, UTC+02:00), R6 WP-95: the contract is named v5,
%   locals are the SCOPED localParams table, robust IRLS runs on every kind, and
%   Section 4 returns the per-observation accuracy with the day (D.Cep).
%   With assistance from Claude (Anthropic), model Claude Opus 5.
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-17 (CEST, UTC+02:00).  Rewritten from the old "how to read an
%   adjustment" walk-through into a build-your-own-case guide (HANDOVER item 1b).
%   Updated 2026-08-31 (CEST, UTC+02:00), v4 SCOPED tables (Level L11,
%   WP-77): the declaration tables carry the scope column (commonParams
%   'common'|'day', localParams 'arc'|'day'|'system'); localsPerArc left the
%   USER surface. With assistance from Claude (Anthropic), model
%   Claude Fable 5.
%   Updated 2026-08-31 (CEST, UTC+02:00), v3 data contract (Level L10,
%   WP-74): USER surface only -- every run is lsqvce.fit on a data struct;
%   the engine names remain as the MAP of what fit dispatches to. With
%   assistance from Claude (Anthropic), model Claude Fable 5.
%
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-124: references to the data
%   contract point at the API reference. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.

%% ------------------------------------------------------------------------
% SECTION 0.  What lsqvce needs from you: a lsqvce.TaskDefinition.
% -------------------------------------------------------------------------
% Three KINDS, one constructor each (you state structure, never counts --
% everything countable is derived):
%
%   (B) IN-MEMORY      lsqvce.TaskDefinition.fromMatrices(A, l, Weights=p)
%                      A (m x u), l (m x 1, o - c), p (m x 1 weights)
%   (C) BATCHED        lsqvce.TaskDefinition.fromBatches(@batchFcn, NumBatches=K, NumUnknowns=u)
%                      (linear form; with X0=x0, u = numel(x0) and the
%                      callback returns the Jacobian at x)
%   (A) STREAMING      lsqvce.TaskDefinition(Days=days), filled with add-calls:
%                      addCommon, addSystem (its Prepare), addLocal and
%                      addStream (its Build and Noise) -- days, arcs,
%                      per-system locals, coloured noise (the GRACE-FO shape)
%
% The batched callback contract (kind C) is
%
%       [A, l, p] = batchFcn(k)          k = 1 .. numBatches
%
%       A   nk x u   design (Jacobian) rows for this batch,
%       l   nk x 1   observed-minus-computed vector (the "o - c"),
%       p   nk x 1   WEIGHTS = diag(P) = 1 ./ sigma_i^2   (or [] for unit).
%
% batchFcn is called TWICE per batch (accumulate, then diagnose) and must
% return identical data both times. Nothing of size n is ever handed to the
% toolbox at once -- that is the whole point of streaming in batches. If a
% batch is expensive to read from disk, cache it inside the closure; if it is
% cheap to regenerate, don't bother.

clear; clc; close all;
rs = RandStream('mt19937ar', 'Seed', 2026);

%% ------------------------------------------------------------------------
% SECTION 1.  Turn your data into a batched declaration.
% -------------------------------------------------------------------------
% In real work each batch is one file / one arc / one day: you read it, form the
% partials, and return [A, l, p]. Here we synthesise a problem so the script is
% self-contained, then slice it -- the slicing is exactly the pattern you copy.

n  = 3000;      % observations (never materialised inside the toolbox)
u  = 8;         % parameters
nb = 30;        % batches

A  = randn(rs, n, u);
xt = randn(rs, u, 1);              % the TRUE parameters (known here only)
sg = 0.5 + rand(rs, n, 1);         % per-observation sigma
p  = 1 ./ sg.^2;                   % WEIGHT = 1/sigma^2  (NOT sigma!)
l  = A*xt + sg .* randn(rs, n, 1); % clean observations

% The slicing pattern: precompute batch edges, then a closure that returns the
% k-th slice. Replace the three deal() arguments with your own file reads.
edges = round(linspace(0, n, nb+1));
batchFcn = @(k) deal(A(edges(k)+1:edges(k+1), :), ...   % A_k : your partials
                     l(edges(k)+1:edges(k+1)),   ...    % l_k : your o - c
                     p(edges(k)+1:edges(k+1)));         % p_k : your weights

fprintf('Your batchFcn: %d observations in %d batches, %d parameters.\n', ...
        n, nb, u);

% SANITY CHECK your own batchFcn before trusting any result: the sizes must be
% consistent and the weights positive. This one line catches most wiring bugs.
[A1, l1, p1] = batchFcn(1);
assert(size(A1,2) == u && numel(l1) == size(A1,1) && all(p1 > 0), ...
       'batchFcn wiring is inconsistent');
fprintf('batchFcn contract OK (A is nk x %d, l is nk x 1, p > 0).\n', u);

task = lsqvce.TaskDefinition.fromBatches(batchFcn, NumBatches=nb, NumUnknowns=u);

%% ------------------------------------------------------------------------
% SECTION 2.  Which KIND, what to ESTIMATE?  A decision tree.
% -------------------------------------------------------------------------
% The constructor you called fixed the KIND of `task`; you say what to
% estimate in Options. Every kind runs on the ONE accumulator (block-arrow: a
% shared common block, per-system local blocks eliminated on the fly):
%
%   kind B / C  (fromMatrices | fromBatches)
%       In memory or batched, one shared parameter block: the batches become
%       the days of one local-free system. Use when: every observation sees
%       the SAME u unknowns and there are no per-system locals. Estimate="none"
%       (a single solve); Robust="huber" switches to M-estimation; X0=x0 makes
%       it Gauss-Newton; task.addPrior(name, L=L, Sigma=s) gives a free network
%       its datum. This section uses it.
%
%   kind A  (days, arcs, systems, named parameters with a scope)
%       Use when: several systems share a common field but each has its own
%       nuisance parameters (arc biases, scale factors), OR the full normal is
%       too big to form but the reduced common one is not. Section 3 shows it.
%
%   kind A + Estimate="noise"
%       The same solve PLUS the joint variance-component estimator: the
%       spectral SHAPE per pool group (Pool="all"|"obsType"|"stream"), one
%       LEVEL per stream, and per-arc factors (ArcFactors) -- passes until the
%       residuals whiten and the factors settle (MaxPasses, Tolerance). Use
%       when you do NOT know the noise model a-priori and want the data to
%       supply it. Watch it pass by pass with Options.OnPass; look at it
%       afterwards with R.plot("psd" | "convergence" | "residuals").
%
%   RULE OF THUMB
%       one parameter block, known noise ............ kind B/C, default Options
%       common + per-system locals, known noise ..... kind A,  Estimate="none"
%       + you must DISCOVER the noise model ......... kind A,  Estimate="noise", Pool="stream", ArcFactors=false
%
% Our Section-1 problem is flat with known-ish weights, so: the default.
R = lsqvce.fit(task);
D = R.diag;

%% ------------------------------------------------------------------------
% SECTION 2b.  Read the verdict (the part that matters).
% -------------------------------------------------------------------------
% Every Results carries the diagnostics on R.diag and the verdict through
% R.report(). Read them in this ORDER -- a failure high up makes everything
% below meaningless.
%
%   1. INTEGRITY   D.sumRRelErr ~ 0   the identity sum(r_i) = m - u holds.
%                  If not, your batchFcn is inconsistent between its two calls.
%   2. VARIANCE    D.varRatio (s0^2)  ~ 1 means the weights match the scatter.
%                  >> 1 : weights too optimistic, or blunders, or model error.
%                  << 1 : weights too pessimistic, or over-parameterised.
%   3. BLUNDERS    D.flagIdx / D.nFlag  which observations data-snooping caught.
%   4. RELIABILITY the smallest blunder each observation could HIDE (MDB) and
%                  how far it could shift a parameter (external reliability).
%
% R.report() rolls these into one line (V.verdict); the numbers are how you
% act on it. (R5 M2: R.source is the sweep's solution struct on every kind.)
V = R.report('Verbose', false);
fprintf('\n--- verdict ---\n');
fprintf('  verdict      : %s\n', V.verdict);
fprintf('  integrity    : sum(r) rel.err = %.2e  (want ~0)\n', D.sumRRelErr);
fprintf('  variance     : s0^2 = %.3f            (want ~1)\n', D.varRatio);
fprintf('  blunders     : %d flagged\n', D.nFlag);

%% ------------------------------------------------------------------------
% SECTION 3.  Common vs local parameters (the block-arrow choice).
% -------------------------------------------------------------------------
% If different systems share SOME unknowns (the common field) but each keeps its
% OWN nuisance parameters, do NOT flatten them into one giant u. DECLARE it:
%
%       COMMON (nc)  : the shared field       -- task.addCommon(names), one
%                      unique name per unknown, estimated once, seen by all.
%       LOCAL  (nl)  : per-system nuisances   -- sys.addLocal(names, Scope=...)
%                      on the SYSTEM, Scope "arc" | "day" | "system"; the
%                      names are TEMPLATES and nl is DERIVED from their
%                      expansion over the calendar. Eliminated on the fly,
%                      never forming the full (nc + sum nl) normal.
%
%       task = lsqvce.TaskDefinition(Days=1:30);
%       task.addCommon(fieldNames);                % one name per unknown, nc derived
%       sys = task.addSystem("sat1", Prepare=@myPrepare);   % once per system and day
%       sys.addLocal(["b1" "b2" "b3" "b4"], Scope="system");
%                                                  % 4 "system" names -> nl = 4,
%                                                  % this system's own ("arc"
%                                                  % gives one block per (day,
%                                                  % arc) and needs Arcs= on the
%                                                  % system, "day" one per day
%                                                  % with arcs)
%       sys.addStream("rr", Build=@myBuildRr, Noise=lsqvce.Noise.arma(ar, sigma));
%                                                  % the one function per stream
%       sys2 = task.addSystem("sat2", ...);        % ... the next system
%       R = lsqvce.fit(task);                      % Estimate="none": one solve
%
% where out = myBuildRr(day, ctx) returns that day's raw design and
% observations: out = ctx.newDay(); out.addArc(A, l, t, Local={names, P})
% for each arc, the local partials by name (LSQVCE.DAYCONTEXT, LSQVCE.DAYDATA).
%
% WHEN to reach for this instead of the flat kinds:
%   * you have MORE THAN ONE system with its own biases/scales, or
%   * the flat normal (nc + sum nl)^2 would not fit, but the reduced nc^2 will.
% The eliminated locals are handled EXACTLY (Helmert reduction), not dropped --
% and the redundancy accounting uses the FULL Qxx including them.
%
% (No runnable block here: it needs your multi-system Build functions. See
% examples/demo_tutorialStreaming.m for a complete runnable TaskDefinition with a
% Prepare and two streams, demo_combine.m for the declaration pattern,
% demo_sensors.m for a complete system/sensor/day hierarchy, and demo_grace.m
% for the full Estimate="noise", Pool="stream", ArcFactors=false loop.)
fprintf('\n[Section 3 is a recipe -- see demo_combine / demo_sensors / demo_grace.]\n');

%% ------------------------------------------------------------------------
% SECTION 4.  Wrong weights? Estimate them: variance components.
% -------------------------------------------------------------------------
% When s0^2 is inflated but the residuals are CLEAN and show NO structure, the
% weights themselves are wrong. Declare one STREAM per observation GROUP and
% let Estimate="noise", Pool="stream", ArcFactors=false estimate one variance
% component per group.
%
% Here: two instrument groups, the second secretly 3x noisier, but we CLAIM
% equal noise -- the normal state of affairs when instrument headers disagree
% with reality.
g      = [ones(n/2,1); 2*ones(n/2,1)];
sgTrue = sg;  sgTrue(g==2) = sgTrue(g==2) * 3;
l5     = A*xt + sgTrue .* randn(rs, n, 1);
t5     = (0:n-1).' * 5;                             % seconds (any monotone works)

vtask = lsqvce.TaskDefinition(Days=1);              % the whole set is one "day"
vtask.addCommon(compose("x%d", 1:u));               % the same u unknowns, one name each
vsys = vtask.addSystem("groups");                   % one system, no locals
% One stream per group, each with its own Build (groupArc, at the end of this
% file): the rows of ONE group as ONE arc, with the CLAIMED accuracy of every
% row, one sigma per row. That accuracy is DATA -- what the preprocessing knows
% about each observation -- so Build passes it WITH the arc, as Sigma=. The
% stream's lsqvce.Noise keeps the TEMPORAL law and its LEVEL only: sigma is 1
% because Sigma= already carries the accuracy, and the stream's
% Estimate="scale" asks for that one LEVEL per stream -- the factor group 2 is
% really off by (LSQVCE.NOISE, LSQVCE.DAYDATA).
vsys.addStream("g1", Build=@(day, ctx) groupArc(ctx, A, l5, t5, sg, g == 1), ...
               Noise=lsqvce.Noise.white(1), Estimate="scale");
vsys.addStream("g2", Build=@(day, ctx) groupArc(ctx, A, l5, t5, sg, g == 2), ...
               Noise=lsqvce.Noise.white(1), Estimate="scale");

Rv = lsqvce.fit(vtask, lsqvce.Options(Estimate="noise", Pool="stream", ArcFactors=false, ...
         Verbose=false));

sigEnd = Rv.history.sigma(end, :);
ratio  = sigEnd(contains(Rv.history.streams, 'g2')) / ...
         sigEnd(contains(Rv.history.streams, 'g1'));
fprintf('\n[Estimate="noise", Pool="stream", ArcFactors=false]  recovered sigma ratio: %.3f   (truth: 3.000)\n', ratio);
fprintf('               the loop rescaled the group covariances and re-solved: sigma0 = %.3f\n', ...
        Rv.sigma0);

% A variance component estimated from r_g observations has a relative standard
% error of about 0.5*sqrt(2/r_g). A "correction" smaller than that is NOISE:
% do not rescale your weights by 8% on the strength of a 2-sigma wobble.
% (At gravity-field sizes -- u in the thousands -- the per-lag redundancies
% are estimated with Hutchinson probes: Options.Probes, the default 100.)

%% ------------------------------------------------------------------------
% SECTION 5.  The decision tree, condensed.
% -------------------------------------------------------------------------
disp(' ');
disp('  THE TASK (its constructor fixes the kind fit runs)');
disp('    one parameter block, known noise ............ TaskDefinition.fromMatrices or .fromBatches');
disp('    common + per-system locals, known noise ..... TaskDefinition(Days=...), addSystem');
disp('    must DISCOVER the noise model ............... addSystem, Estimate="noise", Pool="stream", ArcFactors=false');
disp(' ');
disp('  READING THE RESULT (in this order)');
disp('    integrity fails            -> fix your batchFcn/Build. Nothing else matters.');
disp('    ratio ~ 1                  -> done.');
disp('    ratio >> 1, contamination  -> R.diag.flagIdx, or Robust="huber"');
disp('                                  (any kind, with Estimate="none"; not with X0).');
disp('    ratio >> 1, clean, STRUCTURE -> plot R.diag.w vs time etc. Fix the MODEL.');
disp('    ratio >> 1, clean, no structure -> Estimate="noise", Pool="stream", ArcFactors=false. Fix the WEIGHTS.');
disp('    ratio << 1                 -> over-parameterised?');
disp(' ');

%% ------------------------------------------------------------------ local
function out = groupArc(ctx, A, l, t, sg, inGroup)
%GROUPARC  The Build of Section 4: the rows of ONE group as ONE arc (no
%   locals), with the CLAIMED accuracy of each row as the arc's Sigma.
    out = ctx.newDay();
    out.addArc(A(inGroup, :), l(inGroup), t(inGroup), Sigma=sg(inGroup));
end
