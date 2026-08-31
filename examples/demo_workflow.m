%% DEMO_WORKFLOW  Build your OWN least-squares case, from scratch.
%
%   This is the from-scratch recipe. demo_lsqvce.m shows six FAILURE MODES on a
%   ready-made problem; this script shows how to WIRE YOUR OWN DATA into the
%   `data` struct, let lsqvce.fit pick the right engine from its SHAPE, choose
%   the model with lsqvce.Options, and read what comes back.
%
%   THE USER PATH is one verb on a plain data struct (DATA_CONTRACT.md is the
%   binding spec; examples/demo_fit.m shows it on every data kind):
%       R = lsqvce.fit(data, lsqvce.Options(...));   R.report
%
%   THE MAP (five questions, in order)
%     0.  What does lsqvce need from me?          -> a data struct (Section 0)
%     1.  How do I turn my data into that?        -> Section 1
%     2.  Which KIND / MODEL?                     -> Section 2 (decision tree)
%     3.  Common vs local parameters?             -> Section 3
%     4.  Wrong weights? Estimate them.           -> Section 4
%     5.  How do I read the result?               -> Section 5
%
%   Run it section by section (Ctrl+Enter in the Editor).
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-17 (CEST, UTC+02:00).  Rewritten from the old "how to read an
%   adjustment" walk-through into a build-your-own-case guide (HANDOVER item 1b).
%   Updated 2026-08-31 (CEST, UTC+02:00), v3 data contract (Level L10,
%   WP-74): USER surface only -- every run is lsqvce.fit on a data struct;
%   the engine names remain as the MAP of what fit dispatches to. With
%   assistance from Claude (Anthropic), model Claude Fable 5.

%% ------------------------------------------------------------------------
% SECTION 0.  What lsqvce needs from you: a data struct.
% -------------------------------------------------------------------------
% Three KINDS, told apart by the fields present (v3 contract: you state
% structure, never counts -- everything countable is derived):
%
%   (B) IN-MEMORY      data.A (m x u), data.l (m x 1, o - c), data.p (weights)
%   (C) BATCHED        data.batchFcn, data.numBatches, data.numUnknowns
%                      (linear form; with data.x0, u = numel(x0) and the
%                      callback returns the Jacobian at x)
%   (A) STREAMING      data.commonParams, data.days, data.systems(k) with
%                      your buildDay and covariance -- days, arcs, per-system
%                      locals, coloured noise (the GRACE-FO shape)
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

data = struct('batchFcn', batchFcn, 'numUnknowns', u, 'numBatches', nb);

%% ------------------------------------------------------------------------
% SECTION 2.  Which KIND, which MODEL?  A decision tree.
% -------------------------------------------------------------------------
% fit picks the ENGINE from the SHAPE of `data`; you pick the MODEL in Options:
%
%   kind B / C  (A,l,p | batchFcn)          -> lsqvce.engine.adjust
%       FLAT, one shared parameter block, batches streamed. Use when: every
%       observation sees the SAME u unknowns; u is modest; one covariance, no
%       per-system locals. Model="gls" (the only one these kinds take);
%       Robust="irls" switches to M-estimation. This section uses it.
%
%   kind A, Model="gls"                     -> lsqvce.engine.adjustBlock
%       BLOCK-ARROW (Helmert reduction). A shared COMMON block of nc unknowns,
%       plus PER-SYSTEM LOCAL blocks eliminated on the fly. One GLS solve.
%       Use when: several systems share a common field but each has its own
%       nuisance parameters (arc biases, scale factors), OR the full normal is
%       too big to form but the reduced common one is not. Section 3 shows it.
%
%   kind A, Model="vce"                     -> lsqvce.engine.adjustVCE
%       adjustBlock PLUS the iterative covariance loop: it re-estimates the
%       temporal correlation (the stochastic model) and the per-group variance
%       factor until the residuals whiten and the components agree.
%       Use when: you do NOT know the noise model a-priori and want the data to
%       supply it (separable covariance, coloured noise, unknown sigma ratios).
%
%   RULE OF THUMB
%       one parameter block, known noise ............ kind B/C, default Options
%       common + per-system locals, known noise ..... kind A,  Model="gls"
%       + you must DISCOVER the noise model ......... kind A,  Model="vce"
%
% Our Section-1 problem is flat with known-ish weights, so: the default.
R = lsqvce.fit(data);
D = R.diag;

%% ------------------------------------------------------------------------
% SECTION 2b.  Read the verdict (the part that matters).
% -------------------------------------------------------------------------
% Every Results carries the diagnostics on R.diag and the driver's verdict on
% R.source. Read them in this ORDER -- a failure high up makes everything
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
% R.source.verdict rolls these into one line; the numbers are how you act on it.
fprintf('\n--- verdict ---\n');
fprintf('  verdict      : %s\n', R.source.verdict);
fprintf('  integrity    : sum(r) rel.err = %.2e  (want ~0)\n', D.sumRRelErr);
fprintf('  variance     : s0^2 = %.3f            (want ~1)\n', D.varRatio);
fprintf('  blunders     : %d flagged\n', D.nFlag);

%% ------------------------------------------------------------------------
% SECTION 3.  Common vs local parameters (the block-arrow choice).
% -------------------------------------------------------------------------
% If different systems share SOME unknowns (the common field) but each keeps its
% OWN nuisance parameters, do NOT flatten them into one giant u. DECLARE it:
%
%       COMMON (nc)  : the shared field       -- one row of data.commonParams
%                      per unknown, estimated once, seen by all.
%       LOCAL  (nl)  : per-system nuisances   -- declared on the SYSTEM
%                      (localParams: you number them; or localsPerArc: the
%                      toolbox numbers per-arc blocks), eliminated on the fly,
%                      never forming the full (nc + sum nl) normal.
%
%       data.commonParams = commonTbl;             % nc = rows
%       data.days = 1:30;
%       s.buildDay      = @myBuildDay;             % the one function you write
%       s.covariance.rr = struct('type','arma', ...);   % fieldnames = the types
%       s.localParams   = num2cell((1:4).');       % nl = 4, this system's own
%       data.systems(1) = s;                       % ... (2) for the next system
%       R = lsqvce.fit(data);                      % Model="gls": one solve
%
% where buildDay(day, obsType, common, local, procOpt) returns that day's raw
% design and observations (A, B, l, t -- DATA_CONTRACT.md section 3).
%
% WHEN to reach for this instead of the flat kinds:
%   * you have MORE THAN ONE system with its own biases/scales, or
%   * the flat normal (nc + sum nl)^2 would not fit, but the reduced nc^2 will.
% The eliminated locals are handled EXACTLY (Helmert reduction), not dropped --
% and the redundancy accounting uses the FULL Qxx including them.
%
% (No runnable block here: it needs your multi-system buildDay. See
% examples/demo_combine.m for the declaration pattern, demo_sensors.m for a
% complete system/sensor/day hierarchy, and demo_grace.m for the full
% Model="vce" loop.)
fprintf('\n[Section 3 is a recipe -- see demo_combine / demo_sensors / demo_grace.]\n');

%% ------------------------------------------------------------------------
% SECTION 4.  Wrong weights? Estimate them: variance components.
% -------------------------------------------------------------------------
% When s0^2 is inflated but the residuals are CLEAN and show NO structure, the
% weights themselves are wrong. Declare one covariance field per observation
% GROUP and let Model="vce" estimate one variance component per group.
%
% Here: two instrument groups, the second secretly 3x noisier, but we CLAIM
% equal noise -- the normal state of affairs when instrument headers disagree
% with reality.
g      = [ones(n/2,1); 2*ones(n/2,1)];
sgTrue = sg;  sgTrue(g==2) = sgTrue(g==2) * 3;
l5     = A*xt + sgTrue .* randn(rs, n, 1);
t5     = (0:n-1).' * 5;                             % seconds (any monotone works)

vdata.commonParams = num2cell((1:u).');             % the same u unknowns
vdata.days = 1;                                     % the whole set is one "day"
vs.buildDay = @(day, ot) struct( ...                % rows of ONE group
        'A', A(g == 1 + strcmp(ot,'g2'), :), ...
        'B', zeros(sum(g == 1 + strcmp(ot,'g2')), 0), ...
        'l', l5(g == 1 + strcmp(ot,'g2')), ...
        't', t5(g == 1 + strcmp(ot,'g2')));
vs.covariance.g1 = struct('type', 'diag', 'sigma', sg(g==1));   % the CLAIMED
vs.covariance.g2 = struct('type', 'diag', 'sigma', sg(g==2));   % (mis-stated) noise
vdata.systems = vs;

Rv = lsqvce.fit(vdata, lsqvce.Options(Model="vce", GroupBy="obsType", ...
         Verbose=false, ComponentsPerEpoch=struct('g1',1,'g2',1)));

sigEnd = Rv.history.sigma(end, :);
ratio  = sigEnd(contains(Rv.history.streams, 'g2')) / ...
         sigEnd(contains(Rv.history.streams, 'g1'));
fprintf('\n[Model="vce"]  recovered sigma ratio: %.3f   (truth: 3.000)\n', ratio);
fprintf('               the loop rescaled the group covariances and re-solved: sigma0 = %.3f\n', ...
        Rv.sigma0);

% A variance component estimated from r_g observations has a relative standard
% error of about 0.5*sqrt(2/r_g). A "correction" smaller than that is NOISE:
% do not rescale your weights by 8% on the strength of a 2-sigma wobble.
% (At gravity-field sizes -- u in the thousands -- the exact per-group traces
% dominate; the ENGINE tier keeps a matrix-free Hutchinson backend for that
% regime: doc lsqvce.engine.varComp.)

%% ------------------------------------------------------------------------
% SECTION 5.  The decision tree, condensed.
% -------------------------------------------------------------------------
disp(' ');
disp('  THE DATA STRUCT (fit picks the engine from its shape)');
disp('    one parameter block, known noise ............ data.A/l/p or batchFcn');
disp('    common + per-system locals, known noise ..... data.systems, Model="gls"');
disp('    must DISCOVER the noise model ............... data.systems, Model="vce"');
disp(' ');
disp('  READING THE RESULT (in this order)');
disp('    integrity fails            -> fix your batchFcn/buildDay. Nothing else matters.');
disp('    ratio ~ 1                  -> done.');
disp('    ratio >> 1, contamination  -> R.diag.flagIdx, or Robust="irls" (kinds B/C).');
disp('    ratio >> 1, clean, STRUCTURE -> plot R.diag.w vs time etc. Fix the MODEL.');
disp('    ratio >> 1, clean, no structure -> Model="vce". Fix the WEIGHTS.');
disp('    ratio << 1                 -> over-parameterised?');
disp(' ');
