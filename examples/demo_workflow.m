%% DEMO_WORKFLOW  Build your OWN least-squares case, from scratch.
%
%   This is the from-scratch recipe. demo_lsqvce.m shows six FAILURE MODES on a
%   ready-made problem; this script shows how to WIRE YOUR OWN DATA in, choose
%   the right engine for its size and structure, choose exact vs matrix-free
%   variance-component estimation, and read what comes back.
%
%   THE MAP (five questions, in order)
%     0.  What does lsqvce need from me?          -> a batchFcn (Section 0)
%     1.  How do I turn my data into that?        -> Section 1
%     2.  Which ENGINE do I call?                 -> Section 2 (decision tree)
%     3.  Common vs local parameters?             -> Section 3
%     4.  Exact vs Hutchinson variance estimate?  -> Section 4
%     5.  How do I read the result?               -> Section 5
%
%   Run it section by section (Ctrl+Enter in the Editor). Every name here is the
%   CURRENT one (adjust / adjustBlock / adjustVCE; varComp / vceHutchinson;
%   makeLoader long-form options).
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-17 (CEST, UTC+02:00).  Rewritten from the old "how to read an
%   adjustment" walk-through into a build-your-own-case guide (HANDOVER item 1b).

%% ------------------------------------------------------------------------
% SECTION 0.  What lsqvce needs from you: a batchFcn.
% -------------------------------------------------------------------------
% The ONLY hard requirement of the in-memory front door lsqvce.engine.adjust is a
% function
%
%       [A, l, p] = batchFcn(k)          k = 1 .. nBatch
%
% returning the k-th block of:
%       A   nk x u   design (Jacobian) rows for this block,
%       l   nk x 1   observed-minus-computed vector (the "o - c"),
%       p   nk x 1   WEIGHTS = diag(P) = 1 ./ sigma_i^2   (or [] for unit).
%
% batchFcn is called TWICE per batch (accumulate, then diagnose) and must
% return identical data both times. Nothing of size n is ever handed to the
% toolbox at once -- that is the whole point of streaming in batches. If a batch
% is expensive to read from disk, cache it inside the closure; if it is cheap to
% regenerate, don't bother.
%
% That is the entire contract. Everything below is about producing that handle
% from YOUR data and picking the engine that consumes it.

clear; clc; close all;
rs = RandStream('mt19937ar', 'Seed', 2026);

%% ------------------------------------------------------------------------
% SECTION 1.  Turn your data into a batchFcn.
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
                     l(edges(k)+1:edges(k+1)),    ...   % l_k : your o - c
                     p(edges(k)+1:edges(k+1)));         % p_k : your 1/sigma^2

fprintf('Your batchFcn: %d observations in %d batches, %d parameters.\n', ...
        n, nb, u);

% SANITY CHECK your own batchFcn before trusting any result: the sizes must be
% consistent and the weights positive. This one line catches most wiring bugs.
[A1, l1, p1] = batchFcn(1);
assert(size(A1,2) == u && numel(l1) == size(A1,1) && all(p1 > 0), ...
       'batchFcn contract violated: check A/l/p sizes and that p is 1/sigma^2.');
fprintf('batchFcn contract OK (A is nk x %d, l is nk x 1, p > 0).\n', u);

%% ------------------------------------------------------------------------
% SECTION 2.  Which ENGINE?  A decision tree.
% -------------------------------------------------------------------------
% lsqvce has THREE front doors. Pick by the SHAPE of your problem, not by taste:
%
%   lsqvce.engine.adjust(u, batchFcn, nBatch)
%       FLAT, in-memory normal/QR. One shared parameter block, batches streamed.
%       Use when: every observation sees the SAME u unknowns; u is modest
%       (hundreds to a few thousand); one covariance, no per-system locals.
%       This section uses it.
%
%   lsqvce.engine.adjustBlock(systems, days)
%       BLOCK-ARROW (Helmert reduction). A shared COMMON block of nc unknowns,
%       plus PER-SYSTEM LOCAL blocks eliminated on the fly. One GLS solve.
%       Use when: several systems share a common field but each has its own
%       nuisance parameters (arc biases, scale factors), OR the full normal is
%       too big to form but the reduced common one is not. Section 3 shows it.
%
%   lsqvce.engine.adjustVCE(systems, days, 'GroupBy', ..., 'Shape', struct(...), 'Scale', struct(...))
%       adjustBlock PLUS the iterative covariance loop: it re-estimates the
%       temporal correlation (the stochastic model) and the per-group variance factor until
%       the residuals whiten and the components agree.
%       Use when: you do NOT know the noise model a-priori and want the data to
%       supply it (separable covariance, coloured noise, unknown sigma ratios).
%
%   RULE OF THUMB
%       one parameter block, known noise ............ adjust
%       common + per-system locals, known noise ..... adjustBlock
%       + you must DISCOVER the noise model ......... adjustVCE
%
% Our Section-1 problem is flat with known-ish weights, so: adjust.
[sol, D, V] = lsqvce.engine.adjust(u, batchFcn, nb);

%% ------------------------------------------------------------------------
% SECTION 2b.  Read the verdict (the part that matters).
% -------------------------------------------------------------------------
% Every engine returns a diagnostics struct D and a verdict V. Read them in
% this ORDER -- a failure high up makes everything below meaningless.
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
% V.verdict rolls these into one line; the numbers are how you act on it.
fprintf('\n--- verdict ---\n');
fprintf('  verdict      : %s\n', V.verdict);
fprintf('  integrity    : sum(r) rel.err = %.2e  (want ~0)\n', D.sumRRelErr);
fprintf('  variance     : s0^2 = %.3f            (want ~1)\n', D.varRatio);
fprintf('  blunders     : %d flagged\n', D.nFlag);

%% ------------------------------------------------------------------------
% SECTION 3.  Common vs local parameters (the block-arrow choice).
% -------------------------------------------------------------------------
% If different systems share SOME unknowns (the common field) but each keeps its
% OWN nuisance parameters, do NOT flatten them into one giant u. Split them:
%
%       COMMON (nc)  : the shared field       -- estimated once, seen by all.
%       LOCAL  (nl)  : per-system nuisances   -- eliminated on the fly, never
%                      forming the full (nc + sum nl) normal.
%
% You describe each system with a Loader. The ergonomic way to build one is
% lsqvce.engine.makeLoader, whose options now carry CANONICAL long-form names:
%
%       L = lsqvce.engine.makeLoader(buildDay, ...
%               'observationTypes', {'pos','rr'}, ...   (alias: obsTypes)
%               'covariance',       covStruct,    ...   (alias: cov)
%               'systemId',         1,            ...   (alias: sysID)
%               'numLocal',         4,            ...   (alias: nl)
%               'commonParams',     commonRows,   ...   (alias: common)
%               'localParams',      localRows);         (alias: local)
%
% where buildDay(day, obsType, common, local, procOpt) returns that day's design
% and observations. Then hand the loaders to the block-arrow engine:
%
%       R = lsqvce.engine.adjustBlock({L1, L2}, days);     % one GLS solve
%       R.report();
%
% WHEN to reach for this instead of adjust:
%   * you have MORE THAN ONE system with its own biases/scales, or
%   * the flat normal (nc + sum nl)^2 would not fit, but the reduced nc^2 will.
% The eliminated locals are handled EXACTLY (Helmert reduction), not dropped --
% and the variance-component redundancy uses the FULL Qxx including them, which
% is the subtlety cmb.varComp / lsqvce.engine.vce exist to get right.
%
% (No runnable block here: it needs your multi-system loaders. See
% examples/demo_sensors.m for a complete common+local hierarchy, and
% examples/demo_grace.m for the full adjustVCE loop over makeLoader loaders.)
fprintf('\n[Section 3 is a recipe -- see demo_sensors / demo_grace to run it.]\n');

%% ------------------------------------------------------------------------
% SECTION 4.  Exact vs Hutchinson variance-component estimation.
% -------------------------------------------------------------------------
% When s0^2 is inflated but the residuals are CLEAN and show NO structure, the
% weights themselves are wrong. Estimate them from the data with Helmert VCE:
% one variance component per observation group.
%
% First, a group index per observation (here: two instrument groups, the second
% secretly 3x noisier, but we CLAIM equal weights -- the normal state of affairs
% when instrument headers disagree with reality).
g      = [ones(n/2,1); 2*ones(n/2,1)];
sgTrue = sg;  sgTrue(g==2) = sgTrue(g==2) * 3;
l5     = A*xt + sgTrue .* randn(rs, n, 1);
src5   = @(k) deal(A(edges(k)+1:edges(k+1), :), ...
                   l5(edges(k)+1:edges(k+1)),   ...
                   p(edges(k)+1:edges(k+1)));       % ...but CLAIM the old p
gsrc   = @(k) g(edges(k)+1:edges(k+1));             % group index per batch

% ---- THE EXACT PATH (default). Forms N^-1 N_g per group: O(G*u^3).
%      It also gives per-observation redundancy, snooping, and reliability.
VC = lsqvce.engine.varComp(u, src5, nb, gsrc, 2, 'verbose', false);
fprintf('\n[exact varComp]  recovered sigma ratio: %.3f   (truth: 3.000)\n', ...
        VC.sigma(2)/VC.sigma(1));
fprintf('                 significant? %s   resolution: %s %%\n', ...
        mat2str(VC.significant(:).'), ...
        mat2str(round(100*VC.resolution(:).', 1)));

% ---- THE HUTCHINSON PATH (the large-u method). Never forms N^-1: it estimates
%      each trace with Z random probes against the Cholesky factor, costing
%      O(u^3 + G*Z*u^2). Crossover is u ~ Z, so:
%
%          exact       below a few thousand unknowns   (faster AND gives r_i)
%          hutchinson  above that                       (matrix-free, r_g only)
%
%      Here u = 8, so exact is the right call; we run Hutchinson to show the SAME
%      answer comes back. Z = 4000 keeps the Monte-Carlo error (~1/sqrt(Z))
%      small enough that the agreement is tight. NOTE: no 'checkD' on this path
%      -- that guard is per-observation, and Hutchinson is per-group only.
VCh = lsqvce.engine.varComp(u, src5, nb, gsrc, 2, ...
                     'backend', 'hutchinson', 'probes', 4000, ...
                     'seed', 2026, 'verbose', false);
fprintf('[hutchinson, Z=%d] recovered sigma ratio: %.3f   (exact: %.3f)\n', ...
        VCh.probes, VCh.sigma(2)/VCh.sigma(1), VC.sigma(2)/VC.sigma(1));
fprintf('                 same answer, matrix-free -- reach for it when u is in the thousands.\n');

% A variance component estimated from r_g observations has a relative standard
% error of about 0.5*sqrt(2/r_g). A "correction" smaller than VC.resolution is
% NOISE: do not rescale your weights by 8% on the strength of a 2-sigma wobble.

%% ------------------------------------------------------------------------
% SECTION 5.  The decision tree, condensed.
% -------------------------------------------------------------------------
disp(' ');
disp('  ENGINE');
disp('    one parameter block, known noise ............ lsqvce.engine.adjust');
disp('    common + per-system locals, known noise ..... lsqvce.engine.adjustBlock');
disp('    must DISCOVER the noise model ............... lsqvce.engine.adjustVCE');
disp(' ');
disp('  READING THE RESULT (in this order)');
disp('    integrity fails            -> fix your batchFcn. Nothing else matters.');
disp('    ratio ~ 1                  -> done.');
disp('    ratio >> 1, contamination  -> snoop (test=''mad'') or robust().');
disp('    ratio >> 1, clean, STRUCTURE -> systematic(). Fix the MODEL.');
disp('    ratio >> 1, clean, no structure -> varComp(). Fix the WEIGHTS.');
disp('    ratio << 1                 -> over-parameterised? ftest().');
disp(' ');
disp('  VARIANCE COMPONENTS');
disp('    u below a few thousand -> exact varComp (also gives r_i, snooping).');
disp('    u in the thousands     -> varComp(...,''backend'',''hutchinson'') or');
disp('                              lsqvce.engine.vceHutchinson (matrix-free, r_g only).');
disp(' ');
