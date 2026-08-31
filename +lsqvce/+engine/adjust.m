function [sol, D, V] = adjust(u, batchFcn, nBatch, varargin)
%LSQVCE.ENGINE.ADJUST  Solve, diagnose and judge a batch least-squares problem.
%
%   [sol, D, V] = LSQVCE.ENGINE.ADJUST(u, batchFcn, nBatch) runs the whole pipeline:
%   accumulate -> solve -> second-pass diagnostics -> verdict.
%
%   batchFcn is a handle with signature
%       [A, l, p] = batchFcn(k)      k = 1 .. nBatch
%   where A is nk x u, l is nk x 1, and p is the nk x 1 weight vector (the
%   diagonal of P), or [] for unit weights.
%
%   It is called TWICE per batch -- once to accumulate, once to diagnose --
%   and must return identical data both times. If reading from disk is
%   expensive, cache inside the closure; if the batches are cheap to
%   regenerate, don't bother.
%
%   INPUT
%     u         1x1  double  number of unknowns
%     batchFcn  1x1  fh      [A, l, p] = batchFcn(k), k = 1..nBatch
%     nBatch    1x1  double  number of batches
%
%   OUTPUT
%     sol  1x1 struct  LSQVCE.ENGINE.ACCUMULATOR solve() output, i.e.
%                        .x       u x 1  estimate
%                        .R       u x u  upper-triangular factor
%                        .n .u    counts; .defect datum defect; .rank u - defect
%                        .r       redundancy n - rank
%                        .vPv     weighted residual sum of squares
%                        .s0sq    vPv / r (the A-POSTERIORI variance factor)
%                        .sigma0  the A-PRIORI Sigma0 you passed (default 1)
%                        .QxxCorr u x d datum correction of the cofactor
%                      There is NO .Cxx here: the covariance s0sq * Qxx is
%                      formed by LSQVCE.ENGINE.DIAGNOSE (D.Qxx) -- and by LSQVCE.FIT's
%                      bridge onto a Results, where sigma0 becomes a posteriori
%                      (R4 WP-61: this line used to promise .Cxx and .sigma0 in
%                      the streaming sense; the bridge, not the help, was right).
%     D    1x1 struct  diagnostics: residuals, redundancy numbers, w/tau,
%                      leverage, flagIdx, MDB -- see LSQVCE.ENGINE.DIAGNOSE
%     V    1x1 struct  the level-0..4 verdict -- see LSQVCE.REPORT
%
%   NAME-VALUE  (WP-39: adjust's OWN five options carry canonical
%   UpperCamelCase names with typed validators; any casing still matches
%   (inputParser, case-insensitive -- same behaviour as the
%   arguments-migrated functions). Every OTHER pair is FORWARDED verbatim
%   to LSQVCE.ENGINE.DIAGNOSE -- the deliberate open pass-through that keeps
%   this front door on inputParser; forwarded names below use diagnose's
%   current lowercase surface until its own WP-39 stage.)
%     'Method'      "qr" (default) | "normal" -- enum VALUE is
%                   case-sensitive lowercase (WP-39)
%     'Sparse'      logical, sparse accumulation           (default false)
%     'DatumB'      datum constraint matrix                (default [])
%     'Sigma0'      a-priori sigma0, scalar > 0           (default 1)
%     'alpha'       significance per test                 (default 1e-3)
%     'bonferroni'  correct for n simultaneous tests      (default true)
%     'beta'        MDB type-II error                     (default 0.20)
%     'rWeak'       weak-control threshold                (default 0.10)
%     'test'        'auto' | 'pope' | 'baarda'            (default 'auto')
%     'store'       'all' | 'flagged' | 'none'            (default 'all')
%     'Verbose'     print the report                      (default true)
%
%   EXAMPLE
%       % Two blunders planted at observations 40 and 150.
%       rng(0);  n = 200;  u = 4;
%       A  = randn(n,u);   xt = randn(u,1);
%       sg = 0.5 + rand(n,1);   p = 1 ./ sg.^2;      % p = WEIGHT = 1/sigma^2
%       l  = A*xt + sg.*randn(n,1);
%       l([40 150]) = l([40 150]) + 12*sg([40 150]);
%
%       idx = reshape(1:n, [], 4);                   % 4 batches of 50
%       src = @(k) deal(A(idx(:,k),:), l(idx(:,k)), p(idx(:,k)));
%
%       [sol, D, V] = lsqvce.engine.adjust(u, src, 4);
%
%       % V.verdict      -> 'PASS WITH WARNINGS'
%       % D.flagIdx      -> [40 150]      both blunders found
%       % D.s0           -> 1.675         variance factor inflated by them
%       % D.sumRRelErr   -> 0             the integrity identity holds
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (a) round a2
%   (F-07): canonical UpperCamelCase names + typed validators on the five
%   own options; enum VALUE 'Method' case-sensitive lowercase; the open
%   pass-through to LSQVCE.ENGINE.DIAGNOSE keeps this front door on inputParser
%   BY DESIGN (pinned deviation, §9 WP-39). With assistance from Claude
%   (Anthropic), model Claude Fable 5.
%
%   ACCURACY
%     adjust adds no numerics of its own: sol is Accumulator.solve() and D
%     is diagnose() on the same batches. The estimator is certified by the
%     §6.1.4 authorities (closed forms, Baarda/Pope textbook values, the
%     shDesign oracle); through fit it reproduces byte-identically
%     (dev/tests/test_fit.m, kinds B and C).
%
%   See also LSQVCE.ENGINE.DIAGNOSE, LSQVCE.REPORT, LSQVCE.ENGINE.SNOOP, LSQVCE.FIT.
%
%   Updated 2026-09-01 (CEST, UTC+02:00), R4 WP-61 (Level L3): OUTPUT sol
%   described as it IS (no .Cxx; sigma0 a priori); ACCURACY block. With
%   assistance from Claude (Anthropic), model Claude Fable 5.
%
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.

% WP-39 stage (a) round a2 (F-07), PINNED DEVIATION: adjust KEEPS
% inputParser BY DESIGN -- KeepUnmatched is a documented FEATURE here
% (unmatched options are forwarded to LSQVCE.ENGINE.DIAGNOSE below), which an
% `arguments` block cannot express. The F-07 substance lands anyway:
% canonical UpperCamelCase parameter names (inputParser matches
% case-insensitively, so every existing lowercase call site keeps
% working -- the same behaviour the arguments-migrated functions have),
% and TYPED validators on all five known options. Enum values are
% case-sensitive lowercase (strcmp), matching the migrated trio.
ip = inputParser;
ip.KeepUnmatched = true;
addParameter(ip, 'Method',  'qr', ...
    @(s) (ischar(s) || isstring(s)) && any(strcmp(char(s), {'qr','normal'})));
addParameter(ip, 'Sparse',  false, @(v) islogical(v) && isscalar(v));
addParameter(ip, 'DatumB',  [], @(v) isnumeric(v));
addParameter(ip, 'Sigma0',  1, @(v) isnumeric(v) && isscalar(v) && v > 0);
addParameter(ip, 'Verbose', true, ...
    @(v) (islogical(v) || isnumeric(v)) && isscalar(v));
parse(ip, varargin{:});

accArgs = {'method', ip.Results.Method, 'sparse', ip.Results.Sparse};
if ~isempty(ip.Results.DatumB)
    accArgs = [accArgs, {'datumB', ip.Results.DatumB}];
end

% ---------------------------------------------------------------- pass 1
acc = lsqvce.engine.Accumulator(u, accArgs{:});
for k = 1:nBatch
    [A, l, p] = callBatch(batchFcn, k);
    acc.add(A, l, p);
end
sol = acc.solve('sigma0', ip.Results.Sigma0);

% ---------------------------------------------------------------- pass 2
dArgs = varargin;
dArgs = stripArgs(dArgs, {'Method', 'Sigma0', 'Verbose', 'Sparse', 'DatumB'});  % strcmpi: any casing stripped
D = lsqvce.engine.diagnose(sol, batchFcn, 'nBatch', nBatch, dArgs{:});

% ---------------------------------------------------------------- verdict
V = lsqvce.report(D, 'verbose', ip.Results.Verbose);

end

function [A, l, p] = callBatch(f, k)
%CALLBATCH  Call a batch source f(k), tolerating sources that return only [A,l]
%   (no weights): returns p = [] when the source has no third output.
%   INPUT   f 1x1 fh  batch source;  k 1x1 double  batch index
%   OUTPUT  A nk x u  design;  l nk x 1  observations;  p nk x 1 | []  weights
    try
        [A, l, p] = f(k);
    catch ME
        if strcmp(ME.identifier, 'MATLAB:TooManyOutputs') || ...
           strcmp(ME.identifier, 'MATLAB:maxlhs') || ...
           strcmp(ME.identifier, 'MATLAB:unassignedOutputs')
            [A, l] = f(k);
            p = [];
        else
            rethrow(ME);
        end
    end
end

function a = stripArgs(a, names)
%STRIPARGS  Remove named Name/Value pairs from a varargin cell.
%   INPUT   a 1xN cell  varargin Name/Value pairs; names 1xK cell  names to drop
%   OUTPUT  a 1xM cell  the list with each matching Name and its Value removed
    keep = true(1, numel(a));
    i = 1;
    while i <= numel(a) - 1
        if ischar(a{i}) && any(strcmpi(a{i}, names))
            keep(i) = false;  keep(i+1) = false;
            i = i + 2;
        else
            i = i + 2;
        end
    end
    a = a(keep);
end
