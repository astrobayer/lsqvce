function [sol, D, V] = adjust(u, batchFcn, nBatch, varargin)
%LSQVCE.ADJUST  Solve, diagnose and judge a batch least-squares problem.
%
%   [sol, D, V] = LSQVCE.ADJUST(u, batchFcn, nBatch) runs the whole pipeline:
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
%     sol  1x1 struct  .x (u x 1 estimate), .Cxx, .sigma0, .rank, ...
%     D    1x1 struct  diagnostics: residuals, redundancy numbers, w/tau,
%                      leverage, flagIdx, MDB -- see LSQVCE.DIAGNOSE
%     V    1x1 struct  the level-0..4 verdict -- see LSQVCE.REPORT
%
%   NAME-VALUE  (WP-39: adjust's OWN five options carry canonical
%   UpperCamelCase names with typed validators; any casing still matches
%   (inputParser, case-insensitive -- same behaviour as the
%   arguments-migrated functions). Every OTHER pair is FORWARDED verbatim
%   to LSQVCE.DIAGNOSE -- the deliberate open pass-through that keeps
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
%       [sol, D, V] = lsqvce.adjust(u, src, 4);
%
%       % V.verdict      -> 'PASS WITH WARNINGS'
%       % D.flagIdx      -> [40 150]      both blunders found
%       % D.s0           -> 1.675         variance factor inflated by them
%       % D.sumRRelErr   -> 0             the integrity identity holds
%
%   Updated 2026-08-19 (CEST, UTC+02:00), WP-39 stage (a) round a2
%   (F-07): canonical UpperCamelCase names + typed validators on the five
%   own options; enum VALUE 'Method' case-sensitive lowercase; the open
%   pass-through to LSQVCE.DIAGNOSE keeps this front door on inputParser
%   BY DESIGN (pinned deviation, §9 WP-39). With assistance from Claude
%   (Anthropic), model Claude Fable 5.
%
%   See also LSQVCE.DIAGNOSE, LSQVCE.REPORT, LSQVCE.SNOOP.
%
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.

% WP-39 stage (a) round a2 (F-07), PINNED DEVIATION: adjust KEEPS
% inputParser BY DESIGN -- KeepUnmatched is a documented FEATURE here
% (unmatched options are forwarded to LSQVCE.DIAGNOSE below), which an
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
acc = lsqvce.Accumulator(u, accArgs{:});
for k = 1:nBatch
    [A, l, p] = callBatch(batchFcn, k);
    acc.add(A, l, p);
end
sol = acc.solve('sigma0', ip.Results.Sigma0);

% ---------------------------------------------------------------- pass 2
dArgs = varargin;
dArgs = stripArgs(dArgs, {'Method', 'Sigma0', 'Verbose', 'Sparse', 'DatumB'});  % strcmpi: any casing stripped
D = lsqvce.diagnose(sol, batchFcn, 'nBatch', nBatch, dArgs{:});

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
