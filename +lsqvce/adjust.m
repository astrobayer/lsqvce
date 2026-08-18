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
%   NAME-VALUE  (all forwarded, see Accumulator/diagnose/report)
%     'method'      'qr' (default) | 'normal'
%     'sigma0'      a-priori sigma0                       (default 1)
%     'alpha'       significance per test                 (default 1e-3)
%     'bonferroni'  correct for n simultaneous tests      (default true)
%     'beta'        MDB type-II error                     (default 0.20)
%     'rWeak'       weak-control threshold                (default 0.10)
%     'test'        'auto' | 'pope' | 'baarda'            (default 'auto')
%     'store'       'all' | 'flagged' | 'none'            (default 'all')
%     'verbose'     print the report                      (default true)
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
%   See also LSQVCE.DIAGNOSE, LSQVCE.REPORT, LSQVCE.SNOOP.
%
%
%   Generated with assistance from Claude (Anthropic), Claude Opus 4.8,
%   2026-07-11.

ip = inputParser;
ip.KeepUnmatched = true;
addParameter(ip, 'method',  'qr');
addParameter(ip, 'sparse',  false);
addParameter(ip, 'datumB',  []);
addParameter(ip, 'sigma0',  1);
addParameter(ip, 'verbose', true);
parse(ip, varargin{:});

accArgs = {'method', ip.Results.method, 'sparse', ip.Results.sparse};
if ~isempty(ip.Results.datumB)
    accArgs = [accArgs, {'datumB', ip.Results.datumB}];
end

% ---------------------------------------------------------------- pass 1
acc = lsqvce.Accumulator(u, accArgs{:});
for k = 1:nBatch
    [A, l, p] = callBatch(batchFcn, k);
    acc.add(A, l, p);
end
sol = acc.solve('sigma0', ip.Results.sigma0);

% ---------------------------------------------------------------- pass 2
dArgs = varargin;
dArgs = stripArgs(dArgs, {'method', 'sigma0', 'verbose', 'sparse', 'datumB'});
D = lsqvce.diagnose(sol, batchFcn, 'nBatch', nBatch, dArgs{:});

% ---------------------------------------------------------------- verdict
V = lsqvce.report(D, 'verbose', ip.Results.verbose);

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
