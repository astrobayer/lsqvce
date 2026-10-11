function [sol, D, V] = flatAdjust(u, batchFcn, nBatch, varargin)
%FLATADJUST  The retired flat driver's [sol, D, V] contract, on the block path.
%   [sol, D, V] = flatAdjust(u, batchFcn, nBatch, 'Name', value, ...)
%
%   R5 M5 5b' (WP-87): lsqvce.engine.adjust (the flat batch driver) is
%   retired; the 96 diagnostics-oracle points of the legacy suites keep
%   running against the ONE accumulator through lsqvce.fit on the batch
%   source. This is a TEST helper, not an alias in the package: it maps the
%   flat driver's options onto lsqvce.Options / data.priors and re-shapes
%   the Results into the (sol, D, V) triple the assertions read.
%
%   MAPPING
%     'Method'   -> Options.Backend
%     'DatumB'   -> data.priors (L = B', the exact minimum-constraint datum)
%     'Sigma0'   -> REFUSED (retired: the a-posteriori sigma0 is the report)
%     'Sparse'   -> ignored (retired, measured moot)
%     diagnose's lowercase names ('alpha', 'bonferroni', 'beta', 'rMin',
%     'rWeak', 'test', 'store', 'nWorst') -> Options.Diagnostics.<Name>
%     'Verbose'  -> Options.Verbose
%   OUTPUT (the flat contract)
%     sol.x, .n, .u, .r (= dof), .vPv, .s0sq, .sigma0 (a-priori 1), .rank,
%     .defect, .method, .datumB, .Cxx, .Q (the diagnostic-frame cofactor,
%     lsqvce.engine.BlockCofactor -- WP-86 part 2)
%     D   = R.diag (the same lsqvce.engine.verdict.diagnose on the whitened rows: D.v
%           is sqrt(p).*v -- the whitened residual)
%     V   = R.report(Verbose=false)
%
%   Created 2026-09-03 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M5, WP-87 round 5b'.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-86 part 2: sol carries the
%   block cofactor Q instead of the dense R. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
    p = inputParser;  p.KeepUnmatched = false;
    addParameter(p, 'Method', 'qr');
    addParameter(p, 'DatumB', []);
    addParameter(p, 'Sigma0', 1);
    addParameter(p, 'Sparse', false);
    addParameter(p, 'Verbose', false);
    addParameter(p, 'alpha', 1e-3);   addParameter(p, 'bonferroni', true);
    addParameter(p, 'beta', 0.2);     addParameter(p, 'rMin', 0.01);
    addParameter(p, 'rWeak', 0.1);    addParameter(p, 'test', 'auto');
    addParameter(p, 'store', 'all');  addParameter(p, 'nWorst', 1000);
    addParameter(p, 'robustCheck', true);
    parse(p, varargin{:});
    o = p.Results;
    if o.Sigma0 ~= 1
        error('lsqvce:test:sigma0Retired', ...
              'Sigma0 is retired (the a-posteriori sigma0 is the report); scale the weights instead.');
    end
    data = struct('batchFcn', batchFcn, 'numUnknowns', u, 'numBatches', nBatch);
    if ~isempty(o.DatumB)
        data.priors = struct('L', o.DatumB.', 'sigma', 1, 'name', 'datum');
    end
    opts = lsqvce.Options(Backend=string(lower(char(o.Method))), Verbose=logical(o.Verbose), ...
             Diagnostics=struct('Alpha', o.alpha, 'Bonferroni', o.bonferroni, 'Beta', o.beta, ...
                                'RMin', o.rMin, 'RWeak', o.rWeak, 'Test', string(o.test), ...
                                'Store', string(o.store), 'NWorst', o.nWorst, 'RobustCheck', o.robustCheck));
    R = lsqvce.fit(data, opts);
    D = R.diag;
    V = R.report('Verbose', false);
    d = numel([]);
    if ~isempty(o.DatumB), d = size(o.DatumB, 2); end
    sol = struct('x', R.x(:), 'n', D.n, 'u', u, 'r', D.n - u + d, ...
                 'vPv', R.sol.vtv, 's0sq', R.sigma0^2, 'sigma0', 1, ...
                 'rank', u - d, 'defect', d, 'method', char(o.Method), ...
                 'datumB', o.DatumB, 'Cxx', R.sol.Cxx);
    sol.Q = R.source.Q;                 % the block cofactor (WP-86 part 2)
end
