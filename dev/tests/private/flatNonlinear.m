function [sol, D, V, IT] = flatNonlinear(u, batchFcn, nBatch, x0, varargin)
%FLATNONLINEAR  The retired flat Gauss-Newton driver's [sol, D, V, IT] on the block path.
%   R5 M5 5b' test bridge (see flatAdjust): 'maxIter' -> MaxPasses,
%   'tolX' -> Tolerance, 'tolF'/'lambda0' -> GaussNewton, 'verbose' -> Verbose.
%   IT.converged, IT.iterations, IT.history (the pass record), IT.reason.
%   Created 2026-09-03 with assistance from Claude (Anthropic), model Claude Fable 5.1.
    p = inputParser;
    addParameter(p, 'maxIter', 20);  addParameter(p, 'tolX', 1e-8);
    addParameter(p, 'tolF', 1e-10);  addParameter(p, 'lambda0', 1e-6);
    addParameter(p, 'verbose', false);  addParameter(p, 'Method', 'qr');
    addParameter(p, 'sigma0', 1);    addParameter(p, 'DatumB', []);
    parse(p, varargin{:});  o = p.Results;
    data = struct('batchFcn', batchFcn, 'numUnknowns', u, 'numBatches', nBatch, 'x0', x0(:));
    if ~isempty(o.DatumB), data.priors = struct('L', o.DatumB.', 'sigma', 1, 'name', 'datum'); end
    w = warning('off', 'lsqvce:adjustNonlinear:maxIter');  c1 = onCleanup(@() warning(w)); %#ok<NASGU>
    w2 = warning('off', 'lsqvce:adjustNonlinear:diverged'); c2 = onCleanup(@() warning(w2)); %#ok<NASGU>
    opts = lsqvce.Options(MaxPasses=o.maxIter, Tolerance=o.tolX, Verbose=logical(o.verbose), ...
             Backend=string(lower(char(o.Method))), GaussNewton=struct('TolF', o.tolF, 'Lambda0', o.lambda0));
    R = lsqvce.fit(data, opts);
    D = R.diag;  V = R.report('Verbose', false);
    h = R.history;
    IT = struct('converged', R.converged, 'iterations', numel(h), 'nIter', numel(h), 'history', h, ...
                'reason', R.noise.reason, 'vPv', [h.vPv].');
    sol = struct('x', R.x(:), 'n', D.n, 'u', u, 'r', D.n - u, 's0sq', R.sigma0^2, 'sigma0', 1, 'Cxx', R.sol.Cxx);
end
