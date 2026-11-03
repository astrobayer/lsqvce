function [sol, D, ROB] = flatRobust(u, batchFcn, nBatch, varargin)
%FLATROBUST  The retired flat robust driver's [sol, D, ROB] on the block path.
%   R5 M5 5b' test bridge (see flatAdjust): 'influence' -> Options.Robust,
%   'maxIter'/'tol'/'tuning' -> RobustOpts, 'verbose' -> Verbose.
%   ROB.factor = the final per-observation factors, ROB.iterations = passes,
%   ROB.converged, ROB.scale.
%   Created 2026-09-03 with assistance from Claude (Anthropic), model Claude Fable 5.1.
    p = inputParser;
    addParameter(p, 'influence', 'huber');  addParameter(p, 'maxIter', 30);
    addParameter(p, 'tol', 1e-8);           addParameter(p, 'tuning', []);
    addParameter(p, 'verbose', false);      addParameter(p, 'Method', 'qr');
    addParameter(p, 'sigma0', 1);           addParameter(p, 'test', 'auto');
    parse(p, varargin{:});  o = p.Results;
    data = struct('batchFcn', batchFcn, 'numUnknowns', u, 'numBatches', nBatch);
    opts = lsqvce.Options(Robust=string(lower(char(o.influence))), Verbose=logical(o.verbose), ...
             Backend=string(lower(char(o.Method))), ...
             RobustOpts=struct('MaxIter', o.maxIter, 'Tol', o.tol, 'Tuning', o.tuning), ...
             Diagnostics=struct('Test', string(o.test)));
    R = lsqvce.engine.fitDeclaration(data, opts);   % WP-162: the engine tier takes the struct
    D = R.diag;
    h = R.history;
    ROB = struct('factor', h(end).factor(:), 'iterations', numel(h), 'iter', numel(h), ...
                 'converged', R.converged, 'scale', h(end).scale, 'influence', char(o.influence));
    sol = struct('x', R.x(:), 'n', D.n, 'u', u, 'r', D.n - u, 's0sq', R.sigma0^2, 'sigma0', 1, 'Cxx', R.sol.Cxx);
end
