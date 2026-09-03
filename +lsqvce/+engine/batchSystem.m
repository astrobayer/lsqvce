function L = batchSystem(u, batchFcn, nBatch, opt)
%LSQVCE.ENGINE.BATCHSYSTEM  A batch source as ONE local-free streaming system.
%   L = lsqvce.engine.batchSystem(u, batchFcn, nBatch)
%   L = lsqvce.engine.batchSystem(u, batchFcn, nBatch, 'State', xState, ...)
%
%   THE adapter that lets the in-memory and batched fit kinds ride the block
%   path (R5 M2, WP-81 rounds 2a/2c'; finding H-02: a batch problem is a
%   streaming problem with nl = 0): batch k becomes day k, one arc per day,
%   no local block, rows PRE-WHITENED by sqrt(p) so the declared law is
%   diag(1). Returns the Loader object lsqvce.engine.makeLoader builds for it,
%   ready for adjustBlock / adjustRobust / adjustNonlinear.
%
%   Linear source:      [A, l]    = batchFcn(k)   or   [A, l, p] = batchFcn(k)
%   Nonlinear source:   [A, w, p] = batchFcn(k, x)  -- when 'State' is given,
%   the day builder reads the CURRENT x from xState('x') at every call, so a
%   Gauss-Newton driver re-linearises by writing xState('x') and re-reading.
%
%   INPUT
%     u         1x1 double          number of unknowns (columns of A)
%     batchFcn  1x1 function_handle the batch source (called through
%                                   lsqvce.util.callUser; count learned once)
%     nBatch    1x1 double          number of batches (>= 1)
%   NAME-VALUE
%     'State'   containers.Map | []  when non-empty, batchFcn(k, State('x'))
%     'Name'    char               the system's label (default 'batch')
%   OUTPUT
%     L  1x1 struct   the Loader (see LSQVCE.ENGINE.MAKELOADER)
%
%   ERRORS
%     lsqvce:batchSystem:badSource   batchFcn is not a function handle
%     lsqvce:batchSystem:badMatrix   a batch returned A/l of inconsistent size
%     lsqvce:batchSystem:weightMatrix a batch returned a weight argument that is
%                                    neither m x 1 nor an m x m SPD matrix
%     + whatever batchFcn throws, as lsqvce:callUser:threw with the cause
%
%   ACCURACY  exact: pre-whitening by sqrt(p) is the diagonal-P GLS identity;
%   the per-observation residual the block path reports is the whitened one,
%   sqrt(p) .* v.
%
%   EXAMPLE
%     L = lsqvce.engine.batchSystem(3, @(k) deal(A, l, p), 1);
%     R = lsqvce.engine.adjustBlock({L}, 1);
%
%   See also LSQVCE.FIT, LSQVCE.ENGINE.MAKELOADER, LSQVCE.ENGINE.ADJUSTNONLINEAR
%
%   Created 2026-09-02 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M2, WP-81 round 2c'
%   (moved out of lsqvce.fit so the engine drivers can build the same system).

arguments
    u (1,1) double {mustBeInteger, mustBePositive}
    batchFcn
    nBatch (1,1) double {mustBeInteger, mustBePositive}
    opt.State = []
    opt.Name (1,:) char = 'batch'
end
    if ~isa(batchFcn, 'function_handle')
        error('lsqvce:batchSystem:badSource', 'batchFcn must be a function handle.');
    end
    po = struct('batchFcn', batchFcn, 'state', opt.State, 'nc', u);   % nc carried here: the
                                                                    % non-scoped loader passes no info struct
    % R5 M5 5b': a RE-LINEARISING source (State given) must never be served
    % from the loader's day cache -- the cache would return the day built at
    % the previous x, and Gauss-Newton would take one step and then see the
    % same numbers for ever (found by the trilateration oracle: "converged"
    % to the wrong point). Linear sources keep the cache.
    L = lsqvce.engine.makeLoader(@batchDay, 'Cache', isempty(opt.State), ...
            'ObservationTypes', {'obs'}, ...
            'Covariance', struct('obs', struct('type', 'diag', 'sigma', 1)), ...
            'Arcs', @(day, obsType) 1, ...
            'NumCommon', u, 'NumLocal', 0, 'Days', 1:nBatch, ...
            'SystemId', 1, 'Name', opt.Name, 'ProcessingOptions', po);
end

function D = batchDay(day, ~, ~, ~, po, ~)
%BATCHDAY  One batch as one local-free day with unit-weight rows.
    if isempty(po.state)
        [~, A, l, p] = lsqvce.util.callUser(po.batchFcn, 2, 3, 'batchFcn(k)', day);
    else
        x = po.state('x');
        [~, A, l, p] = lsqvce.util.callUser(po.batchFcn, 2, 3, 'batchFcn(k, x)', day, x);
    end
    l = l(:);  m = numel(l);
    if ~isnumeric(A) || size(A, 1) ~= m || size(A, 2) ~= po.nc
        error('lsqvce:batchSystem:badMatrix', ...
              'batch %d returned A %s and l %d x 1; expected A %d x %d.', ...
              day, mat2str(size(A)), m, m, po.nc);
    end
    if ~isempty(p)
        if isvector(p) && numel(p) == m
            s = sqrt(p(:));
            A = s .* A;  l = s .* l;
        elseif ismatrix(p) && isequal(size(p), [m m])
            % a full weight matrix P = W'W (correlated rows within the batch):
            % whitening by W is the exact GLS identity (R5 M5 5b'); P must be SPD
            [W, flag] = chol(full(p));
            if flag > 0
                error('lsqvce:batchSystem:weightMatrix', ...
                      'batch %d returned a %d x %d weight matrix that is not positive definite.', day, m, m);
            end
            A = W * A;  l = W * l;
        else
            error('lsqvce:batchSystem:weightMatrix', ...
                  ['batch %d returned a %s weight argument; expected m x 1 weights ' ...
                   'or an m x m positive-definite weight matrix (m = %d).'], day, mat2str(size(p)), m);
        end
    end
    D.A = A;  D.B = zeros(m, 0);  D.l = l;  D.t = (1:m).';
    D.arcIdx = {1:m};  D.arcId = 1;
end
