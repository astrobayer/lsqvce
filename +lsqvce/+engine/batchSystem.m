function L = batchSystem(u, batchFcn, nBatch, opt)
%LSQVCE.ENGINE.BATCHSYSTEM  A batch source as ONE streaming system (local-free by default).
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
%     'LocalColumns'  1 x nl double  columns of the source's rows that are this
%                                  system's LOCAL parameters (default [], none):
%                                  the rows are then full width, A = rows(:, 1:u)
%                                  and B = rows(:, LocalColumns) -- how snoop's
%                                  clean re-solve rebuilds a block-arrow system
%                                  from its batches (WP-86 part 3)
%     'SystemId'      1x1          the system's identifier (default 1)
%   OUTPUT
%     L  1x1 struct   the Loader (see LSQVCE.ENGINE.MAKELOADER)
%
%   ERRORS
%     lsqvce:batchSystem:badSource   batchFcn is not a function handle
%     lsqvce:batchSystem:badMatrix   a batch returned A/l of inconsistent size
%     lsqvce:batchSystem:weightMatrix a batch returned a weight argument that is
%                                    neither m x 1 nor an m x m SPD matrix, or
%                                    m x 1 weights that are negative or not finite
%     + whatever batchFcn throws, as lsqvce:callUser:threw with the cause
%
%   ACCURACY  exact: pre-whitening by sqrt(p) is the diagonal-P GLS identity;
%   the per-observation residual the block path reports is the whitened one,
%   sqrt(p) .* v.
%   A row of weight p = 0 keeps its index but is NOT an observation (WP-133):
%   it leaves the degrees of freedom and the verdict calls it 'rejected'.
%
%   EXAMPLE
%     L = lsqvce.engine.batchSystem(3, @(k) deal(A, l, p), 1);
%     R = lsqvce.engine.adjustBlock({L}, 1);
%
%   See also LSQVCE.FIT, LSQVCE.ENGINE.MAKELOADER, LSQVCE.ENGINE.ADJUSTNONLINEAR
%
%   Updated 2026-10-13 (CEST, UTC+02:00), WP-133: m x 1 weights must be >= 0
%   and finite; a zero weight marks its row excluded (out of the dof). With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-11 (CEST, UTC+02:00), R7 WP-86 part 3 (R7-24):
%   'LocalColumns' / 'SystemId': a batch source of full-width rows can stand
%   for a system WITH locals (snoop's block re-solve). With assistance from
%   Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-09-09 (CEST, UTC+02:00), R6 WP-102 (v6): the batch day
%   returns its one arc as one cell (the v6 day contract). With assistance
%   from Claude (Anthropic), model Claude Fable 5.1.
%   Created 2026-09-02 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M2, WP-81 round 2c'
%   (moved out of lsqvce.fit so the engine drivers can build the same system).

arguments
    u (1,1) double {mustBeInteger, mustBePositive}
    batchFcn
    nBatch (1,1) double {mustBeInteger, mustBePositive}
    opt.State = []
    opt.Name (1,:) char = 'batch'
    opt.LocalColumns (1,:) double {mustBeInteger, mustBePositive} = zeros(1, 0)
    opt.SystemId = 1
end
    if ~isa(batchFcn, 'function_handle')
        error('lsqvce:batchSystem:badSource', 'batchFcn must be a function handle.');
    end
    % nc / cols carried here: the non-scoped loader passes no info struct
    po = struct('batchFcn', batchFcn, 'state', opt.State, 'nc', u, 'cols', opt.LocalColumns);
    % R5 M5 5b': a RE-LINEARISING source (State given) must never be served
    % from the loader's day cache -- the cache would return the day built at
    % the previous x, and Gauss-Newton would take one step and then see the
    % same numbers for ever (found by the trilateration oracle: "converged"
    % to the wrong point). Linear sources keep the cache.
    L = lsqvce.engine.makeLoader(@batchDay, 'Cache', isempty(opt.State), ...
            'ObservationTypes', {'obs'}, ...
            'Covariance', struct('obs', struct('type', 'diag', 'sigma', 1)), ...
            'Arcs', @(day, obsType) 1, ...
            'NumCommon', u, 'NumLocal', numel(opt.LocalColumns), 'Days', 1:nBatch, ...
            'SystemId', opt.SystemId, 'Name', opt.Name, 'ProcessingOptions', po);
end

function D = batchDay(day, ~, ~, ~, po, ~)
%BATCHDAY  One batch as one day with unit-weight rows (its locals, if any,
%   sliced from the full-width rows).
    if isempty(po.state)
        [~, A, l, p] = lsqvce.util.callUser(po.batchFcn, 2, 3, 'batchFcn(k)', day);
    else
        x = po.state('x');
        [~, A, l, p] = lsqvce.util.callUser(po.batchFcn, 2, 3, 'batchFcn(k, x)', day, x);
    end
    l = l(:);  m = numel(l);  ex = {};
    width = max([po.nc, po.cols]);                % full-width rows when locals are sliced
    badWidth = size(A, 2) < width || (isempty(po.cols) && size(A, 2) ~= width);
    if ~isnumeric(A) || size(A, 1) ~= m || badWidth
        error('lsqvce:batchSystem:badMatrix', ...
              'batch %d returned A %s and l %d x 1; expected A %d x %d.', ...
              day, mat2str(size(A)), m, m, width);
    end
    if ~isempty(p)
        if isvector(p) && numel(p) == m
            % WP-133 (review 3 finding 17): a negative weight made the rows
            % complex (the QR and normal backends then disagreed: |p| against
            % p); a zero weight deletes the row but kept it in the count --
            % sigma0 0.753 for 0.949 with 10 of 30 rows at weight 0
            if any(~isfinite(p(:))) || any(p(:) < 0)
                error('lsqvce:batchSystem:weightMatrix', ...
                      'batch %d returned weights that are negative or not finite; a weight is >= 0.', day);
            end
            s = sqrt(p(:));
            A = s .* A;  l = s .* l;
            if any(p(:) == 0), ex = {p(:) == 0}; end
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
    B = A(:, po.cols);  A = A(:, 1:po.nc);        % whitened together, then split
    D.A = {A};  D.B = {B};  D.l = {l};  D.t = {(1:m).'};   % v6: one cell = one arc
    if ~isempty(ex), D.excluded = ex; end         % weight 0: kept for its index, not an observation
end
