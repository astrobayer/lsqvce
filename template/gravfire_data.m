function data = gravfire_data(commonData, myOptions, systemsData, opt)
%GRAVFIRE_DATA  Map a GRAVFIRE preprocessing bundle onto the lsqvce `data` struct.
%   data = gravfire_data(commonData, myOptions, systemsData)
%   data = gravfire_data(..., 'DesignFcn', @myDesign)      % your A/B builder
%   data = gravfire_data(..., 'DryRun', true)              % random A/B: pipeline check
%
%   TEMPLATE (R5 M0 close-out, HANDOVER §12.8): the GRAVFIRE bundle
%   (commonData / myOptions / systemsData, as saved in testdata.mat) already
%   contains everything the v5 data contract asks for EXCEPT the design
%   matrices. This function is the declaration half; the design half is the
%   'DesignFcn' you supply (the variational-equation partials of GRAVFIRE).
%   Every rule of the API reference is applied here once, and the file is
%   meant to be READ as the step-by-step guide dev/GRAVFIRE_INTERFACE.md:
%   the section numbers in the comments below are its steps.
%
%   INPUT
%     commonData   1x1 struct   .parameter  R x 14 cell, the common table
%     myOptions    1x1 struct   the GRAVFIRE run options (.days, .observationGroups, ...)
%     systemsData  1xS cell     one bundle per system; bundle{1}.observations is a
%                               (days x groups) cell of per-day structs, .parameter the
%                               EXPANDED local table (one row per day and satellite),
%                               .numberOfPieces (days x groups), .variance (1 x groups)
%   NAME-VALUE
%     'DesignFcn'  function_handle   [A, B] = f(o, day, obsType, common, local, po, info)
%                  o = the (day, group) bundle struct; must return A (m x info.nc) and
%                  B (m x info.nl), rows EPOCH-MAJOR for multi-axis streams, m = EVERY
%                  epoch of the day times its axes -- the flagged epochs included
%                  (their rows are zeroed and weighted out, whatever you return). Default:
%                  a function that ERRORS with gravfire:designMissing (so a forgotten
%                  design step fails loudly, wrapped by lsqvce:makeLoader:buildDay).
%     'DryRun'     logical, default false   replace the design step by DETERMINISTIC
%                  random A/B (seeded per day and observation type) so that the
%                  declaration, the arc cut, Cep, nl/nc enforcement and the whole
%                  fit pipeline can be exercised BEFORE the design step exists.
%                  The numbers of a dry run mean nothing.
%   OUTPUT
%     data  1x1 struct   the v6 lsqvce declaration: commonParams, days, priors,
%           systems(1) with id, name, buildDay, covariance, localParams,
%           priors, processingOptions -- and NO arc schedule: the pieces
%           buildDay returns as cells ARE the arcs of the day (v6, R6 WP-102).  The two priors fields are the DECLARED
%           regularisation (the API reference): data.priors from column 11
%           of the common table, data.systems.priors from the local templates.
%
%   ERRORS
%     gravfire:designMissing        DryRun false and no DesignFcn given (raised at read time)
%     gravfire:unknownQuantity      an observation group's observationQuantity is not
%                                   'kinematicPosition' or 'range-rate'
%     gravfire:designShape          DesignFcn returned A/B with the wrong row/column count
%     gravfire:badConstraint        column 11 is malformed, or the rows of one label
%                                   disagree on sigma or on the activation pass, or
%                                   (local rows) on the a-priori value, column 10
%     (+ every declaration check of lsqvce.fit, raised before any read; the API reference, chapter buildDay)
%
%   EXAMPLE
%     S    = load('testdata.mat');
%     data = gravfire_data(S.commonData, S.myOptions, S.systemsData, 'DryRun', true);
%     R    = lsqvce.fit(data, lsqvce.Options(Verbose=1));   % pipeline check
%     data = gravfire_data(S.commonData, S.myOptions, S.systemsData, 'DesignFcn', @gravfireDesign);
%     R    = lsqvce.fit(data, lsqvce.Options(Estimate="noise", Noise=struct('Grid', 2160), MaxPasses=7));
%
%   See also LSQVCE.FIT, LSQVCE.OPTIONS, the API reference, dev/GRAVFIRE_INTERFACE.md
%
%   Updated 2026-09-09 (CEST, UTC+02:00), R6 WP-102 (v6): THE CELLS ARE THE
%   ARCS -- buildGravfireDay returns one cell per piece (A, B, l, t, Cep),
%   the arc schedule (STEP 5, arcsFcn) is gone: the locals are 'day'-scoped,
%   so no plan needs the piece count before a read, and the toolbox learns
%   the pieces of a day from its cells. With assistance from Claude
%   (Anthropic), model Claude Fable 5.1.
%
%   VERIFIED 2026-09-02 on the user's testdata.mat (R2026a): the declaration
%   passes every declaration check of fit; DryRun + Estimate="none" accumulates 2 days x 16
%   pieces (4 posA + 4 posB + 8 kbr), nc = 117, nl = 144 derived from the 72
%   templates, Cep accepted, in 2.0 s. Estimate="noise", Pool="stream", ArcFactors=false runs (62 s, 2 passes).
%   R5 M3 made the 3-axis pieces fast on the joint estimator (38.6 s for both
%   days at 2160 lags, HANDOVER §12.5 round 3c).
%
%   Updated 2026-09-07 (CEST, UTC+02:00), R6 WP-95: the v5 data contract is
%   named as such, STEP 6 emits ONE prior struct per label with the member
%   rows stacked into L (no 'group' beside L), and STEP 7 documents the
%   declared data.systems.priors instead of the retired priorDay handle.
%   With assistance from Claude (Anthropic), model Claude Opus 5.
%   Updated 2026-09-06 (CEST, UTC+02:00), R6 WP-90: the positions declare
%   the temporal law only and return D.Cep with the day (no 'separable'
%   type, no axis count in the spec). With assistance from Claude
%   (Anthropic), model Claude Fable 5.1.
%   Created 2026-09-02 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M0 close-out (§12.8).
%
%   Updated 2026-10-17 (CEST, UTC+02:00), WP-148: the flagged epochs stay on
%   the grid with weight 0 (D.excluded) instead of being removed -- a removed
%   epoch left a gap a 'psd' law took for one step. DesignFcn returns the rows
%   of every epoch. With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-17 (CEST, UTC+02:00), WP-144 (A4): the members of a
%   labelled LOCAL group must share fromPass and the a-priori value too (a
%   later member's were dropped). With assistance from Claude (Anthropic),
%   model Claude Opus 5.5.
%   Updated 2026-10-12 (CEST, UTC+02:00), WP-124: references to the data
%   contract point at the API reference. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.

arguments
    commonData  (1,1) struct
    myOptions   (1,1) struct
    systemsData (1,:) cell
    opt.DesignFcn = @designMissing
    opt.DryRun (1,1) logical = false
end
    sd   = systemsData{1};
    days = double(myOptions.days(:).');                       % 1 x D, integers (MJD)

    % ---- STEP 1: the common table -- scope token first, your 14 columns after
    data.commonParams = [repmat({'common'}, size(commonData.parameter, 1), 1), ...
                         commonData.parameter];

    % ---- STEP 2: the calendar
    data.days = days;

    % ---- STEP 3: observation types are FIELD NAMES of `covariance`
    groups = myOptions.observationGroups;                     % 1 x G cell of structs
    G      = numel(groups);
    obsTypes = cell(1, G);
    for g = 1:G
        obsTypes{g} = obsTypeOf(groups{g});
    end
    % The covariance spec declares the TEMPORAL law only (white to start; the
    % estimator updates it). The per-epoch error information of the positions
    % -- the 3x3 covariance from iSpp -- comes back from buildDay WITH the day
    % (D.Cep, STEP 10) and the toolbox composes C = D (R_t kron I_3) D (R6
    % WP-90: no 'separable' type, no axis count stated here).
    cov = struct();
    for g = 1:G
        cov.(obsTypes{g}) = struct('type', 'diag', 'sigma', sqrt(sd.variance(g)));
    end

    % ---- STEP 4: the local table -- ONE day's templates, all 'day'-scoped,
    %      times as seconds-of-day (the day itself is NOT in a template)
    localParams = templateFromExpanded(sd.parameter, days(1));

    % ---- STEP 5: the arc schedule -- NONE (v6). The pieces of a (day, group)
    %      come back from buildDay as cells, one per piece, and the toolbox
    %      takes them as the arcs of the day. numberOfPieces is not restated.
    %      (A schedule is required only for 'arc'-scoped template rows.)

    % ---- STEP 6: what buildDay needs, bound once
    po = struct();
    po.days      = days;
    po.obsTypes  = obsTypes;                                  % 1 x G cell
    po.bundles   = sd.observations;                           % (days x groups) cell
    po.designFcn = opt.DesignFcn;
    po.dryRun    = opt.DryRun;
    po.myOptions = myOptions;

    % ---- STEP 7: the system
    s = struct();
    s.id                = 1;                                  % numeric or ONE char (API reference, The data struct)
    s.name              = 'GRACE';                            % the label lives here
    s.buildDay          = @buildGravfireDay;
    s.covariance        = cov;
    s.localParams       = localParams;
    s.processingOptions = po;
    data.systems = s;

    % ---- STEP 6 (R5-D8 round 1): column 11 -> data.priors ------------------
    % Your constraining factor {'Tikhonov', k, sigma[, label]} on a COMMON row
    % j constrains that parameter to x0 = column 10 (the a-priori, [] -> 0)
    % with standard deviation sigma, from pass k on. ONE struct per GROUP:
    % without a label the row is its own group, L = e_j' and sigma fixed
    % (estimate "none" -- one row has redundancy <= 1); a 4th entry is a group
    % LABEL, and its rows become ONE struct whose L stacks their unit rows,
    % with one shared sigma estimated by VCE (estimate "scale", Koch-Kusche).
    % A group is declared by L OR by a name pattern, never both (round 6). In
    % your 15-column table (token in column 1) these are columns 11 and 12.
    data.priors = priorsFromColumn11(commonData.parameter, opt.DryRun);
    % ---- STEP 7 (R5-D8 round 6): column 11 on LOCAL rows -> systems.priors --
    % Regularisation of local parameters is DECLARED, never computed:
    % data.systems(k).priors(l) is one struct per GROUP, given EITHER by
    % 'group' -- a regexp on the NAME column of localParams (column 2 of the
    % 15-column template), from which the toolbox builds L with one row per
    % matched template row per INSTANCE, so a 'day' row contributes one row
    % per day and the group spans the calendar -- OR by 'L', written in the
    % columns of B as lsqvce.layout(data) announces them. Never both. Rows
    % sharing a label (4th entry of column 11, e.g. 'accscale') are ONE group
    % with one sigma and estimate "scale"; without a label each constrained
    % row is its own group with a fixed sigma. How many groups your bundle
    % yields therefore depends on how many template rows carry a constraint
    % and how they are labelled. The helper below emits the 'group' form.
    data.systems.priors = gravfireLocalPriors(data, localParams);
end

% ======================================================================
function D = buildGravfireDay(day, obsType, common, local, po, info)
%BUILDGRAVFIREDAY  One (day, observation type): raw A, B, l, t, Cep -- ONE CELL PER PIECE.
    g  = groupIndex(po.obsTypes, obsType);
    d  = find(po.days == day, 1);
    o  = po.bundles{d, g};

    nAx = o.numberOfAxes;
    % The flagged epochs STAY on the grid with weight 0 (WP-148): removing one
    % leaves a gap inside its piece that a 'psd' law, counting lags in
    % samples, takes for ONE step (lsqvce:covModel:unevenGrid). Their rows are
    % zeroed below and marked D.excluded; the toolbox deletes them from the
    % law's covariance, so their values -- NaN included -- never matter.
    flag  = logical(o.nanIndexOfPiece{1}(:));                 % nEp x 1
    lobs  = o.lobs{1};                                        % nEp x nAx
    lobs(flag, :) = 0;
    tEp   = o.time{1}(:, 2);                                  % seconds of day
    piece = o.indexOfPiece{1}(:);                             % nEp x 1
    nEp   = size(lobs, 1);
    m     = nEp * nAx;
    ex    = reshape(repmat(flag.', nAx, 1), [], 1);           % m x 1, epoch-major

    % STEP 8: epoch-major rows (x1 y1 z1 x2 y2 z2 ...), the whole day first
    l = reshape(lobs.', [], 1);
    t = repelem(tEp(:), nAx);

    % STEP 9: the ERROR INFORMATION -- the per-epoch 3x3 covariance from
    % iSpp as nEpoch x 3 x 3 (the axis count is derived from it, WP-90). A
    % scalar stream with a per-epoch sigma would carry it as nEpoch x 1.
    Cep = [];
    if nAx > 1
        C   = o.iSpp{1};                                      % nEp x 9 (3x3 per epoch)
        C(flag, :) = repmat(mean(C(~flag, :), 1), nnz(flag), 1);   % deleted: any valid 3x3
        Cep = permute(reshape(C.', nAx, nAx, []), [3 1 2]);   % nEp x 3 x 3
    end

    % STEP 10: the design matrices -- yours, or the dry-run stand-in (whole
    % day, m x nc / m x nl; the cut into pieces follows)
    if po.dryRun
        sd = mod(7919 * day + sum(double(obsType)), 2^31 - 1);
        rs = RandStream('mt19937ar', 'Seed', sd);             % deterministic per call
        A = randn(rs, m, info.nc);
        B = randn(rs, m, info.nl);
    else
        [A, B] = po.designFcn(o, day, obsType, common, local, po, info);
        if ~isequal(size(A), [m info.nc]) || ~isequal(size(B), [m info.nl])
            error('gravfire:designShape', ...
                  ['DesignFcn returned A %s / B %s; expected A [%d %d] and B [%d %d] ' ...
                   '(day %d, %s).'], mat2str(size(A)), mat2str(size(B)), ...
                   m, info.nc, m, info.nl, day, obsType);
        end
    end
    A(ex, :) = 0;  B(ex, :) = 0;                              % weight 0: zero rows (API reference)

    % STEP 11: THE CELLS ARE THE ARCS (v6) -- one cell per piece, in piece
    % order; the toolbox takes them as this day's arcs (ids 1:nPiece)
    ids  = unique(piece(:)).';
    nP   = numel(ids);
    D = struct('A', {cell(1, nP)}, 'B', {cell(1, nP)}, 'l', {cell(1, nP)}, 't', {cell(1, nP)}, ...
               'excluded', {cell(1, nP)});
    if nAx > 1, D.Cep = cell(1, nP); end
    for k = 1:nP
        ep   = find(piece == ids(k));                          % this piece's epochs
        rows = (ep(:) - 1) * nAx + (1:nAx);                   % nEp_k x nAx, epoch per row
        rows = reshape(rows.', [], 1);                         % their epoch-major rows
        D.A{k} = A(rows, :);  D.B{k} = B(rows, :);  D.l{k} = l(rows);  D.t{k} = t(rows);
        D.excluded{k} = ex(rows);
        if nAx > 1, D.Cep{k} = Cep(ep, :, :); end
    end
end

% ======================================================================
function [A, B] = designMissing(varargin)
    error('gravfire:designMissing', ...
          ['No DesignFcn given: the design matrices are not part of the ' ...
           'preprocessing bundle. Pass ''DesignFcn'', @yourDesign (or ''DryRun'', true ' ...
           'for a pipeline check).']);
end

function name = obsTypeOf(group)
%OBSTYPEOF  'GraceA' + kinematicPosition -> 'posA'; range-rate -> 'kbr'.
    q = char(group.observationQuantity);
    switch q
        case 'kinematicPosition'
            name = ['pos' upper(group.name(end))];
        case 'range-rate'
            name = 'kbr';
        otherwise
            error('gravfire:unknownQuantity', ...
                  'observation group %s has quantity %s (expected kinematicPosition or range-rate).', ...
                  group.name, q);
    end
end

function g = groupIndex(obsTypes, obsType)
    g = find(strcmp(obsTypes, char(obsType)), 1);
end

function T = templateFromExpanded(P, day0)
%TEMPLATEFROMEXPANDED  One day's rows of the EXPANDED local table, as templates.
%   Selects the rows whose flag column (6) names day0, converts the absolute
%   [MJD sec] spans (columns 2, 3) to seconds-of-day, converts the knot table
%   (column 4: [MJD-fraction value]) to seconds-of-day, and drops the day from
%   the flags. Column 1 becomes the scope token 'day'; the 14 GRAVFIRE columns
%   follow unchanged in position. When GRAVFIRE emits templates itself, this
%   helper disappears.
    flags = P(:, 6);
    sel   = cellfun(@(f) ~isempty(f) && f(1) == day0, flags);
    P     = P(sel, :);
    n     = size(P, 1);
    for k = 1:n
        for c = [2 3]                                         % [MJD sec] -> sec
            v = P{k, c};
            if numel(v) == 2, P{k, c} = v(2); end
        end
        v = P{k, 4};                                          % knots [MJD.frac value]
        if ~isempty(v) && size(v, 2) == 2
            P{k, 4} = [(v(:, 1) - day0) * 86400, v(:, 2)];
        end
        f = P{k, 6};
        if ~isempty(f), P{k, 6} = f(2:end); end               % drop the day
    end
    T = [repmat({'day'}, n, 1), P];
end

%==========================================================================
function priors = priorsFromColumn11(P, ~)
%PRIORSFROMCOLUMN11  Common-row constraints -> data.priors (R5-D8 round 1).
%   P is the ORIGINAL 14-column table (before the token is prepended):
%   column 10 = a-priori value ([] -> 0), column 11 = 'none' | [] |
%   {'Tikhonov', fromPass, sigma[, label]}.
%   ONE struct per GROUP, regrouped for R5-D8 round 6: rows sharing a label
%   become one struct whose L STACKS their unit rows e_j', with one shared
%   sigma and one activation pass, estimated by VCE (estimate "scale",
%   Koch-Kusche -- several rows are what make a sigma estimable); a row
%   without a label is its own struct, L = e_j' and sigma fixed (estimate
%   "none": one row has redundancy <= 1). A group is declared EITHER by
%   'group' -- a name pattern on commonParams -- OR by L, never both
%   (lsqvce:fit:badPrior); the common names are SYNTHESISED here rather than
%   a column of the table, so every group is declared by L and no 'group'
%   field is written. Members of one label must agree on sigma and fromPass.
    priors = struct('L', {}, 'x0', {}, 'sigma', {}, 'name', {}, 'estimate', {}, 'from', {});
    nc = size(P, 1);
    groups = containers.Map('KeyType', 'char', 'ValueType', 'any');  order = {};
    for j = 1:nc
        c = P{j, 11};
        if isNoneConstraint(c), continue, end
        if ~iscell(c) || numel(c) < 3 || ~ischar(c{1})
            error('gravfire:badConstraint', 'row %d: column 11 must be ''none'' or {''Tikhonov'', fromPass, sigma[, label]}.', j);
        end
        if ~strcmpi(c{1}, 'Tikhonov')
            error('gravfire:badConstraint', 'row %d: unknown constraint kind ''%s'' (only Tikhonov is mapped).', j, c{1});
        end
        x0 = 0;
        if ~isempty(P{j, 10}), x0 = double(P{j, 10}); end
        lab = '';
        if numel(c) >= 4 && ~isempty(c{4}), lab = char(string(c{4})); end
        nm  = sprintf('%s_%d_%d', char(string(P{j, 1})), double(P{j, 2}), double(P{j, 3}));
        key = lab;  if isempty(key), key = nm; end
        if ~isKey(groups, key)
            groups(key) = struct('rows', [], 'x0', [], 'sigma', double(c{3}), ...
                                 'from', double(c{2}), 'labelled', ~isempty(lab));
            order{end+1} = key; %#ok<AGROW>
        end
        g = groups(key);
        if g.sigma ~= double(c{3})
            error('gravfire:badConstraint', ...
                  'group ''%s'': members must share one sigma (row %d has %g, the group %g).', ...
                  key, j, double(c{3}), g.sigma);
        end
        if g.from ~= double(c{2})
            error('gravfire:badConstraint', ...
                  'group ''%s'': members must share one fromPass (row %d has %g, group %g).', ...
                  key, j, double(c{2}), g.from);
        end
        g.rows(end+1) = j;  g.x0(end+1, 1) = x0;
        groups(key) = g;
    end
    for k = 1:numel(order)
        g = groups(order{k});
        L = zeros(numel(g.rows), nc);
        for r = 1:numel(g.rows), L(r, g.rows(r)) = 1; end
        est = 'none';  if g.labelled, est = 'scale'; end
        priors(end+1) = struct('L', L, 'x0', g.x0, 'sigma', g.sigma, 'name', order{k}, ...
                               'estimate', est, 'from', g.from); %#ok<AGROW>
    end
    if isempty(priors), priors = []; end
end

%==========================================================================
function P = gravfireLocalPriors(~, localParams)
%GRAVFIRELOCALPRIORS  Column 11 of the LOCAL template rows -> data.systems.priors (R5-D8 round 6).
%   One struct per GROUP in the PATTERN form: rows sharing a label (4th entry
%   of column 11 -> column 12 here) become one struct whose 'group' is the
%   regexp of their names ('^(GraceA_xaccscale|GraceB_xaccscale)$'); rows
%   without a label are one struct each. The toolbox builds L from the names
%   (column 2 of the 15-column template) for every instance -- per day for a
%   'day' row. x0 = column 10 (-> 11), sigma and fromPass from the entry;
%   estimate "scale" for labelled groups (Koch-Kusche), "none" otherwise.
    P = struct('group', {}, 'x0', {}, 'sigma', {}, 'from', {}, 'name', {}, 'estimate', {});
    groups = containers.Map('KeyType', 'char', 'ValueType', 'any');  order = {};
    for j = 1:size(localParams, 1)
        c = localParams{j, 12};
        if isNoneConstraint(c), continue, end
        if ~iscell(c) || numel(c) < 3 || ~ischar(c{1}) || ~strcmpi(c{1}, 'Tikhonov')
            error('gravfire:badConstraint', 'local row %d: column 11 must be ''none'' or {''Tikhonov'', fromPass, sigma[, label]}.', j);
        end
        x0 = 0;  if ~isempty(localParams{j, 11}), x0 = double(localParams{j, 11}); end
        lab = '';  if numel(c) >= 4 && ~isempty(c{4}), lab = char(string(c{4})); end
        nm = char(string(localParams{j, 2}));
        key = lab;  if isempty(key), key = nm; end
        if ~isKey(groups, key)
            groups(key) = struct('names', {{}}, 'x0', x0, 'sigma', double(c{3}), 'from', double(c{2}), ...
                                 'estimate', "none", 'labelled', ~isempty(lab));
            order{end+1} = key; %#ok<AGROW>
        end
        g = groups(key);
        if g.sigma ~= double(c{3})
            error('gravfire:badConstraint', 'group ''%s'': members must share one sigma (row %d has %g, the group %g).', key, j, double(c{3}), g.sigma);
        end
        % WP-144 (A4): a later member's fromPass and a-priori value were
        % DROPPED -- the group kept its first member's (reproduced). fromPass
        % is refused as on the common rows; x0 too, because the pattern form
        % builds L per member and per day in the toolbox's order, so one value
        % per member cannot be written here without restating that order.
        if g.from ~= double(c{2})
            error('gravfire:badConstraint', 'group ''%s'': members must share one fromPass (row %d has %g, the group %g).', key, j, double(c{2}), g.from);
        end
        if g.x0 ~= x0
            error('gravfire:badConstraint', ['group ''%s'': members of a local group must share one a-priori ' ...
                  'value, column 10 (row %d has %g, the group %g); give them separate labels, or ' ...
                  'declare the group by L in the columns lsqvce.layout announces.'], key, j, x0, g.x0);
        end
        g.names{end+1} = regexptranslate('escape', nm);
        groups(key) = g;
    end
    for k = 1:numel(order)
        g = groups(order{k});
        est = "none";  if g.labelled, est = "scale"; end
        P(end+1) = struct('group', ['^(' strjoin(unique(g.names), '|') ')$'], 'x0', g.x0, ...
                          'sigma', g.sigma, 'from', g.from, 'name', order{k}, 'estimate', est); %#ok<AGROW>
    end
    if isempty(P), P = []; end
end
%==========================================================================
function tf = isNoneConstraint(c)
%ISNONECONSTRAINT  'none' in every spelling the tables use: [], 'none', {'none'}.
    if iscell(c) && isscalar(c), c = c{1}; end
    tf = isempty(c) || (ischar(c) && strcmpi(c, 'none')) || (isstring(c) && strcmpi(c, "none"));
end
