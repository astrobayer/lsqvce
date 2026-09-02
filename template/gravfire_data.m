function data = gravfire_data(commonData, myOptions, systemsData, opt)
%GRAVFIRE_DATA  Map a GRAVFIRE preprocessing bundle onto the lsqvce `data` struct.
%   data = gravfire_data(commonData, myOptions, systemsData)
%   data = gravfire_data(..., 'DesignFcn', @myDesign)      % your A/B builder
%   data = gravfire_data(..., 'DryRun', true)              % random A/B: pipeline check
%
%   TEMPLATE (R5 M0 close-out, HANDOVER §12.8): the GRAVFIRE bundle
%   (commonData / myOptions / systemsData, as saved in testdata.mat) already
%   contains everything the v4 data contract asks for EXCEPT the design
%   matrices. This function is the declaration half; the design half is the
%   'DesignFcn' you supply (the variational-equation partials of GRAVFIRE).
%   Every rule of DATA_CONTRACT.md is applied here once, and the file is
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
%                  B (m x info.nl), rows EPOCH-MAJOR for multi-axis streams. Default:
%                  a function that ERRORS with gravfire:designMissing (so a forgotten
%                  design step fails loudly, wrapped by lsqvce:makeLoader:buildDay).
%     'DryRun'     logical, default false   replace the design step by DETERMINISTIC
%                  random A/B (seeded per day and observation type) so that the
%                  declaration, the arc cut, Cep, nl/nc enforcement and the whole
%                  fit pipeline can be exercised BEFORE the design step exists.
%                  The numbers of a dry run mean nothing.
%   OUTPUT
%     data  1x1 struct   the v4 lsqvce declaration: commonParams, days, systems(1)
%           with id, name, buildDay, covariance, arcs, localParams, processingOptions.
%
%   ERRORS
%     gravfire:designMissing        DryRun false and no DesignFcn given (raised at read time)
%     gravfire:unknownQuantity      an observation group's observationQuantity is not
%                                   'kinematicPosition' or 'range-rate'
%     gravfire:designShape          DesignFcn returned A/B with the wrong row/column count
%     (+ every declaration check of lsqvce.fit, raised before any read; see DATA_CONTRACT §3.5)
%
%   EXAMPLE
%     S    = load('testdata.mat');
%     data = gravfire_data(S.commonData, S.myOptions, S.systemsData, 'DryRun', true);
%     R    = lsqvce.fit(data, lsqvce.Options(Model="gls", Verbose=1));   % pipeline check
%     data = gravfire_data(S.commonData, S.myOptions, S.systemsData, 'DesignFcn', @gravfireDesign);
%     R    = lsqvce.fit(data, lsqvce.Options(Model="ellmer", NumLags=2160, MaxOuter=7));
%
%   See also LSQVCE.FIT, LSQVCE.OPTIONS, DATA_CONTRACT.md, dev/GRAVFIRE_INTERFACE.md
%
%   VERIFIED 2026-09-02 on the user's testdata.mat (R2026a): the declaration
%   passes every declaration check of fit; DryRun + Model="gls" accumulates 2 days x 16
%   pieces (4 posA + 4 posB + 8 kbr), nc = 117, nl = 144 derived from the 72
%   templates, Cep accepted, in 2.0 s. Model="vce" runs (62 s, 2 passes).
%   Model="ellmer" on the 3-axis position pieces BLOCKS for minutes (H-15,
%   HANDOVER §12.3): until M1/M3, run ellmer on single-axis streams only.
%
%   Created 2026-09-02 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M0 close-out (§12.8).

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
    cov = struct();
    for g = 1:G
        cov.(obsTypes{g}) = struct('type', 'diag', 'sigma', sqrt(sd.variance(g)));
    end

    % ---- STEP 4: the local table -- ONE day's templates, all 'day'-scoped,
    %      times as seconds-of-day (the day itself is NOT in a template)
    localParams = templateFromExpanded(sd.parameter, days(1));

    % ---- STEP 5: the arc schedule = the pieces, resolvable without reading
    nPieces = sd.numberOfPieces;                              % D x G
    arcsFcn = @(day, obsType) 1:nPieces(day == days, groupIndex(obsTypes, obsType));

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
    s.id                = 1;                                  % numeric or ONE char (contract §2.3)
    s.name              = 'GRACE';                            % the label lives here
    s.buildDay          = @buildGravfireDay;
    s.covariance        = cov;
    s.arcs              = arcsFcn;
    s.localParams       = localParams;
    s.processingOptions = po;
    data.systems = s;
end

% ======================================================================
function D = buildGravfireDay(day, obsType, common, local, po, info)
%BUILDGRAVFIREDAY  One (day, observation type): raw A, B, l, t, pieces, Cep.
    g  = groupIndex(po.obsTypes, obsType);
    d  = find(po.days == day, 1);
    o  = po.bundles{d, g};

    nAx = o.numberOfAxes;
    keep = ~o.nanIndexOfPiece{1}(:);                          % flagged epochs leave
    lobs  = o.lobs{1}(keep, :);                               % nEp x nAx
    tEp   = o.time{1}(keep, 2);                               % seconds of day
    piece = o.indexOfPiece{1}(keep);                          % nEp x 1
    nEp   = size(lobs, 1);
    m     = nEp * nAx;

    % STEP 8: epoch-major rows (x1 y1 z1 x2 y2 z2 ...)
    D.l = reshape(lobs.', [], 1);
    D.t = repelem(tEp(:), nAx);

    % STEP 9: pieces -> arcIdx (row lists) + arcId (their identity)
    pieceRow = repelem(piece(:), nAx);
    ids      = unique(piece(:)).';
    D.arcIdx = arrayfun(@(k) find(pieceRow == k).', ids, 'UniformOutput', false);
    D.arcId  = ids;

    % STEP 10: multi-axis: nAxis + the per-epoch 3x3 covariance from iSpp
    if nAx > 1
        D.nAxis = nAx;
        C = o.iSpp{1}(keep, :);                               % nEp x 9 (3x3 per epoch)
        D.Cep = permute(reshape(C.', nAx, nAx, []), [3 1 2]); % nEp x 3 x 3
    end

    % STEP 11: the design matrices -- yours, or the dry-run stand-in
    if po.dryRun
        sd = mod(7919 * day + sum(double(obsType)), 2^31 - 1);
        rs = RandStream('mt19937ar', 'Seed', sd);             % deterministic per call
        D.A = randn(rs, m, info.nc);
        D.B = randn(rs, m, info.nl);
    else
        [D.A, D.B] = po.designFcn(o, day, obsType, common, local, po, info);
        if ~isequal(size(D.A), [m info.nc]) || ~isequal(size(D.B), [m info.nl])
            error('gravfire:designShape', ...
                  ['DesignFcn returned A %s / B %s; expected A [%d %d] and B [%d %d] ' ...
                   '(day %d, %s).'], mat2str(size(D.A)), mat2str(size(D.B)), ...
                   m, info.nc, m, info.nl, day, obsType);
        end
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
