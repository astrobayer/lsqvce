function [commonData, myOptions, systemsData] = fixture_gravfire_bundle(file)
%FIXTURE_GRAVFIRE_BUNDLE  A small GENERATED GRAVFIRE bundle for template/gravfire_data.m.
%   [commonData, myOptions, systemsData] = fixture_gravfire_bundle()
%   fixture_gravfire_bundle(file)       % also saves the three to a .mat file
%
%   WHY (WP-144, finding A7 of HANDOVER 12.77). template/gravfire_data.m (~400
%   lines: the common and local tables, column 11 -> priors, the epoch-major
%   rows, Cep, the cut into pieces) ran in no test: its only caller is
%   demo_gravfire, and test_examples gives it a path that cannot exist so
%   that every tier runs the same suite. The user's real bundle is
%   measurement data and never enters this public repository. This fixture
%   GENERATES a bundle with the real structure -- the shapes
%   dev/GRAVFIRE_INTERFACE.md describes and the template reads -- so the whole
%   template runs on every tier. Its numbers are random and mean nothing.
%
%   THE BUNDLE (two MJD days 59396-59397, three observation groups):
%     commonData.parameter   6 x 14 cell: {'c'|'s', degree, order, ...}; column
%                            10 the a-priori value, column 11 the constraint --
%                            rows 1-2 one labelled Tikhonov group 'stokes'
%                            (estimate "scale"), row 3 an unlabelled one, the
%                            rest 'none'
%     myOptions              .days (1 x 2) and .observationGroups (1 x 3 cell):
%                            GraceA / GraceB kinematicPosition, GraceX range-rate
%     systemsData{1}         .observations (2 x 3 cell of per-day structs:
%                            numberOfAxes, lobs, time, indexOfPiece,
%                            nanIndexOfPiece, iSpp), .parameter -- the EXPANDED
%                            local table, 2 days x 2 satellites x 3 unknowns, with
%                            [MJD sec] spans and the day first in the flags; the
%                            two accscale rows one labelled group 'accscale' --
%                            .numberOfPieces (2 x 3), .variance (1 x 3)
%   Each (day, group) has 2 pieces of 30 epochs (5 s), epoch 7 of piece 2 flagged.
%
%   INPUT   file  1xK char/string  optional: save commonData, myOptions,
%                                  systemsData there (as demo_gravfire loads)
%   OUTPUT  the three bundle variables
%
%   EXAMPLE
%     [c, m, s] = fixture_gravfire_bundle();
%     data = gravfire_data(c, m, s, 'DryRun', true);
%     R = lsqvce.engine.fitDeclaration(data, Verbose=0);
%
%   See also TEMPLATE/GRAVFIRE_DATA, DEMO_GRAVFIRE, FIXTURE_KBR_L1B
%
%   Created 2026-10-17 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- WP-144 (A7).
    days = [59396 59397];
    rs = RandStream('mt19937ar', 'Seed', 144);

    % ---- the common table: 6 Stokes rows, 14 columns -------------------------
    P = cell(6, 14);
    nm = {'c', 2, 0; 'c', 2, 1; 's', 2, 1; 'c', 2, 2; 's', 2, 2; 'c', 3, 0};
    for j = 1:6
        P(j, 1:3) = nm(j, :);
        P{j, 4} = 'even';  P{j, 10} = [];  P{j, 11} = 'none';  P{j, 13} = 1e-9;
    end
    P{1, 11} = {'Tikhonov', 1, 1e-3, 'stokes'};
    P{2, 11} = {'Tikhonov', 1, 1e-3, 'stokes'};
    P{3, 11} = {'Tikhonov', 1, 1e-2};
    P{1, 10} = 0.5;                                  % an a-priori value on one member
    commonData = struct('parameter', {P});

    % ---- the options ---------------------------------------------------------
    myOptions = struct('days', days, 'observationGroups', {{ ...
        struct('name', 'GraceA', 'observationQuantity', 'kinematicPosition'), ...
        struct('name', 'GraceB', 'observationQuantity', 'kinematicPosition'), ...
        struct('name', 'GraceX', 'observationQuantity', 'range-rate')}});

    % ---- the expanded local table: per day, per satellite, 3 unknowns --------
    unk = {'x0', 'y0', 'accscale'};
    L = cell(0, 14);
    for d = days
        for sat = {'GraceA', 'GraceB'}
            for u = 1:numel(unk)
                r = cell(1, 14);
                r{1} = sprintf('%s_%s', sat{1}, unk{u});
                r{2} = [d 0];  r{3} = [d 86390];  r{4} = [];
                r{5} = 1;      r{6} = [d 1 1 1 0];
                r{10} = [];    r{11} = 'none';    r{13} = 1;
                if strcmp(unk{u}, 'accscale')
                    r{11} = {'Tikhonov', 1, 1e-6, 'accscale'};
                end
                L(end+1, :) = r; %#ok<AGROW>
            end
        end
    end

    % ---- the observations: 2 pieces x 30 epochs per (day, group) -------------
    nPiece = 2;  nPer = 30;  dt = 5;
    nAx = [3 3 1];
    obs = cell(numel(days), 3);
    for di = 1:numel(days)
        for g = 1:3
            nEp = nPiece * nPer;
            t = (0:nEp-1).' * dt + 3600 * di;
            o = struct();
            o.numberOfAxes    = nAx(g);
            o.lobs            = {randn(rs, nEp, nAx(g))};
            o.time            = {[days(di) * ones(nEp, 1), t]};
            o.indexOfPiece    = {repelem((1:nPiece).', nPer)};
            nanIdx = false(nEp, 1);  nanIdx(nPer + 7) = true;
            o.nanIndexOfPiece = {nanIdx};
            if nAx(g) == 3
                o.iSpp = {repmat(reshape(diag([1 2 3] * 1e-4), 1, 9), nEp, 1)};
            end
            obs{di, g} = o;
        end
    end
    systemsData = {struct('observations', {obs}, 'parameter', {L}, ...
                          'numberOfPieces', nPiece * ones(numel(days), 3), ...
                          'variance', [0.01 0.01 1e-12])};

    if nargin >= 1 && strlength(string(file)) > 0
        save(char(file), 'commonData', 'myOptions', 'systemsData');
    end
end
