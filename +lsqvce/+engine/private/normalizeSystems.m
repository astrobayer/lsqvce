function [specs, loaders] = normalizeSystems(systems, days)
%NORMALIZESYSTEMS  Loader objects / spec structs -> BlockAccumulator specs.
%   [specs, loaders] = normalizeSystems(systems, days)
%
%   THE one normalisation of a `systems` argument (cell of Loader objects
%   from lsqvce.engine.makeLoader, or a struct array of raw specs) into the
%   spec array BlockAccumulator.addSystems takes. Shared by
%   lsqvce.engine.adjustBlock and lsqvce.engine.adjustRobust (R5 M2 round
%   2b: moved verbatim out of adjustBlock so there is ONE implementation).
%   Also returns LOADERS: a 1xN cell of the Loader objects (or [] for a raw
%   closure) so a caller can reach the Loader objects after the solve.
%
%   INPUT   systems 1xM cell | 1xM struct;  days 1xD double  fallback calendar
%   OUTPUT  specs 1xM struct (.sysID .nl .dataFcn .dayList .name);  loaders 1xM cell
%   ERRORS  lsqvce:adjustBlock:badCell / :badStruct  (identifiers kept: the
%           user called adjustBlock or adjustRobust with a malformed entry)
%
%   Created 2026-09-02 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M2, WP-81 round 2b.
    % MUST list exactly the fields ONESPEC returns: MATLAB rejects
    % specs(i) = struct(...) when the field sets differ ('Subscripted
    % assignment between dissimilar structures'), and the template is the
    % easy half to forget when a field is added. It was, on 2026-07-19.
    tmpl = struct('sysID',[], 'nl',[], 'dataFcn',[], 'dayList',[], 'name','');
    if iscell(systems)
        n = numel(systems);
        specs   = repmat(tmpl, 1, n);
        loaders = cell(1, n);
        for i = 1:n
            L = systems{i};
            if ~isLoaderObj(L)
                error('lsqvce:adjustBlock:badCell', ...
                    ['Cell entries must be Loader objects from lsqvce.engine.makeLoader. ' ...
                     'For raw closures, use the struct-array form with .sysID/.nl/.loader.']);
            end
            specs(i)   = oneSpec(L.sysID, L.nl, L.closure, ...
                                 loaderDays(L, days, i), loaderName(L));
            loaders{i} = L;
        end
    elseif isstruct(systems)
        need = {'sysID','nl','loader'};
        if ~all(isfield(systems, need))
            error('lsqvce:adjustBlock:badStruct', ...
                'systems struct array needs fields %s.', strjoin(need, ', '));
        end
        n = numel(systems);
        specs   = repmat(tmpl, 1, n);
        loaders = cell(1, n);
        for i = 1:n
            L   = systems(i).loader;
            fcn = loaderFcn(L, i);
            dl  = days;
            if isfield(systems, 'days') && ~isempty(systems(i).days), dl = systems(i).days; end
            if isLoaderObj(L), dl = loaderDays(L, dl, i); loaders{i} = L; end
            specs(i) = oneSpec(systems(i).sysID, systems(i).nl, fcn, dl, ...
                               lsqvce.util.getfielddef(systems(i), 'name', ''));
        end
    else
        error('lsqvce:adjustBlock:badSystems', ...
            'systems must be a struct array or a cell array of Loader objects.');
    end
end

function nm = loaderName(L)
%LOADERNAME  The Loader's display label, or '' if it did not set one.
    nm = '';
    if isstruct(L) && isfield(L, 'name') && ~isempty(L.name)
        nm = char(L.name);
    end
end

function s = oneSpec(sysID, nl, fcn, dayList, name)
%ONESPEC  Build one normalised system-spec struct for the block-arrow engine.
%   INPUT
%     sysID    1x1  system identifier (numeric/char/string)
%     nl       1x1  number of local parameters (nonnegative integer)
%     fcn      1x1 fh  the 4-arg loader closure @(day,obsType,arcID,opt)
%     dayList  1xD  the days this system spans (nonempty)
%     name     1xN  char, OPTIONAL. The loader's human-readable label, carried
%                   through so every message can say "demo-CD-pair" instead of
%                   the mangled key MATLAB needs internally ('x1'). Added
%                   2026-07-19: it was left out as "cosmetic", but the readable
%                   name is what makes the accumulation, residual, stochastic
%                   model and convergence blocks legible, so it is part of the
%                   contract now. Empty means "no name given, use the sysID".
%   OUTPUT
%     s  1x1 struct  fields .sysID .nl .dataFcn .dayList .name
    if isempty(dayList)
        error('lsqvce:adjustBlock:noDays', ...
            'No day list for system %s: pass DAYS, systems.days, or set it on the Loader.', ...
            string(sysID));
    end
    validateattributes(nl, {'numeric'}, {'scalar','integer','nonnegative'}, 'lsqvce.engine.adjustBlock', 'nl');
    if nargin < 5 || isempty(name), name = ''; end
    s = struct('sysID',sysID, 'nl',double(nl), 'dataFcn',fcn, ...
               'dayList',dayList, 'name',char(name));
end

function tf = isLoaderObj(L)
%ISLOADEROBJ  True if L is a makeLoader Loader object (a scalar struct carrying
%   a .closure and .isLoader == true).
%   INPUT  L  any    the candidate
%   OUTPUT tf 1x1 logical
    tf = isstruct(L) && isscalar(L) && isfield(L, 'closure') && ...
         isfield(L, 'isLoader') && isequal(L.isLoader, true);
end

function fcn = loaderFcn(L, i)
%LOADERFCN  Extract the 4-arg loader closure from a Loader object or a bare
%   function handle (errors otherwise).
%   INPUT   L  1x1  a Loader object or a function handle;  i 1x1  its index (msg)
%   OUTPUT  fcn 1x1 fh  the @(day,obsType,arcID,opt) closure
    if isLoaderObj(L)
        fcn = L.closure;
    elseif isa(L, 'function_handle')
        fcn = L;
    else
        error('lsqvce:adjustBlock:badLoader', ...
            'systems(%d).loader must be a Loader object or a 4-arg closure.', i);
    end
end

function dl = loaderDays(L, fallback, i)
%LOADERDAYS  The day list for Loader L: its own .days if present, else fallback.
%   INPUT   L 1x1 Loader;  fallback 1xD  default days;  i 1x1  index (msg)
%   OUTPUT  dl 1xD  the resolved day list (errors if both are empty)
    if isfield(L, 'days') && ~isempty(L.days)
        dl = L.days;
    else
        dl = fallback;
    end
    if isempty(dl)
        error('lsqvce:adjustBlock:noDays', ...
            'Loader %d carries no .days and no DAYS default was given.', i);
    end
end
