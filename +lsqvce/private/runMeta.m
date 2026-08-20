function meta = runMeta(seeds, options)
%RUNMETA  Compose the run-provenance struct carried by LSQVCE.RESULTS (WP-44).
%   meta = runMeta(seeds, options)
%
%   THE ONE COMPOSER (rule 3): every lsqvce.Results receives its .meta from
%   this function -- no call site builds the struct by hand, so the field
%   list can only ever change HERE (add-only across versions, never rename;
%   the plot-hook rule).
%
%   INPUT
%     seeds    1x1 struct   Every RNG seed the run RESOLVED, one field per
%                           seed (a resolved-but-unset seed is stored as its
%                           resolved value, typically []). struct() when the
%                           producer resolves none.
%     options  1x1 struct   The RESOLVED option struct the run actually used
%                           (not the user's raw input). struct() for a bare
%                           Results construction.
%
%   OUTPUT
%     meta     1x1 struct   .toolboxVersion  1xK char  from Contents.m, the
%                                            ONE version source
%                           .matlabRelease   1xK char  version('-release')
%                           .createdUTC      1xK char  ISO 8601, trailing Z
%                           .seeds           1x1 struct  as given
%                           .options         1x1 struct  as given
%
%   ERRORS
%     lsqvce:runMeta:badArg  seeds or options is not a scalar struct
%
%   See also LSQVCE.RESULTS, LSQVCE.REPORT
%
%   Created 2026-08-19 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.

    narginchk(2, 2);
    if ~isstruct(seeds) || ~isscalar(seeds)
        error('lsqvce:runMeta:badArg', 'seeds must be a 1x1 struct.');
    end
    if ~isstruct(options) || ~isscalar(options)
        error('lsqvce:runMeta:badArg', 'options must be a 1x1 struct.');
    end

    meta = struct( ...
        'toolboxVersion', toolboxVersionLine(), ...
        'matlabRelease',  version('-release'), ...
        'createdUTC',     char(datetime('now', 'TimeZone', 'UTC', ...
                               'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z''')), ...
        'seeds',          seeds, ...
        'options',        options);
end

function v = toolboxVersionLine()
% Parse the version from +lsqvce/Contents.m -- the ONE version source.
% This file lives in +lsqvce/private/, so Contents.m is one level up.
    pkgDir = fileparts(fileparts(mfilename('fullpath')));
    txt = fileread(fullfile(pkgDir, 'Contents.m'));
    tok = regexp(txt, '^%\s*Version\s+(\S+)', 'tokens', 'once', 'lineanchors');
    if isempty(tok)
        v = 'unknown';            % Contents.m without a Version line: honest
    else
        v = tok{1};
    end
end
