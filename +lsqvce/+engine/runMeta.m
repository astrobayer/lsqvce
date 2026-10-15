function meta = runMeta(seeds, options)
%LSQVCE.ENGINE.RUNMETA  Compose the run-provenance struct carried by LSQVCE.RESULTS (WP-44).
%   meta = lsqvce.engine.runMeta(seeds, options)
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
%                           .warnings        1xW cell    every warning identifier
%                                            raised through lsqvce.util.log
%                                            during the run, first-fired order,
%                                            each once (R5 M1; fit resets the
%                                            inventory at its start)
%
%   ACCURACY
%     Not applicable: it composes strings (version, release, UTC time) and
%     passes seeds and options through unchanged -- no floating-point
%     result (test_errorDoctrine (testRunMetaContract)). WP-139.
%
%   ERRORS
%     lsqvce:runMeta:badArg  seeds or options is not a scalar struct
%
%   See also LSQVCE.RESULTS, LSQVCE.REPORT
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-140: the H1 line names the unit
%   in full, lsqvce.engine.runMeta, as HANDOVER 6.2.1 asks (it read the bare
%   name). With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-139: the help states its
%   accuracy, "Not applicable" with its reason. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Created 2026-08-19 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.
%   Moved 2026-09-01 (CEST, UTC+02:00) from +lsqvce/private to the public
%   ENGINE tier (R4-D1, WP-64a): used by USER (Results, fit) and ENGINE
%   (adjustBlock, adjustVCE), so it cannot stay package-private.

    narginchk(2, 2);
    if ~isstruct(seeds) || ~isscalar(seeds)
        error('lsqvce:runMeta:badArg', 'seeds must be a 1x1 struct.');
    end
    if ~isstruct(options) || ~isscalar(options)
        error('lsqvce:runMeta:badArg', 'options must be a 1x1 struct.');
    end

    meta = struct( ...
        'toolboxVersion', lsqvce.util.toolboxVersion(), ...
        'matlabRelease',  version('-release'), ...
        'createdUTC',     char(datetime('now', 'TimeZone', 'UTC', ...
                               'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z''')), ...
        'seeds',          seeds, ...
        'options',        options, ...
        'warnings',       {lsqvce.util.log.inventory()});   % braces: a cell VALend
end
