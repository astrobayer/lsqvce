function outPath = showExport(figH, name, showSeconds, exportDir)
%LSQVCE.UTIL.SHOWEXPORT  Show a figure briefly, then export it as a PNG.
%
%   outPath = LSQVCE.UTIL.SHOWEXPORT(figH, name, showSeconds, exportDir)
%
%   The engine-side half of the adjustVCE plot hooks: a user handle only
%   BUILDS and RETURNS a figure; this helper owns the timed show and the PNG
%   export, so every hook behaves consistently and the caller's plotting code
%   stays free of boilerplate.
%
%   WHAT IT DOES
%     1. If figH is a valid figure, brings it up and (when showSeconds > 0)
%        pauses for showSeconds so the user can see it, then continues.
%     2. If exportDir is non-empty, writes <exportDir>/<name>.png (creating the
%        folder if needed) and returns the full path.
%     3. Closes the figure it was given (the hook created it; the engine owns its
%        lifetime), UNLESS showSeconds < 0, which means "leave it open".
%
%   It is defensive: a headless session (no display) or a hook that returned []
%   or an invalid handle is a NO-OP that returns '' rather than erroring -- a
%   plot hook must never be able to break the solve.
%
%   INPUT
%     figH        1x1   handle   a figure handle returned by a plot hook (or []).
%     name        1xk   char     base file name (no extension); also the export
%                                stem. Sanitised for the filesystem.
%     showSeconds 1x1   double   seconds to display before continuing. 0 = do not
%                                pause (just export). < 0 = show and DO NOT close
%                                (pause abs(showSeconds), leave the figure open).
%     exportDir   1xk   char     folder for the PNG. '' or [] -> do not export.
%
%   OUTPUT
%     outPath     1xk   char     full path of the written PNG, or '' if nothing
%                                was exported.
%
%   EXAMPLE
%     f = figure; plot(randn(100,1));
%     p = lsqvce.util.showExport(f, 'postSolve_x1_posC_pass1', 2, '/tmp/figs');
%
%   IDENTIFIERS ALSO RAISED (R5 M1, documented = raised; audit check 20)
%     lsqvce:showExport:exportFailed             [warning] Could not export figure 'X' (X)
%
%   See also LSQVCE.ENGINE.ADJUSTVCE
%
%   Created 2026-07-16 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.

    outPath = '';
    if nargin < 3 || isempty(showSeconds), showSeconds = 0; end
    if nargin < 4, exportDir = ''; end

    % A hook that returned nothing useful -> silent no-op (never break the solve).
    if isempty(figH) || ~allValidFig(figH)
        return
    end
    figH = figH(1);

    keepOpen = showSeconds < 0;
    waitSec  = abs(showSeconds);

    % ---- show ------------------------------------------------------------
    try
        figure(figH);                 % bring to front
        drawnow;
        if waitSec > 0, pause(waitSec); end
    catch ME
        % headless / no display: the show is skipped and SAID (R5 M1); the
        % export below still runs.
        lsqvce.util.log.msg('[showExport: display unavailable (%s); export only]\n', ...
                            ME.identifier);
    end

    % ---- export ----------------------------------------------------------
    if ~isempty(exportDir)
        try
            if exist(exportDir, 'dir') ~= 7, mkdir(exportDir); end
            fname   = [sanitise(name), '.png'];
            outPath = fullfile(exportDir, fname);
            exportOnePng(figH, outPath);
        catch ME
            lsqvce.util.log.warn('lsqvce:showExport:exportFailed', ...
                'Could not export figure ''%s'' (%s).', name, ME.identifier);
            outPath = '';
        end
    end

    % ---- close -----------------------------------------------------------
    if ~keepOpen
        try, close(figH); catch, end
    end
end

% ==========================================================================
function tf = allValidFig(h)
%ALLVALIDFIG  True if h is a graphics handle that ishghandle and is a figure.
    % ishghandle answers false for non-handles, so no exception is possible
    % before the && short-circuits: a CONDITION, not a catch (R5 M1).
    tf = ~isempty(h) && all(ishghandle(h(:))) && all(strcmp(get(h(:), 'type'), 'figure'));
end

function exportOnePng(figH, outPath)
%EXPORTONEPNG  Write a PNG, preferring exportgraphics (R2020a+) then print.
    if exist('exportgraphics', 'file') == 2 || exist('exportgraphics', 'builtin') == 5
        exportgraphics(figH, outPath, 'Resolution', 150);
    else
        print(figH, outPath, '-dpng', '-r150');   % fallback for older MATLAB
    end
end

function s = sanitise(name)
%SANITISE  Make a filesystem-safe stem: keep alnum, dash, underscore, dot; the
%   stream key 'satKey|obsType' has a '|' which becomes '_'.
    s = char(string(name));
    s = regexprep(s, '[^A-Za-z0-9._-]', '_');
    if isempty(s), s = 'figure'; end
end
