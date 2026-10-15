function v = toolboxVersion()
%LSQVCE.UTIL.TOOLBOXVERSION  The toolbox version string, read from +lsqvce/Contents.m.
%   v = lsqvce.util.toolboxVersion()
%   The ONE place the version is read (runMeta writes it into every Results'
%   provenance, lsqvce_setup prints it). The source is the "% Version X.Y  date"
%   line of +lsqvce/Contents.m.
%   INPUT   none
%   OUTPUT  v  char   e.g. '5.2' ('unknown' when the line is missing)
%   ACCURACY
%           Not applicable: the version is read as text from Contents.m.
%           WP-139.
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-139: the help states its
%   accuracy, "Not applicable" with its reason. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%
%   Created 2026-09-04 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
%   model Claude Fable 5.1 -- v5.2 docs level (lifted from runMeta).
    pkgDir = fileparts(fileparts(mfilename('fullpath')));    % +lsqvce/+util -> +lsqvce
    txt = fileread(fullfile(pkgDir, 'Contents.m'));
    tok = regexp(txt, '^%\s*Version\s+(\S+)', 'tokens', 'once', 'lineanchors');
    if isempty(tok)
        v = 'unknown';            % Contents.m without a Version line: honest
    else
        v = tok{1};
    end
end
