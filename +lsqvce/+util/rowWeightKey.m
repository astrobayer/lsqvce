function key = rowWeightKey(sysID, day, obsType, arcId)
%ROWWEIGHTKEY  The one key convention for per-arc row weights.
%   key = lsqvce.util.rowWeightKey(sysID, day, obsType, arcId)
%
%   'sysID|day|obsType|arcId' as char, used by lsqvce.engine.adjustRobust
%   (writer) and lsqvce.engine.makeLoader (reader) for the containers.Map of
%   IRLS weight factors (R5 M2 round 2b). One definition, so the two sides
%   can never disagree on formatting (a numeric sysID prints as %g, a char
%   one verbatim; day and arcId as %g).
%
%   INPUT   sysID double|char;  day 1x1 double;  obsType char/string;
%           arcId 1x1 double|char
%   OUTPUT  key 1xK char
%   ERRORS  none of its own
%   EXAMPLE lsqvce.util.rowWeightKey(1, 59396, 'kbr', 3)  -> '1|59396|kbr|3'
%
%   See also LSQVCE.ENGINE.ADJUSTROBUST, LSQVCE.ENGINE.MAKELOADER
%
%   Created 2026-09-02 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M2, WP-81 round 2b.
    key = sprintf('%s|%s|%s|%s', asText(sysID), asText(day), char(obsType), asText(arcId));
end

function s = asText(v)
    if isnumeric(v) || islogical(v)
        s = sprintf('%g', double(v));
    else
        s = char(v);
    end
end
