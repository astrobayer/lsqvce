function key = rowWeightKey(sysID, day, obsType, arcId)
%LSQVCE.UTIL.ROWWEIGHTKEY  The one key convention for per-arc row weights.
%   key = lsqvce.util.rowWeightKey(sysID, day, obsType, arcId)
%
%   'sysID|day|obsType|arcId' as char, used by lsqvce.engine.adjustRobust
%   (writer) and lsqvce.engine.makeLoader (reader) for the containers.Map of
%   IRLS weight factors (R5 M2 round 2b). One definition, so the two sides
%   can never disagree on formatting: an integer-valued number prints with
%   %d (exact up to 2^53, the format the loader's day cache uses), any other
%   number with %.17g (round-trip exact), a char/string verbatim.
%
%   INPUT   sysID double|char;  day 1x1 double;  obsType char/string;
%           arcId 1x1 double|char
%   OUTPUT  key 1xK char
%   ACCURACY
%           Not applicable: exact string formatting -- integers print with
%           %d, other numbers with %.17g (round-trip exact; R7-06: %g kept
%           six digits and merged day ids). WP-139.
%   ERRORS  none of its own
%   EXAMPLE lsqvce.util.rowWeightKey(1, 59396, 'kbr', 3)     -> '1|59396|kbr|3'
%           lsqvce.util.rowWeightKey(1, 20260901, 'kbr', 3)  -> '1|20260901|kbr|3'
%
%   See also LSQVCE.ENGINE.ADJUSTROBUST, LSQVCE.ENGINE.MAKELOADER
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-140: the H1 line names the unit
%   in full, lsqvce.util.rowWeightKey, as HANDOVER 6.2.1 asks (it read the
%   bare name). With assistance from Claude (Anthropic), model Claude Opus
%   5.5.
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-139: the help states its
%   accuracy, "Not applicable" with its reason. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-115 (finding R7-06): numbers
%   were printed with %g, which keeps six significant digits, so day ids of
%   seven or more digits (YYYYMMDD, Julian Day Numbers) collapsed onto one
%   key and one day's robust factors were applied to another. With
%   assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Created 2026-09-02 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.1 -- R5 Level M2, WP-81 round 2b.
    key = sprintf('%s|%s|%s|%s', asText(sysID), asText(day), char(obsType), asText(arcId));
end

function s = asText(v)
    if isnumeric(v) || islogical(v)
        v = double(v);
        if v == round(v) && abs(v) < 2^53
            s = sprintf('%d', v);
        else
            s = sprintf('%.17g', v);
        end
    else
        s = char(v);
    end
end
