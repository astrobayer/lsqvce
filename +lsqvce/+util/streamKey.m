function [key, satKey] = streamKey(sysID, obsType)
%LSQVCE.UTIL.STREAMKEY  The one key of a stream: 'satKey|obsType', and its satKey.
%   [key, satKey] = lsqvce.util.streamKey(sysID, obsType)
%   [key, satKey] = lsqvce.util.streamKey(sysID)        key = satKey
%
%   A STREAM is one (system, observation type) pair. Its key names the
%   stream in every per-stream map of the joint estimator (levels, laws,
%   estimate policies, robust scales), in R.noise and in the covariance
%   override bySys.(satKey).(obsType). satKey is the system id made a legal
%   struct field name (matlab.lang.makeValidName: 1 -> 'x1'). Defined ONCE
%   here (review R7, R7-43, coding rule R1): it was composed inline in eight
%   places and recovered from the arc key by a regex. A caller that needs the
%   parts carries them; nothing parses a key back apart.
%
%   INPUT   sysID    1x1 double | char | string   the system id
%           obsType  char | string                the observation type
%                                                 (optional)
%   OUTPUT  key      1xK char   'satKey|obsType' (satKey without obsType)
%           satKey   1xS char   the system's struct-field name
%   ACCURACY
%           Not applicable: exact string formatting; test_streamUtils pins every
%           byte (testStreamAndArcKeysAreDefinedOnce). WP-139.
%   ERRORS  none of its own
%   EXAMPLE lsqvce.util.streamKey(1, 'kbr')      -> 'x1|kbr'
%           [~, s] = lsqvce.util.streamKey("C")  -> s = 'C'
%
%   See also LSQVCE.UTIL.ARCKEY, LSQVCE.UTIL.ROWWEIGHTKEY,
%            LSQVCE.ENGINE.ADJUSTVCEJOINT, LSQVCE.ENGINE.JOINTMODELOPT
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-140: the H1 line names the unit
%   in full, lsqvce.util.streamKey, as HANDOVER 6.2.1 asks (it read the bare
%   name). With assistance from Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-139: the help states its
%   accuracy, "Not applicable" with its reason. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Created 2026-10-11 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- R7 WP-115 part 4 (R7-43).
    satKey = matlab.lang.makeValidName(char(string(sysID)));
    if nargin < 2
        key = satKey;
    else
        key = [satKey '|' char(string(obsType))];
    end
end
