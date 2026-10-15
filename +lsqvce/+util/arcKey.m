function key = arcKey(sysID, obsType, day, arc)
%ARCKEY  The one key of an arc in the joint estimator: 's<sys>_<obsType>_d<day>_a<arc>'.
%   key = lsqvce.util.arcKey(sysID, obsType, day, arc)
%
%   The identity the joint estimator pins its arc factors on (WP-34:
%   out.arcKeys, Sigma2_0) and the 'arc' preset of the redundancy / stream
%   contribution group keys. Defined ONCE here (review R7, R7-43): it was
%   formatted in two places and parsed back by a regex in three, which an
%   obsType containing '_d<digits>_a' could defeat. Nothing parses it now:
%   the joint estimator carries each arc's stream key next to it
%   (out.arcStreams).
%
%   INPUT   sysID 1x1 double | char;  obsType char | string;
%           day 1x1 double;  arc 1x1 double | char
%   OUTPUT  key 1xK char
%   ACCURACY
%           Not applicable: exact string formatting; test_streamUtils pins every
%           byte (testStreamAndArcKeysAreDefinedOnce). WP-139.
%   ERRORS  none of its own
%   EXAMPLE lsqvce.util.arcKey(71, 'pos', 12, 3)   -> 's71_pos_d12_a3'
%
%   See also LSQVCE.UTIL.STREAMKEY, LSQVCE.UTIL.ROWWEIGHTKEY,
%            LSQVCE.ENGINE.ADJUSTVCEJOINT
%
%   Updated 2026-10-15 (CEST, UTC+02:00), WP-139: the help states its
%   accuracy, "Not applicable" with its reason. With assistance from Claude
%   (Anthropic), model Claude Opus 5.5.
%   Created 2026-10-11 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Opus 5.5 -- R7 WP-115 part 4 (R7-43).
    key = sprintf('s%s_%s_d%s_a%s', string(sysID), string(obsType), string(day), string(arc));
end
