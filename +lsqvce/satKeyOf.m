function key = satKeyOf(satId)
%LSQVCE.SATKEYOF  Valid struct field name for a satellite identifier.
%   key = lsqvce.satKeyOf(satId)
%
%   The single place that defines the satellite-key convention shared by the
%   GRACE-FO example drivers ITERATESEPARABLE (which WRITES
%   temporalModel.(key).(obsType)) and GRACELOADER (which READS it) and by
%   LSQVCE.RESOLVESATSPECS, so they can never drift.
%
%   INPUT
%     satId   1x1 | 1xk   char | string | numeric   Satellite identifier, e.g.
%                                                    'C', "D", or 1. Required.
%
%   OUTPUT
%     key     1xm         char   A legal MATLAB struct field name: "s_" + satId
%                                 passed through matlab.lang.makeValidName
%                                 ('C' -> 's_C', 1 -> 's_1').
%
%   See also LSQVCE.RESOLVESATSPECS, LSQVCE.NORMALIZESATLIST
%
%   Created 2026-07-15 13:40 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.
    narginchk(1, 1);
    key = char(matlab.lang.makeValidName("s_" + string(satId)));
end
