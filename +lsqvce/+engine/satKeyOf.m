function key = satKeyOf(satId)
%LSQVCE.ENGINE.SATKEYOF  Valid struct field name for a satellite identifier.
%   key = lsqvce.engine.satKeyOf(satId)
%
%   The single place that defines the 's_' + satId satellite-key convention,
%   so every caller that writes a per-satellite field and every caller that
%   reads it form the same name.
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
%   See also MATLAB.LANG.MAKEVALIDNAME
%
%   Updated 2026-10-10 (CEST, UTC+02:00), R7 WP-118 part 1 (R7-34): the
%   references to drivers that exist nowhere removed. With assistance from
%   Claude (Anthropic), model Claude Opus 5.5.
%   Created 2026-07-15 13:40 (CEST, UTC+02:00) with assistance from
%   Claude (Anthropic), model Claude Opus 4.8.
%   Moved 2026-09-01 (CEST, UTC+02:00) into +lsqvce/+engine (R4 decision R4-D1,
%   WP-64a): the ENGINE tier leaves the USER surface; every qualified
%   reference was rewritten mechanically by docs_src/refactor_tiers.py.
    narginchk(1, 1);
    key = char(matlab.lang.makeValidName("s_" + string(satId)));
end
