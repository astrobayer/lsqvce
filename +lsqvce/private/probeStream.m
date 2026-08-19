function rs = probeStream()
%PROBESTREAM  The package's PRIVATE probe RandStream, derived once.
%   rs = probeStream()
%
%   Returns the ONE persistent RandStream all default-seeded probe draws
%   come from (WP-36, fixes F-06). It is created ONCE per MATLAB session,
%   independent of -- and never advancing -- the GLOBAL stream, which is
%   what the probe-consumer help blocks promise: calling any lsqvce
%   estimator with an empty seed must leave rng state untouched.
%
%   INPUT
%     (none)
%
%   OUTPUT
%     rs   1x1 RandStream   'threefry4x64_20', seeded by 'shuffle' at first
%                           use (per-session entropy; counter-based, so it
%                           is cheap and state-robust). The SAME object on
%                           every later call in the session.
%
%   ERRORS
%     (none)
%
%   EXAMPLE
%     s = randi(probeStream(), [0 intmax('int32')], 1, 1);  % pin a seed
%
%   See also PROBEMATRIX, RANDSTREAM
%
%   Created 2026-08-19 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.

    persistent st
    if isempty(st) || ~isvalid(st)
        st = RandStream('threefry4x64_20', 'Seed', 'shuffle');
    end
    rs = st;
end
