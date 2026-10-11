function muteWarn(tc, id)
%MUTEWARN  Silence one warning identifier for the duration of one test.
%   muteWarn(tc, id)
%
%   THE one way a test silences a warning it deliberately provokes (WP-47,
%   external handover 2.5): `warning off <id>` now, restore via
%   tc.addTeardown -- cleanup-based, so a FAILING assertion cannot leave the
%   session's warning state altered for later suites (the source project
%   lost six unrelated suites to exactly that). Do NOT call this where the
%   identifier cannot fire: an unnecessary suppression teaches the next
%   reader that the identifier is unavoidable.
%
%   INPUT
%     tc   1x1 matlab.unittest.TestCase   the running test case
%     id   1xK char/string                the warning identifier to silence
%
%   OUTPUT  none (session warning state altered until teardown)
%
%   ERRORS
%     Input:
%       lsqvce:muteWarn:badId   id is not a nonempty char row / string scalar
%
%   EXAMPLE
%     muteWarn(tc, 'lsqvce:snoop:maxReject');
%     [sol, D, REJ] = lsqvce.engine.verdict.snoop(...);        % provokes it on purpose
%
%   See also RUN_ALL_LSQVCE, GATEVERDICT
%
%   Created 2026-08-20 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5.

    narginchk(2, 2);
    if isstring(id) && isscalar(id), id = char(id); end
    if ~ischar(id) || isempty(id) || size(id, 1) ~= 1
        error('lsqvce:muteWarn:badId', ...
              'id must be a nonempty char row or string scalar.');
    end
    prev = warning('query', id);
    warning('off', id);
    tc.addTeardown(@() warning(prev.state, id));
end
