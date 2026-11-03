function [cfg, L, data] = fixture_gravity()
%FIXTURE_GRAVITY  ENGINE-tier loaders over the demo_gravity_dataLoader data.
%
%   [cfg, L, data] = fixture_gravity()
%
%   The physically structured gravity fixture for the ENGINE-tier tests
%   (test_physOracle): since L10 (WP-74) the example returns the USER-tier
%   declaration only ([cfg, data], v3 contract) and this fixture builds the
%   Loaders from it -- through LSQVCE.ENGINE.DECLARATION, the step LSQVCE.FIT
%   performs internally, where engine calls are legitimate.
%
%   INPUT   (none; examples/ must be on the path -- the callers add it)
%   OUTPUT
%     cfg   1x1 struct  the demo's configuration and planted truth
%     L     1x3 cell    Loaders, one per system (LSQVCE.ENGINE.MAKELOADER)
%     data  1x1 struct  the loader's compiled declaration, passed through
%
%   ERRORS  none of its own; makeLoader construction errors pass through.
%
%   Updated 2026-11-02 (CET, UTC+01:00), WP-161: the Loaders come from
%   lsqvce.engine.declaration (the loader's task compiles to a struct whose
%   generated buildDay reads the commons' plan). With assistance from
%   Claude (Anthropic), model Claude Opus 5.5.
%   Updated 2026-09-09 (CEST, UTC+02:00), R6 WP-102 (v6): 'Arcs' is passed
%   only when the declaration carries a schedule. With assistance from
%   Claude (Anthropic), model Claude Fable 5.1.
%   Created 2026-08-31 (CEST, UTC+02:00) with assistance from Claude
%   (Anthropic), model Claude Fable 5 -- Level L10, WP-74.
    [cfg, data] = demo_gravity_dataLoader();  data = data.declaration();   % WP-162: the compiled struct
    % WP-161: the loader declares a lsqvce.TaskDefinition; its compiled struct
    % reaches the Loaders through the ONE declaration step fit uses (with the
    % commons' plan, which a hand-written makeLoader call does not build)
    [~, L] = lsqvce.engine.declaration(data);
end
